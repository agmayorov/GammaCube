#include "Flux/RadiationBelt.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "Flux/GeoMagnetic.hh"

namespace RadiationBelt {
    namespace {
        std::string Trim(const std::string& s) {
            const size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return "";
            return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
        }

        std::vector<double> ReadAxis(const std::string& text, const std::string& key, const std::string& path) {
            std::istringstream in(text);
            size_t n = 0;
            if (!(in >> n) || n < 2) throw std::runtime_error(path + ": ось " + key + " задана неверно");
            std::vector<double> axis(n);
            for (auto& v : axis)
                if (!(in >> v)) throw std::runtime_error(path + ": в оси " + key + " меньше " + std::to_string(n) + " значений");
            for (size_t i = 1; i < n; ++i)
                if (!(axis[i] > axis[i - 1])) throw std::runtime_error(path + ": ось " + key + " не возрастает");
            return axis;
        }

        double Real(const std::map<std::string, std::string>& p, const std::string& key) {
            const auto it = p.find(key);
            if (it == p.end()) throw std::runtime_error("В Belt_params.txt нет параметра " + key);
            std::istringstream in(it->second);
            double v = 0.;
            if (!(in >> v) || !(in >> std::ws).eof())
                throw std::runtime_error("Belt_params.txt: " + key + " — не число: " + it->second);
            return v;
        }

        double Real(const std::map<std::string, std::string>& p, const std::string& key, const double fallback) {
            return p.count(key) ? Real(p, key) : fallback;
        }

        std::string Text(const std::map<std::string, std::string>& p, const std::string& key, const std::string& fallback) {
            const auto it = p.find(key);
            return it == p.end() ? fallback : it->second;
        }

        double OmniFlux(const std::vector<double>& E, const std::vector<double>& J, const double lo, const double hi) {
            if (!(hi > lo)) return 0.;
            auto at = [&](const double e) {
                if (e <= E.front()) return J.front();
                if (e >= E.back()) return J.back();
                const size_t j = std::upper_bound(E.begin(), E.end(), e) - E.begin();
                const double f0 = J[j - 1], f1 = J[j];
                const double t = (e - E[j - 1]) / (E[j] - E[j - 1]);
                if (f0 > 0. && f1 > 0.) return f0 * std::pow(f1 / f0, std::log(e / E[j - 1]) / std::log(E[j] / E[j - 1]));
                return f0 + t * (f1 - f0);
            };
            const int n = 4000;
            double sum = 0.;
            for (int k = 0; k < n; ++k) {
                const double e0 = lo + (hi - lo) * k / n, e1 = lo + (hi - lo) * (k + 1) / n;
                sum += 0.5 * (at(e0) + at(e1)) * (e1 - e0);
            }
            return 4. * M_PI * sum;
        }

        std::string Fixed(const double v, const int digits) {
            std::ostringstream s;
            s << std::fixed << std::setprecision(digits) << v;
            std::string out = s.str();
            if (out.find('.') != std::string::npos) {
                out.erase(out.find_last_not_of('0') + 1);
                if (out.back() == '.') out.pop_back();
            }
            return out;
        }
    }

    std::map<std::string, std::string> ReadParams(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Не удалось открыть " + path);
        std::map<std::string, std::string> params;
        std::string line;
        while (std::getline(in, line)) {
            line = line.substr(0, line.find('#'));
            const size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            const std::string key = Trim(line.substr(0, colon)), value = Trim(line.substr(colon + 1));
            if (!key.empty() && !value.empty()) params[key] = value;
        }
        return params;
    }

    Model::Model(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("Не удалось открыть модель пояса " + path + " (её пишет Processing/convert_belt_mat.py)");

        std::string line;
        bool data = false;
        while (!data && std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            const size_t colon = line.find(':');
            if (colon == std::string::npos) continue;
            const std::string key = Trim(line.substr(0, colon)), value = line.substr(colon + 1);
            if (key == "model") name = Trim(value);
            else if (key == "particle") particle = Trim(value);
            else if (key == "alpha_eq_deg") alpha = ReadAxis(value, key, path);
            else if (key == "energy_MeV") energy = ReadAxis(value, key, path);
            else if (key == "L") L = ReadAxis(value, key, path);
            else if (key == "J") data = true;
        }
        if (!data || alpha.empty() || energy.empty() || L.empty() || particle.empty())
            throw std::runtime_error(path + ": нет particle, осей или блока J");

        J.resize(alpha.size() * energy.size() * L.size());
        for (auto& v : J)
            if (!(in >> v)) throw std::runtime_error(path + ": в блоке J меньше значений, чем задают оси");
        for (auto& v : J) v = std::max(0., v);
    }

    std::vector<double> Model::MeanIntensity(const double Lvalue, const double BoverB0) const {
        const size_t nE = energy.size();
        std::vector<double> mean(nE, 0.);
        if (!(Lvalue >= L.front() && Lvalue <= L.back()) || !(BoverB0 > 0.)) return mean;

        const size_t il = std::min<size_t>(std::upper_bound(L.begin(), L.end(), Lvalue) - L.begin() - 1, L.size() - 2);
        const double w = (Lvalue - L[il]) / (L[il + 1] - L[il]);
        std::vector<double> onShell(alpha.size() * nE);
        for (size_t ia = 0; ia < alpha.size(); ++ia)
            for (size_t ie = 0; ie < nE; ++ie) {
                const double a = At(ia, ie, il), b = At(ia, ie, il + 1);
                onShell[ia * nE + ie] = a > 0. && b > 0. ? a * std::pow(b / a, w) : a + w * (b - a);
            }

        const int steps = 900;
        const double dAlpha = 0.5 * M_PI / steps;
        for (int k = 0; k < steps; ++k) {
            const double a = (k + 0.5) * dAlpha;
            const double s = std::sin(a) / std::sqrt(BoverB0);
            const double aeq = s >= 1. ? 90. : std::asin(s) * 180. / M_PI;
            if (aeq < alpha.front()) continue;
            const size_t ja = std::min<size_t>(std::upper_bound(alpha.begin(), alpha.end(), aeq) - alpha.begin() - 1,
                                               alpha.size() - 2);
            const double t = std::clamp((aeq - alpha[ja]) / (alpha[ja + 1] - alpha[ja]), 0., 1.);
            const double weight = std::sin(a) * dAlpha;
            for (size_t ie = 0; ie < nE; ++ie)
                mean[ie] += weight * ((1. - t) * onShell[ja * nE + ie] + t * onShell[(ja + 1) * nE + ie]);
        }
        return mean;
    }

    Result Prepare(const std::string& configPath, const double EminMeV, const double EmaxMeV, const Overrides& manual) {
        const auto p = ReadParams(configPath);
        Result r;

        r.model = manual.model.empty() ? Text(p, "model", "") : manual.model;
        if (r.model.empty()) throw std::runtime_error("Не задана модель пояса (model в Belt_params.txt или --belt-model)");
        const std::string modelPath = r.model.find('/') != std::string::npos || r.model.size() > 4 &&
                                      r.model.substr(r.model.size() - 4) == ".txt"
                                          ? r.model
                                          : "../RadiationBelts/" + r.model + ".txt";
        const Model model(modelPath);
        r.model = model.Name();
        r.particle = model.Particle();
        r.modelEmin = model.Energies().front();
        r.modelEmax = model.Energies().back();
        r.dipoleMoment = Real(p, "dipole_moment_nT", GeoMagnetic::kDipoleMomentAP8);
        r.position = Text(p, "position", "point");
        if ((manual.B > 0. || manual.BoverB0 > 0.) && !(manual.L > 0.))
            throw std::runtime_error("--belt-B и --belt-BB0 задаются вместе с --belt-L");
        if (manual.L > 0.) r.position = "point";

        const double lo = std::max(EminMeV, r.modelEmin), hi = std::min(EmaxMeV, r.modelEmax);
        const auto& E = model.Energies();
        std::vector<double> spectrum;

        auto bOverB0 = [&](const double L, const double B) { return B * L * L * L / r.dipoleMoment; };

        if (r.position == "point") {
            r.L = manual.L > 0. ? manual.L : Real(p, "L");
            double B = -1., ratio = -1.;
            if (manual.B > 0. || manual.BoverB0 > 0.) {
                B = manual.B;
                ratio = manual.BoverB0;
            } else {
                B = Real(p, "B_nT", -1.);
                ratio = Real(p, "B_B0", -1.);
            }
            if (B > 0. && ratio > 0.) throw std::runtime_error("Belt: задайте поле одним способом — B (нТ) или B/B0");
            if (!(r.L > 0.) || !(B > 0. || ratio > 0.))
                throw std::runtime_error("Belt, position: point — нужны L > 0 и B_nT > 0 (или B_B0 > 0)");
            const bool byRatio = ratio > 0.;
            r.BoverB0 = byRatio ? ratio : bOverB0(r.L, B);
            r.B = byRatio ? ratio * r.dipoleMoment / (r.L * r.L * r.L) : B;
            spectrum = model.MeanIntensity(r.L, r.BoverB0);
            r.description = byRatio ? "L and B/B0 given" : "L and B given";
            r.tag = r.model + "_L" + Fixed(r.L, 3) + (byRatio ? "_BB0-" + Fixed(r.BoverB0, 3) : "_B" + Fixed(r.B, 0) + "nT");
        } else if (r.position == "geo" || r.position == "orbit") {
            r.year = GeoMagnetic::FractionalYear(Text(p, "date", "2026-01-01"));
            const GeoMagnetic::IGRF field(Text(p, "igrf_path", "../RadiationBelts/IGRF14.shc"), r.year);
            const double maxL = model.Lmax() + 1.;

            if (r.position == "geo") {
                const double lat = Real(p, "latitude_deg"), lon = Real(p, "longitude_deg"), alt = Real(p, "altitude_km");
                const auto shell = GeoMagnetic::McIlwainL(field, GeoMagnetic::GeodeticToECEF(lat, lon, alt),
                                                          r.dipoleMoment, maxL);
                r.L = shell.L;
                r.B = shell.B;
                r.BoverB0 = std::isfinite(r.L) ? bOverB0(r.L, r.B) : 0.;
                spectrum = model.MeanIntensity(r.L, r.BoverB0);
                r.description = "lat = " + Fixed(lat, 3) + " deg, lon = " + Fixed(lon, 3) + " deg, h = " +
                                Fixed(alt, 1) + " km, date = " + Fixed(r.year, 3);
                r.tag = r.model + "_geo_lat" + Fixed(lat, 2) + "_lon" + Fixed(lon, 2) + "_h" + Fixed(alt, 0);
            } else {
                const double alt = Real(p, "altitude_km"), incl = Real(p, "inclination_deg");
                const double days = Real(p, "orbit_days", 1.), node = Real(p, "node_longitude_deg", 0.);
                const int n = static_cast<int>(Real(p, "orbit_points", 2000.));
                if (n < 1 || !(days > 0.)) throw std::runtime_error("Belt_params.txt: orbit_points ≥ 1, orbit_days > 0");
                r.points = n;

                const double radius = GeoMagnetic::kEquatorialRadius + alt;
                const double motion = std::sqrt(GeoMagnetic::kGM / (radius * radius * radius));
                const double i = incl * M_PI / 180., omega = node * M_PI / 180.;
                const double duration = days * 86400.;

                struct Sample {
                    double t, lat, lon, L, B, bb0, omni;
                };
                std::vector<Sample> samples(n);
                std::vector<std::vector<double>> sums;
                std::vector<int> trapped;

                const unsigned workers = std::max(1u, std::min<unsigned>(std::thread::hardware_concurrency(), n));
                sums.assign(workers, std::vector<double>(E.size(), 0.));
                trapped.assign(workers, 0);
                std::vector<std::thread> pool;
                std::vector<std::string> errors(workers);
                for (unsigned w = 0; w < workers; ++w) {
                    pool.emplace_back([&, w] {
                        try {
                            for (int k = static_cast<int>(w); k < n; k += static_cast<int>(workers)) {
                                const double t = (k + 0.5) * duration / n;
                                const double u = motion * t, earth = GeoMagnetic::kEarthRotation * t;
                                const double xi = radius * (std::cos(omega) * std::cos(u) - std::sin(omega) * std::sin(u) * std::cos(i));
                                const double yi = radius * (std::sin(omega) * std::cos(u) + std::cos(omega) * std::sin(u) * std::cos(i));
                                const double zi = radius * std::sin(u) * std::sin(i);
                                const G4ThreeVector x(xi * std::cos(earth) + yi * std::sin(earth),
                                                      -xi * std::sin(earth) + yi * std::cos(earth), zi);

                                const auto shell = GeoMagnetic::McIlwainL(field, x, r.dipoleMoment, maxL);
                                const double bb0 = std::isfinite(shell.L) ? bOverB0(shell.L, shell.B) : 0.;
                                const auto J = model.MeanIntensity(shell.L, bb0);
                                for (size_t e = 0; e < E.size(); ++e) sums[w][e] += J[e];
                                const double omni = OmniFlux(E, J, lo, hi);
                                if (omni > 0.) ++trapped[w];
                                samples[k] = {t, std::asin(zi / radius) * 180. / M_PI,
                                              std::atan2(x.y(), x.x()) * 180. / M_PI, shell.L, shell.B, bb0, omni};
                            }
                        } catch (const std::exception& ex) {
                            errors[w] = ex.what();
                        }
                    });
                }
                for (auto& th : pool) th.join();
                for (const auto& e : errors)
                    if (!e.empty()) throw std::runtime_error(e);

                spectrum.assign(E.size(), 0.);
                int nTrapped = 0;
                for (unsigned w = 0; w < workers; ++w) {
                    for (size_t e = 0; e < E.size(); ++e) spectrum[e] += sums[w][e] / n;
                    nTrapped += trapped[w];
                }
                r.trappedFraction = static_cast<double>(nTrapped) / n;

                r.description = "circular, h = " + Fixed(alt, 1) + " km, i = " + Fixed(incl, 2) + " deg, node lon = " +
                                Fixed(node, 2) + " deg, " + Fixed(days, 3) + " days, " + std::to_string(n) +
                                " points, date = " + Fixed(r.year, 3);
                r.tag = r.model + "_orbit_h" + Fixed(alt, 0) + "_i" + Fixed(incl, 1);
                r.profilePath = "belt_" + r.tag + "_orbit.csv";

                std::ofstream prof(r.profilePath);
                if (!prof) throw std::runtime_error("Не удалось записать " + r.profilePath);
                prof << "t_s,lat_deg,lon_deg,B_nT,L,B_over_B0,omni_flux_cm2s\n";
                prof << std::setprecision(7);
                for (const auto& s : samples)
                    prof << s.t << ',' << s.lat << ',' << s.lon << ',' << s.B << ','
                         << (std::isfinite(s.L) ? s.L : -1.) << ',' << s.bb0 << ',' << s.omni << '\n';
            }
        } else {
            throw std::runtime_error("Belt_params.txt: position = " + r.position + "; допустимо point, geo, orbit");
        }

        r.omniFlux = OmniFlux(E, spectrum, lo, hi);
        if (!(r.omniFlux > 0.)) {
            std::ostringstream msg;
            msg << "Поток " << r.model << " в заданном положении равен нулю в " << lo << "–" << hi << " МэВ";
            if (r.position != "orbit")
                msg << " (L = " << r.L << ", B/B0 = " << r.BoverB0 << ", модель задана для L " << model.Lmin()
                    << "–" << model.Lmax() << ")";
            throw std::runtime_error(msg.str());
        }

        r.spectrumPath = "belt_" + r.tag + ".csv";
        std::ofstream out(r.spectrumPath);
        if (!out) throw std::runtime_error("Не удалось записать " + r.spectrumPath);
        out << "# E, <J>: MeV, cm^-2 s^-1 sr^-1 MeV^-1 (" << r.model << ", " << r.description << ")\n";
        out << std::setprecision(8);
        for (size_t e = 0; e < E.size(); ++e) out << E[e] << ',' << spectrum[e] << '\n';
        return r;
    }
}
