#include "Flux/TableFlux.hh"


namespace {
    std::string trim(const std::string &s) {
        const auto a = s.find_first_not_of(" \t\r\n");
        if (a == std::string::npos) return "";
        const auto b = s.find_last_not_of(" \t\r\n");
        return s.substr(a, b - a + 1);
    }

    std::vector<G4double> extract_numbers(const std::string &line) {
        static const std::regex re(R"(([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?))");
        std::vector<G4double> out;
        for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end; ++it) {
            out.push_back(std::stod((*it)[1].str()));
        }
        return out;
    }

    constexpr G4double NOT_SET = std::numeric_limits<G4double>::quiet_NaN();

    bool is_set(const G4double x) { return std::isfinite(x); }
}

TableFlux::TableFlux() {
    name = "proton";
    GetParams();
    BuildCDF();
}

void TableFlux::GetParams() {
    path.clear();

    Emin = NOT_SET;
    Emax = NOT_SET;

    const std::string filepath = "../Flux_config/Table_Params.txt";
    std::ifstream fin(filepath);
    if (!fin.is_open()) {
        G4Exception("TableFlux::GetParams", "FILE_OPEN_FAIL",
                    JustWarning, ("Cannot open " + filepath + ", using fallbacks.").c_str());
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        auto pos = line.find(':');
        if (pos == std::string::npos) continue;
        const std::string key = trim(line.substr(0, pos));
        const std::string val = trim(line.substr(pos + 1));

        if (key == "path") {
            path = val;
        } else if (key == "Emin") {
            if (!val.empty()) {
                try { Emin = std::stod(val) * MeV; } catch (...) { Emin = NOT_SET; }
            }
        } else if (key == "Emax") {
            if (!val.empty()) {
                try { Emax = std::stod(val) * MeV; } catch (...) { Emax = NOT_SET; }
            }
        }
    }
    fin.close();
}

void TableFlux::BuildCDF() {
    EList.clear();
    CDF.clear();

    if (path.empty()) {
        G4Exception("TableFlux::BuildCDF", "NO_PATH",
                    JustWarning, "CSV path is empty. Using trivial 2-point spectrum.");
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::ifstream in(path.c_str());
    if (!in) {
        G4Exception("TableFlux::BuildCDF", "CSV_OPEN_FAIL",
                    JustWarning, ("Cannot open " + path + ", using trivial spectrum.").c_str());
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::vector<Row> rows;
    rows.reserve(2048);

    std::string line;
    while (std::getline(in, line)) {
        const auto nums = extract_numbers(line);
        if (nums.size() < 2) continue;
        const G4double E_G4 = nums[0] * MeV;
        const G4double flx = nums[1];
        if (!(std::isfinite(E_G4) && std::isfinite(flx))) continue;
        rows.push_back({E_G4, flx});
    }
    in.close();

    if (rows.size() < 2) {
        G4Exception("TableFlux::BuildCDF", "CSV_NO_ROWS",
                    JustWarning, "Not enough data rows (need >=2). Using trivial spectrum.");
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row &a, const Row &b) { return a.E_MeV < b.E_MeV; });

    {
        std::vector<Row> uniq;
        uniq.reserve(rows.size());
        for (const auto &r: rows) {
            if (!uniq.empty() && std::abs(r.E_MeV - uniq.back().E_MeV) <= 1e-12 * MeV) continue;
            uniq.push_back(r);
        }
        rows.swap(uniq);
    }

    G4double dataEmin = rows.front().E_MeV;
    G4double dataEmax = rows.back().E_MeV;

    const bool useBounds =
        is_set(Emin) || is_set(Emax);

    G4double lo = dataEmin;
    G4double hi = dataEmax;

    if (useBounds) {
        if (is_set(Emin)) lo = std::max(Emin, dataEmin);
        if (is_set(Emax)) hi = std::min(Emax, dataEmax);
        if (!(hi > lo)) {
            G4Exception("TableFlux::BuildCDF", "CSV_RANGE_INVALID",
                        JustWarning, "Invalid Emin/Emax vs data range, falling back to data bounds.");
            lo = dataEmin;
            hi = dataEmax;
        }
    } else {
        Emin = dataEmin;
        Emax = dataEmax;
    }

    std::vector<G4double> Es;
    std::vector<G4double> lEs;
    std::vector<G4double> f_sub;

    Es.reserve(rows.size());
    f_sub.reserve(rows.size());

    auto interp_at = [&](const G4double E) -> G4double {
        const auto it = std::upper_bound(rows.begin(), rows.end(), E,
                                   [](const G4double x, const Row &r) { return x < r.E_MeV; });
        if (it == rows.begin()) return rows.front().flux;
        if (it == rows.end()) return rows.back().flux;
        const size_t j = static_cast<size_t>(it - rows.begin());
        const size_t i = j - 1;
        const G4double x0 = rows[i].E_MeV, x1 = rows[j].E_MeV;
        const G4double y0 = rows[i].flux, y1 = rows[j].flux;
        const G4double t = (E - x0) / std::max(x1 - x0, 1e-300);
        return y0 + t * (y1 - y0);
    };

    if (lo > rows.front().E_MeV) {
        Es.push_back(lo);
        f_sub.push_back(std::max(0.0, interp_at(lo)));
    }

    for (const auto &[E, flux]: rows) {
        if (E < lo || E > hi) continue;
        Es.push_back(E);
        f_sub.push_back(std::max(0.0, flux));
    }

    if (hi < rows.back().E_MeV && (Es.empty() || Es.back() < hi - 1e-12 * MeV)) {
        Es.push_back(hi);
        f_sub.push_back(std::max(0.0, interp_at(hi)));
    }

    if (Es.size() < 2) {
        G4Exception("TableFlux::BuildCDF", "CSV_RANGE_TOO_NARROW",
                    JustWarning, "Energy range too narrow (fewer than 2 points). Using trivial spectrum.");
        EList = {1. * MeV, 10. * MeV};
        CDF = {0.0, 1.0};
        return;
    }

    lEs.resize(Es.size());
    for (size_t i = 0; i < Es.size(); ++i) lEs[i] = std::log(Es[i]);

    EList = std::vector<G4double>(Es.begin(), Es.end());
    CDF.assign(EList.size(), 0.0);

    long double acc = 0.0L;
    for (size_t i = 0; i + 1 < EList.size(); ++i) {
        const long double gi = static_cast<long double>(f_sub[i]) * static_cast<long double>(EList[i]);
        const long double gip1 = static_cast<long double>(f_sub[i + 1]) * static_cast<long double>(EList[i + 1]);
        const long double dlnE = lEs[i + 1] - lEs[i];
        const long double wi = 0.5L * (gi + gip1) * dlnE;
        acc += wi;
        CDF[i + 1] = static_cast<G4double>(acc);
    }

    if (acc <= 0.0L || !std::isfinite(static_cast<G4double>(acc))) {
        for (size_t i = 0; i < CDF.size(); ++i) {
            CDF[i] = static_cast<G4double>(i) / static_cast<G4double>(CDF.size() - 1);
        }
        return;
    }

    const G4double norm = static_cast<G4double>(acc);
    for (auto &v: CDF) v /= norm;

    CDF.front() = 0.0;
    CDF.back() = 1.0;
}

G4double TableFlux::SampleEnergy() {
    if (EList.size() < 2) return 1.0 * MeV;

    const G4double u = G4UniformRand();

    const auto it = std::upper_bound(CDF.begin(), CDF.end(), u);
    if (it == CDF.begin()) return EList.front();
    if (it == CDF.end()) return EList.back();

    const size_t j = static_cast<size_t>(std::distance(CDF.begin(), it));
    const size_t i = j - 1;

    const G4double C0 = CDF[i];
    const G4double C1 = CDF[j];
    const G4double t = (u - C0) / std::max(C1 - C0, 1e-12);

    const G4double lnE0 = std::log(EList[i]);
    const G4double lnE1 = std::log(EList[j]);
    const G4double lnEu = lnE0 + t * (lnE1 - lnE0);
    return std::exp(lnEu);
}
