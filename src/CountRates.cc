#include "CountRates.hh"


double fluxPLAW(const double E, const double A, double const alpha, const double E_piv) {
    return A * std::pow(E / E_piv, -alpha);
}

double fluxCOMP(const double E, const double A, const double alpha, const double E_piv, const double E_peak) {
    return A * std::pow(E / E_piv, -alpha) * std::exp((alpha - 2.0) * (E / E_peak));
}

// SEP: read CSV
static std::vector<double> readSepRow(int year, int order, const std::string &csvPath) {
    std::ifstream in(csvPath);
    if (!in.is_open()) {
        throw std::runtime_error("SEP: cannot open coefficients CSV: " + csvPath);
    }
    std::string line;
    std::getline(in, line);

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> cols;
        while (std::getline(ss, cell, ',')) {
            try {
                cols.push_back(std::stod(cell));
            } catch (...) {
                cols.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        if (cols.size() < 3) continue;

        int y = static_cast<int>(cols[0]);
        int k = static_cast<int>(cols[1]);
        if (y == year && k == order) {
            cols.pop_back();
            return cols;
        }
    }
    throw std::runtime_error("SEP: row not found for year=" + std::to_string(year) +
                             " order=" + std::to_string(order));
}

double fluxSEP(const double E, const int year, const int order, const std::string &csvPath) {
    static int cached_year = 0, cached_order = 0;
    static std::string cached_path;
    static std::vector<double> coeffs;

    if (year != cached_year || order != cached_order || csvPath != cached_path || coeffs.empty()) {
        coeffs = readSepRow(year, order, csvPath);
        cached_year = year;
        cached_order = order;
        cached_path = csvPath;
    }

    if (coeffs.size() < 4) return 0.0;
    const double x = std::log10(E);
    double poly = 0.0;
    for (size_t i = 2; i < cached_order + 3; i++) {
        poly += coeffs[i] * std::pow(x, static_cast<int>(i - 2));
    }
    return std::pow(10.0, poly);
}

// --- Uniform ---
double fluxUniform(double) { return 1.0; }


// --- Galactic ---
double J_Proton(const double E_GeV) {
    constexpr double E0 = 1.0;
    constexpr double mp = 0.938272;

    const double ETot = E_GeV + mp;
    double beta_sq = 1.0 - mp * mp / (ETot * ETot);
    if (beta_sq <= 0.0) beta_sq = 1e-12;

    double term1 = 2620.0 / beta_sq * std::pow(E_GeV / E0, 1.1);
    term1 *= std::pow((std::pow(E_GeV / E0, 0.98) + std::pow(0.7, 0.98)) / (1.0 + std::pow(0.7, 0.98)), -4.0);

    const double term2 = 30.0 * std::pow(E_GeV / E0, 2.0) * std::pow((E_GeV / E0 + 8.0) / 9.0, -12.0);

    return term1 + term2;
}


double J_Electron(const double E) {
    constexpr double E0 = 1.0;
    constexpr double me = 0.000511;

    const double ETot = E + me;
    double beta_sq = 1.0 - me * me / (ETot * ETot);
    if (beta_sq <= 0) {
        beta_sq = 1e-6;
    }

    double term1 = 255.0 / beta_sq * std::pow(E / E0, -1.0);
    term1 *= std::pow((E / E0 + 0.63) / 1.63, -2.43);

    const double term2 = 6.4 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 15.0) / 16.0, -26.0);

    return term1 + term2;
}

double J_Positron(const double E) {
    constexpr double E0 = 1.0;
    constexpr double mp = 0.000511;

    const double ETot = E + mp;
    double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) beta = 1e-6;

    double term1 = 25.0 / (beta * beta) * std::pow(E / E0, 0.1);
    term1 *= std::pow((std::pow(E / E0, 1.1) + std::pow(0.2, 1.1)) / (1.0 + std::pow(0.2, 1.1)), -3.31);

    const double term2 = 23.0 * std::pow(E / E0, 0.5) * std::pow((E / E0 + 2.2) / 3.2, -9.5);

    return term1 + term2;
}

double J_Alpha(const double E) {
    constexpr double E0 = 1.0;
    constexpr double malpha = 3.727379;

    const double ETot = E + malpha;
    double beta_sq = 1.0 - malpha * malpha / (ETot * ETot);
    if (beta_sq <= 0) beta_sq = 1e-6;

    double term1 = 163.4 / beta_sq * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.97) + std::pow(0.58, 0.97)) / (1.0 + std::pow(0.58, 0.97)), -4.0);

    return term1;
}

double fluxGalactic(const double E, const double phiMV, const std::string &name) {
    double mass = 0;
    const int Z = name == "alpha" ? 2 : 1;
    const double phiGV = phiMV * 1e-3 * Z;

    const double ELis = E + phiGV;
    if (ELis <= 0) return 0.0;
    if (name == "alpha") {
        mass = 3.727379;
    } else if (name == "e-" or name == "e+") {
        mass = 0.000511;
    } else if (name == "proton") {
        mass = 0.938272;
    }
    const double num = E * (E + 2.0 * mass);
    const double den = ELis * (ELis + 2.0 * mass);
    if (den <= 0) {
        return 0.0;
    }

    double J_LIS = 0;
    if (name == "alpha") {
        J_LIS = J_Alpha(ELis);
    } else if (name == "e-") {
        J_LIS = J_Electron(ELis);
    } else if (name == "proton") {
        J_LIS = J_Proton(ELis);
    } else if (name == "e+") {
        J_LIS = J_Positron(ELis);
    }
    return num / den * J_LIS;
}

// ---------------- Area ----------------

double effectiveArea_cm2(const double R_mm, const double H_mm, const FluxDir dir) {
    const double R_cm = R_mm / 10.0;
    const double H_cm = H_mm / 10.0;

    if (dir == FluxDir::Vertical) {
        return 4.0 * R_cm * R_cm;
    }
    const double val = std::sqrt(R_cm * R_cm + H_cm * H_cm) + 0.5;
    return 2.0 * M_PI * M_PI * val * val;
}

// ---------------- Integral ----------------

static double simpson(const std::function<double(double)> &f, const double a, const double b) {
    const double c = 0.5 * (a + b);
    return (b - a) * (f(a) + 4.0 * f(c) + f(b)) / 6.0;
}

static double adaptiveSimpsonRec(const std::function<double(double)> &f,
                                 const double a, const double b, const double eps,
                                 const double whole, const int depth) {
    const double c = 0.5 * (a + b);
    const double left = simpson(f, a, c);
    const double right = simpson(f, c, b);
    const double delta = left + right - whole;
    if (depth <= 0 || std::fabs(delta) <= 15.0 * eps) {
        return left + right + delta / 15.0;
    }
    return adaptiveSimpsonRec(f, a, c, eps / 2.0, left, depth - 1) +
           adaptiveSimpsonRec(f, c, b, eps / 2.0, right, depth - 1);
}

double integrateAdaptiveSimpson(const std::function<double(double)> &f,
                                const double a, const double b,
                                const double rel_tol, const int max_depth) {
    const double initial = simpson(f, a, b);
    const double eps = rel_tol * std::max(1.0, std::fabs(initial));
    return adaptiveSimpsonRec(f, a, b, eps, initial, max_depth);
}


RateResult computeRate(const FluxType type,
                       const FluxParams &p,
                       EnergyRange eRange,
                       double A_eff_cm2,
                       const int N_histories,
                       const RateCounts &detCounts) {
    std::function<double(double)> f;

    switch (type) {
        case FluxType::PLAW:
            f = [=](const double E) { return fluxPLAW(E, p.A, p.alpha, p.E_piv); };
            break;
        case FluxType::COMP:
            f = [=](const double E) { return fluxCOMP(E, p.A, p.alpha, p.E_piv, p.E_peak); };
            break;
        case FluxType::SEP:
            f = [=](const double E) { return fluxSEP(E, p.sep_year, p.sep_order, p.sep_csv_path); };
            break;
        case FluxType::UNIFORM:
            f = [=](const double E) { return fluxUniform(E); };
            break;
        case FluxType::GALACTIC: {
            eRange.Emin /= 1000.0;
            eRange.Emax /= 1000.0;
            A_eff_cm2 /= 10000.0;
            f = [=](const double E_GeV) { return fluxGalactic(E_GeV, p.phiMV, p.particle); };
            break;
        }
        default:
            throw std::runtime_error("Unknown flux type");
    }

    const double integral = integrateAdaptiveSimpson(f, eRange.Emin, eRange.Emax, 1e-6, 22);
    const double Ndot = A_eff_cm2 * integral;

    RateResult R;
    R.area = A_eff_cm2;
    R.integral = integral;
    R.Ndot = Ndot;
    R.rateCrystal = N_histories > 0 ? (detCounts.crystalOnly + 0.0) * Ndot / N_histories : 0.0;
    const int bothDet = detCounts.crystalOnly + detCounts.crystalAndVeto;
    R.rateBoth = N_histories > 0 ? (bothDet + 0.0) * Ndot / N_histories : 0.0;
    return R;
}
