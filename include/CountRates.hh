#ifndef COUNTRATES_HH
#define COUNTRATES_HH

#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

enum class FluxType { PLAW, COMP, SEP, UNIFORM, GALACTIC };

struct EnergyRange {
    double Emin;
    double Emax;
};

struct FluxParams {
    // PLAW / COMP
    double A = 0.0;
    double alpha = 0.0;
    double E_piv = 1.0;
    double E_peak = 1.0; // only COMP

    // SEP
    int sep_year = 0;
    int sep_order = 0;
    std::string sep_csv_path; // path to CSV with coefficients

    // Galactic
    double phiMV = 600.0;
    std::string particle = "proton";
};

struct RateCounts {
    int crystalOnly = 0; // N_det (Crystal && !Veto)
    int crystalAndVeto = 0; // N_det (Crystal && Veto)
};

struct RateResult {
    double area = 0.0;
    double integral = 0.0; // ∫ flux(E) dE
    double Ndot = 0.0; // A_eff * integral
    double rateCrystal = 0.0; // crystalOnly / (N / Ndot)
    double rateBoth = 0.0; // (crystalOnly+crystalAndVeto) / (N / Ndot)
};

double fluxPLAW(double E, double A, double alpha, double E_piv);

double fluxCOMP(double E, double A, double alpha, double E_piv, double E_peak);

double fluxSEP(double E, int year, int order, const std::string &csvPath);

double fluxUniform(double E);

double fluxGalactic(double E);

double J_proton(double E_GeV);

enum class FluxDir { Vertical, Isotropic };

double effectiveArea_cm2(double R_mm, double H_mm, FluxDir dir);

double integrateAdaptiveSimpson(const std::function<double(double)> &f,
                                double a, double b,
                                double rel_tol = 1e-6, int max_depth = 20);


RateResult computeRate(FluxType type,
                       const FluxParams &p,
                       EnergyRange eRange,
                       double A_eff_cm2,
                       int N_histories,
                       const RateCounts &detCounts);

#endif //COUNTRATES_HH
