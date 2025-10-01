#include "Flux/COMPFlux.hh"

COMPFlux::COMPFlux() {
    name = "gamma";
    GetParams();
    BuildCDF();
}

void COMPFlux::BuildCDF() {
    const double k = (2.0 - alpha) / E_Peak;
    if (k <= 0.0) {
        G4Exception("COMPFlux::BuildCDF", "BAD_PARAM",
                    FatalException, "k <= 0 (alpha >= 2, распределение не определено)");
    }

    auto pdf = [&](double E) {
        return std::pow(E, -alpha) * std::exp(-k * E);
    };

    const int N = 2000;
    energyGrid.resize(N);
    cdfGrid.resize(N);

    double logEmin = std::log(Emin);
    double logEmax = std::log(Emax);

    std::vector<double> integrand(N);
    for (int i = 0; i < N; ++i) {
        double frac = double(i) / (N - 1);
        double Ei = std::exp(logEmin + frac * (logEmax - logEmin));
        energyGrid[i] = Ei;
        integrand[i] = pdf(Ei);
    }

    std::vector<double> cum(N, 0.0);
    for (int i = 1; i < N; ++i) {
        double dE = energyGrid[i] - energyGrid[i - 1];
        cum[i] = cum[i - 1] + 0.5 * (integrand[i] + integrand[i - 1]) * dE;
    }

    double Z = cum.back();
    if (!(Z > 0.0)) {
        G4Exception("COMPFlux::BuildCDF", "BAD_NORM",
                    FatalException, "Нормировка не посчиталась (Z <= 0)");
    }

    for (int i = 0; i < N; ++i) {
        cdfGrid[i] = cum[i] / Z;
    }

    cdfGrid.front() = 0.0;
    cdfGrid.back() = 1.0;
}

void COMPFlux::GetParams() {
    const std::string filepath = "../Flux_config/COMP_params.txt";
    std::ifstream paramFile(filepath);
    if (!paramFile.is_open()) {
        G4Exception("COMPFlux::GetParams", "FILE_OPEN_FAIL",
                    JustWarning, ("Cannot open " + filepath).c_str());
        alpha = -1.0;
        E_Peak = 500.0 * MeV;
        Emin = 10.0 * keV;
        Emax = 10000.0 * MeV;
        paramFile.close();
        return;
    }

    std::string line;
    alpha = MAXFLOAT;
    E_Peak = MAXFLOAT;
    Emin = MAXFLOAT;
    Emax = MAXFLOAT;

    while (std::getline(paramFile, line)) {
        if (line.find("alpha") != std::string::npos) {
            alpha = std::stod(line.substr(line.find(':') + 1));
        } else if (line.find("E_Peak") != std::string::npos) {
            E_Peak = std::stod(line.substr(line.find(':') + 1)) * MeV;
        } else if (line.find("E_min") != std::string::npos) {
            Emin = std::stod(line.substr(line.find(':') + 1)) * MeV;
        } else if (line.find("E_max") != std::string::npos) {
            Emax = std::stod(line.substr(line.find(':') + 1)) * MeV;
        }

        if (alpha != MAXFLOAT && E_Peak != MAXFLOAT && Emin != MAXFLOAT && Emax != MAXFLOAT) {
            paramFile.close();
            return;
        }
    }
    G4Exception("COMPFlux::GetParams", "POOR_CONTENT",
                JustWarning, ("Cannot find all parameters in file " + filepath).c_str());
    alpha = -1.0;
    E_Peak = 1.18511 * MeV;
    Emin = 10.0 * keV;
    Emax = 10000.0 * MeV;
    paramFile.close();
}

double COMPFlux::SampleEnergy() {
    double u = G4UniformRand();

    const auto it = std::lower_bound(cdfGrid.begin(), cdfGrid.end(), u);
    if (it == cdfGrid.begin()) return energyGrid.front();
    if (it == cdfGrid.end()) return energyGrid.back();

    const int idx = std::distance(cdfGrid.begin(), it);
    const double u1 = cdfGrid[idx - 1];
    const double u2 = cdfGrid[idx];
    const double E1 = energyGrid[idx - 1];
    double E2 = energyGrid[idx];

    const double t = (u - u1) / (u2 - u1);
    return E1 + t * (E2 - E1);
}
