#include "Flux/PLAWFlux.hh"


PLAWFlux::PLAWFlux() {
    name = "gamma";
    GetParams();
}


void PLAWFlux::GetParams() {
    const std::string filepath = "../Flux_config/PLAW_params.txt";
    std::ifstream paramFile(filepath);
    if (!paramFile.is_open()) {
        G4Exception("PLAWFlux::GetParams", "FILE_OPEN_FAIL",
                    JustWarning, ("Cannot open " + filepath).c_str());
        alpha = 1.411103;
        Emin = 0.01 * MeV;
        Emax = 100 * MeV;
        paramFile.close();
        return;
    }
    std::string line;
    alpha = MAXFLOAT;
    Emin = MAXFLOAT;
    Emax = MAXFLOAT;
    while (std::getline(paramFile, line)) {
        if (line.find("alpha") != std::string::npos) {
            alpha = std::stod(line.substr(line.find(':') + 1));
        } else if (line.find("E_min") != std::string::npos) {
            Emin = std::stod(line.substr(line.find(':') + 1)) * MeV;
        } else if (line.find("E_max") != std::string::npos) {
            Emax = std::stod(line.substr(line.find(':') + 1)) * MeV;
        }

        if (alpha != MAXFLOAT && Emin != MAXFLOAT && Emax != MAXFLOAT) {
            paramFile.close();
            return;
        }
    }
    G4Exception("PLAWFlux::GetParams", "POOR_CONTENT",
                JustWarning, ("Cannot find value of alpha in file " + filepath).c_str());
    alpha = 1.411103;
    Emin = 0.01 * MeV;
    Emax = 100 * MeV;
    paramFile.close();
}


double PLAWFlux::SampleEnergy() {
    double u = G4UniformRand();
    if (std::abs(alpha - 1.0) < 1e-12) {
        return Emin * std::pow(Emax / Emin, u);
    }
    double EminPow = std::pow(Emin, 1.0 - alpha);
    double EmaxPow = std::pow(Emax, 1.0 - alpha);
    double val = EminPow + u * (EmaxPow - EminPow);
    return std::pow(val, 1.0 / (1.0 - alpha));
}
