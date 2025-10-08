#include "Flux/GalacticFlux.hh"

GalacticFlux::GalacticFlux() {
    GetParams();
    if (!G4ParticleTable::GetParticleTable()->FindParticle(name)) {
        G4Exception("GalacticFlux::GalacticFlux",
                    "InvalidSetup", FatalException,
                    ("Particle " + name + " not found in table").c_str());
    }

    BuildCDF();
}


void GalacticFlux::GetParams() {
    const std::string filepath = "../Flux_config/Galactic_params.txt";
    std::ifstream paramFile(filepath);
    if (!paramFile.is_open()) {
        G4Exception("GalacticFlux::GetParams", "FILE_OPEN_FAIL",
                    JustWarning, ("Cannot open " + filepath).c_str());
        name = "proton";
        phiMV = 600;
        Emin = 0.001 * GeV;
        Emax = 1000 * GeV;
        paramFile.close();
        return;
    }
    std::string line;
    name = "";
    phiMV = 0;
    Emin = MAXFLOAT;
    Emax = MAXFLOAT;
    while (std::getline(paramFile, line)) {
        if (line.find("particle") != std::string::npos) {
            name = line.substr(line.find(':') + 2);
        } else if (line.find("phiMV") != std::string::npos) {
            phiMV = std::stod(line.substr(line.find(':') + 2));
        } else if (line.find("E_min") != std::string::npos) {
            Emin = std::stod(line.substr(line.find(':') + 2)) * MeV / GeV;
        } else if (line.find("E_max") != std::string::npos) {
            Emax = std::stod(line.substr(line.find(':') + 2)) * MeV / GeV;
        }

        if (!name.empty() && phiMV != 0 && Emin != MAXFLOAT && Emax != MAXFLOAT) {
            paramFile.close();
            return;
        }
    }
    G4Exception("GalacticFlux::GetParams", "POOR_CONTENT",
                JustWarning, ("Cannot find correct values in file " + filepath).c_str());
    name = "proton";
    phiMV = 600;
    Emin = 0.001 * GeV;
    Emax = 1000 * GeV;
    paramFile.close();
}

void GalacticFlux::BuildCDF() {
    constexpr G4int NBins = 1000;
    energyGrid.resize(NBins);
    cdfGrid.resize(NBins);

    const G4double logEmin = std::log(Emin);
    const G4double logEmax = std::log(Emax);

    for (int i = 0; i < NBins; i++) {
        const G4double lg = logEmin + (logEmax - logEmin) * i / (NBins - 1);
        energyGrid[i] = std::exp(lg);
    }

    G4double integral = 0.0;
    cdfGrid[0] = 0.0;
    for (int i = 1; i < NBins; i++) {
        const G4double E1 = energyGrid[i - 1];
        const G4double E2 = energyGrid[i];
        const G4double f1 = J_TOA_GeV(E1);
        const G4double f2 = J_TOA_GeV(E2);
        integral += 0.5 * (f1 + f2) * (E2 - E1);
        cdfGrid[i] = integral;
    }

    if (integral <= 0.0 || !std::isfinite(integral)) {
        G4Exception("GalacticFlux::BuildCDF", "BAD_INTEGRAL", FatalException,
                    "Integral of spectrum is non-positive or non-finite.");
    }
    for (int i = 0; i < NBins; i++) {
        cdfGrid[i] /= integral;
    }
    cdfGrid.back() = 1.0;
}


G4double GalacticFlux::J_Proton(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double mp = 0.938272;

    const G4double ETot = E + mp;
    G4double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) {
        beta = 1e-6;
    }
    G4double term1 = 2620.0 / (beta * beta) * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.98) + std::pow(0.7, 0.98)) / (1.0 + std::pow(0.7, 0.98)), -4.0);

    const G4double term2 = 30.0 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 8.0) / 9.0, -12.0);

    return term1 + term2;
}


G4double GalacticFlux::J_Electron(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double me = 0.000511;

    const G4double ETot = E + me;
    G4double beta_sq = 1.0 - me * me / (ETot * ETot);
    if (beta_sq <= 0) {
        beta_sq = 1e-6;
    }

    G4double term1 = 255.0 / beta_sq * std::pow(E / E0, -1.0);
    term1 *= std::pow((E / E0 + 0.63) / 1.63, -2.43);

    const G4double term2 = 6.4 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 15.0) / 16.0, -26.0);

    return term1 + term2;
}


G4double GalacticFlux::J_Positron(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double mp = 0.000511;

    const G4double ETot = E + mp;
    G4double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) beta = 1e-6;

    G4double term1 = 25.0 / (beta * beta) * std::pow(E / E0, 0.1);
    term1 *= std::pow((std::pow(E / E0, 1.1) + std::pow(0.2, 1.1)) / (1.0 + std::pow(0.2, 1.1)), -3.31);

    const G4double term2 = 23.0 * std::pow(E / E0, 0.5) * std::pow((E / E0 + 2.2) / 3.2, -9.5);

    return term1 + term2;
}


G4double GalacticFlux::J_Alpha(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double malpha = 3.727379;

    const G4double ETot = E + malpha;
    G4double beta_sq = 1.0 - malpha * malpha / (ETot * ETot);
    if (beta_sq <= 0) beta_sq = 1e-6;

    G4double term1 = 163.4 / beta_sq * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.97) + std::pow(0.58, 0.97)) / (1.0 + std::pow(0.58, 0.97)), -4.0);

    return term1;
}

G4double GalacticFlux::J_TOA_GeV(const G4double E) const {
    constexpr G4double mp = 0.938272;
    constexpr G4double me = 0.000511;
    constexpr G4double malpha = 3.727379;
    G4double mass = 0;
    const G4int Z = name == "alpha" ? 2 : 1;
    const G4double phiGeV = phiMV * 1e-3 * Z;

    const G4double ELis = E + phiGeV;
    if (ELis <= 0) return 0.0;
    if (name == "alpha") {
        mass = malpha;
    } else if (name == "e-" or name == "e+") {
        mass = me;;
    } else if (name == "proton") {
        mass = mp;;
    }
    const G4double num = E * (E + 2.0 * mass);
    const G4double den = ELis * (ELis + 2.0 * mass);
    if (den <= 0) {
        return 0.0;
    }

    G4double J_LIS = 0;
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

G4double GalacticFlux::SampleEnergy() {
    const G4double u = G4UniformRand(); // равномерное [0,1)

    const auto it = std::lower_bound(cdfGrid.begin(), cdfGrid.end(), u);
    const int idx = std::max(1, static_cast<int>(it - cdfGrid.begin()));

    const G4double u1 = cdfGrid[idx - 1];
    const G4double u2 = cdfGrid[idx];
    const G4double E1 = energyGrid[idx - 1];
    const G4double E2 = energyGrid[idx];

    const G4double t = (u - u1) / (u2 - u1);
    return (E1 + t * (E2 - E1)) * GeV;
}
