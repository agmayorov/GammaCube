#include "Flux/GalacticFlux.hh"

GalacticFlux::GalacticFlux(const G4double EminGeV, const G4double EmaxGeV, const G4double pMV)
    : phiMV(pMV) {
    Emin = EminGeV / GeV;
    Emax = EmaxGeV / GeV;
    name = "proton";

    particleDef = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if (!particleDef) {
        G4Exception("GalacticFlux::GalacticFlux",
                    "InvalidSetup", FatalException,
                    ("Particle " + name + " not found in table").c_str());
    }

    BuildCDF();
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

    for (int i = 0; i < NBins; i++) {
        cdfGrid[i] /= integral;
    }
}


G4double GalacticFlux::J_LIS(const G4double E) const {
    constexpr G4double E0 = 1.0;
    constexpr G4double mp = 0.938272;

    const G4double ETot = E + mp;
    G4double beta = std::sqrt(1.0 - mp * mp / (ETot * ETot));
    if (beta <= 0) beta = 1e-6;

    G4double term1 = 2620.0 / (beta * beta) * std::pow(E / E0, 1.1);
    term1 *= std::pow((std::pow(E / E0, 0.98) + std::pow(0.7, 0.98)) / (1.0 + std::pow(0.7, 0.98)), -4.0);

    const G4double term2 = 30.0 * std::pow(E / E0, 2.0) * std::pow((E / E0 + 8.0) / 9.0, -12.0);

    return term1 + term2;
}

G4double GalacticFlux::J_TOA_GeV(const G4double E) const {
    constexpr G4double mp = 0.938272;
    const G4double phiGeV = phiMV * 1e-3;

    const G4double ELis = E + phiGeV;
    if (ELis <= 0) return 0.0;

    const G4double num = E * (E + 2.0 * mp);
    const G4double den = ELis * (ELis + 2.0 * mp);
    if (den <= 0) {
        return 0.0;
    }

    return num / den * J_LIS(ELis);
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
