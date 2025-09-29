#include "Flux/UniformFlux.hh"


void UniformFlux::PickParticle() {
    if (const G4double r = G4UniformRand(); r < 0.25) {
        name = "gamma";
        Emin = 0.001 * MeV;
        Emax = 500 * MeV;
    } else if (r < 0.50) {
        name = "e-";
        Emin = 0.001 * MeV;
        Emax = 100. * MeV;
    } else if (r < 0.75) {
        name = "proton";
        Emin = 1. * MeV;
        Emax = 10000. * MeV;
    } else {
        name = "alpha";
        Emin = 10. * MeV;
        Emax = 10000. * MeV;
    }
}


double UniformFlux::SampleEnergy() {
    PickParticle();
    return Emin * std::pow(Emax / Emin, G4UniformRand());
}
