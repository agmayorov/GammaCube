#include "Flux/PLAWFlux.hh"


PLAWFlux::PLAWFlux(double minE, double maxE, double a)
    : alpha(a) {
    Emin = minE;
    Emax = maxE;
    name = "gamma";
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
