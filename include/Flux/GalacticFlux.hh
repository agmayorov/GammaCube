#ifndef GALACTICFLUX_HH
#define GALACTICFLUX_HH

#include "Flux/Flux.hh"


class GalacticFlux : public Flux {
public:
    GalacticFlux(G4double, G4double, G4double);

protected:
    G4double SampleEnergy() override;

private:
    void BuildCDF();

    [[nodiscard]] G4double J_LIS(G4double E) const;
    [[nodiscard]] G4double J_TOA_GeV(G4double E) const;

    G4double phiMV;
    G4ParticleDefinition* particleDef{};

    std::vector<G4double> energyGrid;
    std::vector<G4double> cdfGrid;
};

#endif //GALACTICFLUX_HH
