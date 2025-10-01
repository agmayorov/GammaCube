#ifndef GALACTICFLUX_HH
#define GALACTICFLUX_HH

#include "Flux/Flux.hh"


class GalacticFlux : public Flux {
public:
    GalacticFlux();

protected:
    G4double SampleEnergy() override;

private:
    G4double phiMV{};

    std::vector<G4double> energyGrid;
    std::vector<G4double> cdfGrid;

    void GetParams();
    void BuildCDF();

    [[nodiscard]] G4double J_Proton(G4double E) const;
    [[nodiscard]] G4double J_Electron(G4double E) const;
    [[nodiscard]] G4double J_Positron(G4double E) const;
    [[nodiscard]] G4double J_Alpha(G4double E) const;
    [[nodiscard]] G4double J_TOA_GeV(G4double E) const;
};

#endif //GALACTICFLUX_HH
