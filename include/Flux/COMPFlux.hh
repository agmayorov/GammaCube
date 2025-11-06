#ifndef COMPFLUX_HH
#define COMPFLUX_HH

#include "Flux/Flux.hh"


class COMPFlux : public Flux {
public:
    COMPFlux();

private:
    G4double alpha{};
    G4double E_Peak{};

    std::vector<double> energyGrid;
    std::vector<double> cdfGrid;

    void BuildCDF();

    G4double SampleEnergy() override;
};

#endif //COMPFLUX_HH
