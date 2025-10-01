#ifndef COMPFLUX_HH
#define COMPFLUX_HH

#include "Flux/Flux.hh"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>


class COMPFlux : public Flux {
public:
    COMPFlux();

protected:
    G4double SampleEnergy() override;

private:
    G4double alpha{};
    G4double E_Peak{};

    std::vector<double> energyGrid;
    std::vector<double> cdfGrid;

    void GetParams();
    void BuildCDF();
};

#endif //COMPFLUX_HH
