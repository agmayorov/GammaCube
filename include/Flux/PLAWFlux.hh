#ifndef PLAWFLUX_HH
#define PLAWFLUX_HH

#include "Flux.hh"

class PLAWFlux : public Flux {
public:
    PLAWFlux(G4double, G4double, G4double);

protected:
    G4double SampleEnergy() override;

private:
    G4double alpha;
};


#endif //PLAWFLUX_HH
