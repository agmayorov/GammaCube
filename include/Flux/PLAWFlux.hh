#ifndef PLAWFLUX_HH
#define PLAWFLUX_HH

#include "Flux.hh"

class PLAWFlux : public Flux {
public:
    PLAWFlux();

protected:
    G4double SampleEnergy() override;

private:
    G4double alpha{};

    void GetParams();
};


#endif //PLAWFLUX_HH
