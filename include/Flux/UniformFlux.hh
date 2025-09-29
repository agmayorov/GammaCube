#ifndef UNIFORMFLUX_HH
#define UNIFORMFLUX_HH

#include "Flux.hh"

class UniformFlux : public Flux {
public:
    UniformFlux() = default;

private:
    G4double SampleEnergy() override;
    void PickParticle();
};


#endif //UNIFORMFLUX_HH