#ifndef FLUX_HH
#define FLUX_HH

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4ParticleTable.hh"
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

struct ParticleInfo {
    std::string name;
    G4int pdg;
    G4ParticleDefinition *def;
    G4double energy;
};


class Flux {
public:
    virtual ~Flux() = default;

    virtual ParticleInfo GenerateParticle();

protected:
    G4String name;
    G4double Emin{};
    G4double Emax{};

    virtual G4double SampleEnergy() = 0;
};


#endif //FLUX_HH
