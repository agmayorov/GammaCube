#ifndef PRMIARYGENERATIONACTION_HH
#define PRMIARYGENERATIONACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4VVisManager.hh>
#include <G4Event.hh>
#include <G4Circle.hh>
#include <G4EventManager.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <cmath>
#include <fstream>
#include <numeric>

#include "EventAction.hh"
#include "Geometry.hh"
#include "Flux/Flux.hh"
#include "Flux/UniformFlux.hh"
#include "Flux/PLAWFlux.hh"
#include "Flux/SEPFlux.hh"
#include "Flux/GalacticFlux.hh"


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const Geometry *, G4bool, const G4String& );
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event *evt) override;

private:
    G4ParticleGun *particleGun = nullptr;

    G4double radius;
    G4ThreeVector center;
    G4ThreeVector detectorHalfSize;

    G4bool verticalFlux;
    ParticleInfo pInfo{};

    Flux *flux;

    void GenerateOnSphere(G4ThreeVector &pos, G4ThreeVector &dir) const;
};

#endif //PRMIARYGENERATIONACTION_HH
