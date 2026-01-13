#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include <G4VUserActionInitialization.hh>

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "Geometry.hh"

class ActionInitialization : public G4VUserActionInitialization {
public:
    ActionInitialization(G4String, G4String, G4bool, G4bool, G4double, G4double, G4double, G4int, G4double, G4double,
                         G4String);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4String fluxDirection;
    G4String fluxType;
    G4bool useOptics;
    G4bool saveSecondaries;

    G4double eCrystalThreshold;
    G4double eVetoThreshold;

    G4int nBins;
    G4double EminMeV;
    G4double EmaxMeV;
    G4double area;

    G4String fileName;
};

#endif //ACTIONINITIALIZATION_HH
