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
    ActionInitialization(G4String, G4String, G4double, G4double, G4bool);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4ThreeVector center;
    G4ThreeVector modelSize;
    G4String fluxDirection;
    G4String fluxType;
    G4bool lightCollection;

    G4double eCrystalThreshold;
    G4double eVetoThreshold;
};

#endif //ACTIONINITIALIZATION_HH