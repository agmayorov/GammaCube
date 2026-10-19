#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4EventManager.hh>
#include <G4TrackingManager.hh>
#include <G4RunManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4TouchableHistory.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserSteppingAction.hh>
#include <G4OpticalPhoton.hh>

#include "EventAction.hh"


class EventAction;

class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction() = default;
    void UserSteppingAction(const G4Step* step) override;

private:
    static void ApplyPhotonTimeLimit(const G4Step* step, EventAction* ea);
};

#endif //STEPPINGACTION_HH
