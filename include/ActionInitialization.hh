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
    ActionInitialization(Geometry *, G4bool);
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4bool verticalFlux;
    Geometry *geometry = nullptr;
};

#endif //ACTIONINITIALIZATION_HH