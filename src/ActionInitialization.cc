#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(G4bool vFlux, G4String fType)
    : verticalFlux(vFlux), fluxType(std::move(fType)) {
}

void ActionInitialization::BuildForMaster() const {
    RunAction *runAct = new RunAction();
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction *runAct = new RunAction();
    SetUserAction(runAct);

    EventAction *eventAct = new EventAction(runAct->analysisManager, runAct);
    SetUserAction(eventAct);

    PrimaryGeneratorAction *primaryGenerator = new PrimaryGeneratorAction(verticalFlux, fluxType);
    SetUserAction(primaryGenerator);

    SteppingAction *stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
