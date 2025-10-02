#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(Sizes &s, G4ThreeVector mSize, G4ThreeVector c, G4bool vFlux, G4String fType)
    : sizes(s), modelSize(std::move(mSize)), center(std::move(c)), verticalFlux(vFlux), fluxType(std::move(fType)) {
}

void ActionInitialization::BuildForMaster() const {
    RunAction *runAct = new RunAction(sizes);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction *runAct = new RunAction(sizes);
    SetUserAction(runAct);

    EventAction *eventAct = new EventAction(runAct->analysisManager, sizes);
    SetUserAction(eventAct);

    PrimaryGeneratorAction *primaryGenerator = new PrimaryGeneratorAction(center, modelSize, verticalFlux, fluxType);
    SetUserAction(primaryGenerator);

    SteppingAction *stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
