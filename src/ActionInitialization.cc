#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(Geometry *geom, G4bool vFlux) : geometry(geom), verticalFlux(vFlux) {}

void ActionInitialization::BuildForMaster() const {
    RunAction *runAct = new RunAction(geometry->sizes);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction *runAct = new RunAction(geometry->sizes);
    SetUserAction(runAct);

    EventAction *eventAct = new EventAction(runAct->analysisManager, geometry->sizes);
    SetUserAction(eventAct);

    PrimaryGeneratorAction *primaryGenerator = new PrimaryGeneratorAction(geometry, verticalFlux);
    SetUserAction(primaryGenerator);

    SteppingAction *stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
