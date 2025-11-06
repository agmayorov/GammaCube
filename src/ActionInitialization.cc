#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(G4String fDir, G4String fType, const G4double eCrystalThr,
                                           const G4double eVetoThr,
                                           const G4bool lightCollect) : fluxDirection(std::move(fDir)),
                                                                        fluxType(std::move(fType)),
                                                                        eCrystalThreshold(eCrystalThr),
                                                                        eVetoThreshold(eVetoThr),
                                                                        lightCollection(lightCollect) {
}

void ActionInitialization::BuildForMaster() const {
    RunAction *runAct = new RunAction();
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction *runAct = new RunAction();
    SetUserAction(runAct);

    EventAction *eventAct = new EventAction(runAct->analysisManager, runAct, eCrystalThreshold, eVetoThreshold,
                                            lightCollection);
    SetUserAction(eventAct);

    PrimaryGeneratorAction *primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType);
    SetUserAction(primaryGenerator);

    SteppingAction *stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
