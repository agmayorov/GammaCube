#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(G4String fDir, G4String fType, const G4double eCrystalThr,
                                           const G4double eVetoThr, const G4bool useOpt,
                                           const G4bool saveSec) : fluxDirection(std::move(fDir)),
                                                                   fluxType(std::move(fType)),
                                                                   eCrystalThreshold(eCrystalThr),
                                                                   eVetoThreshold(eVetoThr),
                                                                   useOptics(useOpt),
                                                                   saveSecondaries(saveSec) {}

void ActionInitialization::BuildForMaster() const {
    RunAction* runAct = new RunAction();
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction* runAct = new RunAction();
    SetUserAction(runAct);

    EventAction* eventAct = new EventAction(runAct->analysisManager, runAct, eCrystalThreshold, eVetoThreshold,
                                            useOptics, saveSecondaries);
    SetUserAction(eventAct);

    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType);
    SetUserAction(primaryGenerator);

    SteppingAction* stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
