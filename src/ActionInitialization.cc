#include "ActionInitialization.hh"

using namespace Configuration;

ActionInitialization::ActionInitialization(const G4double a, const G4double Emin,
                                           const G4double Emax) : EminMeV(Emin),
                                                                  EmaxMeV(Emax),
                                                                  area(a) {
    if (eCrystalThreshold > EmaxMeV) {
        G4Exception("ActionInitialization", "EnergyRange", FatalException,
                    "The energy threshold for a crystal must be less than the maximum value in a given energy range");
    }
}

void ActionInitialization::BuildForMaster() const {
    RunAction* runAct = new RunAction(area, EminMeV, EmaxMeV);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction* runAct = new RunAction(area, EminMeV, EmaxMeV);
    SetUserAction(runAct);

    EventAction* eventAct = new EventAction(runAct->analysisManager, runAct);
    SetUserAction(eventAct);

    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType, eCrystalThreshold);
    SetUserAction(primaryGenerator);

    if (saveSecondaries || savePhotons) {
        SteppingAction* stepAct = new SteppingAction();
        SetUserAction(stepAct);
    }
}
