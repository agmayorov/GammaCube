#include "ActionInitialization.hh"


ActionInitialization::ActionInitialization(G4String fDir, G4String fType, const G4bool useOpt, const G4bool saveSec,
                                           const G4double eCrystalThr, const G4double eVetoThr, const G4double a,
                                           const G4int nbins, const G4double Emin, const G4double Emax,
                                           G4String file) : fluxDirection(std::move(fDir)),
                                                            fluxType(std::move(fType)),
                                                            useOptics(useOpt),
                                                            saveSecondaries(saveSec),
                                                            eCrystalThreshold(eCrystalThr),
                                                            eVetoThreshold(eVetoThr),
                                                            nBins(nbins),
                                                            EminMeV(Emin),
                                                            EmaxMeV(Emax),
                                                            area(a),
                                                            fileName(std::move(file)) {
    if (eCrystalThreshold > EmaxMeV) {
        G4Exception("ActionInitialization", "EnergyRange", FatalException,
                    "The energy threshold for a crystal must be less than the maximum value in a given energy range");
    }
}

void ActionInitialization::BuildForMaster() const {
    RunAction* runAct = new RunAction(nBins, EminMeV, EmaxMeV, eCrystalThreshold, area, fileName);
    SetUserAction(runAct);
}

void ActionInitialization::Build() const {
    RunAction* runAct = new RunAction(nBins, EminMeV, EmaxMeV, eCrystalThreshold, area, fileName);
    SetUserAction(runAct);

    EventAction* eventAct = new EventAction(runAct->analysisManager, runAct, eCrystalThreshold, eVetoThreshold,
                                            useOptics, saveSecondaries);
    SetUserAction(eventAct);

    PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction(fluxDirection, fluxType, eCrystalThreshold);
    SetUserAction(primaryGenerator);

    SteppingAction* stepAct = new SteppingAction();
    SetUserAction(stepAct);
}
