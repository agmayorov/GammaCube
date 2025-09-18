#ifndef RUNACTION_HH
#define RUNACTION_HH

#include <G4UserRunAction.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>
#include <G4AnalysisManager.hh>
#include <G4Threading.hh>
#include <iomanip>
#include <sstream>

#include "Sizes.hh"
#include "AnalysisManager.hh"


class RunAction : public G4UserRunAction {
public:
    AnalysisManager *analysisManager;

    RunAction(const Sizes &ss);
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *) override;
    void EndOfRunAction(const G4Run *) override;

private:
    Sizes sizes;
};

#endif //RUNACTION_HH
