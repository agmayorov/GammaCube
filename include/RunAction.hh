#ifndef RUNACTION_HH
#define RUNACTION_HH

#include <G4UserRunAction.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <G4Accumulable.hh>
#include <G4AccumulableManager.hh>
#include <G4Run.hh>
#include <G4ios.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>
#include <G4AnalysisManager.hh>
#include <G4Threading.hh>
#include <iomanip>
#include <sstream>

#include "Sizes.hh"
#include "AnalysisManager.hh"

struct ParticleCounts {
    G4int crystalOnly = 0;
    G4int crystalAndVeto = 0;
};

class RunAction : public G4UserRunAction {
public:
    AnalysisManager *analysisManager;

    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *) override;
    void EndOfRunAction(const G4Run *) override;

    void AddCrystalOnly(const G4int v) { crystalOnly += v; }
    void AddCrystalAndVeto(const G4int v) { crystalAndVeto += v; }

    [[nodiscard]] const ParticleCounts& GetCounts() const { return totals; }

private:
    G4Accumulable<G4int> crystalOnly;   // Crystal && !Veto
    G4Accumulable<G4int> crystalAndVeto;   // Crystal && Veto
    ParticleCounts totals;
};

#endif //RUNACTION_HH
