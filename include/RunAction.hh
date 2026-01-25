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
    RunAction(int nbins, double Emin_MeV, double Emax_MeV, double cThreshold, double Agen_cm2,
              const std::string& rootFileName = "GammaCube");
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *) override;
    void EndOfRunAction(const G4Run *) override;

    void AddCrystalOnly(const G4int v) { crystalOnly += v; }
    void AddCrystalAndVeto(const G4int v) { crystalAndVeto += v; }

    void AddGenerated(double E_MeV);
    void AddTriggeredCrystalOnly(double E_MeV);

    [[nodiscard]] const ParticleCounts& GetCounts() const { return totals; }

    [[nodiscard]] const std::vector<double>& GetEffArea() const { return effArea; }

private:
    G4Accumulable<G4int> crystalOnly{0};   // Crystal && !Veto
    G4Accumulable<G4int> crystalAndVeto{0};   // Crystal && Veto
    ParticleCounts totals{};

    int nBins{0};
    double EminMeV{0.0};
    double EmaxMeV{0.0};
    double eCrystalThreshold{0.0};
    double area{0.0};

    double logEmin{0.0};
    double logEmax{0.0};
    double invDlogE{0.0};

    std::vector<G4Accumulable<G4double>> genCounts;
    std::vector<G4Accumulable<G4double>> trigCounts;
    std::vector<G4double> effArea;

    [[nodiscard]] int FindBinLog(double E_MeV) const;
    [[nodiscard]] double BinCenterMeV(int i) const;
    [[nodiscard]] double BinWidthMeV(int i) const;

    void BookAccumulables();
    void FillDerivedHists();
};

#endif //RUNACTION_HH
