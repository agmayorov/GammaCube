#include "RunAction.hh"


RunAction::RunAction() : crystalOnly(0), crystalAndVeto(0) {
    auto formatDouble = [](G4double value) {
        std::ostringstream ss;
        ss << std::setprecision(4) << value;
        return ss.str();
    };

    G4String fileName = "GammaCube.root";
    analysisManager = new AnalysisManager(fileName);

    auto *mgr = G4AccumulableManager::Instance();
    mgr->Register(crystalOnly);
    mgr->Register(crystalAndVeto);
}

RunAction::~RunAction() {
    delete analysisManager;
}

void RunAction::BeginOfRunAction(const G4Run *) {
    analysisManager->Open();
    auto *mgr = G4AccumulableManager::Instance();
    mgr->Reset();
}

void RunAction::EndOfRunAction(const G4Run *) {
    analysisManager->Close();
    auto *mgr = G4AccumulableManager::Instance();
    mgr->Merge();
    if (G4Threading::IsMasterThread()) {
        totals.crystalAndVeto = crystalAndVeto.GetValue();
        totals.crystalOnly = crystalOnly.GetValue();
    }
}
