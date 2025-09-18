#include "RunAction.hh"


RunAction::RunAction(const Sizes &ss) : sizes(ss) {
    auto formatDouble = [](G4double value) {
        std::ostringstream ss;
        ss << std::setprecision(4) << value;
        return ss.str();
    };

    G4String fileName = "GammaCube_S" + formatDouble(sizes.shellThick / mm) +
                        "_V" + formatDouble(sizes.vetoThick) +
                        "_G" + formatDouble(sizes.gapSize / mm) +
                        "_T" + formatDouble(sizes.tapeThick / mm) +
                        "_L" + formatDouble(sizes.lidThick / mm) + ".root";
    analysisManager = new AnalysisManager(fileName);
}

RunAction::~RunAction() {
    delete analysisManager;
}

void RunAction::BeginOfRunAction(const G4Run *) {
    analysisManager->Open();
}

void RunAction::EndOfRunAction(const G4Run *) {
    analysisManager->Close();
}
