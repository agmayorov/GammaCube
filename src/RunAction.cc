#include "RunAction.hh"

using namespace Configuration;

RunAction::RunAction() {
    const G4String fileName = "GammaCube.root";
    analysisManager = new AnalysisManager(fileName);

    auto* mgr = G4AccumulableManager::Instance();
    mgr->Register(crystalOnly);
    mgr->Register(crystalAndVeto);
}

RunAction::RunAction(const double Agen_cm2, const double Emin_MeV,
                     const double Emax_MeV) : EminMeV(Emin_MeV),
                                              EmaxMeV(Emax_MeV),
                                              area(Agen_cm2) {
    analysisManager = new AnalysisManager(outputFile, nBins, EminMeV, EmaxMeV);
    if (nBins < 1) {
        throw std::runtime_error("RunAction: nbins must be >= 1");
    }
    if (EminMeV <= 0.0 || EmaxMeV <= 0.0 || EmaxMeV < EminMeV) {
        throw std::runtime_error("RunAction: invalid Emin/Emax");
    }
    G4double E_low = EmaxMeV > EminMeV ? EminMeV : eCrystalThreshold;
    logEmin = std::log10(E_low);
    logEmax = std::log10(EmaxMeV);
    invDlogE = static_cast<double>(nBins) / (logEmax - logEmin);

    effArea.assign(nBins, 0.0);
    effAreaOpt.assign(nBins, 0.0);

    BookAccumulables();
}

void RunAction::BookAccumulables() {
    auto* mgr = G4AccumulableManager::Instance();

    mgr->Register(crystalOnly);
    mgr->Register(crystalAndVeto);

    mgr->Register(crystalOnlyOpt);
    mgr->Register(crystalAndVetoOpt);

    genCounts.clear();
    trigCounts.clear();
    trigOptCounts.clear();
    genCounts.reserve(EminMeV < EmaxMeV ? nBins : 1);
    trigCounts.reserve(nBins);
    trigOptCounts.reserve(nBins);

    for (int i = 0; i < nBins; ++i) {
        trigCounts.emplace_back(0.0);
        mgr->Register(trigCounts.back());
        trigOptCounts.emplace_back(0.0);
        mgr->Register(trigOptCounts.back());
    }
    for (int i = 0; i < (EminMeV < EmaxMeV ? nBins : 1); ++i) {
        genCounts.emplace_back(0.0);
        mgr->Register(genCounts.back());
    }
}

RunAction::~RunAction() {
    delete analysisManager;
}

void RunAction::BeginOfRunAction(const G4Run*) {
    analysisManager->Open();
    auto* mgr = G4AccumulableManager::Instance();
    mgr->Reset();

    totals = {};
    totalsOpt = {};
    std::fill(effArea.begin(), effArea.end(), 0.0);
    std::fill(effAreaOpt.begin(), effAreaOpt.end(), 0.0);
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto* mgr = G4AccumulableManager::Instance();
    mgr->Merge();
    if (G4Threading::IsMasterThread()) {
        totals.crystalAndVeto = crystalAndVeto.GetValue();
        totals.crystalOnly = crystalOnly.GetValue();
        totalsOpt.crystalAndVeto = crystalAndVetoOpt.GetValue();
        totalsOpt.crystalOnly = crystalOnlyOpt.GetValue();
        if (EminMeV < EmaxMeV) {
            FillDerivedHists();
        }
    }

    analysisManager->Close();
}

int RunAction::FindBinLog(double E_MeV) const {
    const double E_low = EmaxMeV > EminMeV ? EminMeV : eCrystalThreshold;
    if (E_MeV < E_low || E_MeV >= EmaxMeV) return -1;
    const double le = std::log10(E_MeV);
    int idx = static_cast<int>((le - logEmin) * invDlogE);

    if (idx < 0) idx = 0;
    if (idx >= nBins) idx = nBins - 1;
    return idx;
}

double RunAction::BinCenterMeV(int i) const {
    const double ratio = EmaxMeV / EminMeV;

    const double t1 = static_cast<double>(i) / static_cast<double>(nBins);
    const double t2 = static_cast<double>(i + 1) / static_cast<double>(nBins);

    const double e1 = EminMeV * std::pow(ratio, t1);
    const double e2 = EminMeV * std::pow(ratio, t2);

    return std::sqrt(e1 * e2);
}

double RunAction::BinWidthMeV(int i) const {
    const double ratio = EmaxMeV / EminMeV;

    const double t1 = static_cast<double>(i) / static_cast<double>(nBins);
    const double t2 = static_cast<double>(i + 1) / static_cast<double>(nBins);

    const double e1 = EminMeV * std::pow(ratio, t1);
    const double e2 = EminMeV * std::pow(ratio, t2);

    return e2 - e1;
}

void RunAction::AddGenerated(double E_MeV) {
    const int i = EminMeV < EmaxMeV ? FindBinLog(E_MeV) : 0;
    if (i < 0) return;

    genCounts[i] += 1.0;

    if (analysisManager and EminMeV < EmaxMeV) {
        analysisManager->FillGenEnergyHist(E_MeV, 1.0);
    }
}

void RunAction::AddTriggeredCrystalOnly(double E_MeV) {
    const int i = FindBinLog(E_MeV);
    if (i < 0) return;

    trigCounts[i] += 1.0;

    if (analysisManager and EminMeV < EmaxMeV) {
        analysisManager->FillTrigEnergyHist(E_MeV, 1.0);
    }
}

void RunAction::AddTriggeredCrystalOnlyOpt(double E_MeV) {
    const int i = FindBinLog(E_MeV);
    if (i < 0) return;

    trigOptCounts[i] += 1.0;

    if (analysisManager and EminMeV < EmaxMeV) {
        analysisManager->FillTrigOptEnergyHist(E_MeV, 1.0);
    }
}

void RunAction::FillDerivedHists() {
    for (int i = 0; i < nBins; ++i) {
        const double nGen = genCounts[i].GetValue();
        const double nTrig = trigCounts[i].GetValue();
        const double nTrigOpt = trigOptCounts[i].GetValue();

        const double centerE = BinCenterMeV(i);
        double aEff = 0.0;
        double aEffOpt = 0.0;
        double sens = 0.0;
        double sensOpt = 0.0;
        if (nGen > 0.0) {
            aEff = area * (nTrig / nGen);
            aEffOpt = area * (nTrigOpt / nGen);
            sens = aEff;
            sensOpt = aEffOpt;
        }
        effArea[i] = aEff;
        effAreaOpt[i] = aEffOpt;
        if (fluxDirection.find("isotropic") != std::string::npos) {
            analysisManager->FillSensitivityHist(centerE, sens);
            analysisManager->FillSensitivityOptHist(centerE, sensOpt);
        }
        analysisManager->FillEffAreaHist(centerE, aEff);
        analysisManager->FillEffAreaOptHist(centerE, aEffOpt);
    }
}
