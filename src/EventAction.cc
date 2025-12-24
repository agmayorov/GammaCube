#include "EventAction.hh"

using namespace Sizes;

EventAction::EventAction(AnalysisManager* an, RunAction* r, const G4double eCrystalThr, const G4double eVetoThr,
                         const G4bool saveOpt) : analysisManager(an), run(r), eCrystalThreshold(eCrystalThr),
                                                eVetoThreshold(eVetoThr), saveOptics(saveOpt) {
    detMap = {
        {"DetectorSD/EdepHits", 0, "Crystal"},
        {"VetoSD/EdepHits", 1, "Veto"},
        {"BottomVetoSD/EdepHits", 2, "BottomVeto"},
        // {"TyvekOutSD/EdepHits",     3, "TyvekOut"},
        // {"TyvekMidSD/EdepHits",     4, "TyvekMid"},
        // {"TyvekInSD/EdepHits",      5, "TyvekIn"},
        // {"TyvekBottomSD/EdepHits",  6, "TyvekBottom"},
    };
    HCIDs.assign(detMap.size(), -1);
}

void EventAction::BeginOfEventAction(const G4Event*) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;
    hasCrystal = false;
    hasVeto = false;
}

void EventAction::EndOfEventAction(const G4Event* evt) {
    // const int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());
    primBuf.clear();

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    nEdepHits = WriteEdepFromSD_(evt, eventID);

    // analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);

    if (hasCrystal && !hasVeto) run->AddCrystalOnly(1);
    if (hasCrystal && hasVeto) run->AddCrystalAndVeto(1);

    if (saveOptics) WriteSiPMFromSD_(eventID);
}

// ------------------- приватные помощники -------------------

void EventAction::WritePrimaries_(int eventID) {
    for (const auto& p : primBuf) {
        analysisManager->FillPrimaryRow(eventID, p.name, p.E_MeV, p.dir, p.pos_mm);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    int n = 0;
    for (const auto& r : interBuf) {
        // analysisManager->FillInteractionRow(eventID,
        //                                     r.trackID, r.parentID,
        //                                     r.process, r.volumeName, r.pos_mm,
        //                                     r.secIndex, r.secName,
        //                                     r.secE_MeV, r.secDir);
        n++;
    }
    return n;
}

int EventAction::WriteEdepFromSD_(const G4Event* evt, int eventID) {
    auto* hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    auto* sdm = G4SDManager::GetSDMpointer();

    for (size_t i = 0; i < detMap.size(); ++i) {
        if (HCIDs[i] < 0) {
            const auto& hcName = std::get<0>(detMap[i]);
            HCIDs[i] = sdm->GetCollectionID(hcName);
        }
    }

    int nHitsTotal = 0;

    for (size_t i = 0; i < detMap.size(); ++i) {
        const int hcID = HCIDs[i];
        if (hcID < 0) continue;

        auto* hc = dynamic_cast<SDHitCollection*>(hce->GetHC(hcID));
        if (!hc) continue;

        const auto& det_name = std::get<2>(detMap[i]);

        const auto N = hc->GetSize();
        for (unsigned j = 0; j < N; ++j) {
            auto* h = (*hc)[j];
            double edep_MeV = h->edep / MeV;

            if ((det_name == "Veto" or det_name == "BottomVeto") and edep_MeV <= eVetoThreshold) {
                edep_MeV = 0;
            }
            if (det_name == "Crystal" and edep_MeV <= eCrystalThreshold) {
                edep_MeV = 0;
            }
            analysisManager->FillEdepRow(eventID, det_name, edep_MeV);
            if (edep_MeV > 0.0) {
                if (det_name == "Crystal") MarkCrystal();
                else if (det_name == "Veto" or det_name == "BottomVeto") MarkVeto();
            }
        }
        nHitsTotal += static_cast<int>(N);
    }

    return nHitsTotal;
}

void EventAction::WriteSiPMFromSD_(int eventID) {
    auto* sdm = G4SDManager::GetSDMpointer();
    if (!sdm) return;

    auto* sdBase = sdm->FindSensitiveDetector("SiPMOpticalSD", false);
    auto* sipmSD = dynamic_cast<SiPMOpticalSD*>(sdBase);
    if (!sipmSD) return;

    const int npeC = sipmSD->GetNpeCrystal();
    const int npeV = sipmSD->GetNpeVeto();
    const int npeB = sipmSD->GetNpeBottomVeto();

    analysisManager->FillSiPMEventRow(eventID, npeC, npeV, npeB);


    for (const auto& kv : sipmSD->GetPerChannelCrystal()) {
        const int ch = kv.first;
        const int npe = kv.second;
        analysisManager->FillSiPMChannelRow(eventID, "Crystal", ch, npe);
    }

    for (const auto& kv : sipmSD->GetPerChannelVeto()) {
        const int ch = kv.first;
        const int npe = kv.second;
        analysisManager->FillSiPMChannelRow(eventID, "Veto", ch, npe);
    }

    for (const auto& kv : sipmSD->GetPerChannelBottom()) {
        const int ch = kv.first;
        const int npe = kv.second;
        analysisManager->FillSiPMChannelRow(eventID, "BottomVeto", ch, npe);
    }
}
