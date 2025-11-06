#include "EventAction.hh"

using namespace Sizes;

EventAction::EventAction(AnalysisManager *an, RunAction *r, const G4double eCrystalThr, const G4double eVetoThr,
                         const G4bool lightCollection) : analysisManager(an), run(r),
                                                         eCrystalThreshold(eCrystalThr), eVetoThreshold(eVetoThr) {
    detMap = {{"DetectorSD/EdepHits", 0, "Crystal"}};
    if (tunaCanMinSize > 0.) {
        detMap.emplace_back("TunaCanSD/EdepHits", 1, "TunaCan");
    }
    if (vetoThickBottom + vetoThickTop + vetoThickWall > 0.) {
        detMap.emplace_back("VetoSD/EdepHits", 2, "Veto");
    }
    if (tyvekOutMinSize > 0.) {
        detMap.emplace_back("TyvekOutSD/EdepHits", 3, "TyvekOut");
    }
    if (tyvekInMinSize > 0.) {
        detMap.emplace_back("TyvekInSD/EdepHits", 4, "TyvekIn");
    }
    if (tyvekMidMinSize > 0.) {
        detMap.emplace_back("TyvekMidSD/EdepHits", 5, "TyvekMid");
    }
    if (rubberMinSize > 0.) {
        detMap.emplace_back("RubberSD/EdepHits", 6, "Rubber");
    }
    if (AlMinSize > 0.) {
        detMap.emplace_back("AluminiumSD/EdepHits", 7, "Aluminium");
    }
    if (crystalLEDMinSize > 0.) {
        detMap.emplace_back("CrystalLEDSD/EdepHits", 8, "CrystalLED");
    }
    if (vetoLEDMinSize > 0.) {
        detMap.emplace_back("VetoLEDSD/EdepHits", 9, "VetoLED");
    }
    if (vetoBottomLEDMinSize > 0.) {
        detMap.emplace_back("VetoBottomLEDSD/EdepHits", 10, "VetoBottomLED");
    }
    if (lightCollection) {
        optMap.emplace_back("CrystalSensSurfSD/OptHits", 0, "CrystalSensSurf");
        optMap.emplace_back("VetoSensSurfSD/OptHits", 1, "VetoSensSurf");
        optMap.emplace_back("VetoBottomSensSurfSD/OptHits", 2, "VetoBottomSensSurf");
    }
    HCIDs.assign(detMap.size(), -1);
    optHCIDs.assign(optMap.size(), -1);
}

void EventAction::BeginOfEventAction(const G4Event *) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;
    hasCrystal = false;
    hasVeto = false;
}

void EventAction::EndOfEventAction(const G4Event *evt) {
    // const int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());
    primBuf.clear();

    WriteOptFromSD_(evt);

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    nEdepHits = WriteEdepFromSD_(evt, eventID);

    // analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);

    if (hasCrystal && !hasVeto) run->AddCrystalOnly(1);
    if (hasCrystal && hasVeto) run->AddCrystalAndVeto(1);
}

// ------------------- приватные помощники -------------------

void EventAction::WritePrimaries_(int eventID) {
    for (const auto &p: primBuf) {
        analysisManager->FillPrimaryRow(eventID, p.name, p.E_MeV, p.dir, p.pos_mm);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    int n = 0;
    for (const auto &r: interBuf) {
        // analysisManager->FillInteractionRow(eventID,
        //                                     r.trackID, r.parentID,
        //                                     r.process, r.volumeName, r.pos_mm,
        //                                     r.secIndex, r.secName,
        //                                     r.secE_MeV, r.secDir);
        n++;
    }
    return n;
}

int EventAction::WriteEdepFromSD_(const G4Event *evt, int eventID) {
    auto *hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    auto *sdm = G4SDManager::GetSDMpointer();

    for (size_t i = 0; i < detMap.size(); ++i) {
        if (HCIDs[i] < 0) {
            const auto &hcName = std::get<0>(detMap[i]);
            HCIDs[i] = sdm->GetCollectionID(hcName);
        }
    }

    int nHitsTotal = 0;

    for (size_t i = 0; i < detMap.size(); ++i) {
        const int hcID = HCIDs[i];
        if (hcID < 0) continue;

        auto *hc = dynamic_cast<SDHitCollection *>(hce->GetHC(hcID));
        if (!hc) continue;

        const auto &det_name = std::get<2>(detMap[i]);

        const auto N = hc->GetSize();
        for (unsigned j = 0; j < N; ++j) {
            auto *h = (*hc)[j];
            double edep_MeV = h->edep / MeV;

            if (det_name == "Veto" and edep_MeV <= eVetoThreshold) { edep_MeV = 0; }
            if (det_name == "Crystal" and edep_MeV <= eCrystalThreshold) { edep_MeV = 0; }
            analysisManager->FillEdepRow(eventID, det_name, edep_MeV);
            if (edep_MeV > 0.0) {
                if (det_name == "Crystal") MarkCrystal();
                else if (det_name == "Veto") MarkVeto();
            }
        }
        nHitsTotal += static_cast<int>(N);
    }

    return nHitsTotal;
}


int EventAction::WriteOptFromSD_(const G4Event *evt) {
    auto *hce = evt->GetHCofThisEvent();
    if (!hce) return 0;

    auto *sdm = G4SDManager::GetSDMpointer();
    for (size_t i = 0; i < optMap.size(); ++i) {
        if (optHCIDs[i] < 0) {
            const auto &hcName = std::get<0>(optMap[i]);
            optHCIDs[i] = sdm->GetCollectionID(hcName);
        }
    }

    int n = 0;
    for (size_t i = 0; i < optMap.size(); ++i) {
        const int hcID = optHCIDs[i];
        if (hcID < 0) continue;

        auto *hc = dynamic_cast<SDHitCollection *>(hce->GetHC(hcID));
        if (!hc) continue;

        const auto &tag = std::get<2>(optMap[i]);
        const auto N = hc->GetSize();

        for (unsigned j = 0; j < N; ++j) {
            auto *h = (*hc)[j];
            if (!h->isOptical) continue;

            if (tag == "CrystalSensSurf") {
                analysisManager->FillCrystalH2(h->x_loc_mm, h->y_loc_mm);
            } else if (tag == "VetoSensSurf") {
                analysisManager->FillVetoH2(h->x_loc_mm, h->y_loc_mm);
            } else if (tag == "VetoBottomSensSurf") {
                analysisManager->FillVetoBottomH2(h->z_loc_mm, h->phi);
            }
        }
        n += static_cast<int>(N);
    }
    return n;
}
