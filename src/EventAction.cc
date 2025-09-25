#include "EventAction.hh"


EventAction::EventAction(AnalysisManager *an, const Sizes &sizes) : analysisManager(an), detSizes(sizes) {
    detMap = {{"DetectorSD/EdepHits", 0, "NaI"}};
    if (detSizes.shellThick > 0.) {
        detMap.emplace_back("ShellSD/EdepHits", 1, "Shell");
    }
    if (detSizes.tunaCanThick > 0.) {
        detMap.emplace_back("TunaCanSD/EdepHits", 2, "TunaCan");
    }
    if (detSizes.vetoThick > 0.) {
        detMap.emplace_back("VetoSD/EdepHits", 3, "Veto");
    }
    if (detSizes.tyvekThick > 0.) {
        detMap.emplace_back("TyvekOutSD/EdepHits", 4, "TyvekOut");
        detMap.emplace_back("TyvekInSD/EdepHits", 5, "TyvekIn");
    }
    HCIDs.assign(detMap.size(), -1);
}

void EventAction::BeginOfEventAction(const G4Event *) {
    nPrimaries = 0;
    nInteractions = 0;
    nEdepHits = 0;
}

void EventAction::EndOfEventAction(const G4Event *evt) {
    // const int runID = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());
    primBuf.clear();

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    nEdepHits = WriteEdepFromSD_(evt, eventID);

    analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);
}

// ------------------- приватные помощники -------------------

void EventAction::WritePrimaries_(int eventID) {
    for (const auto &p: primBuf) {
        analysisManager->FillPrimaryRow(eventID,
                                        p.index, p.pdg, p.name,
                                        p.E_MeV, p.dir, p.pos_mm, p.t0_ns);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    int n = 0;
    for (const auto &r: interBuf) {
        // analysisManager->FillInteractionRow(eventID,
        //                                     r.trackID, r.parentID,
        //                                     r.process,
        //                                     r.volumeName, r.volumeID,
        //                                     r.pos_mm, r.t_ns,
        //                                     r.secIndex, r.secPDG, r.secName,
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

        const int det_id = std::get<1>(detMap[i]);
        const auto &det_name = std::get<2>(detMap[i]);

        const auto N = hc->GetSize();
        for (unsigned j = 0; j < N; ++j) {
            auto *h = (*hc)[j];
            const double edep_MeV = h->edep / MeV;
            const double tmin_ns = (h->tmin < DBL_MAX / 2 ? h->tmin / ns : -1.0);
            analysisManager->FillEdepRow(eventID, det_id, det_name,
                                         h->volumeID, edep_MeV, tmin_ns);
        }
        nHitsTotal += static_cast<int>(N);
    }

    return nHitsTotal;
}
