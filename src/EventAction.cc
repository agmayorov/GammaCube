#include "EventAction.hh"

using namespace Sizes;
using namespace Configuration;

EventAction::EventAction(AnalysisManager* an, RunAction* r) : analysisManager(an), run(r) {
    detMap = {
        {"DetectorSD/EdepHits", 0, "Crystal"},
        {"VetoSD/EdepHits", 1, "Veto"},
        {"BottomVetoSD/EdepHits", 2, "BottomVeto"},
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
    const int eventID = evt->GetEventID();

    WritePrimaries_(eventID);
    nPrimaries = static_cast<int>(primBuf.size());

    double primaryE_MeV = -1.0;
    if (!primBuf.empty()) {
        primaryE_MeV = primBuf.front().E_MeV;
        if (run) {
            run->AddGenerated(primaryE_MeV);
        }
    }
    primBuf.clear();

    nInteractions = WriteInteractions_(eventID);
    interBuf.clear();

    if (savePhotons) {
        nPhotons = WritePhotons_(eventID);
        photonBuf.clear();
        WritePhotonsCount_(eventID);
        photonCountBuf = {0, 0, 0};
    }


    nEdepHits = WriteEdepFromSD_(evt, eventID);

    if (saveSecondaries) {
        analysisManager->FillEventRow(eventID, nPrimaries, nInteractions, nEdepHits);
    }

    if (run and hasCrystal && !hasVeto) run->AddCrystalOnly(1);
    if (run and hasCrystal && hasVeto) run->AddCrystalAndVeto(1);

    if (primaryE_MeV > 0.0) {
        if (hasCrystal && !hasVeto) {
            if (run) {
                run->AddTriggeredCrystalOnly(primaryE_MeV);
            }
        }
    }

    if (useOptics) {
        WriteSiPMFromSD_(eventID);
        if (run and hasCrystalOpt && !hasVetoOpt) run->AddCrystalOnlyOpt(1);
        if (run and hasCrystalOpt && hasVetoOpt) run->AddCrystalAndVetoOpt(1);

        if (primaryE_MeV > 0.0) {
            if (run and hasCrystalOpt && !hasVetoOpt) run->AddTriggeredCrystalOnlyOpt(primaryE_MeV);
        }
    }
}

void EventAction::WritePrimaries_(int eventID) {
    for (const auto& p : primBuf) {
        analysisManager->FillPrimaryRow(eventID, p.name, p.E_MeV, p.dir, p.pos_mm);
    }
}

int EventAction::WriteInteractions_(int eventID) {
    if (saveSecondaries) {
        for (const auto& r : interBuf) {
            analysisManager->FillInteractionRow(eventID,
                                                r.trackID, r.parentID,
                                                r.process, r.volumeName, r.pos_mm,
                                                r.secIndex, r.secName,
                                                r.secE_MeV, r.secDir);
        }
    }
    return static_cast<int>(interBuf.size());
}

int EventAction::WritePhotonsCount_(int eventID) {
    if (savePhotons) {
        analysisManager->FillPhotonCountRow(eventID, photonCountBuf[0], photonCountBuf[1], photonCountBuf[2]);
    }
    return static_cast<int>(photonCountBuf.size());
}

int EventAction::WritePhotons_(int eventID) {
    for (const auto& photon : photonBuf) {
        analysisManager->FillPhotonRow(eventID, photon.photonID, photon.detName, photon.detCh, photon.energy,
                                       photon.pos_mm.x(), photon.pos_mm.y(), photon.pos_mm.z());
    }
    return static_cast<int>(photonBuf.size());
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
            if (edep_MeV > 0.0) {
                if (det_name == "Crystal") MarkCrystal();
                else if (det_name == "Veto" or det_name == "BottomVeto") MarkVeto();
                analysisManager->FillEdepRow(eventID, det_name, edep_MeV);
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

    int npeC = sipmSD->GetNpeCrystal();
    int npeV = sipmSD->GetNpeVeto();
    int npeB = sipmSD->GetNpeBottomVeto();

    npeC = npeC > oCrystalThreshold ? npeC : 0;
    npeV = npeV > oVetoThreshold ? npeV : 0;
    npeB = npeB > oBottomVetoThreshold ? npeB : 0;

    if (npeC > 0) MarkCrystalOpt();
    if (npeV > 0 or npeB > 0) MarkVetoOpt();

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
