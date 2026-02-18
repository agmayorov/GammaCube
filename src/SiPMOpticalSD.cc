#include "SiPMOpticalSD.hh"

SiPMOpticalSD::SiPMOpticalSD(const G4String& name)
    : G4VSensitiveDetector(name) {}

void SiPMOpticalSD::Initialize(G4HCofThisEvent*) {
    npeCrystal = npeVeto = npeBottom = 0;
    perChCrystal.clear();
    perChVeto.clear();
    perChBottom.clear();
}

G4OpBoundaryProcess* SiPMOpticalSD::GetBoundaryProcess() {
    if (boundary) return boundary;
    auto* pm = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
    if (!pm) return nullptr;

    auto* plist = pm->GetProcessList();
    if (!plist) return nullptr;

    for (int i = 0; i < plist->size(); ++i) {
        auto* p = (*plist)[i];
        if (p && p->GetProcessName() == "OpBoundary") {
            boundary = static_cast<G4OpBoundaryProcess*>(p);
            break;
        }
    }
    return boundary;
}


SiPMGroup SiPMOpticalSD::ClassifyByPVName(const G4VPhysicalVolume* pv) {
    if (!pv) return SiPMGroup::Unknown;
    const std::string name = pv->GetName();

    if (name.find("CrystalSiPM") != std::string::npos) return SiPMGroup::Crystal;
    if (name.find("BottomVetoSiPM") != std::string::npos) return SiPMGroup::Bottom;
    if (name.find("VetoSiPM") != std::string::npos) return SiPMGroup::Veto;

    return SiPMGroup::Unknown;
}

G4bool SiPMOpticalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    if (!step) return false;

    auto* track = step->GetTrack();
    if (!track) return false;

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

    // We only care about geometry boundary crossings
    auto* post = step->GetPostStepPoint();
    auto* pre = step->GetPreStepPoint();
    if (!pre || !post) return false;

    if (post->GetStepStatus() != fGeomBoundary) return false;

    auto* b = GetBoundaryProcess();
    if (!b) return false;

    if (b->GetStatus() != Detection) return false;

    // Determine on which PV we are (prefer PRE, but keep fallback)
    auto* prePV = pre->GetPhysicalVolume();
    auto* postPV = post->GetPhysicalVolume();

    // Prefer pre-step touchable for copy number (normally it's the "window" volume)
    int ch = -1;
    {
        const auto& t = pre->GetTouchableHandle();
        if (t) ch = t->GetCopyNumber(0);
    }
    if (ch < 0) {
        const auto& t = post->GetTouchableHandle();
        if (t) ch = t->GetCopyNumber(0);
    }

    SiPMGroup grp = SiPMGroup::Unknown;

    if (SiPMWindowLV) {
        auto* preLV = prePV ? prePV->GetLogicalVolume() : nullptr;
        auto* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;

        if (preLV == SiPMWindowLV) {
            grp = ClassifyByPVName(prePV);
        } else if (postLV == SiPMWindowLV) {
            grp = ClassifyByPVName(postPV);
        } else {
            // Not on a SiPM window at all (or LV pointers differ); fall back to name check
            grp = ClassifyByPVName(prePV);
            if (grp == SiPMGroup::Unknown) grp = ClassifyByPVName(postPV);
        }
    } else {
        grp = ClassifyByPVName(prePV);
        if (grp == SiPMGroup::Unknown) grp = ClassifyByPVName(postPV);
    }
    G4String detName = "";
    if (grp == SiPMGroup::Crystal) {
        detName = "Crystal";
        ++npeCrystal;
        if (ch >= 0) ++perChCrystal[ch];
    } else if (grp == SiPMGroup::Veto) {
        detName = "Veto";
        ++npeVeto;
        if (ch >= 0) ++perChVeto[ch];
    } else if (grp == SiPMGroup::Bottom) {
        detName = "BottomVeto";
        ++npeBottom;
        if (ch >= 0) ++perChBottom[ch];
    } else {
        // Unknown classification: still kill photon to avoid infinite bouncing after "Detection"
        // but do not count it.
    }

    if (Configuration::savePhotons) {
        if (auto* ea = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction())) {
            PhotonRec rec;
            rec.photonID = track->GetTrackID();
            rec.detName = detName;
            rec.detCh = ch;
            rec.energy = track->GetTotalEnergy() / eV;
            rec.pos_mm = post->GetPosition();
            ea->photonBuf.emplace_back(std::move(rec));
        }
    }
    track->SetTrackStatus(fStopAndKill);
    return true;
}
