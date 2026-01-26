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
    auto* preLV = prePV ? prePV->GetLogicalVolume() : nullptr;
    auto* postLV = postPV ? postPV->GetLogicalVolume() : nullptr;

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

    // Helper: does this step touch a given LV on either side of the boundary?
    auto touchesLV = [&](G4LogicalVolume* lv) -> bool {
        return lv && (preLV == lv || postLV == lv);
    };

    // Classify primarily by LV pointers (stable), fall back to PV name only if needed.
    SiPMGroup grp = SiPMGroup::Unknown;

    // If you have a dedicated LV for the big crystal SiPM window, treat it as Crystal.
    if (touchesLV(crystalSiPMWindowLV)) {
        grp = SiPMGroup::Crystal;
    }
    // If you use a single generic SiPMWindowLV for an array (and names distinguish groups),
    // keep it as a fallback path.
    else if (touchesLV(SiPMWindowLV)) {
        grp = ClassifyByPVName(prePV);
        if (grp == SiPMGroup::Unknown) grp = ClassifyByPVName(postPV);
    } else {
        // Not on expected window LVs; still try name fallback (useful during refactors)
        grp = ClassifyByPVName(prePV);
        if (grp == SiPMGroup::Unknown) grp = ClassifyByPVName(postPV);
    }

    // Count detected photoelectrons
    switch (grp) {
    case SiPMGroup::Crystal:
        ++npeCrystal;
        if (ch >= 0) ++perChCrystal[ch];
        break;
    case SiPMGroup::Veto:
        ++npeVeto;
        if (ch >= 0) ++perChVeto[ch];
        break;
    case SiPMGroup::Bottom:
        ++npeBottom;
        if (ch >= 0) ++perChBottom[ch];
        break;
    case SiPMGroup::Unknown:
    default:
        break;
    }

    track->SetTrackStatus(fStopAndKill);
    return true;
}
