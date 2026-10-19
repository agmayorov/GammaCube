#include "SiPMOpticalSD.hh"

SiPMOpticalSD::SiPMOpticalSD(const G4String& name)
    : G4VSensitiveDetector(name) {}

void SiPMOpticalSD::Initialize(G4HCofThisEvent*) {
    npeCrystal = npeVeto = npeBottom = 0;
    perChCrystal.clear();
    perChVeto.clear();
    perChBottom.clear();
    cellLastFire.clear();
}

// A detected photon lands on a random microcell. The cell fires only if it has
// recharged since its last avalanche (probability 1 - exp(-dt/tau)), which is
// what makes a SiPM saturate. Each avalanche may then set off its neighbours,
// modelled as a geometric cascade with probability Configuration::sipmCrosstalk.
// Photons are not tracked in time order, so the gap between two hits on one cell
// is taken as |dt|: for a pair of hits that gives the same suppression probability,
// only the roles of the two are swapped. Crosstalk cells are not booked as occupied,
// which is harmless at the per-cent occupancies this detector reaches.
G4int SiPMOpticalSD::FireCell(const SiPMGroup grp, const G4int ch, const G4double t) {
    using namespace Configuration;

    if (sipmCells > 0) {
        auto& cells = cellLastFire[static_cast<G4int>(grp) * 1024 + ch];
        const G4int cell = static_cast<G4int>(G4UniformRand() * sipmCells);

        if (const auto it = cells.find(cell); it != cells.end()) {
            const G4double dt = std::abs(t - it->second);
            const G4double recovered = sipmRecoveryTime > 0 ? 1. - std::exp(-dt / sipmRecoveryTime) : 0.;
            if (G4UniformRand() >= recovered) return 0;   // cell is still recharging
            it->second = t;
        } else {
            cells.emplace(cell, t);
        }
    }

    G4int fired = 1;
    if (sipmCrosstalk > 0.) {
        while (G4UniformRand() < sipmCrosstalk) ++fired;
    }

    return fired;
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
    const G4int fired = grp != SiPMGroup::Unknown ? FireCell(grp, ch, track->GetGlobalTime()) : 0;

    G4String detName = "";
    if (grp == SiPMGroup::Crystal) {
        detName = "Crystal";
        npeCrystal += fired;
        if (ch >= 0) perChCrystal[ch] += fired;
    } else if (grp == SiPMGroup::Veto) {
        detName = "Veto";
        npeVeto += fired;
        if (ch >= 0) perChVeto[ch] += fired;
    } else if (grp == SiPMGroup::Bottom) {
        detName = "BottomVeto";
        npeBottom += fired;
        if (ch >= 0) perChBottom[ch] += fired;
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
