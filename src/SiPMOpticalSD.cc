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
    for (int i = 0; i < plist->size(); ++i) {
        auto* p = (*plist)[i];
        if (p && p->GetProcessName() == "OpBoundary") {
            boundary = static_cast<G4OpBoundaryProcess*>(p);
            break;
        }
    }
    return boundary;
}

G4bool SiPMOpticalSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    auto* track = step->GetTrack();
    if (!track) return false;

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

    if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return false;

    auto* b = GetBoundaryProcess();
    if (!b) return false;

    if (b->GetStatus() != Detection) return false;

    const auto& touch = step->GetPreStepPoint()->GetTouchableHandle();
    const int ch = touch->GetCopyNumber(0);

    auto* pv = step->GetPreStepPoint()->GetPhysicalVolume();
    const auto pvName = pv ? pv->GetName() : "";

    if (pvName.find("CrystalSiPMWindowPVP") != std::string::npos) {
        ++npeCrystal;
        ++perChCrystal[ch];
    }
    else if (pvName.find("VetoSiPMWindowPVP") != std::string::npos) {
        ++npeVeto;
        ++perChVeto[ch];
    }
    else if (pvName.find("BottomVetoSiPMWindowPVP") != std::string::npos) {
        ++npeBottom;
        ++perChBottom[ch];
    }

    track->SetTrackStatus(fStopAndKill);
    return true;
}
