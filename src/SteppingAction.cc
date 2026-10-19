#include "SteppingAction.hh"


// Readout gate: an optical photon that is still alive photonTimeLimit after the
// event started depositing energy is dropped, exactly as the SiPM electronics
// would ignore it. The gate opens at the first energy deposit of the event so
// that a delayed decay is timed from its own interaction, not from t = 0.
// Known limitation: with two bursts in one event the gate opens on the first one.
void SteppingAction::ApplyPhotonTimeLimit(const G4Step* step, EventAction* ea) {
    auto* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
        if (ea->eventT0 < 0 and step->GetTotalEnergyDeposit() > 0.)
            ea->eventT0 = step->GetPreStepPoint()->GetGlobalTime();
        return;
    }

    if (track->GetCurrentStepNumber() == 1) ++ea->trackedPhotons;

    // No deposit seen yet (should not happen): time the photon from its own birth
    if (ea->eventT0 < 0) ea->eventT0 = track->GetGlobalTime() - track->GetLocalTime();

    if (track->GetGlobalTime() - ea->eventT0 > Configuration::photonTimeLimit) {
        track->SetTrackStatus(fStopAndKill);
        ++ea->killedPhotons;
    }
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    auto* ea = static_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction());

    if (Configuration::photonTimeLimit > 0 and ea) {
        ApplyPhotonTimeLimit(step, ea);
    }

    if (!Configuration::saveSecondaries and !Configuration::savePhotons) return;

    const auto* post = step->GetPostStepPoint();
    const auto* postProc = post->GetProcessDefinedStep();

    const G4ThreeVector x = post->GetPosition() / mm;
    const G4double t = post->GetGlobalTime() / ns;

    const auto* touch = post->GetTouchable();
    G4String volName = "World";
    int copyNo = -1;
    if (touch && touch->GetVolume()) {
        volName = touch->GetVolume()->GetName();
        copyNo = touch->GetVolume()->GetCopyNo();
    }

    const auto* track = step->GetTrack();
    const G4int trackID = track->GetTrackID();
    const G4int parentID = track->GetParentID();

    auto secs = step->GetSecondaryInCurrentStep();
    if (ea) {
        if (secs && !secs->empty()) {
            for (size_t i = 0; i < secs->size(); ++i) {
                const auto* sc = (*secs)[i];
                const auto* cp = sc->GetCreatorProcess();

                if (sc->GetDefinition()->GetParticleName() == "opticalphoton" and Configuration::savePhotons) {
                    if (volName == "CrystalPVP")
                        ea->photonCountBuf[0] += 1;
                    if (volName == "VetoPVP")
                        ea->photonCountBuf[1] += 1;
                    if (volName == "BottomVetoPVP")
                        ea->photonCountBuf[2] += 1;
                }
                InteractionRec rec;
                rec.trackID = trackID;
                rec.parentID = parentID;
                rec.process = cp ? cp->GetProcessName() : "unknown";
                rec.volumeName = volName;
                rec.volumeID = copyNo;
                rec.pos_mm = x;
                rec.t_ns = t;
                rec.secIndex = static_cast<G4int>(i);
                rec.secPDG = sc->GetDefinition()->GetPDGEncoding();
                rec.secName = sc->GetDefinition()->GetParticleName();
                rec.secE_MeV = sc->GetKineticEnergy() / MeV;
                rec.secDir = sc->GetMomentumDirection();
                ea->interBuf.emplace_back(std::move(rec));
            }
            return;
        }

        if (postProc) {
            const auto& pname = postProc->GetProcessName();
            if (pname != "Transportation") {
                InteractionRec rec;
                rec.trackID = trackID;
                rec.parentID = parentID;
                rec.process = pname;
                rec.volumeName = volName;
                rec.volumeID = copyNo;
                rec.pos_mm = x;
                rec.t_ns = t;
                rec.secIndex = -1;
                rec.secPDG = 0;
                rec.secName = "";
                rec.secE_MeV = 0.0;
                rec.secDir = G4ThreeVector();
                ea->interBuf.emplace_back(std::move(rec));
            }
        }
    }
}
