#include "SensitiveDetector.hh"


SensitiveDetector::SensitiveDetector(const G4String &sdName, G4int detID, const G4String &detName, const G4bool isOpt)
    : G4VSensitiveDetector(sdName), detID(detID), detName(detName), isLED(isOpt) {
    collectionName.insert("EdepHits");
    if (isLED) collectionName.insert("OptHits");
}

void SensitiveDetector::Initialize(G4HCofThisEvent *hce) {
    hits = new SDHitCollection(SensitiveDetectorName, collectionName[0]);
    indexByVol.clear();

    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hits);
    }
    hce->AddHitsCollection(HCID, hits);

    if (isLED) {
        optHC = new SDHitCollection(SensitiveDetectorName, "OptHits");
        auto optID = G4SDManager::GetSDMpointer()->GetCollectionID(optHC);
        hce->AddHitsCollection(optID, optHC);
    } else {
        optHC = nullptr;
    }
}

G4bool SensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *) {
    auto pd = step->GetTrack()->GetDefinition();

    if (isLED and pd == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if (optHC && step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
            auto *post = step->GetPostStepPoint();
            auto touch = post->GetTouchableHandle();

            // глобальная -> локальная точка входа
            const auto T = touch->GetHistory()->GetTopTransform(); // local->global
            const auto x_loc = T.Inverse().TransformPoint(post->GetPosition());

            auto *h = new SDHit();
            h->isOptical = true;
            h->volumeID = detID; // или copyNo, если хотите группировать иначе
            h->x_loc_mm = x_loc.x() / mm;
            h->y_loc_mm = x_loc.y() / mm;
            h->z_loc_mm = x_loc.z() / mm;
            h->phi = std::atan2(x_loc.y(), x_loc.x());
            h->t_ns = post->GetGlobalTime() / ns;
            h->copyNo = touch->GetVolume()->GetCopyNo();

            optHC->insert(h);
            return true;
        }
        return false;
    }


    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return false;

    const G4VTouchable *touch = step->GetPreStepPoint()->GetTouchable();
    const int volumeID = touch->GetVolume()->GetCopyNo();

    const G4double t = step->GetPreStepPoint()->GetGlobalTime();

    SDHit *hit = FindOrCreateHit(volumeID);
    hit->AddEdep(edep);
    hit->UpdateTmin(t);

    return true;
}

SDHit *SensitiveDetector::FindOrCreateHit(G4int volumeID) {
    auto it = indexByVol.find(volumeID);
    if (it != indexByVol.end()) {
        return (*hits)[it->second];
    }

    auto *h = new SDHit(volumeID);
    int idx = hits->insert(h) - 1;
    indexByVol.emplace(volumeID, idx);
    return h;
}