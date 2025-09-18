#include "SensitiveDetector.hh"


SensitiveDetector::SensitiveDetector(const G4String &sdName, G4int detID, const G4String &detName)
    : G4VSensitiveDetector(sdName), detID(detID), detName(detName) {
    collectionName.insert("EdepHits");
}

void SensitiveDetector::Initialize(G4HCofThisEvent *hce) {
    hits = new SDHitCollection(SensitiveDetectorName, collectionName[0]);
    indexByVol.clear();

    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hits);
    }
    hce->AddHitsCollection(HCID, hits);
}

G4bool SensitiveDetector::ProcessHits(G4Step *step, G4TouchableHistory *) {
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