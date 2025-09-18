#include "AnalysisManager.hh"


AnalysisManager::AnalysisManager(const std::string &fName) : fileName(fName) {
    Book();
}


void AnalysisManager::Book() {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFileName(fileName);
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetNtupleActivation(true);

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#endif

    eventNT = analysisManager->CreateNtuple("event", "per-event summary");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("n_primaries");
    analysisManager->CreateNtupleIColumn("n_interactions");
    analysisManager->CreateNtupleIColumn("n_edep_hits");
    analysisManager->FinishNtuple(eventNT);

    primaryNT = analysisManager->CreateNtuple("primary", "per-primary particles");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("primary_index");
    analysisManager->CreateNtupleIColumn("primary_pdg");
    analysisManager->CreateNtupleSColumn("primary_name");
    analysisManager->CreateNtupleDColumn("E_MeV");
    analysisManager->CreateNtupleDColumn("dir_x");
    analysisManager->CreateNtupleDColumn("dir_y");
    analysisManager->CreateNtupleDColumn("dir_z");
    analysisManager->CreateNtupleDColumn("pos_x_mm");
    analysisManager->CreateNtupleDColumn("pos_y_mm");
    analysisManager->CreateNtupleDColumn("pos_z_mm");
    analysisManager->CreateNtupleDColumn("t0_ns");
    analysisManager->FinishNtuple(primaryNT);

    interactionsNT = analysisManager->CreateNtuple("interactions",
                                                   "inelastic/compton/photo/conv vertices and secondaries");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("trackID");
    analysisManager->CreateNtupleIColumn("parentID");
    analysisManager->CreateNtupleSColumn("process");
    analysisManager->CreateNtupleSColumn("volume_name");
    analysisManager->CreateNtupleIColumn("volumeID");
    analysisManager->CreateNtupleDColumn("x_mm");
    analysisManager->CreateNtupleDColumn("y_mm");
    analysisManager->CreateNtupleDColumn("z_mm");
    analysisManager->CreateNtupleDColumn("t_ns");
    analysisManager->CreateNtupleIColumn("sec_index");
    analysisManager->CreateNtupleIColumn("sec_pdg");
    analysisManager->CreateNtupleSColumn("sec_name");
    analysisManager->CreateNtupleDColumn("sec_E_MeV");
    analysisManager->CreateNtupleDColumn("sec_dir_x");
    analysisManager->CreateNtupleDColumn("sec_dir_y");
    analysisManager->CreateNtupleDColumn("sec_dir_z");
    analysisManager->FinishNtuple(interactionsNT);

    edepNT = analysisManager->CreateNtuple("edep", "energy deposition per sensitive channel");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("det_id"); // 0=NaI, 1=veto, ...
    analysisManager->CreateNtupleSColumn("det_name"); // "NaI","veto",...
    analysisManager->CreateNtupleIColumn("volumeID"); // copyNo (или свой канал)
    analysisManager->CreateNtupleDColumn("edep_MeV");
    analysisManager->CreateNtupleDColumn("tmin_ns");
    analysisManager->FinishNtuple(edepNT);
}

void AnalysisManager::Open() {
    G4AnalysisManager::Instance()->OpenFile(fileName);
}

void AnalysisManager::Close() {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

void AnalysisManager::FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(eventNT, 0, eventID);
    analysisManager->FillNtupleIColumn(eventNT, 1, nPrimaries);
    analysisManager->FillNtupleIColumn(eventNT, 2, nInteractions);
    analysisManager->FillNtupleIColumn(eventNT, 3, nEdepHits);
    analysisManager->AddNtupleRow(eventNT);
}

void AnalysisManager::FillPrimaryRow(G4int eventID, G4int primaryIndex,
                                     G4int primaryPDG, const G4String &primaryName,
                                     G4double E_MeV, const G4ThreeVector &dir,
                                     const G4ThreeVector &pos_mm, G4double t0_ns) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(primaryNT, 0, eventID);
    analysisManager->FillNtupleIColumn(primaryNT, 1, primaryIndex);
    // analysisManager->FillNtupleIColumn(primaryNT, 2, primaryPDG);
    analysisManager->FillNtupleSColumn(primaryNT, 3, primaryName);
    analysisManager->FillNtupleDColumn(primaryNT, 4, E_MeV);
    analysisManager->FillNtupleDColumn(primaryNT, 5, dir.x());
    analysisManager->FillNtupleDColumn(primaryNT, 6, dir.y());
    analysisManager->FillNtupleDColumn(primaryNT, 7, dir.z());
    analysisManager->FillNtupleDColumn(primaryNT, 8, pos_mm.x());
    analysisManager->FillNtupleDColumn(primaryNT, 9, pos_mm.y());
    analysisManager->FillNtupleDColumn(primaryNT, 10, pos_mm.z());
    // analysisManager->FillNtupleDColumn(primaryNT, 11, t0_ns);
    analysisManager->AddNtupleRow(primaryNT);
}

void AnalysisManager::FillInteractionRow(G4int eventID,
                                         G4int trackID, G4int parentID,
                                         const G4String &process,
                                         const G4String &volumeName, G4int volumeID,
                                         const G4ThreeVector &x_mm, G4double t_ns,
                                         G4int secIndex, G4int secPDG, const G4String &secName,
                                         G4double secE_MeV, const G4ThreeVector &secDir) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(interactionsNT, 0, eventID);
    analysisManager->FillNtupleIColumn(interactionsNT, 1, trackID);
    analysisManager->FillNtupleIColumn(interactionsNT, 2, parentID);
    analysisManager->FillNtupleSColumn(interactionsNT, 3, process);
    analysisManager->FillNtupleSColumn(interactionsNT, 4, volumeName);
    // analysisManager->FillNtupleIColumn(interactionsNT, 5, volumeID);
    analysisManager->FillNtupleDColumn(interactionsNT, 6, x_mm.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 7, x_mm.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 8, x_mm.z());
    // analysisManager->FillNtupleDColumn(interactionsNT, 9, t_ns);
    analysisManager->FillNtupleIColumn(interactionsNT, 10, secIndex);
    analysisManager->FillNtupleIColumn(interactionsNT, 11, secPDG);
    analysisManager->FillNtupleSColumn(interactionsNT, 12, secName);
    analysisManager->FillNtupleDColumn(interactionsNT, 13, secE_MeV);
    analysisManager->FillNtupleDColumn(interactionsNT, 14, secDir.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 15, secDir.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 16, secDir.z());
    analysisManager->AddNtupleRow(interactionsNT);
}

void AnalysisManager::FillEdepRow(G4int eventID, G4int det_id, const G4String &det_name,
                                  G4int volumeID, G4double edep_MeV, G4double tmin_ns) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(edepNT, 0, eventID);
    // analysisManager->FillNtupleIColumn(edepNT, 1, det_id);
    analysisManager->FillNtupleSColumn(edepNT, 2, det_name);
    // analysisManager->FillNtupleIColumn(edepNT, 3, volumeID);
    analysisManager->FillNtupleDColumn(edepNT, 4, edep_MeV);
    // analysisManager->FillNtupleDColumn(edepNT, 5, tmin_ns);
    analysisManager->AddNtupleRow(edepNT);
}
