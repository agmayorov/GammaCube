#include "AnalysisManager.hh"

using namespace Sizes;

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
    analysisManager->CreateNtupleSColumn("primary_name");
    analysisManager->CreateNtupleDColumn("E_MeV");
    analysisManager->CreateNtupleDColumn("dir_x");
    analysisManager->CreateNtupleDColumn("dir_y");
    analysisManager->CreateNtupleDColumn("dir_z");
    analysisManager->CreateNtupleDColumn("pos_x_mm");
    analysisManager->CreateNtupleDColumn("pos_y_mm");
    analysisManager->CreateNtupleDColumn("pos_z_mm");
    analysisManager->FinishNtuple(primaryNT);

    interactionsNT = analysisManager->CreateNtuple("interactions",
                                                   "inelastic/compton/photo/conv vertices and secondaries");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleIColumn("trackID");
    analysisManager->CreateNtupleIColumn("parentID");
    analysisManager->CreateNtupleSColumn("process");
    analysisManager->CreateNtupleSColumn("volume_name");
    analysisManager->CreateNtupleDColumn("x_mm");
    analysisManager->CreateNtupleDColumn("y_mm");
    analysisManager->CreateNtupleDColumn("z_mm");
    analysisManager->CreateNtupleDColumn("t_ns");
    analysisManager->CreateNtupleIColumn("sec_index");
    analysisManager->CreateNtupleSColumn("sec_name");
    analysisManager->CreateNtupleDColumn("sec_E_MeV");
    analysisManager->CreateNtupleDColumn("sec_dir_x");
    analysisManager->CreateNtupleDColumn("sec_dir_y");
    analysisManager->CreateNtupleDColumn("sec_dir_z");
    analysisManager->FinishNtuple(interactionsNT);

    edepNT = analysisManager->CreateNtuple("edep", "energy deposition per sensitive channel");
    analysisManager->CreateNtupleIColumn("eventID");
    analysisManager->CreateNtupleSColumn("det_name");
    analysisManager->CreateNtupleDColumn("edep_MeV");
    analysisManager->FinishNtuple(edepNT);

    crystalH2 = analysisManager->CreateH2("CrystalH2",
                                          "Crystal;X [mm];Y [mm]",
                                          200, -20., 20.,
                                          200, -20., 20.);

    vetoH2 = analysisManager->CreateH2("VetoH2",
                                       "Veto;X [mm];Y [mm]",
                                       240, -31., 31.,
                                       240, -31., 31.);

    vetoBottomH2 = analysisManager->CreateH2("VetoBottomH2",
                                             "VetoBottom side;#phi [rad]:Z [mm]",
                                             180, 0.0, CLHEP::twopi,
                                             200, -45, 10
    );
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

void AnalysisManager::FillPrimaryRow(G4int eventID, const G4String &primaryName,
                                     G4double E_MeV, const G4ThreeVector &dir,
                                     const G4ThreeVector &pos_mm) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(primaryNT, 0, eventID);
    analysisManager->FillNtupleSColumn(primaryNT, 1, primaryName);
    analysisManager->FillNtupleDColumn(primaryNT, 2, E_MeV);
    analysisManager->FillNtupleDColumn(primaryNT, 3, dir.x());
    analysisManager->FillNtupleDColumn(primaryNT, 4, dir.y());
    analysisManager->FillNtupleDColumn(primaryNT, 5, dir.z());
    analysisManager->FillNtupleDColumn(primaryNT, 6, pos_mm.x());
    analysisManager->FillNtupleDColumn(primaryNT, 7, pos_mm.y());
    analysisManager->FillNtupleDColumn(primaryNT, 8, pos_mm.z());
    analysisManager->AddNtupleRow(primaryNT);
}

void AnalysisManager::FillInteractionRow(G4int eventID,
                                         G4int trackID, G4int parentID,
                                         const G4String &process,
                                         const G4String &volumeName,
                                         const G4ThreeVector &x_mm,
                                         G4int secIndex, const G4String &secName,
                                         G4double secE_MeV, const G4ThreeVector &secDir) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(interactionsNT, 0, eventID);
    analysisManager->FillNtupleIColumn(interactionsNT, 1, trackID);
    analysisManager->FillNtupleIColumn(interactionsNT, 2, parentID);
    analysisManager->FillNtupleSColumn(interactionsNT, 3, process);
    analysisManager->FillNtupleSColumn(interactionsNT, 4, volumeName);
    analysisManager->FillNtupleDColumn(interactionsNT, 5, x_mm.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 6, x_mm.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 7, x_mm.z());
    analysisManager->FillNtupleIColumn(interactionsNT, 8, secIndex);
    analysisManager->FillNtupleSColumn(interactionsNT, 9, secName);
    analysisManager->FillNtupleDColumn(interactionsNT, 10, secE_MeV);
    analysisManager->FillNtupleDColumn(interactionsNT, 11, secDir.x());
    analysisManager->FillNtupleDColumn(interactionsNT, 12, secDir.y());
    analysisManager->FillNtupleDColumn(interactionsNT, 13, secDir.z());
    analysisManager->AddNtupleRow(interactionsNT);
}

void AnalysisManager::FillEdepRow(G4int eventID, const G4String &det_name, G4double edep_MeV) {
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleIColumn(edepNT, 0, eventID);
    analysisManager->FillNtupleSColumn(edepNT, 1, det_name);
    analysisManager->FillNtupleDColumn(edepNT, 2, edep_MeV);
    analysisManager->AddNtupleRow(edepNT);
}

void AnalysisManager::FillCrystalH2(const G4double x_mm, const G4double y_mm) {
    G4AnalysisManager::Instance()->FillH2(crystalH2, x_mm, y_mm);
}

void AnalysisManager::FillVetoH2(const G4double x_mm, const G4double y_mm) {
    G4AnalysisManager::Instance()->FillH2(vetoH2, x_mm, y_mm);
}

void AnalysisManager::FillVetoBottomH2(const G4double z_mm, G4double phi_rad) {
    while (phi_rad < 0) phi_rad += CLHEP::twopi;
    while (phi_rad >= CLHEP::twopi) phi_rad -= CLHEP::twopi;
    G4AnalysisManager::Instance()->FillH2(vetoBottomH2, phi_rad, z_mm);
}
