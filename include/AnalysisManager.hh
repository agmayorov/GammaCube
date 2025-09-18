#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <G4AnalysisManager.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

class AnalysisManager {
public:
    G4String fileName = "GammaDetector";

    explicit AnalysisManager(const std::string &);
    ~AnalysisManager() = default;

    void Open();
    void Close();

    void FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits);

    void FillPrimaryRow(G4int eventID, G4int primaryIndex,
                        G4int primaryPDG, const G4String &primaryName,
                        G4double E_MeV, const G4ThreeVector &dir,
                        const G4ThreeVector &pos_mm, G4double t0_ns);

    void FillInteractionRow(G4int eventID,
                            G4int trackID, G4int parentID,
                            const G4String &process,
                            const G4String &volumeName, G4int volumeID,
                            const G4ThreeVector &x_mm, G4double t_ns,
                            G4int secIndex, G4int secPDG, const G4String &secName,
                            G4double secE_MeV, const G4ThreeVector &secDir);

    void FillEdepRow(G4int eventID,
                     G4int det_id, const G4String &det_name,
                     G4int volumeID, G4double edep_MeV, G4double tmin_ns);


private:
    G4int eventNT{-1};
    G4int primaryNT{-1};
    G4int interactionsNT{-1};
    G4int edepNT{-1};

    void Book();
};


#endif //ANALYSISMANAGER_HH
