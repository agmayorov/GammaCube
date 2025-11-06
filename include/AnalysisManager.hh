#ifndef ANALYSISMANAGER_HH
#define ANALYSISMANAGER_HH

#include <G4AnalysisManager.hh>
#include <G4UnitsTable.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <globals.hh>
#include <Sizes.hh>

class AnalysisManager {
public:
    G4String fileName = "GammaDetector";

    explicit AnalysisManager(const std::string &);
    ~AnalysisManager() = default;

    void Open();
    void Close();

    void FillEventRow(G4int eventID, G4int nPrimaries, G4int nInteractions, G4int nEdepHits);

    void FillPrimaryRow(G4int eventID, const G4String &primaryName,
                        G4double E_MeV, const G4ThreeVector &dir,
                        const G4ThreeVector &pos_mm);

    void FillInteractionRow(G4int eventID,
                            G4int trackID, G4int parentID,
                            const G4String &process,
                            const G4String &volumeName,
                            const G4ThreeVector &x_mm,
                            G4int secIndex, const G4String &secName,
                            G4double secE_MeV, const G4ThreeVector &secDir);

    void FillEdepRow(G4int eventID, const G4String &det_name, G4double edep_MeV);

    void FillCrystalH2(G4double x_mm, G4double y_mm);
    void FillVetoBottomH2(G4double z_mm, G4double phi_rad);
    void FillVetoH2(G4double x_mm, G4double y_mm);

private:
    G4int eventNT{-1};
    G4int primaryNT{-1};
    G4int interactionsNT{-1};
    G4int edepNT{-1};

    G4int crystalH2{-1};
    G4int vetoH2{-1};
    G4int vetoBottomH2{-1};

    void Book();
};


#endif //ANALYSISMANAGER_HH
