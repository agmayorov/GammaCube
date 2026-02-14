#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include <G4UserEventAction.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4Event.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4SystemOfUnits.hh>
#include <cfloat>
#include <vector>

#include "Geometry.hh"
#include "RunAction.hh"
#include "Sizes.hh"
#include "Configuration.hh"
#include "AnalysisManager.hh"
#include "SDHit.hh"
#include "SiPMOpticalSD.hh"

class G4Event;
class Geometry;

struct PrimaryRec {
    int index = 0;
    int pdg = 0;
    G4String name;
    double E_MeV = 0.0;
    G4ThreeVector dir;
    G4ThreeVector pos_mm;  // mm
    double t0_ns = 0.0;    // ns
};

struct InteractionRec {
    int trackID = -1;
    int parentID = -1;
    G4String process;
    G4String volumeName;
    int volumeID = -1;
    G4ThreeVector pos_mm;  // mm
    double t_ns = 0.0;     // ns

    int secIndex = -1;
    int secPDG = 0;
    G4String secName;
    double secE_MeV = 0.0;
    G4ThreeVector secDir;
};

struct PhotonRec {
    G4int photonID = -1;
    G4String detName;
    G4int detCh;
    G4double energy;
    G4ThreeVector pos_mm;  // mm
};

class EventAction : public G4UserEventAction {
public:
    std::vector<PrimaryRec> primBuf;
    std::vector<InteractionRec> interBuf;
    std::vector<G4int> photonCountBuf{0, 0, 0};
    std::vector<PhotonRec> photonBuf;

    EventAction(AnalysisManager *, RunAction *);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event *) override;
    void EndOfEventAction(const G4Event *) override;

private:
    void WritePrimaries_(int eventID);
    int WriteInteractions_(int eventID);
    int WritePhotonsCount_(int eventID);
    int WritePhotons_(int eventID);
    int WriteEdepFromSD_(const G4Event *evt, int eventID);

    void WriteSiPMFromSD_(int eventID);

    void MarkCrystal() { hasCrystal = true; }
    void MarkVeto() { hasVeto = true; }

    void MarkCrystalOpt() { hasCrystalOpt = true; }
    void MarkVetoOpt() { hasVetoOpt = true; }

    AnalysisManager *analysisManager = nullptr;

    std::vector<std::tuple<G4String, int, G4String>> detMap;
    std::vector<std::tuple<G4String, int, G4String>> optMap;
    std::vector<int> optHCIDs;
    std::vector<int> HCIDs;

    int nPrimaries = 0;
    int nInteractions = 0;
    int nPhotons = 0;
    int nEdepHits = 0;

    RunAction* run = nullptr;
    bool hasCrystal = false;
    bool hasVeto = false;
    bool hasCrystalOpt = false;
    bool hasVetoOpt = false;
};

#endif //EVENTACTION_HH
