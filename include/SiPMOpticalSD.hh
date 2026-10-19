#ifndef SIPMOPTICALSD_HH
#define SIPMOPTICALSD_HH

#include <G4VSensitiveDetector.hh>
#include <G4OpBoundaryProcess.hh>
#include <Randomize.hh>
#include <cmath>
#include <unordered_map>

#include <G4Step.hh>
#include <G4Track.hh>
#include <G4OpticalPhoton.hh>
#include <G4TouchableHistory.hh>
#include <G4VProcess.hh>
#include <G4ProcessManager.hh>
#include <G4VPhysicalVolume.hh>

#include "EventAction.hh"
#include "AnalysisManager.hh"

enum class SiPMGroup { Unknown, Crystal, Veto, Bottom };

class SiPMOpticalSD : public G4VSensitiveDetector {
public:
    explicit SiPMOpticalSD(const G4String& name);

    void SetSiPMWindowLV(G4LogicalVolume* lv) { SiPMWindowLV = lv; }

    void Initialize(G4HCofThisEvent*) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

    // getters for EventAction
    int GetNpeCrystal() const { return npeCrystal; }
    int GetNpeVeto() const { return npeVeto; }
    int GetNpeBottomVeto() const { return npeBottom; }

    const std::unordered_map<int,int>& GetPerChannelCrystal() const { return perChCrystal; }
    const std::unordered_map<int,int>& GetPerChannelVeto() const { return perChVeto; }
    const std::unordered_map<int,int>& GetPerChannelBottom() const { return perChBottom; }

private:
    G4OpBoundaryProcess* GetBoundaryProcess();

    // Microcell saturation and optical crosstalk; returns the number of fired cells
    G4int FireCell(SiPMGroup grp, G4int ch, G4double t);

    // (group, channel) -> (microcell -> time of its last avalanche)
    std::unordered_map<G4int, std::unordered_map<G4int, G4double>> cellLastFire;

    G4OpBoundaryProcess* boundary{nullptr};
    G4LogicalVolume* SiPMWindowLV{nullptr};

    int npeCrystal{0};
    int npeVeto{0};
    int npeBottom{0};

    std::unordered_map<int,int> perChCrystal;
    std::unordered_map<int,int> perChVeto;
    std::unordered_map<int,int> perChBottom;

    SiPMGroup ClassifyByPVName(const G4VPhysicalVolume* pv);
};

#endif // SIPMOPTICALSD_HH