#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <G4VisAttributes.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4UnionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4Version.hh>
#include <G4LossTableManager.hh>
#include <G4Element.hh>
#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4VSolid.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>

#include "Sizes.hh"


class Detector {
public:
    G4NistManager *nist;
    G4String detectorType;

    G4LogicalVolume *detContainerLV;
    G4VPhysicalVolume *detContainerPVP;
    G4double detContainerSizeZ;

    G4VisAttributes *visTyvekOut{};
    G4VisAttributes *visTyvekMid{};
    G4VisAttributes *visTyvekIn{};
    G4VisAttributes *visCrystal{};
    G4VisAttributes *visVeto{};
    G4VisAttributes *visLED{};
    G4VisAttributes *visAl{};
    G4VisAttributes *visWire{};
    G4VisAttributes *visBoard{};
    G4VisAttributes *visRubber{};

    Detector(G4LogicalVolume *, G4VPhysicalVolume *, G4double, G4NistManager *, G4double, const G4String &, G4bool);
    ~Detector() = default;

    void DefineMaterials();
    void DefineVisual();

    void Construct();

    [[nodiscard]] G4String GetDetectorType() const;
    [[nodiscard]] std::vector<G4LogicalVolume *>GetSensitiveLV() const;

private:
    G4MultiUnion *crystalLED;
    G4MultiUnion *elongatedLED;
    G4MultiUnion *vetoLED;
    G4MultiUnion *vetoBottomLED;
    G4MultiUnion *vetoLEDWire;
    G4MultiUnion *vetoLEDWireInsulation;
    G4MultiUnion *vetoLEDFictive;
    G4MultiUnion *vetoWireFictive;
    G4MultiUnion *vetoBottomLEDFictive;

    G4double viewDeg;

    G4bool lightCollection;

    G4double zCorrection;

    G4double additionalLength;
    G4double sensSurfThick;

    G4ThreeVector crystalSize;
    G4ThreeVector vetoSize;
    G4ThreeVector tyvekInSize;
    G4ThreeVector tyvekMidSize;
    G4ThreeVector tyvekOutSize;

    G4Material *vetoMat{};
    G4Material *CrystalMat{};
    G4Material *tyvekMat{};
    G4Material *LEDMat{};
    G4Material *AlMat{};
    G4Material *rubberMat{};
    G4Material *wireMat{};
    G4Material *boardMat{};
    G4Material *galacticMat{};

    G4LogicalVolume *crystalLV{};
    G4LogicalVolume *tyvekInLV{};
    G4LogicalVolume *AlLV{};
    G4LogicalVolume *rubberLV{};
    G4LogicalVolume *tyvekMidLV{};
    G4LogicalVolume *vetoLV{};
    G4LogicalVolume *tyvekOutLV{};
    G4LogicalVolume *crystalLEDLV{};
    G4LogicalVolume *vetoLEDLV{};
    G4LogicalVolume *vetoBottomLEDLV{};

    G4LogicalVolume *crystalSensSurfLV{};
    G4LogicalVolume *vetoSensSurfLV{};
    G4LogicalVolume *vetoBottomSensSurfLV{};

    G4VPhysicalVolume *crystalPVP{};
    G4VPhysicalVolume *tyvekInPVP{};
    G4VPhysicalVolume *AlPVP{};
    G4VPhysicalVolume *rubberPVP{};
    G4VPhysicalVolume *tyvekMidPVP{};
    G4VPhysicalVolume *vetoPVP{};
    G4VPhysicalVolume *tyvekOutPVP{};
    G4VPhysicalVolume *crystalLEDPVP{};
    G4VPhysicalVolume *vetoLEDPVP{};
    G4VPhysicalVolume *vetoBottomLEDPVP{};

    G4VPhysicalVolume *crystalSensSurfPVP{};
    G4VPhysicalVolume *vetoSensSurfPVP{};
    G4VPhysicalVolume *vetoBottomSensSurfPVP{};

    void ConstructCrystal();
    void ConstructAl();
    void ConstructVeto();

    void ConstructCrystalLED();
    void ConstructVetoLED();
    void ConstructVetoBottomLED();

    void AddOptics();
};


#endif //DETECTOR_HH
