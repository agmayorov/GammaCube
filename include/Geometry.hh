#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4Box.hh>
#include <G4GeometryManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisAttributes.hh>
#include <G4ProductionCuts.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SDManager.hh>
#include <utility>

#include "SensitiveDetector.hh"
#include "Detector.hh"
#include "Sizes.hh"


class Geometry : public G4VUserDetectorConstruction {
public:
    G4NistManager *nist;
    G4Material *worldMat;
    G4double worldHalfSize;
    G4Box *worldBox;
    G4LogicalVolume *worldLV;
    G4VPhysicalVolume *worldPVP;

    G4double viewDeg;
    G4bool doubleLED;

    Detector *detector;
    G4String detectorType;
    G4ThreeVector detContainerSize;
    G4ThreeVector modelSize;

    G4Tubs *detContainer;
    G4LogicalVolume *detContainerLV;
    G4ThreeVector detContainerPos;

    Sizes sizes;

    G4LogicalVolume *crystalLV;
    G4LogicalVolume *vetoLV;
    G4LogicalVolume *tyvekOutLV;
    G4LogicalVolume *tyvekInLV;

    G4Region *vetoTopRegion;
    G4Region *vetoBottomRegion;
    G4Region *vetoLeftRegion;
    G4Region *vetoRightRegion;
    G4Region *vetoFrontRegion;
    G4Region *vetoBackRegion;

    Geometry(G4String, const Sizes &, G4double, G4bool);

    G4VPhysicalVolume *Construct() override;
    void ConstructSDandField() override;

private:
    G4RotationMatrix *zeroRot;

    G4UnionSolid *shell;
    G4LogicalVolume *shellLV;
    G4UnionSolid *tunaCan;
    G4LogicalVolume *tunaCanLV;
    G4double temperature;

    G4UserLimits *vetoStepLimit;
    G4Region *vetoRegion;
    G4ProductionCuts *vetoCuts;

    G4VisAttributes *shellVisAttr;
    G4VisAttributes *tunaCanVisAttr;
    G4VisAttributes *LEDVisAttr;

    void ConstructDetector();
    void ConstructTunaCan();
    void ConstructLED();
    void SetStepLimits();
};

#endif //GEOMETRY_HH
