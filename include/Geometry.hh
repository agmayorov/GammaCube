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
#include <G4VSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SDManager.hh>
#include <utility>

#include "SensitiveDetector.hh"
#include "Detector.hh"
#include "Sizes.hh"


class Geometry : public G4VUserDetectorConstruction {
public:
    Geometry(G4String);

    G4VPhysicalVolume *Construct() override;
    void ConstructSDandField() override;

private:
    G4NistManager *nist;
    G4Material *worldMat;
    G4double worldHalfSize;
    G4Box *worldBox;
    G4LogicalVolume *worldLV;
    G4VPhysicalVolume *worldPVP;

    G4double viewDeg;

    Detector *detector;
    G4String detectorType;
    G4ThreeVector detContainerSize;

    G4Tubs *detContainer;
    G4LogicalVolume *detContainerLV;
    G4ThreeVector detContainerPos;

    G4LogicalVolume *crystalLV;
    G4LogicalVolume *tyvekInLV;
    G4LogicalVolume *AlLV;
    G4LogicalVolume *rubberLV;
    G4LogicalVolume *tyvekMidLV;
    G4LogicalVolume *vetoLV;
    G4LogicalVolume *tyvekOutLV;
    G4LogicalVolume *tunaCanLV;
    G4LogicalVolume *vetoLEDLV;
    G4LogicalVolume *crystalLEDLV;
    G4LogicalVolume *vetoBottomLEDLV;

    G4Region *vetoTopRegion;
    G4Region *vetoBottomRegion;
    G4Region *vetoLeftRegion;
    G4Region *vetoRightRegion;
    G4Region *vetoFrontRegion;
    G4Region *vetoBackRegion;

    G4RotationMatrix *zeroRot;

    G4VSolid *tunaCan;

    G4UserLimits *vetoStepLimit;
    G4Region *vetoRegion;
    G4ProductionCuts *vetoCuts;

    G4VisAttributes *tunaCanVisAttr;

    void CheckSizes();

    void ConstructDetector();
    void ConstructTunaCan();

    void SetStepLimits();
};

#endif //GEOMETRY_HH
