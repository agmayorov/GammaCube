#include "geometry/TunaCan.hh"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVPlacement.hh>

#include "Configuration.hh"
#include "geometry/Sizes.hh"

using namespace Sizes;
using namespace Configuration;


TunaCan::TunaCan(G4LogicalVolume* mother, G4NistManager* nistManager, const G4Transform3D& placement,
                 const G4double innerRadius)
    : motherLV(mother), nist(nistManager), placement(placement), innerRadius(innerRadius) {
    tunaCanVisAttr = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
    tunaCanVisAttr->SetForceSolid(true);
}


void TunaCan::Place(G4LogicalVolume* lv, const G4ThreeVector& pos, const G4String& name) const {
    new G4PVPlacement(placement * G4Translate3D(pos), lv, name, motherLV, false, 0, true);
}


void TunaCan::Construct() {
    if (tunaCanThickWall <= 0. and tunaCanThickTop <= 0.) return;

    G4Material* tunaCanMat = nist->FindOrBuildMaterial("G4_Al");
    G4Material* plateMat = nist->FindOrBuildMaterial("G4_Ti");

    // TunaCan
    G4VSolid* tunaCanWall = new G4Tubs("tunaCanWall", innerRadius, modelRadius, modelHeight / 2, 0, viewDeg);
    G4ThreeVector tunaCanWallPos = G4ThreeVector(0, 0, 0);
    G4LogicalVolume* tunaCanWallLV = new G4LogicalVolume(tunaCanWall, tunaCanMat, "tunaCanWallLV");
    Place(tunaCanWallLV, tunaCanWallPos, "TunaCanWallPVPL");

    G4VSolid* tunaCanTop = new G4Tubs("TunaCanTop", 0, innerRadius, tunaCanThickTop / 2., 0, viewDeg);
    const G4ThreeVector tunaCanTopPos = G4ThreeVector(0, 0, (modelHeight - tunaCanThickTop) / 2.0);
    G4LogicalVolume* tunaCanTopLV = new G4LogicalVolume(tunaCanTop, tunaCanMat, "tunaCanTopLV");
    Place(tunaCanTopLV, tunaCanTopPos, "TunaCanTopPVPL");

    // Plate square part
    G4VSolid* plateIncomplete = new G4Box("PlateIncomplete", plateSize / 2., plateSize / 2. + plateCornerSize,
                                          plateThick / 2.);
    G4VSolid* plateHole = new G4Tubs("PlateHole", 0., plateOuterHoleRadius, plateThick + 5 * mm, 0., 360);
    const G4ThreeVector platePos = G4ThreeVector(0, 0, -(modelHeight + plateThick) / 2.0);
    G4VSolid* plate = new G4SubtractionSolid("Plate", plateIncomplete, plateHole);
    G4LogicalVolume* plateLV = new G4LogicalVolume(plate, plateMat, "PlateLV");
    Place(plateLV, platePos, "PlatePVPL");

    // Plate outer part
    G4VSolid* plateStrip = new G4Box("PlateStrip", plateCornerSize / 2., plateSize / 2., plateThick / 2.);
    const G4ThreeVector plateStripLeftPos = G4ThreeVector((plateSize + plateCornerSize) / 2., 0,
                                                          -(modelHeight + plateThick) / 2.0);
    const G4ThreeVector plateStripRightPos = G4ThreeVector(-(plateSize + plateCornerSize) / 2., 0,
                                                           -(modelHeight + plateThick) / 2.0);
    G4LogicalVolume* plateStripLeftLV = new G4LogicalVolume(plateStrip, plateMat, "PlateStripLeftLV");
    G4LogicalVolume* plateStripRightLV = new G4LogicalVolume(plateStrip, plateMat, "PlateStripRightLV");
    Place(plateStripLeftLV, plateStripLeftPos, "PlateStripLeftPVPL");
    Place(plateStripRightLV, plateStripRightPos, "PlateStripRightPVPL");

    // Plate center tube
    G4VSolid* plateCenter = new G4Tubs("PlateCenter", plateInnerHoleRadius, plateOuterHoleRadius, plateCenterThick / 2,
                                       0, viewDeg);
    const G4ThreeVector plateCenterPos = G4ThreeVector(0., 0., -(modelHeight + plateCenterThick) / 2.0);
    G4LogicalVolume* plateCenterLV = new G4LogicalVolume(plateCenter, plateMat, "PlateCenterLV");
    Place(plateCenterLV, plateCenterPos, "PlateCenterPVPL");

    // Plate center tube cap
    G4VSolid* plateCenterCap = new G4Tubs("PlateCenterCap", plateBottomHoleRadius, plateOuterHoleRadius, plateThick / 2,
                                          0, viewDeg);
    const G4ThreeVector plateCenterCapPos = G4ThreeVector(0., 0., -(modelHeight + plateThick) / 2.0 - plateCenterThick);
    G4LogicalVolume* plateCenterCapLV = new G4LogicalVolume(plateCenterCap, plateMat, "PlateCenterCapLV");
    Place(plateCenterCapLV, plateCenterCapPos, "PlateCenterCapPVPL");

    // Plate bottom cap
    G4VSolid* bottomPartWall = new G4Tubs("BottomPartWall", bottomCapInnerRadius, plateInnerHoleRadius,
                                          bottomCapHeight / 2, 0, viewDeg);
    const G4ThreeVector bottomPartWallPos = G4ThreeVector(
                                                          0., 0.,
                                                          -(modelHeight + bottomCapHeight) / 2.0 - plateCenterThick -
                                                          plateThick);
    G4LogicalVolume* bottomPartWallLV = new G4LogicalVolume(bottomPartWall, tunaCanMat, "BottomPartWallLV");
    Place(bottomPartWallLV, bottomPartWallPos, "BottomPartWallPVPL");

    G4VSolid* bottomPartCap = new G4Tubs("BottomPartCap", 0, bottomCapInnerRadius, bottomCapThick / 2, 0, viewDeg);
    const G4ThreeVector bottomPartCapPos = G4ThreeVector(
                                                         0., 0.,
                                                         -(modelHeight - bottomCapThick) / 2.0 - bottomCapHeight -
                                                         plateCenterThick - plateThick);
    G4LogicalVolume* bottomPartCapLV = new G4LogicalVolume(bottomPartCap, tunaCanMat, "BottomPartCapLV");
    Place(bottomPartCapLV, bottomPartCapPos, "BottomPartCapPVPL");

    tunaCanWallLV->SetVisAttributes(tunaCanVisAttr);
    tunaCanTopLV->SetVisAttributes(tunaCanVisAttr);

    plateLV->SetVisAttributes(tunaCanVisAttr);
    plateStripRightLV->SetVisAttributes(tunaCanVisAttr);
    plateStripLeftLV->SetVisAttributes(tunaCanVisAttr);

    plateCenterLV->SetVisAttributes(tunaCanVisAttr);
    plateCenterCapLV->SetVisAttributes(tunaCanVisAttr);

    bottomPartWallLV->SetVisAttributes(tunaCanVisAttr);
    bottomPartCapLV->SetVisAttributes(tunaCanVisAttr);
}
