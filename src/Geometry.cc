#include "Geometry.hh"

using namespace Sizes;


Geometry::Geometry(G4String detType, const G4bool useOpt, const G4bool lightCollect) : useOptics(useOpt),
    lightCollection(lightCollect), detectorType(std::move(detType)) {
    std::vector<G4String> detectorList = {"NaI", "CsI"};
    if (std::find(detectorList.begin(), detectorList.end(), detectorType) == detectorList.end()) {
        G4Exception("Geometry::ConstructDetector", "DetectorType", FatalException,
                    ("Detector not found: " + detectorType + ".\nAvailable detectors: NaI, CsI").c_str());
    }

    nist = G4NistManager::Instance();
    zeroRot = new G4RotationMatrix(0, 0, 0);

    viewDeg = 360 * deg;

    detContainerSize = G4ThreeVector(0 * mm,
                                     modelRadius - tunaCanThickWall,
                                     (modelHeight - bottomCapThick - tunaCanThickTop) / 2.0);
    detContainerPos = G4ThreeVector(0, 0, -(plateCenterThick + tunaCanThickTop) / 2);

    worldHalfSize = std::max({modelRadius * 2, modelHeight}) * 2;

    tunaCanVisAttr = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
    tunaCanVisAttr->SetForceSolid(true);
    detContVisAttr = new G4VisAttributes(G4VisAttributes::GetInvisible());
    detContVisAttr->SetForceSolid(false);
}

void Geometry::ConstructTunaCan() {
    if (tunaCanThickWall <= 0. and tunaCanThickTop <= 0.) return;

    G4Material* tunaCanMat = nist->FindOrBuildMaterial("G4_Al");
    G4Material* plateMat = nist->FindOrBuildMaterial("G4_Ti");

    // TunaCan
    G4VSolid* tunaCanWall = new G4Tubs("tunaCanWall", detContainerSize.y(), modelRadius, modelHeight / 2, 0, viewDeg);
    G4ThreeVector tunaCanWallPos = G4ThreeVector(0, 0, 0);
    G4LogicalVolume* tunaCanWallLV = new G4LogicalVolume(tunaCanWall, tunaCanMat, "tunaCanWallLV");
    new G4PVPlacement(zeroRot, tunaCanWallPos, tunaCanWallLV, "TunaCanWallPVPL", worldLV, false, 0, true);

    G4VSolid* tunaCanTop = new G4Tubs("TunaCanTop", 0, detContainerSize.y(), tunaCanThickTop / 2., 0, viewDeg);
    const G4ThreeVector tunaCanTopPos = G4ThreeVector(0, 0, (modelHeight - tunaCanThickTop) / 2.0);
    G4LogicalVolume* tunaCanTopLV = new G4LogicalVolume(tunaCanTop, tunaCanMat, "tunaCanTopLV");
    new G4PVPlacement(zeroRot, tunaCanTopPos, tunaCanTopLV, "TunaCanTopPVPL", worldLV, false, 0, true);

    // Plate square part
    G4VSolid* plateIncomplete = new G4Box("PlateIncomplete", plateSize / 2., plateSize / 2. + plateCornerSize,
                                          plateThick / 2.);
    G4VSolid* plateHole = new G4Tubs("PlateHole", 0., plateOuterHoleRadius, plateThick + 5 * mm, 0., 360);
    const G4ThreeVector platePos = G4ThreeVector(0, 0, -(modelHeight + plateThick) / 2.0);
    G4VSolid* plate = new G4SubtractionSolid("Plate", plateIncomplete, plateHole);
    G4LogicalVolume* plateLV = new G4LogicalVolume(plate, plateMat, "PlateLV");
    new G4PVPlacement(zeroRot, platePos, plateLV, "PlatePVPL", worldLV, false, 0, true);

    // Plate outer part
    G4VSolid* plateStrip = new G4Box("PlateStrip", plateCornerSize / 2., plateSize / 2., plateThick / 2.);
    const G4ThreeVector plateStripLeftPos = G4ThreeVector((plateSize + plateCornerSize) / 2., 0,
                                                          -(modelHeight + plateThick) / 2.0);
    const G4ThreeVector plateStripRightPos = G4ThreeVector(-(plateSize + plateCornerSize) / 2., 0,
                                                           -(modelHeight + plateThick) / 2.0);
    G4LogicalVolume* plateStripLeftLV = new G4LogicalVolume(plateStrip, plateMat, "PlateStripLeftLV");
    G4LogicalVolume* plateStripRightLV = new G4LogicalVolume(plateStrip, plateMat, "PlateStripRightLV");
    new G4PVPlacement(zeroRot, plateStripLeftPos, plateStripLeftLV, "PlateStripLeftPVPL", worldLV, false, 0, true);
    new G4PVPlacement(zeroRot, plateStripRightPos, plateStripRightLV, "PlateStripRightPVPL", worldLV, false, 0, true);

    // Plate center tube
    G4VSolid* plateCenter = new G4Tubs("PlateCenter", plateInnerHoleRadius, plateOuterHoleRadius, plateCenterThick / 2,
                                       0, viewDeg);
    const G4ThreeVector plateCenterPos = G4ThreeVector(0., 0., -(modelHeight + plateCenterThick) / 2.0);
    G4LogicalVolume* plateCenterLV = new G4LogicalVolume(plateCenter, plateMat, "PlateCenterLV");
    new G4PVPlacement(zeroRot, plateCenterPos, plateCenterLV, "PlateCenterPVPL", worldLV, false, 0, true);

    // Plate center tube cap
    G4VSolid* plateCenterCap = new G4Tubs("PlateCenterCap", plateBottomHoleRadius, plateOuterHoleRadius, plateThick / 2,
                                          0, viewDeg);
    const G4ThreeVector plateCenterCapPos = G4ThreeVector(0., 0., -(modelHeight + plateThick) / 2.0 - plateCenterThick);
    G4LogicalVolume* plateCenterCapLV = new G4LogicalVolume(plateCenterCap, plateMat, "PlateCenterCapLV");
    new G4PVPlacement(zeroRot, plateCenterCapPos, plateCenterCapLV, "PlateCenterCapPVPL", worldLV, false, 0, true);

    // Plate bottom cap
    G4VSolid* bottomPartWall = new G4Tubs("BottomPartWall", bottomCapInnerRadius, plateInnerHoleRadius,
                                          bottomCapHeight / 2, 0, viewDeg);
    const G4ThreeVector bottomPartWallPos = G4ThreeVector(
        0., 0., -(modelHeight + bottomCapHeight) / 2.0 - plateCenterThick - plateThick);
    G4LogicalVolume* bottomPartWallLV = new G4LogicalVolume(bottomPartWall, plateMat, "BottomPartWallLV");
    new G4PVPlacement(zeroRot, bottomPartWallPos, bottomPartWallLV, "BottomPartWallPVPL", worldLV, false, 0, true);

    G4VSolid* bottomPartCap = new G4Tubs("BottomPartCap", 0, bottomCapInnerRadius, bottomCapThick / 2, 0, viewDeg);
    const G4ThreeVector bottomPartCapPos = G4ThreeVector(
        0., 0., -(modelHeight - bottomCapThick) / 2.0 - bottomCapHeight - plateCenterThick - plateThick);
    G4LogicalVolume* bottomPartCapLV = new G4LogicalVolume(bottomPartCap, plateMat, "BottomPartCapLV");
    new G4PVPlacement(zeroRot, bottomPartCapPos, bottomPartCapLV, "BottomPartCapPVPL", worldLV, false, 0, true);

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


void Geometry::ConstructDetector() {
    G4VSolid* detContTopTube = new G4Tubs("DetContTopTube", 0, detContainerSize.y(),
                                          (modelHeight - tunaCanThickTop + plateCenterThick) / 2, 0,
                                          360 * deg);
    G4VSolid* detContMidTube = new G4Tubs("DetContMidTube", 0, plateBottomHoleRadius, plateThick / 2., 0,
                                          360 * deg);
    G4VSolid* detContBottomTube = new G4Tubs("DetContBottomTube", 0, bottomCapInnerRadius,
                                             (bottomCapHeight - bottomCapThick) / 2., 0,
                                             360 * deg);
    G4ThreeVector detContMidPos = G4ThreeVector(
        0, 0, -(modelHeight - tunaCanThickTop + plateCenterThick + plateThick) / 2);
    G4ThreeVector detContBottomPos = G4ThreeVector(
        0, 0, -(modelHeight - tunaCanThickTop + plateCenterThick + bottomCapHeight - bottomCapThick) / 2 - plateThick);
    G4VSolid* detContIncomplete = new G4UnionSolid("DetContIncomplete", detContTopTube, detContMidTube, zeroRot,
                                                   detContMidPos);

    detContainer = new G4UnionSolid("DetectorContainer", detContIncomplete, detContBottomTube, zeroRot,
                                    detContBottomPos);
    detContainerLV = new G4LogicalVolume(detContainer, worldMat, "DetectorContainerLV");
    detContainerPVPL = new G4PVPlacement(zeroRot, detContainerPos, detContainerLV, "DetectorContainerPVPL", worldLV,
                                         false, 0, true);
    detContainerLV->SetVisAttributes(detContVisAttr);

    detector = new Detector(detContainerLV,  nist, viewDeg, detectorType);
    detector->Construct();
    std::vector<G4LogicalVolume *> sensitiveLV = detector->GetSensitiveLV();
    crystalLV = sensitiveLV.at(0);
    vetoLV = sensitiveLV.at(1);
    bottomVetoLV = sensitiveLV.at(2);
    tyvekOutLV = sensitiveLV.at(3);
    tyvekMidLV = sensitiveLV.at(4);
    tyvekInLV = sensitiveLV.at(5);
    tyvekBottomLV = sensitiveLV.at(6);
}


G4VPhysicalVolume* Geometry::Construct() {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    worldMat = nist->FindOrBuildMaterial("G4_Galactic");

    auto* mptVac = new G4MaterialPropertiesTable();
    const G4int N = 2;
    G4double E[N] = {2.0 * eV, 3.4 * eV};
    G4double n1[N] = {1.0, 1.0};
    mptVac->AddProperty("RINDEX", E, n1, N);
    worldMat->SetMaterialPropertiesTable(mptVac);

    worldBox = new G4Box("World", worldHalfSize, worldHalfSize, worldHalfSize);
    worldLV = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldPVP = new G4PVPlacement(zeroRot, G4ThreeVector(0, 0, 0), worldLV, "WorldPVPL", nullptr, false, 0, false);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    ConstructDetector();
    ConstructTunaCan();

    return worldPVP;
}

void Geometry::SetStepLimits() {
    vetoStepLimit = new G4UserLimits(0.01 * mm, DBL_MAX, DBL_MAX, 0., 0.);
    vetoCuts = new G4ProductionCuts();
    vetoCuts->SetProductionCut(0.001 * mm, "e-");
    vetoCuts->SetProductionCut(0.001 * mm, "e+");
    vetoCuts->SetProductionCut(0.001 * mm, "gamma");

    vetoLV->SetUserLimits(vetoStepLimit);
    vetoRegion = new G4Region("VetoRegion");
    vetoLV->SetRegion(vetoRegion);
    vetoRegion->AddRootLogicalVolume(vetoLV);
    vetoRegion->SetProductionCuts(vetoCuts);
}


void Geometry::ConstructSDandField() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    auto *detectorSD = new SensitiveDetector("DetectorSD", 0, "Crystal");
    sdManager->AddNewDetector(detectorSD);
    crystalLV->SetSensitiveDetector(detectorSD);

    auto *vetoSD = new SensitiveDetector("VetoSD", 1, "Veto");
    sdManager->AddNewDetector(vetoSD);
    vetoLV->SetSensitiveDetector(vetoSD);

    auto *bottomVetoSD = new SensitiveDetector("BottomVetoSD", 2, "BottomVeto");
    sdManager->AddNewDetector(bottomVetoSD);
    bottomVetoLV->SetSensitiveDetector(bottomVetoSD);

    auto *tyvekOutSD = new SensitiveDetector("TyvekOutSD", 3, "TyvekOut");
    sdManager->AddNewDetector(tyvekOutSD);
    tyvekOutLV->SetSensitiveDetector(tyvekOutSD);

    auto *tyvekMidSD = new SensitiveDetector("TyvekMidSD", 4, "TyvekMid");
    sdManager->AddNewDetector(tyvekMidSD);
    tyvekMidLV->SetSensitiveDetector(tyvekMidSD);

    auto *tyvekInSD = new SensitiveDetector("TyvekInSD", 5, "TyvekIn");
    sdManager->AddNewDetector(tyvekInSD);
    tyvekInLV->SetSensitiveDetector(tyvekInSD);

    auto *tyvekBottomSD = new SensitiveDetector("TyvekBottomSD", 6, "TyvekBottom");
    sdManager->AddNewDetector(tyvekBottomSD);
    tyvekBottomLV->SetSensitiveDetector(tyvekBottomSD);

    if (useOptics) {
        auto* sipmSD = new SiPMOpticalSD("SiPMOpticalSD");
        sdManager->AddNewDetector(sipmSD);
        auto* sipmWindowLV = detector->GetSiPMWindowLV();
        sipmWindowLV->SetSensitiveDetector(sipmSD);
    }
}
