#include "Geometry.hh"

Geometry::Geometry(G4String detType, const Sizes &ss, G4double temp, G4bool dLED)
    : doubleLED(dLED), detectorType(std::move(detType)), sizes(ss), temperature(temp) {
    std::vector<G4String> detectorList = {"NaI"};
    if (std::find(detectorList.begin(), detectorList.end(), detectorType) == detectorList.end()) {
        G4Exception("Geometry::ConstructDetector", "DetectorType", FatalException,
                    ("Detector not found: " + detectorType + ".\nAvailable detectors: NaI").c_str());
    }

    nist = G4NistManager::Instance();
    zeroRot = new G4RotationMatrix(0, 0, 0);

    modelSize = G4ThreeVector(0 * mm, 30 * mm + sizes.lidThick, 22.5 * mm + sizes.lidThick / 2.);
    detContainerSize = G4ThreeVector(0 * mm, 30 * mm, 22.5 * mm);
    detContainerPos = G4ThreeVector(0, 0, 0);

    worldHalfSize = std::max({
                        modelSize.x(),
                        modelSize.y(),
                        modelSize.z()
                    }) * 2.;

    viewDeg = 360 * deg;
    sizes.lidThick = 2 * mm;

    lidVisAttr = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
    lidVisAttr->SetForceSolid(true);
}

void Geometry::ConstructLid() {
    if (sizes.lidThick <= 0.) return;

    G4Material *lidMat = nist->FindOrBuildMaterial("G4_Al");
    G4ThreeVector lidSize = G4ThreeVector(modelSize.y(), modelSize.y() + sizes.lidThick,
                                          modelSize.z() + sizes.lidThick / 2.);
    G4Tubs *lidTube = new G4Tubs("LidTube", lidSize.x(), lidSize.y(), lidSize.z(), 0, viewDeg);
    G4Tubs *lidCap = new G4Tubs("LidCap", 0, modelSize.y(), sizes.lidThick / 2., 0, viewDeg);
    G4ThreeVector lidCapPos = G4ThreeVector(0, 0, modelSize.z());
    lid = new G4UnionSolid("Lid", lidTube, lidCap, zeroRot, lidCapPos);
    lidLV = new G4LogicalVolume(lid, lidMat, "LidLV");
    G4ThreeVector lidPos = G4ThreeVector(detectorPos.x(), detectorPos.y(),
                                         detectorPos.z() + sizes.lidThick / 2);
    new G4PVPlacement(zeroRot, lidPos, lidLV, "LidPVPL", worldLV, false, 0, true);
    lidLV->SetVisAttributes(lidVisAttr);
}


void Geometry::ConstructDetector() {
    detContainer = new G4Tubs("DetectorContainer", detContainerSize.x(), detContainerSize.y(), detContainerSize.z(), 0,
                              viewDeg);
    detContainerLV = new G4LogicalVolume(detContainer, worldMat, "DetectorContainerLV");
    new G4PVPlacement(zeroRot, detContainerPos, detContainerLV, "DetectorContainerPVPL", worldLV, false, 0, true);
    detContainerLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    if (detectorType == "NaI") {
        detector = new NaI(detContainerLV, detContainerSize, nist, sizes, viewDeg, doubleLED);
    }
    detector->Construct();
    std::vector<G4LogicalVolume *> sensitiveLV = detector->GetSensitiveLV();
    detectorLV = sensitiveLV.at(0);
    vetoLV = sensitiveLV.at(1);
    tapeOutLV = sensitiveLV.at(2);
    tapeInLV = sensitiveLV.at(3);
    shellLV = sensitiveLV.at(4);
}


G4VPhysicalVolume *Geometry::Construct() {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    worldMat = nist->FindOrBuildMaterial("G4_Galactic");

    worldBox = new G4Box("World", worldHalfSize, worldHalfSize, worldHalfSize);
    worldLV = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldPVP = new G4PVPlacement(zeroRot, G4ThreeVector(0, 0, 0), worldLV, "WorldPVPL", 0, false, 0, false);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    ConstructDetector();
    ConstructLid();

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
    G4SDManager *sdManager = G4SDManager::GetSDMpointer();

    auto *detectorSD = new SensitiveDetector("DetectorSD", 0, "NaI");
    sdManager->AddNewDetector(detectorSD);
    detectorLV->SetSensitiveDetector(detectorSD);

    if (sizes.shellThick > 0) {
        auto *shellSD = new SensitiveDetector("ShellSD", 1, "Shell");
        sdManager->AddNewDetector(shellSD);
        shellLV->SetSensitiveDetector(shellSD);
    }

    if (sizes.lidThick > 0) {
        auto *lidSD = new SensitiveDetector("LidSD", 2, "Lid");
        sdManager->AddNewDetector(lidSD);
        lidLV->SetSensitiveDetector(lidSD);
    }

    if (sizes.vetoThick > 0) {
        auto *vetoSD = new SensitiveDetector("VetoSD", 3, "Veto");
        sdManager->AddNewDetector(vetoSD);
        vetoLV->SetSensitiveDetector(vetoSD);
    }

    if (sizes.tapeThick > 0.) {
        auto *tapeOutSD = new SensitiveDetector("TapeOutSD", 4, "TapeOut");
        sdManager->AddNewDetector(tapeOutSD);
        tapeOutLV->SetSensitiveDetector(tapeOutSD);

        auto *tapeInSD = new SensitiveDetector("TapeInSD", 5, "TapeIn");
        sdManager->AddNewDetector(tapeInSD);
        tapeInLV->SetSensitiveDetector(tapeInSD);
    }
}
