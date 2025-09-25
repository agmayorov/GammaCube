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

    viewDeg = 180 * deg;
    sizes.tunaCanThick = 2 * mm;

    modelSize = G4ThreeVector(0 * mm, 30 * mm + sizes.tunaCanThick, 22.5 * mm + sizes.tunaCanThick / 2.);
    detContainerSize = G4ThreeVector(0 * mm, 30 * mm, 22.5 * mm);
    detContainerPos = G4ThreeVector(0, 0, 0);

    worldHalfSize = std::max({
                        modelSize.x(),
                        modelSize.y(),
                        modelSize.z()
                    }) * 2.;

    tunaCanVisAttr = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
    tunaCanVisAttr->SetForceSolid(true);
}

void Geometry::ConstructTunaCan() {
    if (sizes.tunaCanThick <= 0.) return;

    G4Material *tunaCanMat = nist->FindOrBuildMaterial("G4_Al");
    G4ThreeVector tunaCanSize = G4ThreeVector(detContainerSize.y(), modelSize.y(), modelSize.z());
    G4Tubs *tunaCanTube = new G4Tubs("TunaCanTube", tunaCanSize.x(), tunaCanSize.y(), tunaCanSize.z(), 0, viewDeg);
    G4Tubs *tunaCanCap = new G4Tubs("TunaCanCap", 0, detContainerSize.y(), sizes.tunaCanThick / 2., 0, viewDeg);
    G4ThreeVector tunaCanCapPos = G4ThreeVector(0, 0, detContainerSize.z());
    tunaCan = new G4UnionSolid("TunaCan", tunaCanTube, tunaCanCap, zeroRot, tunaCanCapPos);
    tunaCanLV = new G4LogicalVolume(tunaCan, tunaCanMat, "TunaCanLV");
    G4ThreeVector tunaCanPos = G4ThreeVector(detectorPos.x(), detectorPos.y(),
                                             detectorPos.z() + sizes.gapSize / 2 + sizes.tunaCanThick / 2.);
    new G4PVPlacement(zeroRot, tunaCanPos, tunaCanLV, "TunaCanPVPL", worldLV, false, 0, true);
    tunaCanLV->SetVisAttributes(tunaCanVisAttr);
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
    tyvekOutLV = sensitiveLV.at(2);
    tyvekInLV = sensitiveLV.at(3);
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
    G4SDManager *sdManager = G4SDManager::GetSDMpointer();

    auto *detectorSD = new SensitiveDetector("DetectorSD", 0, "NaI");
    sdManager->AddNewDetector(detectorSD);
    detectorLV->SetSensitiveDetector(detectorSD);

    if (sizes.shellThick > 0) {
        auto *shellSD = new SensitiveDetector("ShellSD", 1, "Shell");
        sdManager->AddNewDetector(shellSD);
        shellLV->SetSensitiveDetector(shellSD);
    }

    if (sizes.tunaCanThick > 0) {
        auto *tunaCanSD = new SensitiveDetector("TunaCanSD", 2, "TunaCan");
        sdManager->AddNewDetector(tunaCanSD);
        tunaCanLV->SetSensitiveDetector(tunaCanSD);
    }

    if (sizes.vetoThick > 0) {
        auto *vetoSD = new SensitiveDetector("VetoSD", 3, "Veto");
        sdManager->AddNewDetector(vetoSD);
        vetoLV->SetSensitiveDetector(vetoSD);
    }

    if (sizes.tyvekThick > 0.) {
        auto *tyvekOutSD = new SensitiveDetector("TyvekOutSD", 4, "TyvekOut");
        sdManager->AddNewDetector(tyvekOutSD);
        tyvekOutLV->SetSensitiveDetector(tyvekOutSD);

        auto *tyvekInSD = new SensitiveDetector("TyvekInSD", 5, "TyvekIn");
        sdManager->AddNewDetector(tyvekInSD);
        tyvekInLV->SetSensitiveDetector(tyvekInSD);
    }
}
