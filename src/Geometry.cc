#include "Geometry.hh"

using namespace Sizes;


void Geometry::CheckSizes() {
    // G4double sumRadius = tunaCanThickWall + gapSizeWall + tyvekOutThickWall + vetoThickWall +
    //                      tyvekMidThickWall + AlThickWall + tyvekInThickWall;
    // G4double sumHeight = tunaCanThickBottom + tunaCanThickTop + gapSizeBottom + gapSizeTop +
    //                      tyvekOutThickBottom + tyvekOutThickTop + vetoThickBottom + vetoThickTop +
    //                      tyvekMidThickBottom + tyvekMidThickTop + rubberHeight +
    //                      AlCapThickBottom + AlThickTop + tyvekInThickBottom + tyvekInThickTop;
}


Geometry::Geometry(G4String detType)
    : detectorType(std::move(detType)) {
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
                                     (modelHeight - tunaCanThickTop - tunaCanThickBottom) / 2.0);
    detContainerPos = G4ThreeVector(0, 0, 0);

    worldHalfSize = std::max({modelRadius * 2, modelHeight}) * 2;

    tunaCanVisAttr = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
    tunaCanVisAttr->SetForceSolid(true);
}

void Geometry::ConstructTunaCan() {
    if (tunaCanThickWall <= 0. and tunaCanThickTop <= 0. and tunaCanThickBottom <= 0.) return;

    G4Material *tunaCanMat = nist->FindOrBuildMaterial("G4_Al");

    G4VSolid *tunaCanIncomplete = nullptr;
    G4VSolid *tunaCanWall = nullptr;
    if (tunaCanThickWall > 0) {
        tunaCanWall = new G4Tubs("tunaCanWall", detContainerSize.y(), modelRadius, modelHeight / 2, 0, viewDeg);
    }

    G4VSolid *tunaCanTop = nullptr;
    const G4ThreeVector tunaCanTopPos = G4ThreeVector(0, 0, (modelHeight - tunaCanThickTop) / 2.0);
    if (tunaCanThickTop > 0) {
        tunaCanTop = new G4Tubs("TunaCanTop", 0, detContainerSize.y(), tunaCanThickTop / 2., 0, viewDeg);
    }
    if (tunaCanThickWall > 0 and tunaCanThickTop > 0) {
        tunaCanIncomplete = new G4UnionSolid("TunaCanIncomplete", tunaCanWall, tunaCanTop, zeroRot,
                                             tunaCanTopPos);
    } else if (tunaCanThickWall > 0 and tunaCanThickTop <= 0) {
        tunaCanIncomplete = tunaCanWall;
    } else if (tunaCanThickWall <= 0 and tunaCanThickTop > 0) {
        tunaCanIncomplete = tunaCanTop;
    }

    G4Tubs *tunaCanBottom = nullptr;
    const G4ThreeVector tunaCanBottomPos = G4ThreeVector(0, 0, -(modelHeight - tunaCanThickBottom) / 2.0);
    if (tunaCanThickBottom > 0) {
        tunaCanBottom = new G4Tubs("TunaCanBottom", 0, detContainerSize.y(), tunaCanThickBottom / 2., 0,
                                   viewDeg);
    }
    if (tunaCanIncomplete != nullptr and tunaCanThickBottom > 0) {
        tunaCan = new G4UnionSolid("TunaCan", tunaCanIncomplete, tunaCanBottom, zeroRot, tunaCanBottomPos);
    } else if (tunaCanIncomplete == nullptr and tunaCanThickBottom > 0) {
        tunaCan = tunaCanBottom;
    } else if (tunaCanIncomplete != nullptr and tunaCanThickBottom <= 0) {
        tunaCan = tunaCanIncomplete;
    }

    tunaCanLV = new G4LogicalVolume(tunaCan, tunaCanMat, "TunaCanLV");
    G4ThreeVector tunaCanPos = G4ThreeVector(0, 0, modelHeight / 2. - detContainerSize.z() - tunaCanThickBottom);
    new G4PVPlacement(zeroRot, tunaCanPos, tunaCanLV, "TunaCanPVPL", worldLV, false, 0, true);
    tunaCanLV->SetVisAttributes(tunaCanVisAttr);
}


void Geometry::ConstructDetector() {
    detContainer = new G4Tubs("DetectorContainer", detContainerSize.x(), detContainerSize.y(), detContainerSize.z(), 0,
                              360 * deg);
    detContainerLV = new G4LogicalVolume(detContainer, worldMat, "DetectorContainerLV");
    new G4PVPlacement(zeroRot, detContainerPos, detContainerLV, "DetectorContainerPVPL", worldLV, false, 0, true);
    detContainerLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    detector = new Detector(detContainerLV, detContainerSize, nist, viewDeg, detectorType);
    detector->Construct();
    std::vector<G4LogicalVolume *> sensitiveLV = detector->GetSensitiveLV();
    crystalLV = sensitiveLV.at(0);
    tyvekInLV = sensitiveLV.at(1);
    AlLV = sensitiveLV.at(2);
    rubberLV = sensitiveLV.at(3);
    tyvekMidLV = sensitiveLV.at(4);
    vetoLV = sensitiveLV.at(5);
    tyvekOutLV = sensitiveLV.at(6);
    crystalLEDLV = sensitiveLV.at(7);
    vetoLEDLV = sensitiveLV.at(8);
    vetoBottomLEDLV = sensitiveLV.at(9);
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

    auto *detectorSD = new SensitiveDetector("DetectorSD", 0, "Crystal");
    sdManager->AddNewDetector(detectorSD);
    crystalLV->SetSensitiveDetector(detectorSD);

    if (tunaCanMinSize > 0) {
        auto *tunaCanSD = new SensitiveDetector("TunaCanSD", 1, "TunaCan");
        sdManager->AddNewDetector(tunaCanSD);
        tunaCanLV->SetSensitiveDetector(tunaCanSD);
    }

    if (vetoMinSize > 0) {
        auto *vetoSD = new SensitiveDetector("VetoSD", 2, "Veto");
        sdManager->AddNewDetector(vetoSD);
        vetoLV->SetSensitiveDetector(vetoSD);
    }

    if (tyvekOutMinSize > 0.) {
        auto *tyvekOutSD = new SensitiveDetector("TyvekOutSD", 3, "TyvekOut");
        sdManager->AddNewDetector(tyvekOutSD);
        tyvekOutLV->SetSensitiveDetector(tyvekOutSD);
    }

    if (tyvekInMinSize > 0.) {
        auto *tyvekInSD = new SensitiveDetector("TyvekInSD", 4, "TyvekIn");
        sdManager->AddNewDetector(tyvekInSD);
        tyvekInLV->SetSensitiveDetector(tyvekInSD);
    }

    if (tyvekMidMinSize > 0.) {
        auto *tyvekMidSD = new SensitiveDetector("TyvekMidSD", 5, "TyvekMid");
        sdManager->AddNewDetector(tyvekMidSD);
        tyvekMidLV->SetSensitiveDetector(tyvekMidSD);
    }

    if (rubberMinSize > 0.) {
        auto *rubberSD = new SensitiveDetector("RubberSD", 6, "Rubber");
        sdManager->AddNewDetector(rubberSD);
        rubberLV->SetSensitiveDetector(rubberSD);
    }

    if (AlMinSize > 0.) {
        auto *AlSD = new SensitiveDetector("AluminiumSD", 7, "Aluminium");
        sdManager->AddNewDetector(AlSD);
        AlLV->SetSensitiveDetector(AlSD);
    }

    if (crystalLEDMinSize > 0.) {
        auto *crystalLEDSD = new SensitiveDetector("CrystalLEDSD", 8, "CrystalLED");
        sdManager->AddNewDetector(crystalLEDSD);
        crystalLEDLV->SetSensitiveDetector(crystalLEDSD);
    }

    if (crystalLEDMinSize > 0.) {
        auto *vetoLEDSD = new SensitiveDetector("VetoLEDSD", 9, "VetoLED");
        sdManager->AddNewDetector(vetoLEDSD);
        vetoLEDLV->SetSensitiveDetector(vetoLEDSD);
    }

    if (crystalLEDMinSize > 0.) {
        auto *vetoBottomLEDSD = new SensitiveDetector("VetoBottomLEDSD", 10, "VetoBottomLED");
        sdManager->AddNewDetector(vetoBottomLEDSD);
        vetoBottomLEDLV->SetSensitiveDetector(vetoBottomLEDSD);
    }
}
