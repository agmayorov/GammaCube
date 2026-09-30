#include "geometry/Engine.hh"

using namespace EngineSizes;

Engine::Engine(G4LogicalVolume* world, G4NistManager* nistManager) : worldLV(world), nist(nistManager)
{
    DefineMaterial();
    DefineVisual();

    auto* engineContainer = new G4Box("EngineContainer", engineContainer::halfX, engineContainer::halfY,
                                      engineContainer::halfZ);
    engineContainerLV = new G4LogicalVolume(engineContainer, vacuumMat, "EngineContainerLV");
    engineContainerLV->SetVisAttributes(G4VisAttributes::GetInvisible());
}

void Engine::DefineVisual()
{
    visBoardTop = new G4VisAttributes(G4Color(0.0, 1.0, 0.0));
    visBoardTop->SetForceSolid(true);
    visBoardBottom = new G4VisAttributes(G4Color(1.0, 0.0, 0.0));
    visBoardBottom->SetForceSolid(true);
    visMagnetic = new G4VisAttributes(G4Color(1.0, 0.5, 0.0));
    visMagnetic->SetForceSolid(true);
    visCapacitor = new G4VisAttributes(G4Color(0.4, 0.4, 0.4));
    visCapacitor->SetForceSolid(true);
    visAl = new G4VisAttributes(G4Color(1, 1, 1));
    visAl->SetForceSolid(true);
    visElectronics = new G4VisAttributes(G4Color(0, 0, 0));
    visElectronics->SetForceSolid(true);
}

void Engine::DefineMaterial()
{
    auto* elH = nist->FindOrBuildElement("H");
    auto* elC = nist->FindOrBuildElement("C");

    auto* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    vacuumMat = nist->FindOrBuildMaterial("G4_Galactic");

    {
        auto* Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, 2);

        Epoxy->AddElement(elH, 2);
        Epoxy->AddElement(elC, 2);

        boardMat = new G4Material("PCB_FR4", 1.86 * g / cm3, 2);

        boardMat->AddMaterial(Epoxy, 0.472);
        boardMat->AddMaterial(SiO2, 0.528);
    }
    electronicsMat = nist->FindOrBuildMaterial("G4_Al");
    magneticMat = nist->FindOrBuildMaterial("G4_Cu");
    AlMat = nist->FindOrBuildMaterial("G4_Al");
    capacitorMat = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
}

void Engine::ConstructBoardBottom()
{
    using namespace EngineSizes::BoardBottom;

    auto* boardSolid = new G4Box("BoardBottom", halfX, halfY, halfZ);
    auto* cornerCutoutSolid = new G4Box("BoardBottomCornerCutout", cornerCutoutHalfY*2, cornerCutoutHalfX*2, cutoutHalfZ*2);

    auto* cornerCutout1 = new G4SubtractionSolid("BoardBottomCornerCutout1", boardSolid, cornerCutoutSolid, nullptr,
                                                 G4ThreeVector(
                                                     -halfX,
                                                     -halfY,
                                                     0)
    );

    auto* cornerCutout2 = new G4SubtractionSolid("BoardBottomCornerCutout2", cornerCutout1, cornerCutoutSolid, nullptr,
                                                 G4ThreeVector(
                                                     halfX,
                                                     -halfY,
                                                     0)
    );

    auto* cornerCutout3 = new G4SubtractionSolid("BoardBottomCornerCutout3", cornerCutout2, cornerCutoutSolid, nullptr,
                                                 G4ThreeVector(
                                                     -halfX,
                                                     halfY,
                                                     0)
    );

    auto* boardWithCorners = new G4SubtractionSolid("BoardBottomWithCornerCutouts", cornerCutout3, cornerCutoutSolid,
                                                    nullptr,
                                                    G4ThreeVector(
                                                        halfX,
                                                        halfY,
                                                        0)
    );

    auto* edgeCutoutX = new G4Box("BoardBottomEdgeCutoutX", edgeCutoutHalfY, edgeCutoutHalfX, cutoutHalfZ*2);
    auto* edgeCutoutY = new G4Box("BoardBottomEdgeCutoutY", edgeCutoutHalfX, edgeCutoutHalfY, cutoutHalfZ*2);

    auto* edgeCutout1 = new G4SubtractionSolid("BoardBottomEdgeCutout1", boardWithCorners, edgeCutoutX, nullptr,
                                               G4ThreeVector(
                                                   -halfX + edgeCutoutHalfY,
                                                   0,
                                                   0)
    );

    auto* edgeCutout2 = new G4SubtractionSolid("BoardBottomEdgeCutout2", edgeCutout1, edgeCutoutX, nullptr,
                                               G4ThreeVector(
                                                   halfX - edgeCutoutHalfY,
                                                   0,
                                                   0)
    );

    auto* edgeCutout3 = new G4SubtractionSolid("BoardBottomEdgeCutout3", edgeCutout2, edgeCutoutY, nullptr,
                                               G4ThreeVector(
                                                   0,
                                                   -halfY + edgeCutoutHalfY,
                                                   0)
    );

    auto* boardFinalSolid = new G4SubtractionSolid("BoardBottomFinal", edgeCutout3, edgeCutoutY, nullptr,
                                                   G4ThreeVector(
                                                       0,
                                                       halfY - edgeCutoutHalfY,
                                                       0)
    );

    boardBottomLV = new G4LogicalVolume(boardFinalSolid, boardMat, "BoardBottomLV");
    boardBottomLV->SetVisAttributes(visBoardBottom);
}

void Engine::ConstructBoardTop()
{
    using namespace EngineSizes::BoardTop;
    auto* boardTopSolid = new G4Box("BoardTop", halfX, halfY, halfZ);

    boardTopLV = new G4LogicalVolume(boardTopSolid, boardMat, "BoardTopLV");
    boardTopLV->SetVisAttributes(visBoardTop);
}

void Engine::ConstructCapacitor()
{
    using namespace EngineSizes::Capacitor;

    auto* capacitorSolid = new G4Box("Capacitor", halfX, halfY, halfZ);

    capacitorLV = new G4LogicalVolume(capacitorSolid, capacitorMat, "CapacitorLV");
    capacitorLV->SetVisAttributes(visCapacitor);
}

void Engine::ConstructElectronics()
{
    using namespace EngineSizes::Electronics;

    auto* electronicsSolid = new G4Box("Electronics", halfX, halfY, halfZ);

    electronicsLV = new G4LogicalVolume(electronicsSolid, electronicsMat, "ElectronicsLV");
    electronicsLV->SetVisAttributes(visElectronics);
}
void Engine::ConstructFieldCoilCentral()
{
    using namespace EngineSizes::FieldCoilCentral;

    auto* fieldCoilCentralSolid = new G4Tubs("FieldCoilCentral", 3.0*mm, radius-1*mm, halfHeight, 0.0, 360.0 * deg);
    auto* fieldCoilCentralSolidCu = new G4Tubs("FieldCoilCentralCu", radius-1*mm, radius, halfHeight, 0.0, 360.0 * deg);

    fieldCoilCentralLV = new G4LogicalVolume(fieldCoilCentralSolid, AlMat, "FieldCoilCentralLV");
    fieldCoilCentralLV->SetVisAttributes(visAl);
    fieldCoilCentralCuLV = new G4LogicalVolume(fieldCoilCentralSolidCu, magneticMat, "FieldCoilCentralCuLV");
    fieldCoilCentralCuLV->SetVisAttributes(visMagnetic);
}

void Engine::ConstructFieldCoil()
{
    using namespace EngineSizes::FieldCoil;

    auto* fieldCoilSolid = new G4Tubs("FieldCoil", innerRadius, outerRadius, halfHeight, 0.0, 360.0 * deg);
    auto* fieldCoilAir = new G4Tubs("FieldCoilAir", innerRadius+wire, outerRadius-wire, halfHeight-wire, 0.0, 360.0 * deg);

    fieldCoilLV = new G4LogicalVolume(fieldCoilSolid, magneticMat, "FieldCoilLV");
    fieldCoilLV->SetVisAttributes(visMagnetic);
    auto* fieldCoilAirLV = new G4LogicalVolume(fieldCoilAir, vacuumMat, "FieldCoilAirLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, 0.0), fieldCoilAirLV, "FieldCoilAirPVP",fieldCoilLV, false,0, true);
}

void Engine::ConstructEngine()
{
    ConstructBoardTop();
    ConstructBoardBottom();
    ConstructCapacitor();
    ConstructFieldCoilCentral();
    ConstructFieldCoil();
    ConstructElectronics();
    auto* boardBottom1PV = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, -engineContainer::halfZ + BoardBottom::halfZ),
                                             boardBottomLV,
                                             "BoardBottom1PV",
                                             engineContainerLV,
                                             false,
                                             0,
                                             true
    );

    const G4double boardBottom2Z = EngineSizes::BoardBottom::halfZ + 4.0 * mm + EngineSizes::BoardBottom::halfZ -engineContainer::halfZ + BoardBottom::halfZ;

    auto* boardBottom2PV = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, boardBottom2Z),
                                             boardBottomLV,
                                             "BoardBottom2PV",
                                             engineContainerLV,
                                             false,
                                             1,
                                             true
    );

    const G4double fieldCoilCentralZ = boardBottom2Z + BoardBottom::halfZ + FieldCoilCentral::halfHeight;

    auto* fieldCoilCentralPV = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, fieldCoilCentralZ),
                                                 fieldCoilCentralLV,
                                                 "FieldCoilCentralPV",
                                                 engineContainerLV,
                                                 false,
                                                 0,
                                                 true
                                                 );

    auto* fieldCoilCentralCuPV = new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, fieldCoilCentralZ),
                                                 fieldCoilCentralCuLV,
                                                 "FieldCoilCentralCuPV",
                                                 engineContainerLV,
                                                 false,
                                                 0,
                                                 true
    );

    auto* capacitor1PV = new G4PVPlacement( nullptr,G4ThreeVector(-4.5 * mm, -29.0 * mm, 29.1 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        capacitorLV,
        "Capacitor1PV",
        engineContainerLV,
        false,
        0,
        true
    );
    auto* capacitor2Rotation = new G4RotationMatrix();
    capacitor2Rotation->rotateZ(90.0 * deg);

    auto* capacitor2PV = new G4PVPlacement( capacitor2Rotation, G4ThreeVector(29.0 * mm, -18.5 * mm, 29.1 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        capacitorLV,
        "Capacitor2PV",
        engineContainerLV,
        false,
        1,
        true
    );

    auto* capacitor3Rotation = new G4RotationMatrix();
    capacitor3Rotation->rotateZ(90.0 * deg);

    auto* capacitor3PV = new G4PVPlacement(capacitor3Rotation, G4ThreeVector(-29.0 * mm, 18.5 * mm, 29.1 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        capacitorLV,
        "Capacitor3PV",
        engineContainerLV,
        false,
        2,
        true
    );

    auto* capacitor4PV = new G4PVPlacement( nullptr,G4ThreeVector(4.5 * mm, 29.0 * mm, 29.1 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        capacitorLV,
        "Capacitor1PV",
        engineContainerLV,
        false,
        0,
        true
    );

    auto* boardTop1PV = new G4PVPlacement(nullptr,  G4ThreeVector(-37.135 * mm, -21.0 * mm, 29.0 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        boardTopLV,
        "BoardTopPV",
        engineContainerLV,
        false,
        0,
        true
    );

    auto* boardTop2PV = new G4PVPlacement(nullptr,  G4ThreeVector(37.135 * mm, 21.0 * mm, 29.0 * mm -engineContainer::halfZ + BoardBottom::halfZ),
       boardTopLV,
       "BoardTopPV",
       engineContainerLV,
       false,
       0,
       true
   );
    auto* fieldCoilRotation = new G4RotationMatrix();
    fieldCoilRotation->rotateY(90.0 * deg);

    auto* fieldCoilPV = new G4PVPlacement(fieldCoilRotation,G4ThreeVector(-32.635 * mm, -21.0 * mm, 24.25 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        fieldCoilLV,
        "FieldCoilPV",
        engineContainerLV,
        false,
        0,
        true
    );

    auto* electronicsPV = new G4PVPlacement(nullptr,G4ThreeVector(34.485 * mm, 21.0 * mm, 24.25 * mm -engineContainer::halfZ + BoardBottom::halfZ),
        electronicsLV,
        "ElectronicsPV",
        engineContainerLV,
        false,
        0,
        true
    );

}

