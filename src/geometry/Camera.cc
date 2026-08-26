#include "geometry/Camera.hh"

using namespace CameraSizes;

Camera::Camera(G4LogicalVolume* world, G4NistManager* nistManager) : worldLV(world), nist(nistManager) {
    DefineMaterial();
    DefineVisual();
}


void Camera::DefineVisual() {
    visLense = new G4VisAttributes(G4Color(G4Color(0.0, 0.86, 1.0)));
    visLense->SetForceSolid(true);
    visLenseShell = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visLenseShell->SetForceSolid(true);
    visRetainer = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visRetainer->SetForceSolid(true);
    visLenseHolder = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visLenseHolder->SetForceSolid(true);
    visElectronics = new G4VisAttributes(G4Color(G4Color(1.0, 0.5, 0.0)));
    visElectronics->SetForceSolid(true);
    visMatrix = new G4VisAttributes(G4Color(G4Color(1.0, 0.0, 1.0)));
    visMatrix->SetForceSolid(true);
    visElectronicsHolder = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visElectronicsHolder->SetForceSolid(true);
    visCameraPlate = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visCameraPlate->SetForceSolid(true);
}

void Camera::DefineMaterial() {
    auto* elH = nist->FindOrBuildElement("H");
    auto* elC = nist->FindOrBuildElement("C");

    G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    vacuumMat = nist->FindOrBuildMaterial("G4_Galactic");
    lenseMat = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
    lenseShellMat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    retainerMat = nist->FindOrBuildMaterial("G4_Al");
    lenseHolderMat = nist->FindOrBuildMaterial("G4_Al");
    electronicsHolderMat = nist->FindOrBuildMaterial("G4_Al");
    {
        auto* Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, 2);
        Epoxy->AddElement(elH, 2);
        Epoxy->AddElement(elC, 2);

        electronicsMat = new G4Material("FiberglassLaminate", 1.86 * g / cm3, 2);
        electronicsMat->AddMaterial(Epoxy, 0.472);
        electronicsMat->AddMaterial(SiO2, 0.528);
    }
    matrixMat = electronicsMat;
    cameraPlateMat = nist->FindOrBuildMaterial("G4_Al");
}


void Camera::ConstructLense() {
    lenseTop = new G4Sphere("LenseTop", Lenses::curvitureRadiusTop - Lenses::thicknessTop, Lenses::curvitureRadiusTop,
                            0, 360 * deg, 0, asin(Lenses::radiusTop / Lenses::curvitureRadiusTop));
    lenseTopLV = new G4LogicalVolume(lenseTop, lenseMat, "LenseTopLV");
    lenseTopLV->SetVisAttributes(visLense);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, CameraBox::halfZ - Lenses::offsetZTop - Lenses::curvitureRadiusTop -
                                    Lenses::lenseOffsetZ), lenseTopLV, "LenseTopPVP", cameraLV, false, 0, true);

    lenseBottom = new G4Tubs("LenseBottom", 0, Lenses::radiusBottom, Lenses::thicknessBottom, 0, 360 * deg);
    lenseBottomLV = new G4LogicalVolume(lenseBottom, lenseMat, "LenseBottomLV");
    lenseBottomLV->SetVisAttributes(visLense);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, CameraBox::halfZ - 2 * (LenseShellTop::height + LenseShellBody::wholeHeight +
                                        LenseShellBottom::height) + Lenses::offsetZBottom + Lenses::thicknessBottom +
                                    LenseShellTop::delta - Lenses::lenseOffsetZ), lenseBottomLV, "LenseBottomPVP",
                      cameraLV, false, 0, true);

    // Shell Top
    G4VSolid* topPlasticBase = new G4Tubs("TopPlasticBase", Lenses::radiusTop, LenseShellTop::radius,
                                          LenseShellTop::height, 0, 360 * deg);

    G4VSolid* sub1 = new G4Cons("TopPlasticSub1", 0, Lenses::radiusTop, 0,
                                Lenses::radiusTop + (LenseShellTop::hole2Radius - Lenses::radiusTop) * (1 +
                                    LenseShellTop::hole3Deep / LenseShellTop::hole2Deep),
                                LenseShellTop::hole2Deep + LenseShellTop::hole3Deep, 0, 360 * deg);
    G4VSolid* base1 = new G4SubtractionSolid("TopPlasticBase1", topPlasticBase, sub1, nullptr,
                                             G4ThreeVector(0, 0, -LenseShellTop::height +
                                                           LenseShellTop::hole1Deep * 2 + LenseShellTop::hole2Deep +
                                                           LenseShellTop::hole3Deep));

    G4VSolid* sub2 = new G4Tubs("TopPlasticSub2", 0, LenseShellTop::hole3Radius, LenseShellTop::height, 0,
                                360 * deg);
    G4VSolid* base2 = new G4SubtractionSolid("TopPlasticBase2", base1, sub2, nullptr,
                                             G4ThreeVector(0, 0, 2 * (LenseShellTop::hole1Deep +
                                                               LenseShellTop::hole2Deep)));

    G4VSolid* sub3 = new G4Cons("TopPlasticSub3", 0, LenseShellTop::hole3Radius, 0,
                                LenseShellTop::hole3Radius + (LenseShellTop::hole4Radius -
                                    LenseShellTop::hole3Radius) *
                                (1 + LenseShellTop::hole4Deep / 5 * LenseShellTop::hole3Deep),
                                LenseShellTop::hole3Deep + LenseShellTop::hole4Deep / 5, 0, 360 * deg);
    G4VSolid* base3 = new G4SubtractionSolid("TopPlasticBase3", base2, sub3, nullptr,
                                             G4ThreeVector(0, 0, -LenseShellTop::height + 2 * (
                                                               LenseShellTop::hole1Deep
                                                               +
                                                               LenseShellTop::hole2Deep) +
                                                           LenseShellTop::hole3Deep));

    G4VSolid* sub4 = new G4Tubs("TopPlasticSub4", 0, LenseShellTop::hole5Radius, LenseShellTop::height, 0,
                                360 * deg);
    G4VSolid* base4 = new G4SubtractionSolid("TopPlasticBase3", base3, sub4, nullptr,
                                             G4ThreeVector(0, 0, 2 * (LenseShellTop::hole1Deep +
                                                               LenseShellTop::hole2Deep +
                                                               LenseShellTop::hole3Deep)));

    G4VSolid* sub5 = new G4Cons("TopPlasticSub5", 0, LenseShellTop::hole5Radius, 0,
                                LenseShellTop::hole5Radius + (LenseShellTop::hole6Radius -
                                    LenseShellTop::hole5Radius) *
                                (1 + LenseShellTop::hole4Deep / LenseShellTop::hole6Deep),
                                LenseShellTop::hole6Deep + LenseShellTop::hole4Deep, 0, 360 * deg);
    G4VSolid* base5 = new G4SubtractionSolid("TopPlasticBase4", base4, sub5, nullptr,
                                             G4ThreeVector(0, 0, LenseShellTop::height -
                                                           LenseShellTop::hole6Deep));

    G4VSolid* sub6 = new G4Tubs("TopPlasticSub6", LenseShellTop::notchRadius, LenseShellTop::radius * 2,
                                LenseShellTop::height, 0, 360 * deg);
    G4VSolid* shellTop = new G4SubtractionSolid("LenseShellTop", base5, sub6, nullptr,
                                                G4ThreeVector(0, 0, -2 * (LenseShellTop::height -
                                                                  LenseShellTop::notchHeight) + LenseShellTop::delta));

    // Shell Body
    G4VSolid* bodyBase0Wall = new G4Tubs("BodyBase0Wall", LenseShellBody::radius1 - LenseShellBody::thickness * 2,
                                         LenseShellBody::radius1, LenseShellBody::height1, 0, 360 * deg);
    G4VSolid* bodyBase0Top = new G4Tubs("BodyBase0Top", Lenses::radiusTop,
                                        LenseShellBody::radius1 - LenseShellBody::thickness, LenseShellBody::thickness,
                                        0, 360 * deg);
    G4VSolid* bodyBase0Bottom = new G4Tubs("BodyBase0Bottom", LenseShellBody::radius2,
                                           LenseShellBody::radius1 - LenseShellBody::thickness,
                                           LenseShellBody::thickness, 0, 360 * deg);
    G4VSolid* bodyBase0Inc = new G4UnionSolid("BodyBase0Inc", bodyBase0Wall, bodyBase0Top, nullptr,
                                              G4ThreeVector(0, 0, LenseShellBody::height1 - LenseShellBody::thickness));
    G4VSolid* bodyBase0 = new G4UnionSolid("BodyBase0", bodyBase0Inc, bodyBase0Bottom, nullptr,
                                           G4ThreeVector(0, 0, LenseShellBody::thickness - LenseShellBody::height1));

    G4VSolid* bodyBase1Inc = new G4Tubs("BodyBase1Inc", LenseShellBody::radius3,
                                        LenseShellBody::radius2,
                                        LenseShellBody::height2 + LenseShellBody::thickness, 0, 360 * deg);
    G4VSolid* bodyBase1Sub = new G4Tubs("BodyBase1Sub", 0, LenseShellBody::radius2 - LenseShellBody::thickness * 2,
                                        LenseShellBody::height2, 0, 360 * deg);
    G4VSolid* bodyBase1 = new G4SubtractionSolid("BodyBase1", bodyBase1Inc, bodyBase1Sub, nullptr,
                                                 G4ThreeVector(0, 0, LenseShellBody::thickness));

    G4VSolid* bodyBaseInc1 = new G4UnionSolid("BodyBaseInc1", bodyBase0, bodyBase1, nullptr,
                                              G4ThreeVector(0, 0, -LenseShellBody::height1 - LenseShellBody::height2 +
                                                            LenseShellBody::thickness));

    G4VSolid* bodyBase2 = new G4Tubs("BodyBase2", LenseShellBody::radius3 - LenseShellBody::thickness * 2,
                                     LenseShellBody::radius3, LenseShellBody::height3 + LenseShellBody::thickness,
                                     0, 360 * deg);


    G4VSolid* bodyBaseInc2 = new G4UnionSolid("BodyBaseInc2", bodyBaseInc1, bodyBase2, nullptr,
                                              G4ThreeVector(0, 0, -LenseShellBody::height1 - LenseShellBody::height2 * 2
                                                            - LenseShellBody::height3 + LenseShellBody::thickness));

    G4VSolid* bodyBase3Inc = new G4Tubs("BodyBase3Inc", LenseShellBody::radius3 - LenseShellBody::thickness * 2,
                                        LenseShellBody::radius4, LenseShellBody::height4, 0, 360 * deg);
    G4VSolid* bodyBase3Sub = new G4Tubs("BodyBase3Sub", 0, LenseShellBody::radius4 - LenseShellBody::thickness * 2,
                                        LenseShellBody::height4, 0, 360 * deg);
    G4VSolid* bodyBase3 = new G4SubtractionSolid("BodyBase3", bodyBase3Inc, bodyBase3Sub, nullptr,
                                                 G4ThreeVector(0, 0, -LenseShellBody::thickness));

    G4double tempOffsetZ = -LenseShellBody::height1 - 2 * (LenseShellBody::height2 + LenseShellBody::height3) -
        LenseShellBody::height4;

    G4VSolid* bodyBaseInc3 = new G4UnionSolid("BodyBaseInc3", bodyBaseInc2, bodyBase3, nullptr,
                                              G4ThreeVector(0, 0, tempOffsetZ));

    G4VSolid* bodyBase4Inc = new G4Tubs("BodyBase4Inc", LenseShellBody::radius4 - LenseShellBody::thickness * 2,
                                        LenseShellBody::radius5, LenseShellBody::height5, 0, 360 * deg);
    G4VSolid* bodyBase4Sub = new G4Tubs("BodyBase4Sub", 0, LenseShellBody::radius5 - LenseShellBody::thickness * 2,
                                        LenseShellBody::height5, 0, 360 * deg);
    G4VSolid* bodyBase4 = new G4SubtractionSolid("BodyBase4", bodyBase4Inc, bodyBase4Sub, nullptr,
                                                 G4ThreeVector(0, 0, -LenseShellBody::thickness));


    G4VSolid* bodyBaseInc4 = new G4UnionSolid("BodyBaseInc4", bodyBaseInc3, bodyBase4, nullptr,
                                              G4ThreeVector(0, 0, tempOffsetZ - LenseShellBody::height5 -
                                                            LenseShellBody::height4));

    G4VSolid* bodyBase5Wall = new G4Tubs("BodyBase5Wall", LenseShellBody::radius6 - LenseShellBody::thickness * 2,
                                         LenseShellBody::radius6, LenseShellBody::height6, 0, 360 * deg);
    G4VSolid* bodyBase5Top = new G4Tubs("BodyBase5Top", LenseShellBody::radius5 - LenseShellBody::thickness * 2,
                                        LenseShellBody::radius6 - LenseShellBody::thickness, LenseShellBody::thickness,
                                        0, 360 * deg);
    G4VSolid* bodyBase5Bottom = new G4Tubs("BodyBase5Bottom", LenseShellBody::radius7 - LenseShellBody::thickness * 2,
                                           LenseShellBody::radius6 - LenseShellBody::thickness,
                                           LenseShellBody::thickness, 0, 360 * deg);
    G4VSolid* bodyBase5Inc = new G4UnionSolid("BodyBase5Inc", bodyBase5Wall, bodyBase5Top, nullptr,
                                              G4ThreeVector(0, 0, LenseShellBody::height6 - LenseShellBody::thickness));
    G4VSolid* bodyBase5 = new G4UnionSolid("BodyBase5", bodyBase5Inc, bodyBase5Bottom, nullptr,
                                           G4ThreeVector(0, 0, LenseShellBody::thickness - LenseShellBody::height6));
    tempOffsetZ = tempOffsetZ - LenseShellBody::height4 - LenseShellBody::height5;
    G4VSolid* bodyBaseInc5 = new G4UnionSolid("BodyBaseInc5", bodyBaseInc4, bodyBase5, nullptr,
                                              G4ThreeVector(0, 0, tempOffsetZ - LenseShellBody::height6 -
                                                            LenseShellBody::height5));

    G4VSolid* bodyBase6 = new G4Tubs("BodyBase6", LenseShellBody::radius7 - LenseShellBody::thickness * 2,
                                     LenseShellBody::radius7, LenseShellBody::height7 + LenseShellBody::thickness,
                                     0, 360 * deg);
    tempOffsetZ = tempOffsetZ - LenseShellBody::height5 - LenseShellBody::height6;
    G4VSolid* bodyBaseInc6 = new G4UnionSolid("BodyBaseInc6", bodyBaseInc5, bodyBase6, nullptr,
                                              G4ThreeVector(0, 0, tempOffsetZ - LenseShellBody::height7 -
                                                            LenseShellBody::height6));

    G4VSolid* bodyBase7Wall = new G4Tubs("BodyBase7Wall", LenseShellBody::radius8 - LenseShellBody::thickness * 2,
                                         LenseShellBody::radius8, LenseShellBody::height8, 0, 360 * deg);
    G4VSolid* bodyBase7Top = new G4Tubs("BodyBase7Top", LenseShellBody::radius7 - LenseShellBody::thickness * 2,
                                        LenseShellBody::radius8 - LenseShellBody::thickness, LenseShellBody::thickness,
                                        0, 360 * deg);
    G4VSolid* bodyBase7Bottom = new G4Tubs("BodyBase7Bottom", LenseShellBottom::radius - LenseShellBody::thickness * 2,
                                           LenseShellBody::radius8 - LenseShellBody::thickness,
                                           LenseShellBody::thickness, 0, 360 * deg);
    G4VSolid* bodyBase7Inc = new G4UnionSolid("BodyBase7Inc", bodyBase7Wall, bodyBase7Top, nullptr,
                                              G4ThreeVector(0, 0, LenseShellBody::height8 - LenseShellBody::thickness));
    G4VSolid* bodyBase7 = new G4UnionSolid("BodyBase7", bodyBase7Inc, bodyBase7Bottom, nullptr,
                                           G4ThreeVector(0, 0, LenseShellBody::thickness - LenseShellBody::height8));
    tempOffsetZ = tempOffsetZ - LenseShellBody::height7 - LenseShellBody::height6;
    G4VSolid* shellBody = new G4UnionSolid("LenseShellBody", bodyBaseInc6, bodyBase7, nullptr,
                                           G4ThreeVector(0, 0, tempOffsetZ - LenseShellBody::height7 -
                                                         LenseShellBody::height8));

    G4VSolid* shellInc = new G4UnionSolid("LenseShellInc", shellTop, shellBody, nullptr,
                                          G4ThreeVector(0, 0, -LenseShellTop::height - LenseShellBody::height1 +
                                                        LenseShellTop::delta));

    // Shell Bottom
    G4VSolid* bottomPlasticBase = new G4Tubs("BottomPlasticBase", Lenses::radiusBottom, LenseShellBottom::radius,
                                             LenseShellBottom::height, 0, 360 * deg);

    G4VSolid* bottomSub1 = new G4Cons("BottomPlasticSub1", 0,
                                      Lenses::radiusBottom + (LenseShellBottom::hole2Radius - Lenses::radiusBottom) * (1
                                          + LenseShellBottom::hole3Deep / LenseShellBottom::hole2Deep), 0,
                                      Lenses::radiusBottom, LenseShellBottom::hole2Deep + LenseShellBottom::hole3Deep,
                                      0, 360 * deg);
    G4VSolid* bottomBase1 = new G4SubtractionSolid("BottomPlasticBase1", bottomPlasticBase, bottomSub1, nullptr,
                                                   G4ThreeVector(0, 0, LenseShellBottom::height - 2 *
                                                                 LenseShellBottom::hole1Deep -
                                                                 LenseShellBottom::hole2Deep -
                                                                 LenseShellBottom::hole3Deep));

    G4VSolid* bottomSub2 = new G4Tubs("BottomPlasticSub2", 0, LenseShellBottom::hole3Radius, LenseShellBottom::height,
                                      0, 360 * deg);
    G4VSolid* bottomBase2 = new G4SubtractionSolid("BottomPlasticBase2", bottomBase1, bottomSub2, nullptr,
                                                   G4ThreeVector(0, 0, -2 * (LenseShellBottom::hole1Deep +
                                                                     LenseShellBottom::hole2Deep)));

    G4VSolid* bottomSub3 = new G4Cons("BottomPlasticSub3", 0,
                                      LenseShellBottom::hole3Radius + (LenseShellBottom::hole4Radius -
                                          LenseShellBottom::hole3Radius) * 2, 0, LenseShellBottom::hole3Radius,
                                      LenseShellBottom::hole3Deep * 2, 0, 360 * deg);
    G4VSolid* shellBottom = new G4SubtractionSolid("LenseShellBottom", bottomBase2, bottomSub3, nullptr,
                                                   G4ThreeVector(0, 0, -LenseShellBottom::height));

    G4VSolid* shell = new G4UnionSolid("LenseShell", shellInc, shellBottom, nullptr,
                                       G4ThreeVector(0, 0, -LenseShellTop::height - LenseShellBottom::height -
                                                     LenseShellBody::wholeHeight * 2 + LenseShellTop::delta));

    lenseShellLV = new G4LogicalVolume(shell, lenseShellMat, "LenseShellLV");
    lenseShellLV->SetVisAttributes(visLenseShell);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, CameraBox::halfZ - LenseShellTop::height - Lenses::lenseOffsetZ),
                      lenseShellLV, "LensePVP", cameraLV, false, 0, true);

    retainer = new G4Tubs("Retainer", LenseShellBody::radius5, LenseShellBody::retainerRadius,
                          LenseShellBody::retainerHeight, 0, 35 * pi / 18);
    retainerLV = new G4LogicalVolume(retainer, retainerMat, "RetainerLV");
    retainerLV->SetVisAttributes(visLenseShell);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, CameraBox::halfZ - (LenseShellTop::height + LenseShellBody::wholeHeight) * 2 -
                                    Lenses::lenseOffsetZ + LenseShellTop::delta + (LenseShellBody::height8 +
                                        LenseShellBody::height7 + LenseShellBody::height6) * 2 +
                                    LenseShellBody::retainerHeight), retainerLV, "RetainerPVP", cameraLV, false, 0,
                      true);
}

void Camera::ConstructLenseHolder() {
    G4VSolid* holderBase = new G4Box("HolderBase", LenseHolder::length, LenseHolder::length, LenseHolder::height);

    G4VSolid* holderSub1 = new G4Tubs("HolderSub1", 0, LenseHolder::radius, LenseHolder::height * 2, 0, 360 * deg);
    G4VSolid* holderBase0 = new G4SubtractionSolid("HolderBase0", holderBase, holderSub1, nullptr,
                                                   G4ThreeVector(0, 0, 0));

    G4VSolid* holderLeg = new G4Tubs("HolderLeg", 0, LenseHolder::legRadiusOuter,
                                     LenseHolder::legHeight + LenseHolder::height / 2, 0, 360 * deg);
    G4VSolid* holderBase1 = new G4UnionSolid("HolderBase1", holderBase0, holderLeg, nullptr,
                                             G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                           LenseHolder::widht - LenseHolder::holeOffsetY,
                                                           -LenseHolder::legHeight - LenseHolder::height));
    G4VSolid* holderBase2 = new G4UnionSolid("HolderBase2", holderBase1, holderLeg, nullptr,
                                             G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                           LenseHolder::widht - LenseHolder::holeOffsetY,
                                                           -LenseHolder::legHeight - LenseHolder::height));
    G4VSolid* holderBase3 = new G4UnionSolid("HolderBase2", holderBase2, holderLeg, nullptr,
                                             G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                           -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                           -LenseHolder::legHeight - LenseHolder::height));
    G4VSolid* holderBase4 = new G4UnionSolid("HolderBase4", holderBase3, holderLeg, nullptr,
                                             G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                           -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                           -LenseHolder::legHeight - LenseHolder::height));

    G4VSolid* holderSub2 = new G4Tubs("HolderSub2", 0, LenseHolder::legRadiusInner, (LenseHolder::height +
                                          LenseHolder::legHeight) * 2, 0, 360 * deg);
    G4VSolid* holderBase5 = new G4SubtractionSolid("HolderBase5", holderBase4, holderSub2, nullptr,
                                                   G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                                 LenseHolder::widht - LenseHolder::holeOffsetY,
                                                                 -LenseHolder::legHeight - LenseHolder::height));
    G4VSolid* holderBase6 = new G4SubtractionSolid("HolderBase6", holderBase5, holderSub2, nullptr,
                                                   G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                                 LenseHolder::widht - LenseHolder::holeOffsetY,
                                                                 -LenseHolder::legHeight - LenseHolder::height));
    G4VSolid* holderBase7 = new G4SubtractionSolid("HolderBase5", holderBase6, holderSub2, nullptr,
                                                   G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                                 -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                                 -LenseHolder::legHeight - LenseHolder::height));
    lenseHolder = new G4SubtractionSolid("LenseHolder", holderBase7, holderSub2, nullptr,
                                         G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                       -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                       -LenseHolder::legHeight - LenseHolder::height));

    lenseHolderLV = new G4LogicalVolume(lenseHolder, lenseHolderMat, "LenseHolderLV");
    lenseHolderLV->SetVisAttributes(visLenseHolder);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, CameraBox::halfZ - LenseHolder::height), lenseHolderLV,
                      "LenseHolderPVP", cameraLV, false, 0, true);
}

void Camera::ConstructCameraPlate() {
    G4VSolid* plateBase0 = new G4Box("CameraPlateBase", CameraPlate::halfX, CameraPlate::halfY, CameraPlate::halfZ);

    G4VSolid* sub1 = new G4Box("CameraPlateSub1", CameraPlate::ledge * 2, CameraPlate::ledge * 2,
                               CameraPlate::halfZ * 2);
    G4VSolid* plateBase1 = new G4SubtractionSolid("CameraPlateBase1", plateBase0, sub1, nullptr,
                                                  G4ThreeVector(CameraPlate::halfX, CameraPlate::halfY));
    G4VSolid* plateBase2 = new G4SubtractionSolid("CameraPlateBase2", plateBase1, sub1, nullptr,
                                                  G4ThreeVector(CameraPlate::halfX, -CameraPlate::halfY));
    G4VSolid* plateBase3 = new G4SubtractionSolid("CameraPlateBase3", plateBase2, sub1, nullptr,
                                                  G4ThreeVector(-CameraPlate::halfX, -CameraPlate::halfY));
    G4VSolid* plateBase4 = new G4SubtractionSolid("CameraPlateBase4", plateBase3, sub1, nullptr,
                                                  G4ThreeVector(-CameraPlate::halfX, CameraPlate::halfY));

    G4VSolid* sub2 = new G4Tubs("CameraPlateSub2", 0, CameraPlate::holeRadius, CameraPlate::halfZ * 2, 0, 360 * deg);
    G4VSolid* plateBase5 = new G4SubtractionSolid("CameraPlateBase5", plateBase4, sub2, nullptr, G4ThreeVector());

    G4VSolid* sub3 = new G4Tubs("CameraPlateSub3", 0, CameraPlate::recessRadius, CameraPlate::recessDeep * 2, 0,
                                360 * deg);
    G4VSolid* plateBase6 = new G4SubtractionSolid("CameraPlateBase6", plateBase5, sub3, nullptr,
                                                  G4ThreeVector(0, 0, -CameraPlate::halfZ));

    G4VSolid* sub4 = new G4Tubs("CameraPlateSub4", 0, LenseHolder::legRadiusOuter * 1,
                                CameraPlate::lenseHolderHolesDeep * 2, 0, 360 * deg);
    G4VSolid* plateBase7 = new G4SubtractionSolid("CameraPlateBase7", plateBase6, sub4, nullptr,
                                                  G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                                LenseHolder::widht - LenseHolder::holeOffsetY,
                                                                CameraPlate::halfZ));
    G4VSolid* plateBase8 = new G4SubtractionSolid("CameraPlateBase8", plateBase7, sub4, nullptr,
                                                  G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                                LenseHolder::widht - LenseHolder::holeOffsetY,
                                                                CameraPlate::halfZ));
    G4VSolid* plateBase9 = new G4SubtractionSolid("CameraPlateBase9", plateBase8, sub4, nullptr,
                                                  G4ThreeVector(LenseHolder::length - LenseHolder::holeOffsetX,
                                                                -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                                CameraPlate::halfZ));
    G4VSolid* plateBase10 = new G4SubtractionSolid("CameraPlateBase10", plateBase9, sub4, nullptr,
                                                   G4ThreeVector(-(LenseHolder::length - LenseHolder::holeOffsetX),
                                                                 -(LenseHolder::widht - LenseHolder::holeOffsetY),
                                                                 CameraPlate::halfZ));

    G4VSolid* sub5 = new G4Tubs("CameraPlateSub5", 0, ElectronicsHolder::legRadiusOuter * 1,
                                CameraPlate::electronicsHolderHolesDeep * 2, 0, 360 * deg);
    G4VSolid* plateBase11 = new G4SubtractionSolid("CameraPlateBase11", plateBase10, sub5, nullptr,
                                                   G4ThreeVector(ElectronicsHolder::length -
                                                                 ElectronicsHolder::holeOffsetX,
                                                                 ElectronicsHolder::widht -
                                                                 ElectronicsHolder::holeOffsetY,
                                                                 -CameraPlate::halfZ));
    G4VSolid* plateBase12 = new G4SubtractionSolid("CameraPlateBase12", plateBase11, sub5, nullptr,
                                                   G4ThreeVector(-(ElectronicsHolder::length -
                                                                     ElectronicsHolder::holeOffsetX),
                                                                 ElectronicsHolder::widht -
                                                                 ElectronicsHolder::holeOffsetY,
                                                                 -CameraPlate::halfZ));
    G4VSolid* plateBase13 = new G4SubtractionSolid("CameraPlateBase13", plateBase12, sub5, nullptr,
                                                   G4ThreeVector(ElectronicsHolder::length -
                                                                 ElectronicsHolder::holeOffsetX,
                                                                 -(ElectronicsHolder::widht -
                                                                     ElectronicsHolder::holeOffsetY),
                                                                 -CameraPlate::halfZ));
    G4VSolid* plateBase14 = new G4SubtractionSolid("CameraPlateBase14", plateBase13, sub5, nullptr,
                                                   G4ThreeVector(-(ElectronicsHolder::length -
                                                                     ElectronicsHolder::holeOffsetX),
                                                                 -(ElectronicsHolder::widht -
                                                                     ElectronicsHolder::holeOffsetY),
                                                                 -CameraPlate::halfZ));

    cameraPlateLV = new G4LogicalVolume(plateBase14, cameraPlateMat, "CameraPlateLV");
    cameraPlateLV->SetVisAttributes(visCameraPlate);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -CameraBox::halfZ + CameraPlate::halfZ), cameraPlateLV,
                      "CameraPlatePVP", cameraLV, false, 0, true);
}


void Camera::ConstructElectronicsHolder() {
    G4VSolid* holderBase0 = new G4Box("EHolderBase", ElectronicsHolder::length, ElectronicsHolder::length,
                                      ElectronicsHolder::height);

    G4VSolid* holderLeg = new G4Tubs("EHolderLeg", 0, ElectronicsHolder::legRadiusOuter,
                                     ElectronicsHolder::legHeight + ElectronicsHolder::height / 2, 0, 360 * deg);
    G4VSolid* holderBase1 = new G4UnionSolid("EHolderBase1", holderBase0, holderLeg, nullptr,
                                             G4ThreeVector(ElectronicsHolder::length - ElectronicsHolder::holeOffsetX,
                                                           ElectronicsHolder::widht - ElectronicsHolder::holeOffsetY,
                                                           ElectronicsHolder::legHeight + ElectronicsHolder::height /
                                                           2));
    G4VSolid* holderBase2 = new G4UnionSolid("EHolderBase2", holderBase1, holderLeg, nullptr,
                                             G4ThreeVector(-(ElectronicsHolder::length -
                                                               ElectronicsHolder::holeOffsetX),
                                                           ElectronicsHolder::widht - ElectronicsHolder::holeOffsetY,
                                                           ElectronicsHolder::legHeight + ElectronicsHolder::height /
                                                           2));
    G4VSolid* holderBase3 = new G4UnionSolid("EHolderBase2", holderBase2, holderLeg, nullptr,
                                             G4ThreeVector(ElectronicsHolder::length - ElectronicsHolder::holeOffsetX,
                                                           -(ElectronicsHolder::widht - ElectronicsHolder::holeOffsetY),
                                                           ElectronicsHolder::legHeight + ElectronicsHolder::height /
                                                           2));
    G4VSolid* holderBase4 = new G4UnionSolid("EHolderBase4", holderBase3, holderLeg, nullptr,
                                             G4ThreeVector(-(ElectronicsHolder::length -
                                                               ElectronicsHolder::holeOffsetX),
                                                           -(ElectronicsHolder::widht - ElectronicsHolder::holeOffsetY),
                                                           ElectronicsHolder::legHeight + ElectronicsHolder::height /
                                                           2));

    G4VSolid* holderSub2 = new G4Tubs("EHolderSub2", 0, ElectronicsHolder::legRadiusInner, (ElectronicsHolder::height +
                                          ElectronicsHolder::legHeight) * 2, 0, 360 * deg);
    G4VSolid* holderBase5 = new G4SubtractionSolid("EHolderBase5", holderBase4, holderSub2, nullptr,
                                                   G4ThreeVector(ElectronicsHolder::length -
                                                                 ElectronicsHolder::holeOffsetX,
                                                                 ElectronicsHolder::widht -
                                                                 ElectronicsHolder::holeOffsetY,
                                                                 ElectronicsHolder::legHeight +
                                                                 ElectronicsHolder::height));
    G4VSolid* holderBase6 = new G4SubtractionSolid("EHolderBase6", holderBase5, holderSub2, nullptr,
                                                   G4ThreeVector(-(ElectronicsHolder::length -
                                                                     ElectronicsHolder::holeOffsetX),
                                                                 ElectronicsHolder::widht -
                                                                 ElectronicsHolder::holeOffsetY,
                                                                 ElectronicsHolder::legHeight +
                                                                 ElectronicsHolder::height));
    G4VSolid* holderBase7 = new G4SubtractionSolid("EHolderBase5", holderBase6, holderSub2, nullptr,
                                                   G4ThreeVector(ElectronicsHolder::length -
                                                                 ElectronicsHolder::holeOffsetX,
                                                                 -(ElectronicsHolder::widht -
                                                                     ElectronicsHolder::holeOffsetY),
                                                                 ElectronicsHolder::legHeight +
                                                                 ElectronicsHolder::height));
    electronicsHolder = new G4SubtractionSolid("ElectronicsHolder", holderBase7, holderSub2, nullptr,
                                               G4ThreeVector(-(ElectronicsHolder::length -
                                                                 ElectronicsHolder::holeOffsetX),
                                                             -(ElectronicsHolder::widht -
                                                                 ElectronicsHolder::holeOffsetY),
                                                             ElectronicsHolder::legHeight + ElectronicsHolder::height));

    electronicsHolderLV = new G4LogicalVolume(electronicsHolder, electronicsHolderMat, "ElectronicsHolderLV");
    electronicsHolderLV->SetVisAttributes(visElectronicsHolder);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -CameraBox::halfZ - CameraBox::bottomHalfZ * 2 + ElectronicsHolder::height),
                      electronicsHolderLV, "ElectronicsHolderPVP", cameraLV, false, 0, true);
}


void Camera::ConstructElectronics() {
    G4VSolid* eBase = new G4Tubs("EBase", Electronics::innerRadius, Electronics::outerRadius, Electronics::height, 0,
                                 360 * deg);

    G4double zPlanes[2] = {-Electronics::pentagonHeight * 2, Electronics::pentagonHeight * 2};
    G4double rInner[2] = {0.0, 0.0};
    G4double rOuter[2] = {Electronics::pentagonRadius, Electronics::pentagonRadius};
    auto* pentagon = new G4Polyhedra("Pentagon", 0.0 * deg, 360.0 * deg, 5, 2, zPlanes, rInner, rOuter);
    G4VSolid* eBase0 = new G4UnionSolid("EBase0", eBase, pentagon, nullptr, G4ThreeVector(0, 0, Electronics::height));

    G4VSolid* sub1 = new G4Tubs("ESub1", 0, Electronics::innerRadius,
                                (Electronics::pentagonHeight + Electronics::height) * 2, 0, 360 * deg);
    G4VSolid* eBase1 = new G4SubtractionSolid("EBase1", eBase0, sub1);

    electronicsLV = new G4LogicalVolume(eBase1, electronicsMat, "ElectronicsLV");
    electronicsLV->SetVisAttributes(visElectronics);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -CameraBox::halfZ + Electronics::height - Electronics::offsetZ),
                      electronicsLV, "ElectronicsPVP", cameraLV, false, 0, true);

    matrix = new G4Tubs("Matrix", 0, Electronics::innerRadius, Electronics::matrixHeight, 0, 360 * deg);

    matrixLV = new G4LogicalVolume(matrix, matrixMat, "MatrixLV");
    matrixLV->SetVisAttributes(visMatrix);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -CameraBox::halfZ + Electronics::height * 2 + Electronics::matrixHeight -
                                    Electronics::offsetZ), matrixLV, "ElectronicsPVP", cameraLV, false, 0, true);
}

void Camera::ConstructCamera() {
    G4VSolid* cameraBase = new G4Box("CameraBase", CameraBox::halfX, CameraBox::halfY, CameraBox::halfZ);

    G4VSolid* sub = new G4Box("CameraSub", CameraBox::ledge * 2, CameraBox::ledge * 2, CameraBox::halfZ * 2);
    G4VSolid* cameraBase1 = new G4SubtractionSolid("CameraBase1", cameraBase, sub, nullptr,
                                                   G4ThreeVector(CameraBox::halfX, CameraBox::halfY));
    G4VSolid* cameraBase2 = new G4SubtractionSolid("CameraBase2", cameraBase1, sub, nullptr,
                                                   G4ThreeVector(CameraBox::halfX, -CameraBox::halfY));
    G4VSolid* cameraBase3 = new G4SubtractionSolid("CameraBase3", cameraBase2, sub, nullptr,
                                                   G4ThreeVector(-CameraBox::halfX, -CameraBox::halfY));
    G4VSolid* cameraBase4 = new G4SubtractionSolid("CameraBase4", cameraBase3, sub, nullptr,
                                                   G4ThreeVector(-CameraBox::halfX, CameraBox::halfY));

    G4VSolid* bottomCube = new G4Box("BottomCube", CameraBox::bottomHalfX, CameraBox::bottomHalfY,
                                     CameraBox::bottomHalfZ + CameraBox::halfZ / 2);
    G4VSolid* camera = new G4UnionSolid("Camera", cameraBase4, bottomCube, nullptr,
                                        G4ThreeVector(0, 0, -CameraBox::halfZ / 2 - CameraBox::bottomHalfZ));

    cameraLV = new G4LogicalVolume(camera, vacuumMat, "CameraLV");
    cameraLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    ConstructLense();
    ConstructLenseHolder();
    ConstructCameraPlate();
    ConstructElectronicsHolder();
    ConstructElectronics();
}
