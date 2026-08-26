#include "geometry/CubeSat_GeoScan_3U.hh"

using namespace CubeSatSizes;

CubeSat_GeoScan_3U::CubeSat_GeoScan_3U(G4LogicalVolume* world, G4NistManager* nistManager) : worldLV(world),
    nist(nistManager) {
    DefineVisual();
    DefineMaterial();

    // cubeSat = new G4Box("CubeSat", CubeSat::halfX, CubeSat::halfY,
    //                     CubeSat::halfZ + CameraSizes::CameraBox::halfZ - Frame::bracketHeight);

    cubeSat = new G4Box("CubeSat", CubeSat::halfX, CubeSat::halfY, CubeSat::halfZ);
    cubeSatLV = new G4LogicalVolume(cubeSat, vacuumMat, "CubeSatLV");
    cubeSatLV->SetVisAttributes(G4VisAttributes::GetInvisible());
}


void CubeSat_GeoScan_3U::DefineVisual() {
    visFrame = new G4VisAttributes(G4Color(G4Color(0.4, 0.4, 0.4)));
    visFrame->SetForceSolid(true);
    visHolder = new G4VisAttributes(G4Color(G4Color(0.6, 0.6, 0.6)));
    visHolder->SetForceSolid(true);
    visSolarPanel = new G4VisAttributes(G4Color(G4Color(0.25, 0.25, 0.25)));
    visSolarPanel->SetForceSolid(true);
    visCover = new G4VisAttributes(G4Color(G4Color(1.0, 1.0, 1.0)));
    visCover->SetForceSolid(true);
    visBoard = new G4VisAttributes(G4Color(G4Color(0.0, 1.0, 0.0)));
    visBoard->SetForceSolid(true);
    visMechanics = new G4VisAttributes(G4Color(G4Color(1.0, 1.0, 1.0)));
    visMechanics->SetForceSolid(true);
}

void CubeSat_GeoScan_3U::DefineMaterial() {
    auto* elH = nist->FindOrBuildElement("H");
    auto* elC = nist->FindOrBuildElement("C");

    G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    frameMat = nist->FindOrBuildMaterial("G4_Al");
    lintelMat = nist->FindOrBuildMaterial("G4_Al");
    holderMat = nist->FindOrBuildMaterial("G4_Al");
    coverMat = nist->FindOrBuildMaterial("G4_Al");
    solarPanelMat = nist->FindOrBuildMaterial("G4_Ge");
    {
        auto* Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, 2);
        Epoxy->AddElement(elH, 2);
        Epoxy->AddElement(elC, 2);

        boardMat = new G4Material("FiberglassLaminate", 1.86 * g / cm3, 2);
        boardMat->AddMaterial(Epoxy, 0.472);
        boardMat->AddMaterial(SiO2, 0.528);
    }
    mechanicsMat = nist->FindOrBuildMaterial("G4_Fe");
    vacuumMat = nist->FindOrBuildMaterial("G4_Galactic");
}

void CubeSat_GeoScan_3U::ConstructCubeSat() {
    // auto camera = Camera(worldLV, nist);
    // camera.ConstructCamera();
    // cameraLV = camera.GetCameraLV();
    // new G4PVPlacement(nullptr,
    //                   G4ThreeVector(0, 0, CubeSat::halfZ - Frame::bracketHeight), cameraLV, "CameraPVP", cubeSatLV,
    //                   false, 0, true);

    ConstructBridge();
    ConstructLintel();
    ConstructFixator();
    ConstructUnitFrame();
    ConstructSolarPanelHolder();
    ConstructSolarPanel();
    ConstructMechanics();
    ConstructCover();
    ConstructServiceSystem();
    ConstructFrame();
}


void CubeSat_GeoScan_3U::ConstructBridge() {
    G4VSolid* base = new G4Box("BridgeBase", Unit::halfX, UnitFrame::halfY, Frame::gap);
    G4VSolid* sub0 = new G4Box("BridgeSub0", Unit::halfX - 2 * Frame::thickness, UnitFrame::halfY, 2 * Frame::gap);
    G4VSolid* sub1 = new G4Box("BridgeSub1", Unit::halfX - 2 * UnitFrame::halfX, 2 * UnitFrame::halfY, 2 * Frame::gap);

    G4VSolid* base0 = new G4SubtractionSolid("BridgeBase0", base, sub0, nullptr,
                                             G4ThreeVector(0, 2 * UnitFrame::thickness, 0));
    bridge = new G4SubtractionSolid("BridgeBase0", base0, sub1, nullptr,
                                    G4ThreeVector(0, 0, 0));
    bridgeLV = new G4LogicalVolume(bridge, frameMat, "BridgeLV");
    bridgeLV->SetVisAttributes(visFrame);
}

void CubeSat_GeoScan_3U::ConstructLintel() {
    G4VSolid* lintelBase = new G4Box("LintelBase", Lintel::thickness, Lintel::length, Lintel::width);
    G4VSolid* lintelBracket = new G4Box("LintelBracket", Lintel::bracketX + Lintel::thickness, Lintel::bracketY,
                                        Lintel::width);
    G4VSolid* lintelLedge = new G4Box("LintelLedge", Lintel::ledgeHeight + Lintel::thickness,
                                      Lintel::length - 2 * Lintel::bracketY, Lintel::ledgeThickness);
    G4VSolid* lintelFixLedgeX = new G4Box("LintelFixLedgeX", Lintel::bracketX + Lintel::thickness,
                                          Lintel::fixLedgeHeight, Lintel::fixLedgeXThickness);
    G4VSolid* lintelFixLedgeY = new G4Box("LintelFixLedgeY", Lintel::fixLedgeZThickness,
                                          Lintel::fixLedgeHeight, Lintel::width);
    G4VSolid* bracket = new G4Box("Bracket", Bridge::fixLedgeLength, Bridge::fixLedgeHeight, Bridge::fixLedgeWidth);

    G4VSolid* lintel0 = new G4UnionSolid("Lintel0", lintelBase, lintelLedge, nullptr,
                                         G4ThreeVector(Lintel::ledgeHeight, 0,
                                                       -Lintel::width + Lintel::ledgeThickness));
    G4VSolid* lintel1 = new G4UnionSolid("Lintel1", lintel0, lintelFixLedgeY, nullptr,
                                         G4ThreeVector(-Lintel::thickness + Lintel::fixLedgeZThickness,
                                                       -Lintel::length + Lintel::fixLedgeHeight, 0));
    G4VSolid* lintel2 = new G4UnionSolid("Lintel2", lintel1, lintelFixLedgeY, nullptr,
                                         G4ThreeVector(-Lintel::thickness + Lintel::fixLedgeZThickness,
                                                       Lintel::length - Lintel::fixLedgeHeight, 0));
    G4VSolid* lintel3 = new G4UnionSolid("Lintel3", lintel2, lintelFixLedgeX, nullptr,
                                         G4ThreeVector(Lintel::bracketX,
                                                       -Lintel::length + Lintel::fixLedgeHeight,
                                                       Lintel::fixLedgeXThickness - Lintel::width));
    G4VSolid* lintel4 = new G4UnionSolid("Lintel4", lintel3, lintelFixLedgeX, nullptr,
                                         G4ThreeVector(Lintel::bracketX,
                                                       Lintel::length - Lintel::fixLedgeHeight,
                                                       Lintel::fixLedgeXThickness - Lintel::width));
    G4VSolid* lintel5 = new G4UnionSolid("Lintel5", lintel4, lintelBracket, nullptr,
                                         G4ThreeVector(Lintel::bracketX,
                                                       Lintel::length - Lintel::bracketY - Lintel::fixLedgeHeight * 2,
                                                       0));
    G4VSolid* lintel6 = new G4UnionSolid("Lintel6", lintel5, lintelBracket, nullptr,
                                         G4ThreeVector(Lintel::bracketX,
                                                       -Lintel::length + Lintel::bracketY + Lintel::fixLedgeHeight * 2,
                                                       0));
    G4VSolid* lintel7 = new G4UnionSolid("Lintel7", lintel6, bracket, nullptr,
                                         G4ThreeVector(Lintel::bracketX * 2 - Bridge::fixLedgeLength +
                                                       Lintel::thickness,
                                                       -Lintel::length - Lintel::gap * 2 + Bridge::fixLedgeHeight,
                                                       Lintel::width - Bridge::fixLedgeWidth));
    lintel = new G4UnionSolid("Lintel", lintel7, bracket, nullptr,
                              G4ThreeVector(Lintel::bracketX * 2 - Bridge::fixLedgeLength + Lintel::thickness,
                                            Lintel::length + Lintel::gap * 2 - Bridge::fixLedgeHeight,
                                            Lintel::width - Bridge::fixLedgeWidth));

    lintelLV = new G4LogicalVolume(lintel, lintelMat, "LintelLV");
    lintelLV->SetVisAttributes(visFrame);
}


void CubeSat_GeoScan_3U::ConstructUnitFrame() {
    G4VSolid* base = new G4Box("Base", Unit::halfX, UnitFrame::halfY, UnitFrame::halfZ);

    G4VSolid* sub0 = new G4Box("Sub0", Unit::halfX - 2 * UnitFrame::thickness, UnitFrame::halfY,
                               UnitFrame::halfZ - 2 * Bridge::width);
    G4VSolid* base0 = new G4SubtractionSolid("Base0", base, sub0, nullptr,
                                             G4ThreeVector(0, 2 * UnitFrame::thickness, 0));

    G4VSolid* sub1 = new G4Box("Sub1", Unit::halfX - 2 * Bridge::bracketX, UnitFrame::halfY,
                               UnitFrame::halfZ - 2 * Bridge::ledgeThickness);
    G4VSolid* base1 = new G4SubtractionSolid("Base1", base0, sub1, nullptr,
                                             G4ThreeVector(0, 2 * (Bridge::thickness + UnitFrame::thickness -
                                                               UnitFrame::ledgeThickness), 0));

    G4VSolid* sub2 = new G4Box("Sub2", Unit::halfX - 2 * (UnitFrame::halfX + UnitFrame::ledgeHeight),
                               UnitFrame::halfY + 2 * UnitFrame::thickness, UnitFrame::halfZ - 2 * Bridge::width);
    G4VSolid* base2 = new G4SubtractionSolid("Base2", base1, sub2, nullptr, G4ThreeVector(0, 0, 0));

    G4VSolid* sub3 = new G4Box("Sub3", Unit::halfX - 2 * Bridge::fixLedgeX, UnitFrame::halfY, UnitFrame::halfZ * 2);
    G4VSolid* base3 = new G4SubtractionSolid("Base3", base2, sub3, nullptr,
                                             G4ThreeVector(0, 2 * (Bridge::bracketY + Bridge::thickness +
                                                               UnitFrame::thickness - UnitFrame::ledgeThickness), 0));

    G4VSolid* sub4 = new G4Box("Sub4", Unit::halfX - 2 * UnitFrame::halfX, UnitFrame::halfY, 2 * UnitFrame::halfZ);
    G4VSolid* base4 = new G4SubtractionSolid("Base4", base3, sub4, nullptr,
                                             G4ThreeVector(0, 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness -
                                                               UnitFrame::halfY), 0));

    G4VSolid* ledgeY = new G4Box("LedgeX", UnitFrame::ledgeThickness, UnitFrame::ledgeHeight * 2, UnitFrame::ledgeY);
    G4VSolid* base5 = new G4UnionSolid("Base5", base4, ledgeY, nullptr,
                                       G4ThreeVector(Unit::halfX - 2 * UnitFrame::thickness + UnitFrame::ledgeThickness,
                                                     UnitFrame::halfY, 0));
    unitFrame = new G4UnionSolid("Base6", base5, ledgeY, nullptr,
                                 G4ThreeVector(-(Unit::halfX - 2 * UnitFrame::thickness +
                                                   UnitFrame::ledgeThickness), UnitFrame::halfY, 0));

    unitFrameLV = new G4LogicalVolume(unitFrame, frameMat, "UnitFrameLV");
    unitFrameLV->SetVisAttributes(visFrame);
}

void CubeSat_GeoScan_3U::ConstructFixator() {
    G4VSolid* fixatorBase = new G4Box("FixatorBase", Fixator::length, UnitFrame::halfY - UnitFrame::thickness,
                                      Fixator::width);

    G4VSolid* sub1 = new G4Box("FixatorSub1", Fixator::holeX * 2, UnitFrame::halfY, Fixator::holeY);
    G4VSolid* base1 = new G4SubtractionSolid("FixatorBase1", fixatorBase, sub1, nullptr,
                                             G4ThreeVector(Fixator::length, 0,
                                                           Fixator::width - 2 * Bridge::width - Fixator::holeY));

    G4VSolid* sub2 = new G4Box("FixatorSub2",
                               (Fixator::length - (UnitFrame::halfX + UnitFrame::ledgeHeight - UnitFrame::thickness)) *
                               2, UnitFrame::halfY - UnitFrame::thickness, (Fixator::width - Bridge::width) * 2);
    G4VSolid* base2 = new G4SubtractionSolid("FixatorBase2", base1, sub2, nullptr,
                                             G4ThreeVector(-Fixator::length, 2 * Fixator::thickness, -Fixator::width));

    G4VSolid* sub3 = new G4Box("FixatorSub3", Fixator::cutX * 2, UnitFrame::halfY, Fixator::cutY * 2);
    G4VSolid* base3 = new G4SubtractionSolid("FixatorBase3", base2, sub3, nullptr,
                                             G4ThreeVector(-Fixator::length, 0, Fixator::width));

    G4VSolid* sub4 = new G4Box("FixatorSub4", UnitFrame::halfX * 2, UnitFrame::halfY, Bridge::width * 2);
    G4VSolid* base4 = new G4SubtractionSolid("FixatorBase4", base3, sub4, nullptr,
                                             G4ThreeVector(Fixator::length, 0, Fixator::width));

    G4VSolid* base5 = new G4SubtractionSolid("FixatorBase5", base4, fixatorBase, nullptr,
                                             G4ThreeVector(0, 2 * (UnitFrame::halfY - UnitFrame::thickness +
                                                               UnitFrame::ledgeThickness - Bridge::ledgeThickness -
                                                               Bridge::bracketY), Fixator::width));

    G4VSolid* sub5 = new G4Box("FixatorSub5", Fixator::recessX, Fixator::recessDeep * 2, Fixator::recessY);
    fixator = new G4SubtractionSolid("Fixator", base5, sub5, nullptr,
                                     G4ThreeVector(-Fixator::length + Fixator::recessX + Fixator::recessGapX * 2,
                                                   -UnitFrame::halfY + UnitFrame::thickness,
                                                   -Fixator::width + Fixator::recessY + Fixator::recessGapY * 2));

    fixatorLV = new G4LogicalVolume(fixator, frameMat, "FixatorLV");
    fixatorLV->SetVisAttributes(visFrame);
}

void CubeSat_GeoScan_3U::ConstructSolarPanelHolder() {
    G4VSolid* holderBase = new G4Box("HolderBase", Holder::halfX, Holder::halfY, Holder::halfZ);

    G4VSolid* ledge = new G4Box("HolderLedge", Holder::ledgeX, Holder::ledgeY, Holder::ledgeZ * 2 + Holder::halfZ);
    G4VSolid* base0 = new G4UnionSolid("HolderBase0", holderBase, ledge, nullptr,
                                       G4ThreeVector(0, Holder::halfY - Holder::ledgeY, 0));

    G4VSolid* sub0 = new G4Box("HolderSub0", Holder::ledgeX - Holder::ledgeThickness * 2, Holder::ledgeY * 2,
                               (Holder::ledgeZ - Holder::ledgeThickness) * 2 + Holder::halfZ);
    G4VSolid* base1 = new G4SubtractionSolid("HolderBase1", base0, sub0, nullptr,
                                             G4ThreeVector(0, Holder::halfY + (Holder::ledgeY -
                                                               Holder::recess1Deep) * 2, 0));

    G4VSolid* sub1 = new G4Box("HolderSub1", Holder::halfX - Holder::ledgeThickness * 2, Holder::ledgeY * 2,
                               Holder::halfZ - Holder::recess1Ledge * 2);
    G4VSolid* base2 = new G4SubtractionSolid("HolderBase2", base1, sub1, nullptr,
                                             G4ThreeVector(0, Holder::halfY + (Holder::ledgeY -
                                                               Holder::recess1Deep) * 2, 0));

    G4VSolid* sub2 = new G4Box("HolderSub2", Holder::recess2X, Holder::recess2Y * 2,
                               Holder::recess2Z - Holder::recess1Radius);
    G4VSolid* base3 = new G4SubtractionSolid("HolderBase3", base2, sub2, nullptr,
                                             G4ThreeVector(0, Holder::halfY - Holder::recess1Deep * 2,
                                                           0));

    G4VSolid* sub3 = new G4Box("HolderSub3", Holder::recess2X - Holder::recess1Radius, Holder::recess2Y * 2,
                               Holder::recess2Z);
    G4VSolid* base4 = new G4SubtractionSolid("HolderBase4", base3, sub3, nullptr,
                                             G4ThreeVector(0, Holder::halfY - Holder::recess1Deep * 2,
                                                           0));

    auto* rotMat = new G4RotationMatrix(0, 180 * deg, 0);
    auto* rotMat0 = new G4RotationMatrix(0, 90 * deg, -5 * deg);
    auto* rotMat1 = new G4RotationMatrix(0, 90 * deg, 85 * deg);
    auto* rotMat2 = new G4RotationMatrix(0, 90 * deg, 265 * deg);
    auto* rotMat3 = new G4RotationMatrix(0, 90 * deg, 175 * deg);
    G4VSolid* sub4 = new G4Tubs("HolderSub4", 0, Holder::recess1Radius, Holder::recess2Y * 2, 0 * deg, 100 * deg);
    G4VSolid* base5 = new G4SubtractionSolid("HolderBase5", base4, sub4, rotMat0,
                                             G4ThreeVector(Holder::recess2X - Holder::recess1Radius,
                                                           Holder::halfY - Holder::recess1Deep * 2,
                                                           Holder::recess2Z - Holder::recess1Radius));
    G4VSolid* base6 = new G4SubtractionSolid("HolderBase6", base5, sub4, rotMat1,
                                             G4ThreeVector(-(Holder::recess2X - Holder::recess1Radius),
                                                           Holder::halfY - Holder::recess1Deep * 2,
                                                           Holder::recess2Z - Holder::recess1Radius));
    G4VSolid* base7 = new G4SubtractionSolid("HolderBase7", base6, sub4, rotMat2,
                                             G4ThreeVector(Holder::recess2X - Holder::recess1Radius,
                                                           Holder::halfY - Holder::recess1Deep * 2,
                                                           -(Holder::recess2Z - Holder::recess1Radius)));
    G4VSolid* base8 = new G4SubtractionSolid("HolderBase8", base7, sub4, rotMat3,
                                             G4ThreeVector(-(Holder::recess2X - Holder::recess1Radius),
                                                           Holder::halfY - Holder::recess1Deep * 2,
                                                           -(Holder::recess2Z - Holder::recess1Radius)));

    G4VSolid* sub5 = new G4Box("HolderSub5", Holder::halfX - Holder::ledgeBackThickness * 2, Holder::recessBackDeep * 2,
                               Holder::halfZ - Holder::recessBackLedgeY * 2);
    G4VSolid* base9 = new G4SubtractionSolid("HolderBase9", base8, sub5, nullptr,
                                             G4ThreeVector(0, -Holder::halfY, 0));

    G4VSolid* sub6 = new G4Box("HolderSub6", Holder::halfX - Holder::recessBackLedgeX * 2, Holder::recessBackDeep * 2,
                               Holder::halfZ - Holder::ledgeBackThickness * 2);
    G4VSolid* base10 = new G4SubtractionSolid("HolderBase10", base9, sub6, nullptr,
                                              G4ThreeVector(0, -Holder::halfY, 0));

    G4VSolid* sub7 = new G4Box("HolderSub7", Holder::recessBackGap, Holder::recessBackDeep * 2,
                               Holder::halfZ * 2);
    G4VSolid* base11 = new G4SubtractionSolid("HolderBase11", base10, sub7, nullptr,
                                              G4ThreeVector(0, -Holder::halfY, 0));

    G4VSolid* sub8 = new G4Box("HolderSub8", Holder::holeX1, Holder::halfY, Holder::holeZ1);
    G4VSolid* base12 = new G4SubtractionSolid("HolderBase12", base11, sub8, nullptr,
                                              G4ThreeVector(0, 0, -Holder::halfZ + Holder::holeGapZ * 2 +
                                                            Holder::holeZ1));

    G4VSolid* sub9 = new G4Box("HolderSub9", Holder::holeX2, Holder::halfY, Holder::holeZ2 * 2);
    G4VSolid* base13 = new G4SubtractionSolid("HolderBase13", base12, sub9, nullptr,
                                              G4ThreeVector(0, 0, -Holder::halfZ + 2 * (Holder::holeGapZ +
                                                                Holder::holeZ) - Holder::holeZ2 * 2));

    G4VSolid* sub10 = new G4Trd("HolderSub10", Holder::holeX1, Holder::holeX2, Holder::holeX1, Holder::holeX2,
                                Holder::holeZ - Holder::holeZ1 - Holder::holeZ2);
    G4VSolid* base14 = new G4SubtractionSolid("HolderBase14", base13, sub10, nullptr,
                                              G4ThreeVector(0, 0, -Holder::halfZ + 2 * (Holder::holeGapZ +
                                                                Holder::holeZ1) + Holder::holeZ - Holder::holeZ1 -
                                                            Holder::holeZ2));
    // Strip
    G4VSolid* stripBaseX = new G4Box("StripBaseX", Holder::stripX - Holder::stripRadius, Holder::recess2Y * 2,
                                     Holder::stripThickness);
    G4VSolid* stripBaseY = new G4Box("StripBaseY", Holder::stripThickness, Holder::recess2Y * 2,
                                     Holder::stripY - Holder::stripRadius);
    G4VSolid* rounding = new G4Tubs("StripRounding", 0, Holder::stripRadius, Holder::recess2Y * 2, 0 * deg, 100 * deg);
    G4VSolid* gap = new G4Box("StripGap", Holder::stripGap, Holder::recess2Y * 3, Holder::stripThickness * 2);

    G4VSolid* stripBase0 = new G4UnionSolid("StripBase0", stripBaseX, rounding, rotMat0,
                                            G4ThreeVector(Holder::stripX - Holder::stripRadius, 0,
                                                          -Holder::stripThickness));
    G4VSolid* stripBase1 = new G4UnionSolid("StripBase1", stripBase0, rounding, rotMat1,
                                            G4ThreeVector(-Holder::stripX + Holder::stripRadius, 0,
                                                          -Holder::stripThickness));
    G4VSolid* stripBase2 = new G4UnionSolid("StripBase2", stripBase1, stripBaseY, nullptr,
                                            G4ThreeVector(-Holder::stripX + Holder::stripThickness, 0,
                                                          -Holder::stripThickness - Holder::stripY +
                                                          Holder::stripRadius));
    G4VSolid* stripBase3 = new G4UnionSolid("StripBase3", stripBase2, stripBaseY, nullptr,
                                            G4ThreeVector(Holder::stripX - Holder::stripThickness, 0,
                                                          -Holder::stripThickness - Holder::stripY +
                                                          Holder::stripRadius));
    G4VSolid* stripBase4 = new G4UnionSolid("StripBase4", stripBase3, stripBase1, rotMat,
                                            G4ThreeVector(0, 0, 2 * (-Holder::stripY + Holder::stripThickness)));
    G4VSolid* strip = new G4SubtractionSolid("Strip", stripBase4, gap, nullptr,
                                             G4ThreeVector(0, 0, 2 * (-Holder::stripY + Holder::stripThickness)));

    G4VSolid* holder = new G4UnionSolid("SolarPanelHolder", base14, strip, nullptr,
                                        G4ThreeVector(0, Holder::halfY - 2 * (Holder::recess1Deep + Holder::recess2Y),
                                                      Holder::stripY - Holder::stripThickness));

    solarPanelHolderLV = new G4LogicalVolume(holder, holderMat, "SolarPanelHolderLV");
    solarPanelHolderLV->SetVisAttributes(visHolder);
}

void CubeSat_GeoScan_3U::ConstructSolarPanel() {
    // Panel
    G4VSolid* panelBaseX = new G4Box("PanelBaseX", Holder::halfX - (Holder::ledgeThickness + SolarPanel::gap) * 2,
                                     SolarPanel::thickness,
                                     Holder::halfZ - (Holder::recess1Ledge + SolarPanel::gap) * 2);
    G4VSolid* panelBaseY = new G4Box("PanelBaseY", Holder::ledgeX - (Holder::ledgeThickness + SolarPanel::gap) * 2,
                                     SolarPanel::thickness,
                                     Holder::halfZ - (Holder::recess1Ledge + SolarPanel::gap - Holder::ledgeZ) * 2);

    G4VSolid* panel = new G4UnionSolid("Panel", panelBaseX, panelBaseY, nullptr, G4ThreeVector(0, 0, 0));
    auto* panelLV = new G4LogicalVolume(panel, solarPanelMat, "PanelLV");
    panelLV->SetVisAttributes(visSolarPanel);


    // Union
    solarPanel = new G4Box("SolarPanelX", Holder::halfX, Holder::halfY, Holder::halfZ + Holder::ledgeZ * 2);
    solarPanelLV = new G4LogicalVolume(solarPanel, vacuumMat, "SolarPanelLV");

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), solarPanelHolderLV, "HolderPVP", solarPanelLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, Holder::halfY - SolarPanel::thickness - SolarPanel::deep * 2, 0),
                      panelLV, "PanelPVP", solarPanelLV, false, 0, true);
    solarPanelLV->SetVisAttributes(G4VisAttributes::GetInvisible());
}

void CubeSat_GeoScan_3U::ConstructMechanics() {
    G4VSolid* ssBase = new G4Box("SSBase", Mechanics::halfX, Mechanics::halfY, Mechanics::halfZ);
    G4VSolid* ssLedge = new G4Box("SSLedge", Mechanics::ledgeX, Mechanics::ledgeY, Mechanics::ledgeZ * 2);

    G4VSolid* sub0 = new G4Box("SSSub0", Mechanics::cut1X * 2, Mechanics::cut1Y * 2, Mechanics::cut1Z * 2);
    G4VSolid* sub1 = new G4Box("SSSub1", Mechanics::cut2X * 2, Mechanics::cut2Y * 2, Mechanics::cut1Z * 2);

    G4VSolid* base0 = new G4SubtractionSolid("SSBase0", ssBase, sub0, nullptr,
                                             G4ThreeVector(Mechanics::halfX, -Mechanics::halfY,
                                                           -Mechanics::halfZ));
    G4VSolid* base1 = new G4SubtractionSolid("SSBase1", base0, sub0, nullptr,
                                             G4ThreeVector(-Mechanics::halfX, -Mechanics::halfY,
                                                           -Mechanics::halfZ));

    G4VSolid* base2 = new G4SubtractionSolid("SSBase2", base1, sub1, nullptr,
                                             G4ThreeVector(Mechanics::halfX, Mechanics::halfY,
                                                           -Mechanics::halfZ));
    G4VSolid* base3 = new G4SubtractionSolid("SSBase3", base2, sub1, nullptr,
                                             G4ThreeVector(-Mechanics::halfX, Mechanics::halfY,
                                                           -Mechanics::halfZ));

    G4VSolid* sub2 = new G4Box("SSSub2", Mechanics::halfX * 2, Mechanics::channalWidth,
                               Mechanics::cut1Z * 2);
    G4VSolid* base4 = new G4SubtractionSolid("SSBase4", base3, sub2, nullptr,
                                             G4ThreeVector(0, Mechanics::halfY - Mechanics::cut2Y * 2 +
                                                           Mechanics::channalWidth, -Mechanics::halfZ));

    G4VSolid* sub3 = new G4Box("SSSub3", Mechanics::halfX * 2,
                               Mechanics::cut2Y * 2 - Mechanics::channalWidth, Mechanics::cut2Z * 2);
    G4VSolid* base5 = new G4SubtractionSolid("SSBase5", base4, sub3, nullptr,
                                             G4ThreeVector(0, Mechanics::halfY, -Mechanics::halfZ));

    mechanics = new G4UnionSolid("Mechanics", base5, ssLedge, nullptr,
                                 G4ThreeVector(0, 0, Mechanics::halfZ));

    mechanicsLV = new G4LogicalVolume(mechanics, mechanicsMat, "MechanicsLV");
    mechanicsLV->SetVisAttributes(visMechanics);
}


void CubeSat_GeoScan_3U::ConstructCover() {
    G4VSolid* coverBase = new G4Box("CoverBase", Cover::halfX, Cover::halfY, Cover::halfZ);

    G4VSolid* sub0 = new G4Box("CoverSub0", Cover::ledgeX * 2, Cover::ledgeY * 2, Cover::halfZ);
    G4VSolid* coverBase0 = new G4SubtractionSolid("CoverBase0", coverBase, sub0, nullptr,
                                                  G4ThreeVector(Cover::halfX, Cover::halfY, Cover::stripThickness * 2));
    G4VSolid* coverBase1 = new G4SubtractionSolid("CoverBase1", coverBase0, sub0, nullptr,
                                                  G4ThreeVector(-Cover::halfX, Cover::halfY,
                                                                Cover::stripThickness * 2));
    G4VSolid* coverBase2 = new G4SubtractionSolid("CoverBase2", coverBase1, sub0, nullptr,
                                                  G4ThreeVector(-Cover::halfX, -Cover::halfY,
                                                                Cover::stripThickness * 2));
    G4VSolid* coverBase3 = new G4SubtractionSolid("CoverBase3", coverBase2, sub0, nullptr,
                                                  G4ThreeVector(Cover::halfX, -Cover::halfY,
                                                                Cover::stripThickness * 2));

    G4VSolid* sub1 = new G4Box("CoverSub1", Cover::holeX, Cover::holeY, Cover::halfZ * 2);
    G4VSolid* coverBase4 = new G4SubtractionSolid("CoverBase4", coverBase3, sub1, nullptr,
                                                  G4ThreeVector(Cover::halfX - Cover::holeX - Cover::stripThickness *
                                                                2, 0, 0));
    G4VSolid* coverBase5 = new G4SubtractionSolid("CoverBase5", coverBase4, sub1, nullptr,
                                                  G4ThreeVector(-(Cover::halfX - Cover::holeX - Cover::stripThickness *
                                                                    2), 0, 0));
    G4VSolid* coverBase6 = new G4SubtractionSolid("CoverBase6", coverBase5, sub1, nullptr,
                                                  G4ThreeVector(Cover::halfX - Cover::stripThickness * 2, 0,
                                                                Cover::stripHeigth * 2 + Cover::halfZ));
    G4VSolid* coverBase7 = new G4SubtractionSolid("CoverBase7", coverBase6, sub1, nullptr,
                                                  G4ThreeVector(-(Cover::halfX - Cover::stripThickness * 2), 0,
                                                                Cover::stripHeigth * 2 + Cover::halfZ));

    G4VSolid* sub2 = new G4Box("CoverSub2", Cover::bottomRecessLength, Cover::bottomRecessWidth, Cover::halfZ);
    G4VSolid* coverBase8 = new G4SubtractionSolid("CoverBase8", coverBase7, sub2, nullptr,
                                                  G4ThreeVector(0, 0, 2 * (-Cover::halfZ + Cover::bottomRecessDeep)));

    G4VSolid* sub3 = new G4Box("CoverSub3", Cover::cut1X, Cover::cut1Y * 2, Cover::halfZ);
    G4VSolid* coverBase9 = new G4SubtractionSolid("CoverBase9", coverBase8, sub3, nullptr,
                                                  G4ThreeVector(0, -Cover::halfY,
                                                                2 * (-Cover::halfZ + Cover::bottomRecessDeep +
                                                                    Cover::cut1Z)));
    G4VSolid* coverBase10 = new G4SubtractionSolid("CoverBase10", coverBase9, sub3, nullptr,
                                                   G4ThreeVector(0, Cover::halfY,
                                                                 2 * (-Cover::halfZ + Cover::bottomRecessDeep)));

    G4VSolid* sub4 = new G4Box("CoverSub4", Cover::cut2X, Cover::cut2Y, Cover::halfZ);
    cover = new G4SubtractionSolid("Cover", coverBase10, sub4, nullptr,
                                   G4ThreeVector(0, 0, 2 * (-Cover::halfZ + Cover::bottomRecessDeep + Cover::cut2Z)));

    coverLV = new G4LogicalVolume(cover, coverMat, "CoverLV");
    coverLV->SetVisAttributes(visCover);
}


void CubeSat_GeoScan_3U::ConstructServiceSystem() {
    serviceSystem = new G4Box("ServiceSystem", Cover::halfX, Cover::halfY, Cover::halfZ + Boards::half1Z);
    serviceSystemLV = new G4LogicalVolume(serviceSystem, vacuumMat, "ServiceSystemLV");
    serviceSystemLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    G4VSolid* board1 = new G4Box("Board1", Boards::half1X, Boards::half1Y, Boards::half1Z);
    G4VSolid* board2 = new G4Box("Board2", Boards::half2X, Boards::half2Y, Boards::half2Z);
    auto* board1LV = new G4LogicalVolume(board1, boardMat, "Board1LV");
    auto* board2LV = new G4LogicalVolume(board2, boardMat, "Board2LV");
    board1LV->SetVisAttributes(visBoard);
    board2LV->SetVisAttributes(visBoard);

    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, Boards::half1Z), coverLV, "CoverPVP", serviceSystemLV, false, 0,
                      true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -Cover::halfZ), board1LV, "Board1PVP",
                      serviceSystemLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(2 * (Boards::gap1X - Boards::gap2X), 0, -Cover::halfZ + Boards::half1Z +
                                             Cover::bottomRecessDeep * 2 - Boards::gapZ * 2 + Boards::half2Z), board2LV,
                      "Board2PVP", serviceSystemLV, false, 0, true);
}

void CubeSat_GeoScan_3U::ConstructFrame() {
    auto* rotMat = new G4RotationMatrix(0, 180 * deg, 180 * deg);
    auto* rotMatY = new G4RotationMatrix(0, 180 * deg, 0);
    auto* rotMatZ = new G4RotationMatrix(0, 0, 180 * deg);
    auto* rotMat90Z = new G4RotationMatrix(0, 0, 90 * deg);
    auto* rotMat270Z = new G4RotationMatrix(0, 0, 270 * deg);

    // Frame
    new G4PVPlacement(nullptr, G4ThreeVector(0, -Unit::halfY + UnitFrame::halfY, -displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 0, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, Unit::halfY - UnitFrame::halfY, -displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 1, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, -Unit::halfY + UnitFrame::halfY,
                                             2 * (UnitFrame::halfZ + Frame::gap) - displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 2, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, Unit::halfY - UnitFrame::halfY,
                                             2 * (UnitFrame::halfZ + Frame::gap) - displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 3, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, -Unit::halfY + UnitFrame::halfY,
                                             -2 * (UnitFrame::halfZ + Frame::gap) - displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 4, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, Unit::halfY - UnitFrame::halfY,
                                             -2 * (UnitFrame::halfZ + Frame::gap) - displacement), unitFrameLV,
                      "UnitFramePVP", cubeSatLV, false, 5, true);

    // Bridges
    new G4PVPlacement(nullptr, G4ThreeVector(0, -Unit::halfY + UnitFrame::halfY,
                                             UnitFrame::halfZ + Frame::gap - displacement), bridgeLV, "BridgePVP",
                      cubeSatLV, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, -Unit::halfY + UnitFrame::halfY,
                                             -(UnitFrame::halfZ + Frame::gap) - displacement), bridgeLV, "BridgePVP",
                      cubeSatLV, false, 1, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, Unit::halfY - UnitFrame::halfY,
                                             -(UnitFrame::halfZ + Frame::gap) - displacement), bridgeLV, "BridgePVP",
                      cubeSatLV, false, 2, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, Unit::halfY - UnitFrame::halfY,
                                             UnitFrame::halfZ + Frame::gap - displacement), bridgeLV, "BridgePVP",
                      cubeSatLV, false, 3, true);

    // Brackets
    G4VSolid* bracket = new G4Box("Bracket", UnitFrame::halfX, UnitFrame::halfY, Frame::bracketHeight);
    auto* bracketLV = new G4LogicalVolume(bracket, frameMat, "BracketLV");
    bracketLV->SetVisAttributes(visFrame);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(Unit::halfX - UnitFrame::halfX, Unit::halfY - UnitFrame::halfY,
                                    Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight - displacement), bracketLV,
                      "BracketPVP", cubeSatLV, false, 0, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(Unit::halfX - UnitFrame::halfX, -(Unit::halfY - UnitFrame::halfY),
                                    Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight - displacement), bracketLV,
                      "BracketPVP", cubeSatLV, false, 1, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(-(Unit::halfX - UnitFrame::halfX), Unit::halfY - UnitFrame::halfY,
                                    Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight - displacement), bracketLV,
                      "BracketPVP", cubeSatLV, false, 2, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(-(Unit::halfX - UnitFrame::halfX), -(Unit::halfY - UnitFrame::halfY),
                                    Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight - displacement), bracketLV,
                      "BracketPVP", cubeSatLV, false, 3, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(Unit::halfX - UnitFrame::halfX, Unit::halfY - UnitFrame::halfY,
                                    -(Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight) - displacement),
                      bracketLV, "BracketPVP", cubeSatLV, false, 0, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(Unit::halfX - UnitFrame::halfX, -(Unit::halfY - UnitFrame::halfY),
                                    -(Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight) - displacement),
                      bracketLV, "BracketPVP", cubeSatLV, false, 1, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(-(Unit::halfX - UnitFrame::halfX), Unit::halfY - UnitFrame::halfY,
                                    -(Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight) - displacement),
                      bracketLV, "BracketPVP", cubeSatLV, false, 2, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(-(Unit::halfX - UnitFrame::halfX), -(Unit::halfY - UnitFrame::halfY),
                                    -(Unit::halfZ * 3 + Frame::gap * 2 + Frame::bracketHeight) - displacement),
                      bracketLV, "BracketPVP", cubeSatLV, false, 3, true);


    // Lintel
    new G4PVPlacement(nullptr,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, Lintel::width - UnitFrame::halfZ - displacement), lintelLV, "LintelPVP",
                      cubeSatLV, false, 0, true);
    new G4PVPlacement(rotMatZ,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, Lintel::width - UnitFrame::halfZ - displacement), lintelLV, "LintelPVP",
                      cubeSatLV, false, 1, true);
    new G4PVPlacement(rotMat,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, -(Lintel::width - UnitFrame::halfZ) - displacement), lintelLV, "LintelPVP",
                      cubeSatLV, false, 2, true);
    new G4PVPlacement(rotMatY,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, -(Lintel::width - UnitFrame::halfZ) - displacement), lintelLV, "LintelPVP",
                      cubeSatLV, false, 3, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, Lintel::width - UnitFrame::halfZ + (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 4, true);
    new G4PVPlacement(rotMatZ,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, Lintel::width - UnitFrame::halfZ + (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 5, true);
    new G4PVPlacement(rotMat,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, -(Lintel::width - UnitFrame::halfZ) + (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 6, true);
    new G4PVPlacement(rotMatY,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, -(Lintel::width - UnitFrame::halfZ) + (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 7, true);
    new G4PVPlacement(nullptr,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, Lintel::width - UnitFrame::halfZ - (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 8, true);
    new G4PVPlacement(rotMatZ,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, Lintel::width - UnitFrame::halfZ - (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 9, true);
    new G4PVPlacement(rotMat,
                      G4ThreeVector(-(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX),
                                    0, -(Lintel::width - UnitFrame::halfZ) - (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 10, true);
    new G4PVPlacement(rotMatY,
                      G4ThreeVector(2 * Bridge::fixLedgeX - (2 * Lintel::bracketX + Lintel::thickness) - Unit::halfX,
                                    0, -(Lintel::width - UnitFrame::halfZ) - (Unit::halfZ + Frame::gap) * 2 -
                                    displacement), lintelLV, "LintelPVP", cubeSatLV, false, 11, true);

    // Fixators
    new G4PVPlacement(nullptr,
                      G4ThreeVector(Unit::halfX - UnitFrame::thickness * 2 - Fixator::length,
                                    Unit::halfY - UnitFrame::thickness - UnitFrame::halfY,
                                    Unit::halfZ * 3 + Frame::gap * 2 - Fixator::width - displacement), fixatorLV,
                      "FixatorPVP", cubeSatLV, false, 0, true);
    new G4PVPlacement(rotMatZ,
                      G4ThreeVector(-(Unit::halfX - UnitFrame::thickness * 2 - Fixator::length),
                                    -(Unit::halfY - UnitFrame::thickness - UnitFrame::halfY),
                                    Unit::halfZ * 3 + Frame::gap * 2 - Fixator::width - displacement), fixatorLV,
                      "FixatorPVP", cubeSatLV, false, 1, true);

    // SolarPanels
    new G4PVPlacement(nullptr, G4ThreeVector(0, Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                             Holder::halfY, -displacement), solarPanelLV, "SolarPanelPVP", cubeSatLV,
                      false, 0, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, -(Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                 Holder::halfY), -displacement), solarPanelLV, "SolarPanelPVP",
                      cubeSatLV, false, 1, true);
    new G4PVPlacement(rotMat270Z, G4ThreeVector(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                Holder::halfY, 0, -displacement), solarPanelLV, "SolarPanelPVP", cubeSatLV,
                      false, 2, true);
    new G4PVPlacement(rotMat90Z, G4ThreeVector(-(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                   Holder::halfY), 0, -displacement), solarPanelLV, "SolarPanelPVP",
                      cubeSatLV, false, 3, true);

    new G4PVPlacement(nullptr, G4ThreeVector(0, Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                             Holder::halfY, (Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 4, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, -(Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                 Holder::halfY), (Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 5, true);
    new G4PVPlacement(rotMat270Z, G4ThreeVector(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                Holder::halfY, 0, (Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 6, true);
    new G4PVPlacement(rotMat90Z, G4ThreeVector(-(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                   Holder::halfY), 0, (Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 7, true);

    new G4PVPlacement(nullptr, G4ThreeVector(0, Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                             Holder::halfY, -(Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 8, true);
    new G4PVPlacement(rotMatZ, G4ThreeVector(0, -(Unit::halfY - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                 Holder::halfY), -(Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 9, true);
    new G4PVPlacement(rotMat270Z, G4ThreeVector(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                Holder::halfY, 0, -(Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 10, true);
    new G4PVPlacement(rotMat90Z, G4ThreeVector(-(Unit::halfX - 2 * (UnitFrame::thickness - UnitFrame::ledgeThickness) +
                                                   Holder::halfY), 0, -(Unit::halfZ + Frame::gap) * 2 - displacement),
                      solarPanelLV, "SolarPanelPVP", cubeSatLV, false, 11, true);

    // Mechanics
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -Unit::halfZ - (Frame::gap + Bridge::width) * 2 - Mechanics::halfZ -
                                    displacement), mechanicsLV, "MechanicsPVP", cubeSatLV, false, 0, true);

    // ServiceSystem
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, -Unit::halfZ - (Frame::gap + Bridge::width) * 2 + Mechanics::serviceSystemGap
                                    * 2 - Cover::halfZ - Boards::half1Z - displacement), serviceSystemLV,
                      "ServiceSystemPVP", cubeSatLV, false, 0, true);


    // G4VPhysicalVolume* cubeSatPVP = new G4PVPlacement(nullptr, G4ThreeVector(), cubeSatLV, "CubeSatPVP", worldLV, false,
    //                                                   0, true);
    // G4GDMLParser parser;
    // parser.Write("CubeSat_Geometry.gdml", cubeSatPVP);
}
