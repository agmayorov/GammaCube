#include "Detector.hh"


using namespace Sizes;

Detector::Detector(G4LogicalVolume *detContLV, const G4ThreeVector &detContSize, G4NistManager *nistMan, G4double vDeg,
                   const G4String &detType) {
    detectorType = detType;
    detContainerLV = detContLV;
    detContainerSize = detContSize;

    nist = nistMan;
    viewDeg = vDeg;

    additionalLength = 0 * mm;

    crystalSize = G4ThreeVector(0,
                                modelRadius - tunaCanThickWall - gapSizeWall - tyvekOutThickWall - vetoThickWall -
                                tyvekMidThickWall - AlThickWall - AlCapThickWall - tyvekInThickWall,
                                modelHeight - tunaCanThickTop - tunaCanThickBottom - gapSizeTop - gapSizeBottom -
                                tyvekOutThickTop - tyvekOutThickBottom - vetoThickTop - vetoThickBottom - rubberHeight -
                                tyvekMidThickTop - tyvekMidThickBottom - AlThickTop - AlCapThickBottom -
                                tyvekInThickTop - boardSpace);

    zCorrection = detContainerSize.z() - (crystalSize.z() / 2 + tyvekInThickTop + AlThickTop + rubberHeight +
                                          tyvekMidThickTop + vetoThickTop + tyvekOutThickTop + gapSizeTop);

    DefineMaterials();
    DefineVisual();
}

inline void AddProp(G4MaterialPropertiesTable *mpt,
                    const char *key,
                    const G4double *e, const G4double *v, size_t n,
                    bool spline = true,
                    bool createNewKey = true) {
#if G4VERSION_NUMBER >= 1100
    mpt->AddProperty(key,
                     std::vector<G4double>(e, e + n),
                     std::vector<G4double>(v, v + n),
                     createNewKey,
                     spline);
#else
    mpt->AddProperty(key, e, v, (int) n)->SetSpline(spline);
#endif
}


void Detector::DefineMaterials() {
    auto *elH = nist->FindOrBuildElement("H");
    auto *elC = nist->FindOrBuildElement("C");
    auto *elNa = nist->FindOrBuildElement("Na");
    auto *elI = nist->FindOrBuildElement("I");
    auto *elTl = nist->FindOrBuildElement("Tl");

    if (detectorType == "NaI") {
        const G4double rhoCrystal = 3.67 * g / cm3;
        CrystalMat = new G4Material("CrystalMat", rhoCrystal, 3, kStateSolid);

        const G4double wTl = 1.0e-3;
        const G4double wNa_noTl = 0.153;
        const G4double wI_noTl = 0.847;
        const G4double scale = 1.0 - wTl;

        CrystalMat->AddElement(elNa, wNa_noTl * scale);
        CrystalMat->AddElement(elI, wI_noTl * scale);
        CrystalMat->AddElement(elTl, wTl);

        const G4int n = 8;
        G4double eph[n] = {
            2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV
        };

        G4double rindex[n];
        G4double abslen[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.85;
            abslen[i] = 200. * cm;
        }

        G4double emitCrystal[n] = {0.02, 0.10, 0.50, 1.00, 0.85, 0.45, 0.12, 0.02};

        auto *mptCrystal = new G4MaterialPropertiesTable();
        AddProp(mptCrystal, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptCrystal, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptCrystal, "SCINTILLATIONCOMPONENT1", eph, emitCrystal, n, true, true);

        mptCrystal->AddConstProperty("SCINTILLATIONYIELD", 38000. / MeV);
        mptCrystal->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptCrystal->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 230. * ns);

        CrystalMat->SetMaterialPropertiesTable(mptCrystal);

        CrystalMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
    } else if (detectorType == "CsI") {
        CrystalMat = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");

        const G4int n = 8;
        G4double eph[n] = {
            2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV,
            2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV
        };

        G4double rindex[n];
        G4double abslen[n];
        for (int i = 0; i < n; i++) {
            rindex[i] = 1.80;
            abslen[i] = 200. * cm;
        }

        G4double emitCsI[n] = {0.01, 0.05, 0.40, 1.00, 0.80, 0.30, 0.05, 0.0};

        auto *mptCsI = new G4MaterialPropertiesTable();
        AddProp(mptCsI, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptCsI, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptCsI, "SCINTILLATIONCOMPONENT1", eph, emitCsI, n, true, true);

        mptCsI->AddConstProperty("SCINTILLATIONYIELD", 54000. / MeV);
        mptCsI->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptCsI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1000. * ns);

        CrystalMat->SetMaterialPropertiesTable(mptCsI);
        CrystalMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
    }

    {
        vetoMat = new G4Material("VetoMat", 1.032 * g / cm3, 2, kStateSolid);
        vetoMat->AddElement(elC, 9);
        vetoMat->AddElement(elH, 10);

        const G4int n = 8;
        G4double eph[n] = {2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV};

        G4double rindex[n];
        G4double abslen[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.58;
            abslen[i] = 380. * cm;
        }

        G4double emitPVT[n] = {0.01, 0.20, 0.60, 1.00, 0.65, 0.25, 0.05, 0.01};

        auto *mptVeto = new G4MaterialPropertiesTable();
        AddProp(mptVeto, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptVeto, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptVeto, "SCINTILLATIONCOMPONENT1", eph, emitPVT, n, true, true);

        mptVeto->AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
        mptVeto->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptVeto->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns);

        vetoMat->SetMaterialPropertiesTable(mptVeto);

        vetoMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }

    {
        const G4double rhoTyvek = 0.38 * g / cm3;

        tyvekInMat = new G4Material("tyvekInMat", rhoTyvek, 2, kStateSolid);
        tyvekOutMat = new G4Material("tyvekOutMat", rhoTyvek, 2, kStateSolid);
        tyvekInMat->AddElement(elC, 1);
        tyvekInMat->AddElement(elH, 2);
        tyvekOutMat->AddElement(elC, 1);
        tyvekOutMat->AddElement(elH, 2);

        const G4int n = 8;
        G4double eph[n] = {2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV};

        G4double rindex[n];
        G4double abslen[n];
        G4double rayl[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.50;
            abslen[i] = 1000. * m;
            rayl[i] = 10. * um;
        }

        auto *mptTyvek = new G4MaterialPropertiesTable();
        AddProp(mptTyvek, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptTyvek, "ABSLENGTH", eph, abslen, n, false, true);
        AddProp(mptTyvek, "RAYLEIGH", eph, rayl, n, false, true);

        tyvekInMat->SetMaterialPropertiesTable(mptTyvek);

        auto *mptTyvekOut = new G4MaterialPropertiesTable(*mptTyvek);
        tyvekOutMat->SetMaterialPropertiesTable(mptTyvekOut);
    }

    tyvekMidMat = tyvekOutMat;

    LEDMat = nist->FindOrBuildMaterial("G4_Galactic");

    AlMat = nist->FindOrBuildMaterial("G4_Al");

    rubberMat = nist->FindOrBuildMaterial("G4_RUBBER_NEOPRENE");

    wireMat = nist->FindOrBuildMaterial("G4_Cu");

    boardMat = nist->FindOrBuildMaterial("G4_Si");
}


G4String Detector::GetDetectorType() const {
    return detectorType;
}


void Detector::DefineVisual() {
    visCrystal = new G4VisAttributes(G4Color(1.0, 1.0, 1.0));
    visCrystal->SetForceSolid(true);
    visTyvekOut = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    visTyvekOut->SetForceSolid(true);
    visTyvekMid = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    visTyvekMid->SetForceSolid(true);
    visTyvekIn = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    visTyvekIn->SetForceSolid(true);
    visVeto = new G4VisAttributes(G4Color(1.0, 0.0, 0.0));
    visVeto->SetForceSolid(true);
    visLED = new G4VisAttributes(G4Color(0.0, 1.0, 1.0));
    visLED->SetForceSolid(true);
    visAl = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visAl->SetForceSolid(true);
    visWire = new G4VisAttributes(G4Color(0.508, 0.301, 0.171));
    visWire->SetForceSolid(true);
    visBoard = new G4VisAttributes(G4Color(0.121, 0.765, 0.175));
    visBoard->SetForceSolid(true);
    visRubber = new G4VisAttributes(G4Color(0., 0., 0.));
    visRubber->SetForceSolid(true);
}


void Detector::Construct() {
    ConstructCrystal();
    ConstructAl();
    ConstructVeto();
}


void Detector::ConstructCrystal() {
    G4ThreeVector pos = G4ThreeVector(0, 0, zCorrection);

    // Construct Crystal

    G4Tubs *Crystal = new G4Tubs("Crystal", crystalSize.x(), crystalSize.y(), crystalSize.z() / 2., 0, viewDeg);
    crystalLV = new G4LogicalVolume(Crystal, CrystalMat, "CrystalLV");
    new G4PVPlacement(nullptr, pos, crystalLV, "CrystalPVPL", detContainerLV, false, 0, true);
    crystalLV->SetVisAttributes(visCrystal);

    ConstructCrystalLED();

    // Construct Board
    G4Box *board = new G4Box("Board", boardLength / 2., boardWidth / 2., boardHeight / 2.);
    G4LogicalVolume *boardLV = new G4LogicalVolume(board, boardMat, "BoardLV");
    G4ThreeVector boardPos = G4ThreeVector(0, 0, -(crystalSize.z() - boardHeight) / 2. - boardSpace + zCorrection);
    new G4PVPlacement(nullptr, boardPos, boardLV, "BoardPVPL", detContainerLV, false, 0, true);
    boardLV->SetVisAttributes(visBoard);

    // Construct TyvekIn
    tyvekInSize = G4ThreeVector(crystalSize.y(),
                                crystalSize.y() + tyvekInThickWall,
                                crystalSize.z() + tyvekInThickTop + tyvekInThickBottom);


    if (tyvekInThickWall <= 0. and tyvekInThickTop <= 0. and tyvekInThickBottom <= 0.) return;

    G4VSolid *tyvekInIncomplete = nullptr;
    G4VSolid *tyvekInWall = nullptr;
    if (tyvekInThickWall > 0) {
        tyvekInWall = new G4Tubs("TyvekInWall", tyvekInSize.x(), tyvekInSize.y(), tyvekInSize.z() / 2, 0, viewDeg);
    }

    G4VSolid *tyvekInTop = nullptr;
    const G4ThreeVector tyvekInTopPos = G4ThreeVector(0, 0, (crystalSize.z() + tyvekInThickTop) / 2.0);
    if (tyvekInThickTop > 0) {
        tyvekInTop = new G4Tubs("TyvekInTop", 0, tyvekInSize.x(), tyvekInThickTop / 2., 0, viewDeg);
    }
    if (tyvekInThickWall > 0 and tyvekInThickTop > 0) {
        tyvekInIncomplete = new G4UnionSolid("TyvekInIncomplete", tyvekInWall, tyvekInTop, nullptr,
                                             tyvekInTopPos);
    } else if (tyvekInThickWall > 0 and tyvekInThickTop <= 0) {
        tyvekInIncomplete = tyvekInWall;
    } else if (tyvekInThickWall <= 0 and tyvekInThickTop > 0) {
        tyvekInIncomplete = tyvekInTop;
    }

    G4VSolid *tyvekInBottom = nullptr;
    const G4ThreeVector tyvekInBottomPos = G4ThreeVector(0, 0, -(crystalSize.z() + tyvekInThickBottom) / 2.0);
    if (tyvekInThickBottom > 0) {
        if (crystalLEDCount > 1) {
            G4Tubs *tyvekInBottomIncomplete = new G4Tubs("TyvekInBottomIncomplete", 0, tyvekInSize.x(),
                                                         tyvekInThickBottom / 2., 0, viewDeg);
            tyvekInBottom = new G4SubtractionSolid("TyvekInBottom", tyvekInBottomIncomplete, elongatedLED);
        } else {
            tyvekInBottom = new G4Tubs("TyvekInBottom", crystalLEDRadius, tyvekInSize.x(),
                                       tyvekInThickBottom / 2., 0, viewDeg);
        }
    }

    G4VSolid *tyvekIn = nullptr;
    if (tyvekInIncomplete != nullptr and tyvekInThickBottom > 0) {
        tyvekIn = new G4UnionSolid("TyvekIn", tyvekInIncomplete, tyvekInBottom, nullptr, tyvekInBottomPos);
    } else if (tyvekInIncomplete == nullptr and tyvekInThickBottom > 0) {
        tyvekIn = tyvekInBottom;
    } else if (tyvekInIncomplete != nullptr and tyvekInThickBottom <= 0) {
        tyvekIn = tyvekInIncomplete;
    }

    tyvekInLV = new G4LogicalVolume(tyvekIn, tyvekInMat, "TyvekInLV");
    G4ThreeVector tyvekInPos = G4ThreeVector(0, 0, zCorrection);
    new G4PVPlacement(nullptr, tyvekInPos, tyvekInLV, "TyvekInPVPL", detContainerLV, false, 0, true);
    tyvekInLV->SetVisAttributes(visTyvekIn);

    // Crystal wire with insulation
    G4double wireLength = AlCapThickBottom + tyvekMidThickBottom + vetoThickBottom + tyvekOutThickBottom +
                          additionalLength;
    G4VSolid *wire = new G4Tubs("CrystalWire", 0, wireRadius, wireLength / 2, 0, 360 * deg);
    G4LogicalVolume *wireLV = new G4LogicalVolume(wire, wireMat, "CrystalWireLV");
    G4ThreeVector wirePos = G4ThreeVector(0, 0, -(crystalSize.z() + wireLength) / 2.0 - boardSpace + zCorrection);
    new G4PVPlacement(nullptr, wirePos, wireLV, "CrystalWirePVPL", detContainerLV, false, 0, true);
    wireLV->SetVisAttributes(visWire);

    G4VSolid *wireInsulation = new G4Tubs("CrystalWireInsulation", wireRadius, wireRadius + wireInsulationThick,
                                          wireLength / 2, 0, 360 * deg);
    G4LogicalVolume *wireInsulationLV = new G4LogicalVolume(wireInsulation, tyvekMidMat, "CrystalWireInsulationLV");
    new G4PVPlacement(nullptr, wirePos, wireInsulationLV, "CrystalWireInsulationPVPL", detContainerLV, false, 0, true);
    wireInsulationLV->SetVisAttributes(visTyvekMid);
}


void Detector::ConstructCrystalLED() {
    if (crystalLEDRadius <= 0) return;

    G4double pinLength = boardSpace - crystalLEDHeight - boardHeight;
    G4Tubs *pin = new G4Tubs("CrystalPin", 0, wireRadius, pinLength / 2., 0, 360 * deg);

    const G4double radiusLEDXY = crystalLEDRadius / 2;
    std::vector<G4ThreeVector> pinXY;
    if (pinCount == 1) {
        pinXY = {G4ThreeVector(0, 0, 0)};
    } else {
        const G4double dphi = 2. * pi / pinCount;
        pinXY.reserve(pinCount);
        for (int i = 0; i < pinCount; ++i) {
            const G4double a = i * dphi;
            pinXY.emplace_back(radiusLEDXY * std::cos(a), radiusLEDXY * std::sin(a), 0);
        }
    }

    G4MultiUnion *pins = new G4MultiUnion("PinLEDMulti");
    for (auto &p: pinXY) {
        G4Transform3D tr(G4RotationMatrix(), G4ThreeVector(p.x(), p.y(), 0.0));
        pins->AddNode(*pin, tr);
    }
    pins->Voxelize();


    const G4double radiusXY = crystalSize.y() / 2.;

    std::vector<G4ThreeVector> LEDXY;
    if (crystalLEDCount == 1) {
        LEDXY = {G4ThreeVector(0, 0, 0)};
    } else {
        const G4double dphi = 2. * pi / crystalLEDCount;
        LEDXY.reserve(crystalLEDCount);
        for (int i = 0; i < crystalLEDCount; ++i) {
            const G4double a = i * dphi;
            LEDXY.emplace_back(radiusXY * std::cos(a), radiusXY * std::sin(a), 0);
        }
    }
    if (LEDXY.empty()) return;

    G4Tubs *LED = new G4Tubs("LEDCyl", 0, crystalLEDRadius, crystalLEDHeight / 2., 0, 360 * deg);
    G4Tubs *eLED = new G4Tubs("ElongatedLEDCyl", 0, crystalLEDRadius, tyvekInThickBottom + 1 * mm, 0, 360 * deg);

    crystalLED = new G4MultiUnion("CrystalLEDMulti");
    G4MultiUnion *pinsAll = new G4MultiUnion("PinsLEDMulti");
    elongatedLED = new G4MultiUnion("ElongatedLEDCylMulti");
    for (auto &p: LEDXY) {
        G4Transform3D tr(G4RotationMatrix(), G4ThreeVector(p.x(), p.y(), 0.0));
        crystalLED->AddNode(*LED, tr);
        pinsAll->AddNode(*pins, tr);
        elongatedLED->AddNode(*eLED, tr);
    }
    crystalLED->Voxelize();
    pinsAll->Voxelize();
    elongatedLED->Voxelize();

    crystalLEDLV = new G4LogicalVolume(crystalLED, LEDMat, "LEDMultiLV");
    const G4double zOff = -(crystalSize.z() / 2. + crystalLEDHeight / 2.) + zCorrection;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOff), crystalLEDLV, "LEDMultiPVPL", detContainerLV, false, 0, true);
    crystalLEDLV->SetVisAttributes(visLED);

    G4LogicalVolume *pinsLV = new G4LogicalVolume(pinsAll, wireMat, "PinsLEDMultiLV");
    const G4double pinOffZ = -(crystalSize.z() + pinLength) / 2. - crystalLEDHeight + zCorrection;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, pinOffZ), pinsLV, "PinsLEDMultiPVPL", detContainerLV, false, 0,
                      true);
    pinsLV->SetVisAttributes(visWire);
}


void Detector::ConstructAl() {
    if (AlCapThickBottom <= 0 and AlCapThickWall <= 0 and AlThickTop <= 0 and AlThickWall <= 0) return;

    // Aluminium container for insulation
    G4ThreeVector AlContSize = G4ThreeVector(crystalSize.y() + tyvekInThickWall,
                                             crystalSize.y() + tyvekInThickWall + AlThickWall,
                                             crystalSize.z() + tyvekInThickTop + AlThickTop);

    G4VSolid *AlCont = nullptr;
    G4VSolid *AlContWall = nullptr;
    if (AlThickWall > 0) {
        AlContWall = new G4Tubs("AlContWall", AlContSize.x(), AlContSize.y(), AlContSize.z() / 2, 0, viewDeg);
    }

    G4VSolid *AlContTop = nullptr;
    if (AlThickTop > 0) {
        AlContTop = new G4Tubs("AlContTop", 0, AlContSize.x(), AlThickTop / 2, 0, viewDeg);
    }
    if (AlContWall != nullptr and AlContTop != nullptr) {
        AlCont = new G4UnionSolid("AlCont", AlContWall, AlContTop, nullptr,
                                  G4ThreeVector(0, 0, (AlContSize.z() - AlThickTop) / 2.));
    } else if (AlContWall == nullptr and AlContTop != nullptr) {
        AlCont = AlContTop;
    } else {
        AlCont = AlContWall;
    }

    // Aluminium container cap for insulation
    G4ThreeVector AlCapSize = G4ThreeVector(crystalSize.y() + tyvekInThickWall + AlThickWall,
                                            crystalSize.y() + tyvekInThickWall + AlThickWall + AlCapThickWall,
                                            boardSpace + AlCapThickBottom);

    G4VSolid *AlCap = nullptr;
    G4VSolid *AlCapWall = nullptr;
    if (AlCapThickWall > 0) {
        AlCapWall = new G4Tubs("AlCapWall", AlCapSize.x(), AlCapSize.y(), AlCapSize.z() / 2, 0, viewDeg);
    }

    G4VSolid *AlCapBottom = nullptr;
    if (AlCapThickBottom > 0) {
        AlCapBottom = new G4Tubs("AlCapBottom", wireRadius + wireInsulationThick, AlCapSize.x(), AlCapThickBottom / 2,
                                 0, viewDeg);
    }
    if (AlCapWall != nullptr and AlCapBottom != nullptr) {
        AlCap = new G4UnionSolid("AlCap", AlCapWall, AlCapBottom, nullptr,
                                 G4ThreeVector(0, 0, -(AlCapSize.z() - AlCapThickBottom) / 2.));
    } else if (AlCapWall == nullptr and AlCapBottom != nullptr) {
        AlCap = AlCapBottom;
    } else {
        AlCap = AlCapWall;
    }

    G4VSolid *Al = nullptr;
    if (AlCont != nullptr and AlCap != nullptr) {
        G4ThreeVector posAlCont = G4ThreeVector(0, 0, (AlContSize.z() + AlCapSize.z()) / 2.);
        Al = new G4UnionSolid("Al", AlCap, AlCont, nullptr, posAlCont);
    } else if (AlCont == nullptr and AlCap != nullptr) {
        Al = AlCap;
    } else {
        Al = AlCont;
    }

    AlLV = new G4LogicalVolume(Al, AlMat, "AlLV");
    G4ThreeVector posAl = G4ThreeVector(0, 0, -(crystalSize.z() + boardSpace + AlCapThickBottom) / 2. + zCorrection);
    new G4PVPlacement(nullptr, posAl, AlLV, "AlPVPL", detContainerLV, false, 0, true);
    AlLV->SetVisAttributes(visAl);

    // Rubber
    G4VSolid *rubber = new G4Tubs("Rubber", 0, rubberRadius, rubberHeight / 2., 0, viewDeg);
    rubberLV = new G4LogicalVolume(rubber, rubberMat, "RubberLV");
    G4ThreeVector posRubber = G4ThreeVector(0, 0,
                                            (crystalSize.z() + rubberHeight) / 2. + tyvekInThickTop + AlThickTop +
                                            zCorrection);
    new G4PVPlacement(nullptr, posRubber, rubberLV, "RubberPVPL", detContainerLV, false, 0, true);
    rubberLV->SetVisAttributes(visRubber);

    // TyvekMid
    tyvekMidSize = G4ThreeVector(AlCapSize.y(),
                                 AlCapSize.y() + tyvekMidThickWall,
                                 AlCapSize.z() + AlContSize.z() + rubberHeight + tyvekMidThickBottom +
                                 tyvekMidThickTop);

    G4VSolid *tyvekMidIncomplete = nullptr;
    G4VSolid *tyvekMidWall = nullptr;
    if (tyvekMidThickWall > 0) {
        tyvekMidWall = new G4Tubs("TyvekMidWall", tyvekMidSize.x(), tyvekMidSize.y(), tyvekMidSize.z() / 2, 0, viewDeg);
    }

    G4VSolid *tyvekMidTop = nullptr;
    const G4ThreeVector tyvekMidTopPos = G4ThreeVector(0, 0, (tyvekMidSize.z() - tyvekMidThickTop) / 2.0);
    if (tyvekMidThickTop > 0) {
        tyvekMidTop = new G4Tubs("TyvekMidTop", 0, tyvekMidSize.x(), tyvekMidThickTop / 2., 0, viewDeg);
    }
    if (tyvekMidThickWall > 0 and tyvekMidThickTop > 0) {
        tyvekMidIncomplete = new G4UnionSolid("TyvekMidIncomplete", tyvekMidWall, tyvekMidTop, nullptr,
                                              tyvekMidTopPos);
    } else if (tyvekMidThickWall > 0 and tyvekMidThickTop <= 0) {
        tyvekMidIncomplete = tyvekMidWall;
    } else if (tyvekMidThickWall <= 0 and tyvekMidThickTop > 0) {
        tyvekMidIncomplete = tyvekMidTop;
    }

    G4VSolid *tyvekMidBottom = nullptr;
    const G4ThreeVector tyvekMidBottomPos = G4ThreeVector(0, 0, (tyvekMidThickBottom - tyvekMidSize.z()) / 2.0);
    if (tyvekMidThickBottom > 0) {
        tyvekMidBottom = new G4Tubs("TyvekMidBottom", wireRadius + wireInsulationThick, tyvekMidSize.x(),
                                    tyvekMidThickBottom / 2., 0, viewDeg);
    }

    G4VSolid *tyvekMid = nullptr;
    if (tyvekMidIncomplete != nullptr and tyvekMidThickBottom > 0) {
        tyvekMid = new G4UnionSolid("TyvekMid", tyvekMidIncomplete, tyvekMidBottom, nullptr, tyvekMidBottomPos);
    } else if (tyvekMidIncomplete == nullptr and tyvekMidThickBottom > 0) {
        tyvekMid = tyvekMidBottom;
    } else if (tyvekMidIncomplete != nullptr and tyvekMidThickBottom <= 0) {
        tyvekMid = tyvekMidIncomplete;
    }

    tyvekMidLV = new G4LogicalVolume(tyvekMid, tyvekMidMat, "TyvekMidLV");
    G4ThreeVector tyvekMidPos = G4ThreeVector(0, 0,
                                              (crystalSize.z() - tyvekMidSize.z()) / 2.0 + tyvekInThickTop + AlThickTop
                                              + rubberHeight + tyvekMidThickTop + zCorrection);
    new G4PVPlacement(nullptr, tyvekMidPos, tyvekMidLV, "TyvekMidPVPL", detContainerLV, false, 0, true);
    tyvekMidLV->SetVisAttributes(visTyvekMid);
}


void Detector::ConstructVeto() {
    vetoSize = G4ThreeVector(tyvekMidSize.y(),
                             tyvekMidSize.y() + vetoThickWall,
                             tyvekMidSize.z());
    ConstructVetoLED();

    G4VSolid *vetoIncomplete = nullptr;
    G4VSolid *vetoWall = nullptr;
    if (vetoThickWall > 0) {
        vetoWall = new G4Tubs("VetoWall", vetoSize.x(), vetoSize.y(), vetoSize.z() / 2, 0, viewDeg);
    }

    G4VSolid *vetoTop = nullptr;
    G4ThreeVector vetoTopPos = G4ThreeVector(0, 0, (vetoSize.z() + vetoThickTop) / 2.0);
    if (vetoThickTop > 0) {
        if (vetoChamferHeigh > 0) {
            G4VSolid *vetoTopTube = new G4Tubs("VetoTopTube", 0, vetoSize.y(), (vetoThickTop - vetoChamferHeigh) / 2, 0,
                                               viewDeg);
            G4VSolid *vetoTopCons = new G4Cons("VetoTopCons", 0, vetoSize.y(), 0, vetoSize.y() - vetoChamferHeigh,
                                               vetoChamferHeigh / 2, 0, viewDeg);

            G4ThreeVector vetoTopConsPos = G4ThreeVector(0, 0, vetoThickTop / 2.);
            vetoTopPos = G4ThreeVector(0, 0, (vetoSize.z() + vetoThickTop - vetoChamferHeigh) / 2.0);
            vetoTop = new G4UnionSolid("VetoTop", vetoTopTube, vetoTopCons, nullptr, vetoTopConsPos);
        } else {
            vetoTop = new G4Tubs("VetoTop", 0, vetoSize.y(), vetoThickTop / 2., 0, viewDeg);
        }
    }
    if (vetoThickWall > 0 and vetoThickTop > 0) {
        vetoIncomplete = new G4UnionSolid("VetoIncomplete", vetoWall, vetoTop, nullptr, vetoTopPos);
    } else if (vetoThickWall > 0 and vetoThickTop <= 0) {
        vetoIncomplete = vetoWall;
    } else if (vetoThickWall <= 0 and vetoThickTop > 0) {
        vetoIncomplete = vetoTop;
    }

    G4VSolid *vetoBottom = nullptr;
    const G4ThreeVector vetoBottomPos = G4ThreeVector(0, 0, -(vetoSize.z() + vetoThickBottom) / 2.0);
    if (vetoThickBottom > 0) {
        if (vetoLEDCount > 0) {
            G4VSolid *vetoBottomIncomplete = new G4Tubs("VetoBottomIncomplete", wireRadius + wireInsulationThick,
                                                        vetoSize.y(), vetoThickBottom / 2., 0, viewDeg);

            G4ThreeVector vetoLEDPos = G4ThreeVector(0, 0, vetoThickBottom - vetoLEDHeight / 2.0);

            G4ThreeVector vetoLEDWirePos = G4ThreeVector(0, 0, 0);
            G4VSolid *vetoBottomIncomplete1 = new G4SubtractionSolid("VetoBottomIncompleteWithWire",
                                                                     vetoBottomIncomplete, vetoWireFictive,
                                                                     nullptr, vetoLEDWirePos);

            vetoBottom = new G4SubtractionSolid("VetoBottom",
                                                vetoBottomIncomplete1, vetoLEDFictive, nullptr,
                                                vetoLEDPos);
        } else {
            vetoBottom = new G4Tubs("VetoBottom", wireRadius + wireInsulationThick, vetoSize.y(),
                                    vetoThickBottom / 2., 0, viewDeg);
        }
    }

    G4VSolid *veto = nullptr;
    if (vetoIncomplete != nullptr and vetoThickBottom > 0) {
        veto = new G4UnionSolid("Veto", vetoIncomplete, vetoBottom, nullptr, vetoBottomPos);
    } else if (vetoIncomplete == nullptr and vetoThickBottom > 0) {
        veto = vetoBottom;
    } else if (vetoIncomplete != nullptr and vetoThickBottom <= 0) {
        veto = vetoIncomplete;
    }

    vetoLV = new G4LogicalVolume(veto, vetoMat, "VetoLV");
    G4ThreeVector vetoPos = G4ThreeVector(0, 0,
                                          (crystalSize.z() - tyvekMidSize.z()) / 2.0 + tyvekInThickTop + AlThickTop
                                          + rubberHeight + tyvekMidThickTop + zCorrection);
    if (vetoThickWall <= 0) {
        vetoPos = G4ThreeVector(0, 0,
                                -(crystalSize.z() + vetoThickBottom) / 2.0 - tyvekInThickBottom - AlCapThickBottom -
                                tyvekMidThickBottom - boardSpace + zCorrection);
    }
    new G4PVPlacement(nullptr, vetoPos, vetoLV, "VetoPVPL", detContainerLV, false, 0, true);
    vetoLV->SetVisAttributes(visVeto);


    // TyvekOut
    tyvekOutSize = G4ThreeVector(vetoSize.y(),
                                 vetoSize.y() + tyvekOutThickWall,
                                 vetoSize.z() + vetoThickBottom + vetoThickTop);

    ConstructVetoBottomLED();

    G4VSolid *tyvekOutIncomplete = nullptr;
    G4VSolid *tyvekOutWall = nullptr;
    if (tyvekOutThickWall > 0) {
        if (vetoBottomLEDCount > 0) {
            G4VSolid *tyvekOutWallIncomplete = new G4Tubs("TyvekOutWallIncomplete", tyvekOutSize.x(), tyvekOutSize.y(),
                                                          (tyvekOutSize.z() - vetoChamferHeigh) / 2, 0, viewDeg);
            G4ThreeVector vetoBottomLEDPos = G4ThreeVector(
                0, 0, -(tyvekOutSize.z() - vetoChamferHeigh - vetoThickBottom) / 2.);
            tyvekOutWall = new G4SubtractionSolid("TyvekOutWall", tyvekOutWallIncomplete, vetoBottomLEDFictive, nullptr,
                                                  vetoBottomLEDPos);
        } else {
            tyvekOutWall = new G4Tubs("TyvekOutWall", tyvekOutSize.x(), tyvekOutSize.y(),
                                      (tyvekOutSize.z() - vetoChamferHeigh) / 2, 0, viewDeg);
        }
    }

    G4VSolid *tyvekOutTop = nullptr;
    G4ThreeVector tyvekOutTopPos = G4ThreeVector(0, 0, (tyvekOutSize.z() + tyvekOutThickTop) / 2.0);
    if (tyvekOutThickTop > 0) {
        if (vetoChamferHeigh > 0) {
            G4VSolid *tyvekOutTopConsOut = new G4Cons("TyvekOutTopConsOut", tyvekOutSize.x(), tyvekOutSize.y(),
                                                      tyvekOutSize.y() - 2 * tyvekOutThickWall - vetoChamferHeigh,
                                                      tyvekOutSize.y() - tyvekOutThickWall - vetoChamferHeigh,
                                                      (vetoChamferHeigh + tyvekOutThickTop) / 2, 0, viewDeg);
            G4VSolid *tyvekOutTopConsIn = new G4Cons("TyvekOutTopConsIn", 0, tyvekOutSize.x() - vetoChamferHeigh, 0,
                                                     tyvekOutSize.y() - 2 * tyvekOutThickWall - vetoChamferHeigh,
                                                     tyvekOutThickTop / 2, 0, viewDeg);

            tyvekOutTopPos = G4ThreeVector(0, 0, (tyvekOutSize.z() + tyvekOutThickWall) / 2.0);

            G4ThreeVector tyvekOutTopConsPos = G4ThreeVector(
                0, 0, (vetoChamferHeigh + tyvekOutThickWall - tyvekOutThickTop) / 2.);
            tyvekOutTop = new G4UnionSolid("TyvekOutTop", tyvekOutTopConsOut, tyvekOutTopConsIn, nullptr,
                                           tyvekOutTopConsPos);
        } else {
            tyvekOutTop = new G4Tubs("TyvekOutTop", 0, tyvekOutSize.y(), tyvekOutThickTop / 2., 0, viewDeg);
        }
    }
    if (tyvekOutThickWall > 0 and tyvekOutThickTop > 0) {
        tyvekOutIncomplete = new G4UnionSolid("TyvekOutIncomplete", tyvekOutWall, tyvekOutTop, nullptr, tyvekOutTopPos);
    } else if (tyvekOutThickWall > 0 and tyvekOutThickTop <= 0) {
        tyvekOutIncomplete = tyvekOutWall;
    } else if (tyvekOutThickWall <= 0 and tyvekOutThickTop > 0) {
        tyvekOutIncomplete = tyvekOutTop;
    }

    G4VSolid *tyvekOutBottom = nullptr;
    const G4ThreeVector tyvekOutBottomPos = G4ThreeVector(
        0, 0, -(tyvekOutSize.z() - vetoChamferHeigh + tyvekOutThickBottom) / 2.0);
    if (tyvekOutThickBottom > 0) {
        if (vetoLEDCount > 0) {
            G4VSolid *tyvekOutBottomIncomplete = new G4Tubs("TyvekOutBottomIncomplete",
                                                            wireRadius + wireInsulationThick,
                                                            tyvekOutSize.y(), tyvekOutThickBottom / 2., 0, viewDeg);

            G4ThreeVector tyvekOutLEDWirePos = G4ThreeVector(0, 0, 0);
            tyvekOutBottom = new G4SubtractionSolid("TyvekOutBottom",
                                                    tyvekOutBottomIncomplete, vetoWireFictive, nullptr,
                                                    tyvekOutLEDWirePos);
        } else {
            tyvekOutBottom = new G4Tubs("TyvekOutBottom", wireRadius + wireInsulationThick, tyvekOutSize.y(),
                                        tyvekOutThickBottom / 2., 0, viewDeg);
        }
    }

    G4VSolid *tyvekOut = nullptr;
    if (tyvekOutIncomplete != nullptr and tyvekOutThickBottom > 0) {
        tyvekOut = new G4UnionSolid("TyvekOut", tyvekOutIncomplete, tyvekOutBottom, nullptr, tyvekOutBottomPos);
    } else if (tyvekOutIncomplete == nullptr and tyvekOutThickBottom > 0) {
        tyvekOut = tyvekOutBottom;
    } else if (tyvekOutIncomplete != nullptr and tyvekOutThickBottom <= 0) {
        tyvekOut = tyvekOutIncomplete;
    }

    tyvekOutLV = new G4LogicalVolume(tyvekOut, tyvekOutMat, "TyvekOutLV");
    G4ThreeVector tyvekOutPos = G4ThreeVector(0, 0,
                                              -(vetoThickBottom + tyvekMidThickBottom +
                                                AlCapThickBottom + boardSpace - tyvekInThickTop - AlThickTop -
                                                rubberHeight - tyvekMidThickTop - vetoThickTop + vetoChamferHeigh) /
                                              2. + zCorrection);

    new G4PVPlacement(nullptr, tyvekOutPos, tyvekOutLV, "TyvekOutPVPL", detContainerLV, false, 0, true);
    tyvekOutLV->SetVisAttributes(visTyvekOut);
}


void Detector::ConstructVetoLED() {
    if (vetoLEDRadius <= 0) return;

    const G4double radiusXY = vetoSize.y() - vetoThickWall / 2.0;

    std::vector<G4ThreeVector> LEDXY;
    if (vetoLEDCount == 1) {
        LEDXY = {G4ThreeVector(0, 0, 0)};
    } else {
        const G4double dphi = 2. * pi / vetoLEDCount;
        LEDXY.reserve(vetoLEDCount);
        for (int i = 0; i < vetoLEDCount; ++i) {
            const G4double a = i * dphi;
            LEDXY.emplace_back(radiusXY * std::cos(a), radiusXY * std::sin(a), 0);
        }
    }
    if (LEDXY.empty()) return;

    G4Tubs *LED = new G4Tubs("vetoLEDCyl", 0, vetoLEDRadius, vetoLEDHeight / 2., 0, 360 * deg);
    G4Tubs *LEDFictive = new G4Tubs("LEDCylFictive", 0, vetoLEDRadius, vetoLEDHeight + vetoThickBottom, 0,
                                    360 * deg);

    G4double vetoWireLength = vetoThickBottom + tyvekOutThickBottom - vetoLEDHeight + additionalLength;
    G4Tubs *vetoWire = new G4Tubs("VetoWire", 0, wireRadius, vetoWireLength / 2, 0,
                                  360 * deg);
    G4Tubs *vetoWireInsulation = new G4Tubs("VetoWireInsulation", wireRadius, wireRadius + wireInsulationThick,
                                            vetoWireLength / 2, 0, 360 * deg);

    G4Tubs *wireFictive = new G4Tubs("VetoWireFictive", 0, wireRadius + wireInsulationThick,
                                     vetoWireLength / 2, 0, 360 * deg);

    vetoLED = new G4MultiUnion("VetoLEDMulti");
    vetoLEDWire = new G4MultiUnion("VetoLEDWireMulti");
    vetoLEDWireInsulation = new G4MultiUnion("VetoLEDWireInsulationMulti");
    vetoLEDFictive = new G4MultiUnion("VetoLEDFictiveMulti");
    vetoWireFictive = new G4MultiUnion("VetoWireFictiveMulti");
    for (auto &p: LEDXY) {
        G4Transform3D tr(G4RotationMatrix(), G4ThreeVector(p.x(), p.y(), 0.0));
        vetoLED->AddNode(*LED, tr);
        vetoLEDWire->AddNode(*vetoWire, tr);
        vetoLEDWireInsulation->AddNode(*vetoWireInsulation, tr);
        vetoLEDFictive->AddNode(*LEDFictive, tr);
        vetoWireFictive->AddNode(*wireFictive, tr);
    }
    vetoLED->Voxelize();
    vetoLEDWire->Voxelize();
    vetoLEDWireInsulation->Voxelize();
    vetoLEDFictive->Voxelize();
    vetoWireFictive->Voxelize();

    vetoLEDLV = new G4LogicalVolume(vetoLED, LEDMat, "VetoLEDMultiLV");
    const G4double zOff = -(crystalSize.z() / 2. + boardSpace + AlCapThickBottom + tyvekMidThickBottom + vetoLEDHeight /
                            2.) + zCorrection;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOff), vetoLEDLV, "VetoLEDMultiPVPL", detContainerLV, false, 0,
                      true);
    vetoLEDLV->SetVisAttributes(visLED);

    G4LogicalVolume *vetoLEDWireLV = new G4LogicalVolume(vetoLEDWire, wireMat, "VetoLEWireMultiLV");
    const G4double zOffWire = -(crystalSize.z() / 2. + boardSpace + AlCapThickBottom + tyvekMidThickBottom +
                                vetoLEDHeight + vetoWireLength / 2.) + zCorrection;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOffWire), vetoLEDWireLV, "VetoLEDWireMultiPVPL", detContainerLV,
                      false, 0, true);
    vetoLEDWireLV->SetVisAttributes(visWire);

    G4LogicalVolume *vetoLEDWireInsulationLV = new G4LogicalVolume(vetoLEDWireInsulation, tyvekMidMat,
                                                                   "VetoLEWireInsulationMultiLV");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOffWire), vetoLEDWireInsulationLV, "VetoLEDWireInsulationMultiPVPL",
                      detContainerLV, false, 0, true);
    vetoLEDWireInsulationLV->SetVisAttributes(visTyvekMid);
}

void Detector::ConstructVetoBottomLED() {
    if (vetoBottomLEDRadius <= 0) return;

    const G4double radiusXY = vetoSize.y() + vetoBottomLEDHeight / 2.0;

    std::vector<G4ThreeVector> LEDXY;
    std::vector<G4double> LEDAlpha;
    if (vetoBottomLEDCount == 1) {
        LEDXY = {G4ThreeVector(0, 0, 0)};
        LEDAlpha = {0};
    } else {
        const G4double dphi = 2. * pi / vetoBottomLEDCount;
        LEDXY.reserve(vetoBottomLEDCount);
        LEDAlpha.reserve(vetoBottomLEDCount);
        for (int i = 0; i < vetoBottomLEDCount; ++i) {
            const G4double a = (i + 0.5) * dphi;
            LEDXY.emplace_back(radiusXY * std::cos(a), radiusXY * std::sin(a), 0);
            LEDAlpha.emplace_back(a);
        }
    }
    if (LEDXY.empty()) return;

    G4Tubs *LED = new G4Tubs("VetoBottomLEDCyl", 0, vetoBottomLEDRadius, vetoBottomLEDHeight / 2., 0, 360 * deg);
    G4Tubs *LEDFictive = new G4Tubs("BottomLEDCylFictive", 0, vetoBottomLEDRadius,
                                    (tyvekOutThickWall + vetoBottomLEDHeight) * 2, 0, 360 * deg);

    vetoBottomLED = new G4MultiUnion("VetoBottomLEDMulti");
    vetoBottomLEDFictive = new G4MultiUnion("VetoBottomLEDFictiveMulti");
    for (int i = 0; i < vetoBottomLEDCount; ++i) {
        const G4ThreeVector &p = LEDXY[i];
        G4Transform3D tr(G4RotationMatrix(0, 90 * deg, 90 * deg - LEDAlpha[i]), G4ThreeVector(p.x(), p.y(), 0.0));
        vetoBottomLED->AddNode(*LED, tr);
        vetoBottomLEDFictive->AddNode(*LEDFictive, tr);
    }
    vetoBottomLED->Voxelize();
    vetoBottomLEDFictive->Voxelize();

    vetoBottomLEDLV = new G4LogicalVolume(vetoBottomLED, LEDMat, "VetoBottomLEDMultiLV");
    const G4double zOff = -(crystalSize.z() / 2. + boardSpace + AlCapThickBottom + tyvekMidThickBottom +
                            vetoThickBottom / 2.) + zCorrection;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zOff), vetoBottomLEDLV, "VetoBottomLEDMultiPVPL", detContainerLV,
                      false, 0, true);
    vetoBottomLEDLV->SetVisAttributes(visLED);
}

std::vector<G4LogicalVolume *> Detector::GetSensitiveLV() const {
    return {
        crystalLV,
        tyvekInLV,
        AlLV,
        rubberLV,
        tyvekMidLV,
        vetoLV,
        tyvekOutLV,
        crystalLEDLV,
        vetoLEDLV,
        vetoBottomLEDLV
    };
}
