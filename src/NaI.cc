#include "NaI.hh"


NaI::NaI(G4LogicalVolume *detContLV, const G4ThreeVector &detContSize, G4NistManager *nistMan, const Sizes &sizes,
         G4double vDeg, G4bool dLED) {
    detectorType = "NaI";
    detContainerLV = detContLV;
    detContainerSize = detContSize;
    nist = nistMan;
    vetoThick = sizes.vetoThick;
    tyvekThick = sizes.tyvekThick;
    shellThick = sizes.shellThick;
    gapSize = sizes.gapSize;
    LEDSize = sizes.LEDSize;
    viewDeg = vDeg;
    doubleLED = dLED;

    init();
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

void NaI::DefineMaterials() {
    auto *elH = nist->FindOrBuildElement("H");
    auto *elC = nist->FindOrBuildElement("C");
    auto *elNa = nist->FindOrBuildElement("Na");
    auto *elI = nist->FindOrBuildElement("I");
    auto *elTl = nist->FindOrBuildElement("Tl");

    {
        const G4double rhoNaI = 3.67 * g / cm3;
        NaIMat = new G4Material("NaIMat", rhoNaI, 3, kStateSolid);

        const G4double wTl = 1.0e-3;
        const G4double wNa_noTl = 0.153;
        const G4double wI_noTl = 0.847;
        const G4double scale = 1.0 - wTl;

        NaIMat->AddElement(elNa, wNa_noTl * scale);
        NaIMat->AddElement(elI, wI_noTl * scale);
        NaIMat->AddElement(elTl, wTl);

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

        G4double emitNaI[n] = {0.02, 0.10, 0.50, 1.00, 0.85, 0.45, 0.12, 0.02};

        auto *mptNaI = new G4MaterialPropertiesTable();
        AddProp(mptNaI, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptNaI, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptNaI, "SCINTILLATIONCOMPONENT1", eph, emitNaI, n, true, true);

        mptNaI->AddConstProperty("SCINTILLATIONYIELD", 38000. / MeV);
        mptNaI->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptNaI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 230. * ns);

        NaIMat->SetMaterialPropertiesTable(mptNaI);

        NaIMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
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

    LEDMat = nist->FindOrBuildMaterial("G4_Galactic");
    shellMat = nist->FindOrBuildMaterial("G4_Al");
}


void NaI::Construct() {
    G4ThreeVector pos = G4ThreeVector(0, 0, -gapSize / 2);

    NaISize = G4ThreeVector(0,
                            detContainerSize.y() - vetoThick - tyvekThick * 2. - gapSize - shellThick,
                            detContainerSize.z() - vetoThick - tyvekThick * 2 - gapSize / 2 - shellThick);

    G4Tubs *NaIBox = new G4Tubs("NaI", NaISize.x(), NaISize.y(), NaISize.z(), 0, viewDeg);
    NaILV = new G4LogicalVolume(NaIBox, NaIMat, "NaILV");
    new G4PVPlacement(nullptr, pos, NaILV, "NaIPVPL", detContainerLV, false, 0, true);
    NaILV->SetVisAttributes(visWhite);

    G4ThreeVector holePos1 = G4ThreeVector(-NaISize.y() / 2, 0, 0);
    G4ThreeVector holePos2 = G4ThreeVector(NaISize.y() / 2, 0, 0);
    if (LEDSize > 0) {
        hole = new G4Tubs("Hole", 0, LEDSize, vetoThick / 2. + tyvekThick + 2 * cm + shellThick / 2, 0, 360 * deg);
    }

    if (tyvekThick > 0.) {
        G4ThreeVector tyvekOutSize = G4ThreeVector(NaISize.y() + tyvekThick + vetoThick,
                                                  NaISize.y() + 2. * tyvekThick + vetoThick,
                                                  NaISize.z() + 2. * tyvekThick + vetoThick);
        G4Tubs *tyvekOutTube = new G4Tubs("TyvekOutTube", tyvekOutSize.x(), tyvekOutSize.y(), tyvekOutSize.z(), 0,
                                         viewDeg);
        G4Tubs *tyvekOutTop = new G4Tubs("TyvekOutTop", 0, tyvekOutSize.x(), tyvekThick / 2., 0,
                                        viewDeg);
        G4VSolid *tyvekOutBottom;
        if (doubleLED && LEDSize > 0) {
            G4SubtractionSolid *tyvekOutBottom1 = new G4SubtractionSolid("TyvekOutBottomIncomplete", tyvekOutTop, hole,
                                                                        nullptr, holePos1);
            tyvekOutBottom = new G4SubtractionSolid("TyvekOutBottom", tyvekOutBottom1, hole, nullptr, holePos2);
        } else {
            tyvekOutBottom = new G4Tubs("TyvekOutBottom", LEDSize, tyvekOutSize.x(), tyvekThick / 2., 0, viewDeg);
        }
        G4ThreeVector tyvekOutBottomPos = G4ThreeVector(0, 0, -tyvekOutSize.z() + tyvekThick / 2.);
        G4UnionSolid *tyvekOutIncomplete = new G4UnionSolid("TyvekOutIncomplete", tyvekOutTube, tyvekOutBottom, nullptr,
                                                           tyvekOutBottomPos);
        G4ThreeVector tyvekOutTopPos = G4ThreeVector(0, 0, tyvekOutSize.z() - tyvekThick / 2.);
        G4UnionSolid *tyvekOut = new G4UnionSolid("TyvekOut", tyvekOutIncomplete, tyvekOutTop, nullptr, tyvekOutTopPos);
        tyvekOutLV = new G4LogicalVolume(tyvekOut, tyvekOutMat, "TyvekOutLV");
        new G4PVPlacement(nullptr, pos, tyvekOutLV, "TyvekOutPVPL", detContainerLV, false, 0, true);
        tyvekOutLV->SetVisAttributes(visBlue);


        G4ThreeVector tyvekInSize = G4ThreeVector(NaISize.y(),
                                                 NaISize.y() + tyvekThick,
                                                 NaISize.z() + tyvekThick);
        G4Tubs *tyvekInTube = new G4Tubs("TyvekInTube", tyvekInSize.x(), tyvekInSize.y(), tyvekInSize.z(), 0, viewDeg);

        G4Tubs *tyvekInTop = new G4Tubs("TyvekInTop", 0, tyvekInSize.x(), tyvekThick / 2., 0,
                                       viewDeg);
        G4VSolid *tyvekInBottom;
        if (doubleLED && LEDSize > 0) {
            G4SubtractionSolid *tyvekInBottom1 = new G4SubtractionSolid("TyvekInBottomIncomplete", tyvekInTop, hole,
                                                                       nullptr, holePos1);
            tyvekInBottom = new G4SubtractionSolid("TyvekInBottom", tyvekInBottom1, hole, nullptr, holePos2);
        } else {
            tyvekInBottom = new G4Tubs("TyvekInBottom", LEDSize, tyvekInSize.x(), tyvekThick / 2., 0, viewDeg);
        }
        G4ThreeVector tyvekInBottomPos = G4ThreeVector(0, 0, -tyvekInSize.z() + tyvekThick / 2.);
        G4UnionSolid *tyvekInIncomplete = new G4UnionSolid("TyvekInIncomplete", tyvekInTube, tyvekInBottom, nullptr,
                                                          tyvekInBottomPos);
        G4ThreeVector tyvekInTopPos = G4ThreeVector(0, 0, tyvekInSize.z() - tyvekThick / 2.);
        G4UnionSolid *tyvekIn = new G4UnionSolid("TyvekIn", tyvekInIncomplete, tyvekInTop, nullptr, tyvekInTopPos);
        tyvekInLV = new G4LogicalVolume(tyvekIn, tyvekInMat, "TyvekInLV");
        new G4PVPlacement(nullptr, pos, tyvekInLV, "TyvekInPVPL", detContainerLV, false, 0, true);
        tyvekInLV->SetVisAttributes(visBlue);
    }

    if (vetoThick > 0.) {
        G4ThreeVector vetoSize = G4ThreeVector(NaISize.y() + tyvekThick,
                                               NaISize.y() + tyvekThick + vetoThick,
                                               NaISize.z() + tyvekThick + vetoThick);
        G4Tubs *vetoTube = new G4Tubs("VetoTube", vetoSize.x(), vetoSize.y(), vetoSize.z(), 0, viewDeg);

        G4Tubs *vetoTop = new G4Tubs("VetoTop", 0, vetoSize.x(), vetoThick / 2., 0, viewDeg);
        G4VSolid *vetoBottom;
        if (doubleLED && LEDSize > 0) {
            G4SubtractionSolid *vetoBottom1 = new G4SubtractionSolid("VetoBottomIncomplete", vetoTop, hole,
                                                                     nullptr, holePos1);
            vetoBottom = new G4SubtractionSolid("VetoBottom", vetoBottom1, hole, nullptr, holePos2);
        } else {
            vetoBottom = new G4Tubs("VetoBottom", LEDSize, vetoSize.x(), vetoThick / 2., 0, viewDeg);
        }
        G4ThreeVector vetoBottomPos = G4ThreeVector(0, 0, -vetoSize.z() + vetoThick / 2.);
        G4UnionSolid *vetoIncomplete = new G4UnionSolid("VetoIncomplete", vetoTube, vetoBottom, nullptr, vetoBottomPos);
        G4ThreeVector vetoTopPos = G4ThreeVector(0, 0, vetoSize.z() - vetoThick / 2.);
        G4UnionSolid *veto = new G4UnionSolid("Veto", vetoIncomplete, vetoTop, nullptr, vetoTopPos);
        vetoLV = new G4LogicalVolume(veto, vetoMat, "VetoLV");
        new G4PVPlacement(nullptr, pos, vetoLV, "VetoPVPL", detContainerLV, false, 0, true);
        vetoLV->SetVisAttributes(visRed);
    }
    ConstructShell();
    ConstructLED();
}


void NaI::ConstructShell() {
    if (shellThick <= 0.) return;

    G4ThreeVector shellSize = G4ThreeVector(detContainerSize.y() - shellThick,
                                            detContainerSize.y(),
                                            detContainerSize.z());
    G4Tubs *shellTube = new G4Tubs("ShellTube", shellSize.x(), shellSize.y(), shellSize.z(), 0, viewDeg);
    G4Tubs *shellTop = new G4Tubs("ShellTop", 0, detContainerSize.y(), shellThick / 2., 0, viewDeg);
    G4ThreeVector shellTopPos = G4ThreeVector(0, 0, detContainerSize.z() - shellThick / 2.);
    G4VSolid *shellBottom;
    if (doubleLED && LEDSize > 0) {
        G4ThreeVector NaISize = G4ThreeVector(0,
                                              detContainerSize.y() - vetoThick - tyvekThick * 2. - gapSize,
                                              detContainerSize.z() - vetoThick - tyvekThick * 2 - gapSize / 2);

        G4ThreeVector holePos1 = G4ThreeVector(-NaISize.y() / 2, 0, 0);
        G4ThreeVector holePos2 = G4ThreeVector(NaISize.y() / 2, 0, 0);
        G4SubtractionSolid *shellBottom1 = new G4SubtractionSolid("ShellBottomIncomplete", shellTop, hole,
                                                                  nullptr, holePos1);
        shellBottom = new G4SubtractionSolid("ShellBottom", shellBottom1, hole, nullptr, holePos2);
    } else {
        shellBottom = new G4Tubs("ShellBottom", LEDSize, shellSize.x(), shellThick / 2., 0, viewDeg);
    }
    G4ThreeVector shellBottomPos = G4ThreeVector(0, 0, -detContainerSize.z() + shellThick / 2.);
    G4UnionSolid *shellIncomplete = new G4UnionSolid("ShellIncomplete", shellTube, shellTop, nullptr, shellTopPos);
    G4UnionSolid *shell = new G4UnionSolid("Shell", shellIncomplete, shellBottom, nullptr, shellBottomPos);
    shellLV = new G4LogicalVolume(shell, shellMat, "ShellLV");
    G4ThreeVector shellPos = G4ThreeVector(0, 0, 0);
    new G4PVPlacement(nullptr, shellPos, shellLV, "ShellPVPL", detContainerLV, false, 0, true);
    shellLV->SetVisAttributes(visGrey);
}


void NaI::ConstructLED() {
    G4Tubs *LED = new G4Tubs("LED", 0, LEDSize, vetoThick / 2. + tyvekThick + shellThick / 2, 0,
                             viewDeg);
    if (doubleLED) {
        G4LogicalVolume *LEDLV1 = new G4LogicalVolume(LED, LEDMat, "LEDLV1");
        G4ThreeVector LEDLV1Pos = G4ThreeVector(NaISize.y() / 2, 0,
                                                -(gapSize / 2 + NaISize.z() + vetoThick / 2. +
                                                  tyvekThick + shellThick / 2));
        new G4PVPlacement(nullptr, LEDLV1Pos, LEDLV1, "LEDLV1PVPL", detContainerLV, false, 0, true);

        G4LogicalVolume *LEDLV2 = new G4LogicalVolume(LED, LEDMat, "LEDLV2");
        G4ThreeVector LEDLV2Pos = G4ThreeVector(-NaISize.y() / 2, 0,
                                                -(gapSize / 2 + NaISize.z() + vetoThick / 2. + tyvekThick +
                                                    shellThick / 2));
        new G4PVPlacement(nullptr, LEDLV2Pos, LEDLV2, "LEDLV2PVPL", detContainerLV, false, 0, true);
        LEDLV1->SetVisAttributes(visCyan);
        LEDLV2->SetVisAttributes(visCyan);
    } else {
        G4LogicalVolume *LEDLV = new G4LogicalVolume(LED, LEDMat, "LEDLV");
        G4ThreeVector LEDLVPos = G4ThreeVector(
            0, 0, -(gapSize / 2 + NaISize.z() + vetoThick / 2. + tyvekThick + shellThick / 2));
        new G4PVPlacement(nullptr, LEDLVPos, LEDLV, "LEDLVPVPL", detContainerLV, false, 0, true);
        LEDLV->SetVisAttributes(visCyan);
    }
}


std::vector<G4LogicalVolume *> NaI::GetSensitiveLV() const {
    return {
        NaILV,
        vetoLV,
        tyvekOutLV,
        tyvekInLV,
        shellLV
    };
}
