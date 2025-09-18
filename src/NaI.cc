#include "NaI.hh"


NaI::NaI(G4LogicalVolume *detContLV, const G4ThreeVector &detContSize, G4NistManager *nistMan, const Sizes &sizes,
         G4double vDeg, G4bool dLED) {
    detectorType = "NaI";
    detContainerLV = detContLV;
    detContainerSize = detContSize;
    nist = nistMan;
    vetoThick = sizes.vetoThick;
    tapeThick = sizes.tapeThick;
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
    // ========= ЭЛЕМЕНТЫ =========
    auto *elH = nist->FindOrBuildElement("H");
    auto *elC = nist->FindOrBuildElement("C");
    auto *elNa = nist->FindOrBuildElement("Na");
    auto *elI = nist->FindOrBuildElement("I");
    auto *elTl = nist->FindOrBuildElement("Tl"); // для легирования NaI(Tl)

    // ========= NaI(Tl) (кристалл) =========
    {
        const G4double rhoNaI = 3.67 * g / cm3;
        NaIMat = new G4Material("NaIMat", rhoNaI, 3, kStateSolid);

        const G4double wTl = 1.0e-3; // 0.1 wt%
        const G4double wNa_noTl = 0.153; // массовые доли Na и I в стех. NaI
        const G4double wI_noTl = 0.847;
        const G4double scale = 1.0 - wTl;

        NaIMat->AddElement(elNa, wNa_noTl * scale);
        NaIMat->AddElement(elI, wI_noTl * scale);
        NaIMat->AddElement(elTl, wTl);

        // Оптика NaI(Tl)
        const G4int n = 8;
        G4double eph[n] = {
            2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV
        }; // ~620–365 nm

        G4double rindex[n];
        G4double abslen[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.85;
            abslen[i] = 200. * cm; // усреднённая длина поглощения
        }

        // Спектр излучения ~415 nm (пик), нормированный профиль
        G4double emitNaI[n] = {0.02, 0.10, 0.50, 1.00, 0.85, 0.45, 0.12, 0.02};

        auto *mptNaI = new G4MaterialPropertiesTable();
        AddProp(mptNaI, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptNaI, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptNaI, "SCINTILLATIONCOMPONENT1", eph, emitNaI, n, true, true);

        mptNaI->AddConstProperty("SCINTILLATIONYIELD", 38000. / MeV);
        mptNaI->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptNaI->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 230. * ns);

        NaIMat->SetMaterialPropertiesTable(mptNaI);

        // Для неорганики Birks обычно не применяют
        NaIMat->GetIonisation()->SetBirksConstant(0. * mm / MeV);
    }

    // ========= PVT (поливинилтолуол) — антисовпадение =========
    {
        // Плотность типично ~1.032 g/cm3
        vetoMat = new G4Material("VetoMat", 1.032 * g / cm3, 2, kStateSolid);
        // Мономер C9H10
        vetoMat->AddElement(elC, 9);
        vetoMat->AddElement(elH, 10);

        const G4int n = 8;
        G4double eph[n] = {2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV};

        // Оптика PVT-подобного пластика (EJ-200/BC-408 близко по параметрам)
        G4double rindex[n];
        G4double abslen[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.58;
            abslen[i] = 380. * cm; // поглощение ~3.8 m
        }

        // Спектр излучения ~425 nm (пик)
        G4double emitPVT[n] = {0.01, 0.20, 0.60, 1.00, 0.65, 0.25, 0.05, 0.01};

        auto *mptVeto = new G4MaterialPropertiesTable();
        AddProp(mptVeto, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptVeto, "ABSLENGTH", eph, abslen, n, true, true);
        AddProp(mptVeto, "SCINTILLATIONCOMPONENT1", eph, emitPVT, n, true, true);

        mptVeto->AddConstProperty("SCINTILLATIONYIELD", 10000. / MeV);
        mptVeto->AddConstProperty("RESOLUTIONSCALE", 1.0);
        mptVeto->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns);

        vetoMat->SetMaterialPropertiesTable(mptVeto);

        // Birks для органики (классическое значение порядка 0.126 mm/MeV)
        vetoMat->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }

    // ========= Tyvek (объёмно, а не как поверхность) =========
    // Приблизим Tyvek как вспененный HDPE (полиэтилен высокой плотности) с низкой эффективной плотностью.
    // Состав (CH2)n, плотность ~0.38 g/cm3. Зададим сильное Рэлеевское рассеяние,
    // чтобы слой вёл себя диффузно (при желании можно ещё добавить оптическую поверхность с LUT).
    {
        const G4double rhoTyvek = 0.38 * g / cm3;

        tapeInMat = new G4Material("tapeInMat", rhoTyvek, 2, kStateSolid);
        tapeOutMat = new G4Material("tapeOutMat", rhoTyvek, 2, kStateSolid);
        tapeInMat->AddElement(elC, 1);
        tapeInMat->AddElement(elH, 2);
        tapeOutMat->AddElement(elC, 1);
        tapeOutMat->AddElement(elH, 2);

        const G4int n = 8;
        G4double eph[n] = {2.0 * eV, 2.2 * eV, 2.4 * eV, 2.6 * eV, 2.8 * eV, 3.0 * eV, 3.2 * eV, 3.4 * eV};

        G4double rindex[n];
        G4double abslen[n];
        G4double rayl[n];
        for (G4int i = 0; i < n; ++i) {
            rindex[i] = 1.50;
            abslen[i] = 1000. * m; // по сути «прозрачный» в объёме
            rayl[i] = 10. * um; // сильное Рэлеевское рассеяние => диффузность
        }

        auto *mptTyvek = new G4MaterialPropertiesTable();
        AddProp(mptTyvek, "RINDEX", eph, rindex, n, true, true);
        AddProp(mptTyvek, "ABSLENGTH", eph, abslen, n, false, true);
        AddProp(mptTyvek, "RAYLEIGH", eph, rayl, n, false, true);

        tapeInMat->SetMaterialPropertiesTable(mptTyvek);

        auto *mptTyvekOut = new G4MaterialPropertiesTable(*mptTyvek);
        tapeOutMat->SetMaterialPropertiesTable(mptTyvekOut);
    }

    LEDMat = nist->FindOrBuildMaterial("G4_Galactic");
    shellMat = nist->FindOrBuildMaterial("G4_Al");
}


void NaI::Construct() {
    G4ThreeVector pos = G4ThreeVector(0, 0, -gapSize / 2);

    NaISize = G4ThreeVector(0,
                            detContainerSize.y() - vetoThick - tapeThick * 2. - gapSize - shellThick,
                            detContainerSize.z() - vetoThick - tapeThick * 2 - gapSize / 2 - shellThick);

    G4Tubs *NaIBox = new G4Tubs("NaI", NaISize.x(), NaISize.y(), NaISize.z(), 0, viewDeg);
    NaILV = new G4LogicalVolume(NaIBox, NaIMat, "NaILV");
    new G4PVPlacement(nullptr, pos, NaILV, "NaIPVPL", detContainerLV, false, 0, true);
    NaILV->SetVisAttributes(visWhite);

    G4ThreeVector holePos1 = G4ThreeVector(-NaISize.y() / 2, 0, 0);
    G4ThreeVector holePos2 = G4ThreeVector(NaISize.y() / 2, 0, 0);
    if (LEDSize > 0) {
        hole = new G4Tubs("Hole", 0, LEDSize, vetoThick / 2. + tapeThick + 2 * cm + shellThick / 2, 0, 360 * deg);
    }

    if (tapeThick > 0.) {
        G4ThreeVector tapeOutSize = G4ThreeVector(NaISize.y() + tapeThick + vetoThick,
                                                  NaISize.y() + 2. * tapeThick + vetoThick,
                                                  NaISize.z() + 2. * tapeThick + vetoThick);
        G4Tubs *tapeOutTube = new G4Tubs("TapeOutTube", tapeOutSize.x(), tapeOutSize.y(), tapeOutSize.z(), 0,
                                         viewDeg);
        G4Tubs *tapeOutTop = new G4Tubs("TapeOutTop", 0, tapeOutSize.x(), tapeThick / 2., 0,
                                        viewDeg);
        G4VSolid *tapeOutBottom;
        if (doubleLED && LEDSize > 0) {
            G4SubtractionSolid *tapeOutBottom1 = new G4SubtractionSolid("TapeOutBottomIncomplete", tapeOutTop, hole,
                                                                        nullptr, holePos1);
            tapeOutBottom = new G4SubtractionSolid("TapeOutBottom", tapeOutBottom1, hole, nullptr, holePos2);
        } else {
            tapeOutBottom = new G4Tubs("TapeOutBottom", LEDSize, tapeOutSize.x(), tapeThick / 2., 0, viewDeg);
        }
        G4ThreeVector tapeOutBottomPos = G4ThreeVector(0, 0, -tapeOutSize.z() + tapeThick / 2.);
        G4UnionSolid *tapeOutIncomplete = new G4UnionSolid("TapeOutIncomplete", tapeOutTube, tapeOutBottom, nullptr,
                                                           tapeOutBottomPos);
        G4ThreeVector tapeOutTopPos = G4ThreeVector(0, 0, tapeOutSize.z() - tapeThick / 2.);
        G4UnionSolid *tapeOut = new G4UnionSolid("TapeOut", tapeOutIncomplete, tapeOutTop, nullptr, tapeOutTopPos);
        tapeOutLV = new G4LogicalVolume(tapeOut, tapeOutMat, "TapeOutLV");
        new G4PVPlacement(nullptr, pos, tapeOutLV, "TapeOutPVPL", detContainerLV, false, 0, true);
        tapeOutLV->SetVisAttributes(visBlue);


        G4ThreeVector tapeInSize = G4ThreeVector(NaISize.y(),
                                                 NaISize.y() + tapeThick,
                                                 NaISize.z() + tapeThick);
        G4Tubs *tapeInTube = new G4Tubs("TapeInTube", tapeInSize.x(), tapeInSize.y(), tapeInSize.z(), 0, viewDeg);

        G4Tubs *tapeInTop = new G4Tubs("TapeInTop", 0, tapeInSize.x(), tapeThick / 2., 0,
                                       viewDeg);
        G4VSolid *tapeInBottom;
        if (doubleLED && LEDSize > 0) {
            G4SubtractionSolid *tapeInBottom1 = new G4SubtractionSolid("TapeInBottomIncomplete", tapeInTop, hole,
                                                                       nullptr, holePos1);
            tapeInBottom = new G4SubtractionSolid("TapeInBottom", tapeInBottom1, hole, nullptr, holePos2);
        } else {
            tapeInBottom = new G4Tubs("TapeInBottom", LEDSize, tapeInSize.x(), tapeThick / 2., 0, viewDeg);
        }
        G4ThreeVector tapeInBottomPos = G4ThreeVector(0, 0, -tapeInSize.z() + tapeThick / 2.);
        G4UnionSolid *tapeInIncomplete = new G4UnionSolid("TapeInIncomplete", tapeInTube, tapeInBottom, nullptr,
                                                          tapeInBottomPos);
        G4ThreeVector tapeInTopPos = G4ThreeVector(0, 0, tapeInSize.z() - tapeThick / 2.);
        G4UnionSolid *tapeIn = new G4UnionSolid("TapeIn", tapeInIncomplete, tapeInTop, nullptr, tapeInTopPos);
        tapeInLV = new G4LogicalVolume(tapeIn, tapeInMat, "TapeInLV");
        new G4PVPlacement(nullptr, pos, tapeInLV, "TapeInPVPL", detContainerLV, false, 0, true);
        tapeInLV->SetVisAttributes(visBlue);
    }

    if (vetoThick > 0.) {
        G4ThreeVector vetoSize = G4ThreeVector(NaISize.y() + tapeThick,
                                               NaISize.y() + tapeThick + vetoThick,
                                               NaISize.z() + tapeThick + vetoThick);
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
                                              detContainerSize.y() - vetoThick - tapeThick * 2. - gapSize,
                                              detContainerSize.z() - vetoThick - tapeThick * 2 - gapSize / 2);

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
    G4Tubs *LED = new G4Tubs("LED", 0, LEDSize, vetoThick / 2. + tapeThick + shellThick / 2, 0,
                             viewDeg);
    if (doubleLED) {
        G4LogicalVolume *LEDLV1 = new G4LogicalVolume(LED, LEDMat, "LEDLV1");
        G4ThreeVector LEDLV1Pos = G4ThreeVector(NaISize.y() / 2, 0,
                                                -(gapSize / 2 + NaISize.z() + vetoThick / 2. +
                                                  tapeThick + shellThick / 2));
        new G4PVPlacement(nullptr, LEDLV1Pos, LEDLV1, "LEDLV1PVPL", detContainerLV, false, 0, true);

        G4LogicalVolume *LEDLV2 = new G4LogicalVolume(LED, LEDMat, "LEDLV2");
        G4ThreeVector LEDLV2Pos = G4ThreeVector(-NaISize.y() / 2, 0,
                                                -(gapSize / 2 + NaISize.z() + vetoThick / 2. + tapeThick +
                                                    shellThick / 2));
        new G4PVPlacement(nullptr, LEDLV2Pos, LEDLV2, "LEDLV2PVPL", detContainerLV, false, 0, true);
        LEDLV1->SetVisAttributes(visCyan);
        LEDLV2->SetVisAttributes(visCyan);
    } else {
        G4LogicalVolume *LEDLV = new G4LogicalVolume(LED, LEDMat, "LEDLV");
        G4ThreeVector LEDLVPos = G4ThreeVector(
            0, 0, -(gapSize / 2 + NaISize.z() + vetoThick / 2. + tapeThick + shellThick / 2));
        new G4PVPlacement(nullptr, LEDLVPos, LEDLV, "LEDLVPVPL", detContainerLV, false, 0, true);
        LEDLV->SetVisAttributes(visCyan);
    }
}


std::vector<G4LogicalVolume *> NaI::GetSensitiveLV() const {
    return {
        NaILV,
        vetoLV,
        tapeOutLV,
        tapeInLV,
        shellLV
    };
}
