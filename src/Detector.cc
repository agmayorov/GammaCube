#include "Detector.hh"

using namespace Sizes;


Detector::Detector(G4LogicalVolume* detContLV, G4NistManager* nistMan, G4double vDeg, const G4int yScale,
                   const G4String& detType) {
    detectorType = detType;
    detContainerLV = detContLV;

    detContainerTopSize = G4ThreeVector(0, modelRadius - tunaCanThickWall,
                                        modelHeight - tunaCanThickTop + plateCenterThick);

    nist = nistMan;
    viewDeg = vDeg;
    yieldScale = yScale;

    crystalSize = G4ThreeVector(0, crystalRadius, crystalHeight);

    DefineMaterials();
    DefineVisual();
}

void Detector::DefineMaterials() {
    auto* elH = nist->FindOrBuildElement("H");
    auto* elC = nist->FindOrBuildElement("C");
    auto* elNa = nist->FindOrBuildElement("Na");
    auto* elCs = nist->FindOrBuildElement("Cs");
    auto* elI = nist->FindOrBuildElement("I");
    auto* elTl = nist->FindOrBuildElement("Tl");
    auto* elN = nist->FindOrBuildElement("N");
    auto* elO = nist->FindOrBuildElement("O");

    G4Material* SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    // --------- helper lambdas
    const std::string base = "../OpticalParameters/";

    auto loadRIndex = [&](const std::string& prefix) -> Utils::Table {
        // rindex is dimensionless
        return Utils::ReadCSV(base + prefix + "_refractive_index.csv", /*valueScale=*/1.0, /*clampNonNegative=*/false);
    };

    auto loadAbsLengthMM = [&](const std::string& prefix) -> Utils::Table {
        // ABSLENGTH values are stored in mm in your files
        return Utils::ReadCSV(base + prefix + "_absorption_length.csv", /*valueScale=*/mm, /*clampNonNegative=*/true);
    };

    auto loadEmission2 = [&](const std::string& prefix,
                             const std::string& suffix = "_normalised_emission_intensity.csv")
        -> Utils::EmissionTables {
        auto e = Utils::ReadEmissionCSV(base + prefix + suffix, /*valueScale=*/1.0, /*clampNonNegative=*/true);

        Utils::NormalizeMaxToOne(e.c1);
        Utils::NormalizeMaxToOne(e.c2);

        return e;
    };

    auto loadConsts = [&](const std::string& prefix) -> Utils::ConstMap {
        return Utils::ReadConstFile(base + prefix + "_optical_consts.txt");
    };

    auto applyYieldScaleIfPresent = [&](Utils::ConstMap& c) {
        auto it = c.find("SCINTILLATIONYIELD");
        if (it != c.end() && yieldScale > 0) {
            it->second /= static_cast<G4double>(yieldScale);
        }
    };

    // Crystal (NaI or CsI)
    {
        const std::string matPrefix = detectorType;

        // --- Create material by composition (as you had)
        if (detectorType == "NaI") {
            const G4double rhoCrystal = 3.67 * g / cm3;
            CrystalMat = new G4Material("CrystalMat", rhoCrystal, 3, kStateSolid);

            const G4double wTl = 0.065;
            const G4double wNa_noTl = 0.153;
            const G4double wI_noTl = 0.847;
            const G4double scale = 1.0 - wTl;

            CrystalMat->AddElement(elNa, wNa_noTl * scale);
            CrystalMat->AddElement(elI, wI_noTl * scale);
            CrystalMat->AddElement(elTl, wTl);
        } else if (detectorType == "CsI") {
            const G4double rhoCrystal = 4.51 * g / cm3;
            CrystalMat = new G4Material("CrystalMat", rhoCrystal, 3, kStateSolid);

            const G4double wTl = 0.0008;
            const G4double wCs_noTl = 0.511549;
            const G4double wI_noTl = 0.488451;
            const G4double scale = 1.0 - wTl;

            CrystalMat->AddElement(elCs, wCs_noTl * scale);
            CrystalMat->AddElement(elI, wI_noTl * scale);
            CrystalMat->AddElement(elTl, wTl);
        } else {
            throw std::runtime_error("Unknown detectorType (expected NaI or CsI): " + std::string(detectorType));
        }

        // --- Load optical tables
        const auto rindex = loadRIndex(matPrefix);
        const auto absl = loadAbsLengthMM(matPrefix);
        auto emission = loadEmission2(matPrefix);

        // --- Load scint consts
        auto c = loadConsts(matPrefix);
        applyYieldScaleIfPresent(c);

        // --- Build MPT
        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", rindex.E, rindex.V, rindex.E.size());
        mpt->AddProperty("ABSLENGTH", absl.E, absl.V, absl.E.size());

        Utils::ApplyScintillation(CrystalMat, mpt, c, emission.c1, emission.c2, true);

        Utils::ApplyBirksIfPresent(CrystalMat, c);
    }

    // Veto (plastic scintillator)
    {
        vetoMat = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

        const std::string p = "Veto";

        const auto rindex = loadRIndex(p);
        const auto absl = loadAbsLengthMM(p);
        auto emission = loadEmission2(p);

        auto c = loadConsts(p);
        applyYieldScaleIfPresent(c);

        auto* mpt = new G4MaterialPropertiesTable();
        mpt->AddProperty("RINDEX", rindex.E, rindex.V, rindex.E.size());
        mpt->AddProperty("ABSLENGTH", absl.E, absl.V, absl.E.size());

        Utils::ApplyScintillation(vetoMat, mpt, c, emission.c1, emission.c2, true);

        Utils::ApplyBirksIfPresent(vetoMat, c);
    }

    // Tyvek material
    {
        const G4double rhoTyvek = 0.38 * g / cm3;
        tyvekMat = new G4Material("Tyvek", rhoTyvek, 2, kStateSolid);
        tyvekMat->AddElement(elC, 2);
        tyvekMat->AddElement(elH, 4);
    }

    // World / Vacuum
    {
        galacticMat = nist->FindOrBuildMaterial("G4_Galactic");

        auto tN = Utils::MakeConstantTable(1.0 * eV, 4.0 * eV, 1.0);
        auto tA = Utils::MakeConstantTable(1.0 * eV, 4.0 * eV, 1e6 * m);
        Utils::ApplyMaterialTable(galacticMat, tN, &tA);
    }

    // SiPM bulk (silicon)
    SiPMMat = nist->FindOrBuildMaterial("G4_Si");

    // SiPM encapsulant (clear epoxy proxy): Bisphenol A diglycidyl ether (BADGE / DGEBA)
    {
        const G4double rho = 1.16 * g / cm3;

        SiPMEncapsulantMat = new G4Material("SiPMEncapsulant_DGEBA", rho, 3, kStateSolid);

        // Stoichiometry from molecular formula C21H24O4
        SiPMEncapsulantMat->AddElement(elC, 21);
        SiPMEncapsulantMat->AddElement(elH, 24);
        SiPMEncapsulantMat->AddElement(elO, 4);

        const auto rindex = loadRIndex("SiPM_Encapsulant");
        const auto absl = loadAbsLengthMM("SiPM_Encapsulant");

        Utils::ApplyMaterialTable(SiPMEncapsulantMat, rindex, &absl);
    }


    // SiPM frame / Aluminum
    SiPMFrameMat = nist->FindOrBuildMaterial("G4_Al");
    AlMat = nist->FindOrBuildMaterial("G4_Al");

    // Rubber
    rubberMat = nist->FindOrBuildMaterial("G4_RUBBER_NEOPRENE");

    // Board
    {
        G4Material* Epoxy = new G4Material("Epoxy", 1.2 * g / cm3, 2);
        Epoxy->AddElement(elH, 2);
        Epoxy->AddElement(elC, 2);

        boardMat = new G4Material("FiberglassLaminate", 1.86 * g / cm3, 2);
        boardMat->AddMaterial(Epoxy, 0.472);
        boardMat->AddMaterial(SiO2, 0.528);
    }

    // Crystal glass interface plate
    {
        glassMat = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
        const auto rindex = loadRIndex("Glass");
        const auto absl = loadAbsLengthMM("Glass");
        Utils::ApplyMaterialTable(glassMat, rindex, &absl);
    }

    // Optical coupling layer
    {
        // DGEBA (C19H20O4), 0.56 by weight
        auto* matDGEBA = new G4Material("Epotek301_DGEBA", 1.16 * g / cm3, 3, kStateSolid);
        matDGEBA->AddElement(elC, 19);
        matDGEBA->AddElement(elH, 20);
        matDGEBA->AddElement(elO, 4);

        // 1,4-Butanediol Diglycidyl Ether (C10H18O4)
        auto* matBDDE = new G4Material("Epotek301_BDDE", 1.10 * g / cm3, 3, kStateSolid);
        matBDDE->AddElement(elC, 10);
        matBDDE->AddElement(elH, 18);
        matBDDE->AddElement(elO, 4);

        // 1,6-Hexanediamine 2,2,4-trimethyl- (C9H22N2)
        auto* matAmine = new G4Material("Epotek301_Amine", 0.865 * g / cm3, 3, kStateSolid);
        matAmine->AddElement(elC, 9);
        matAmine->AddElement(elH, 22);
        matAmine->AddElement(elN, 2);

        // ЕPO-TEK 301-1
        const G4double rhoOpticLayer = 1.19 * g / cm3;
        opticLayerMat = new G4Material("OpticLayerMat_Epotek301_1", rhoOpticLayer, 3, kStateSolid);

        opticLayerMat->AddMaterial(matDGEBA, 0.56);
        opticLayerMat->AddMaterial(matBDDE, 0.24);
        opticLayerMat->AddMaterial(matAmine, 0.20);

        const auto rindex = loadRIndex("OpticLayer");
        const auto absl = loadAbsLengthMM("OpticLayer");
        Utils::ApplyMaterialTable(opticLayerMat, rindex, &absl);
    }

    // Payload
    payloadMat = nist->FindOrBuildMaterial("G4_Si");
}


G4String Detector::GetDetectorType() const {
    return detectorType;
}


void Detector::DefineVisual() {
    visCrystal = new G4VisAttributes(G4Color(1.0, 0.337, 0.0));
    visCrystal->SetForceSolid(true);
    visTyvekOut = new G4VisAttributes(G4Color(0.27, 0.0, 0.65));
    visTyvekOut->SetForceSolid(true);
    visTyvekMid = new G4VisAttributes(G4Color(0.27, 0.0, 0.65));
    visTyvekMid->SetForceSolid(true);
    visTyvekIn = new G4VisAttributes(G4Color(0.27, 0.0, 0.65));
    visTyvekIn->SetForceSolid(true);
    visTyvekBottom = new G4VisAttributes(G4Color(0.27, 0.0, 0.65));
    visTyvekBottom->SetForceSolid(true);
    visVeto = new G4VisAttributes(G4Color(1.0, 0.615, 0.0));
    visVeto->SetForceSolid(true);
    visSpring = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visSpring->SetForceSolid(true);
    visSiPM = new G4VisAttributes(G4Color(1.0, 0.5, 1.0));
    visSiPM->SetForceSolid(true);
    visSiPMFrame = new G4VisAttributes(G4Color(0.0, 0.0, 0.0));
    visSiPMFrame->SetForceSolid(true);
    visAl = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visAl->SetForceSolid(true);
    visShell = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visShell->SetForceSolid(true);
    visHolder = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visHolder->SetForceSolid(true);
    visGlass = new G4VisAttributes(G4Color(0.0, 0.876, 0.96, 0.5));
    visGlass->SetForceSolid(true);
    visOpticLayer = new G4VisAttributes(G4Color(0.635, 0.95, 0.0));
    visOpticLayer->SetForceSolid(true);
    visBottomVetoShell = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visBottomVetoShell->SetForceSolid(true);
    visGasket = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visGasket->SetForceSolid(true);
    visBoard = new G4VisAttributes(G4Color(1.0, 1.0, 0.0));
    visBoard->SetForceSolid(true);
    visPayload = new G4VisAttributes(G4Color(0.5, 0.5, 0.0));
    visPayload->SetForceSolid(true);
    visRubber = new G4VisAttributes(G4Color(0.0, 0.0, 0.0));
    visRubber->SetForceSolid(true);
}


void Detector::Construct() {
    ConstructVeto();
    ConstructShell();
    ConstructBottomVeto();
    ConstructCrystal();
    ConstructSiPM();
    ConstructCrystalSiPM();
    ConstructBottomVetoSiPM();
    ConstructVetoSiPM();
    ConstructOpticalSurfaces();
}

void Detector::ConstructVeto() {
    // Tyvek Out
    tyvekOutSize = G4ThreeVector(vetoRadius, vetoRadius + tyvekOutThickWall, vetoHeight / 2);
    G4ThreeVector tyvekOutPos = G4ThreeVector(0, 0, (detContainerTopSize.z() - vetoHeight -
                                                  vetoChamferHeight - vetoTopRoundedRadius) / 2 - tyvekOutThickTop);

    G4VSolid* tyvekOutWall = new G4Tubs("TyvekOutWall", tyvekOutSize.x(), tyvekOutSize.y(),
                                        tyvekOutSize.z() - (vetoChamferHeight + vetoTopRoundedRadius) / 2, 0, viewDeg);
    G4ThreeVector tyvekOutTopPos = G4ThreeVector(0, 0, tyvekOutSize.z() + tyvekOutThickTop / 2.0);
    G4VSolid* tyvekOutTop;
    if (vetoTopRoundedRadius > 0) {
        G4double tyvekOutTopRoundingRadius = vetoTopRoundedRadius + tyvekOutThickWall;
        G4VSolid* tyvekOutTopRoundingBase = new G4Torus("TyvekOutTopRoundingBase", vetoTopRoundedRadius,
                                                        tyvekOutTopRoundingRadius,
                                                        tyvekOutSize.y() - tyvekOutTopRoundingRadius, 0, viewDeg);
        G4VSolid* cutCube = new G4Box("CutCube", tyvekOutSize.y(), tyvekOutSize.y(), tyvekOutSize.y());
        G4VSolid* cutTube = new G4Tubs("CutTube", 0, tyvekOutSize.y() - tyvekOutTopRoundingRadius, tyvekOutSize.z(), 0,
                                       viewDeg);
        G4VSolid* tyvekOutTopRoundingInc = new G4SubtractionSolid("TyvekOutTopRoundingIncomplete",
                                                                  tyvekOutTopRoundingBase, cutCube, nullptr,
                                                                  G4ThreeVector(0, 0, -tyvekOutSize.y()));
        G4VSolid* tyvekOutTopRounding = new G4SubtractionSolid("TyvekOutTopRounding", tyvekOutTopRoundingInc, cutTube);
        G4VSolid* tyvekOutTopTube = new G4Tubs("TyvekOutTopTube", 0, tyvekOutSize.y() - tyvekOutTopRoundingRadius,
                                               tyvekOutThickTop / 2, 0, viewDeg);
        tyvekOutTop = new G4UnionSolid("TyvekOutTop", tyvekOutTopRounding, tyvekOutTopTube, nullptr,
                                       G4ThreeVector(0, 0, vetoTopRoundedRadius + tyvekOutThickTop / 2.));
        tyvekOutTopPos = G4ThreeVector(0, 0, tyvekOutSize.z() - (vetoChamferHeight + vetoTopRoundedRadius) / 2);
    } else if (vetoChamferHeight != 0) {
        G4VSolid* tyvekOutTopConsOut = new G4Cons("TyvekOutTopConsOut", tyvekOutSize.x(), tyvekOutSize.y(),
                                                  tyvekOutSize.y() - 2 * tyvekOutThickWall - vetoChamferHeight,
                                                  tyvekOutSize.y() - tyvekOutThickWall - vetoChamferHeight,
                                                  (vetoChamferHeight + tyvekOutThickTop) / 2, 0, viewDeg);
        G4VSolid* tyvekOutTopConsIn = new G4Cons("TyvekOutTopConsIn", 0, tyvekOutSize.x() - vetoChamferHeight, 0,
                                                 tyvekOutSize.y() - 2 * tyvekOutThickWall - vetoChamferHeight,
                                                 tyvekOutThickTop / 2, 0, viewDeg);
        G4ThreeVector tyvekOutTopConsPos = G4ThreeVector(0, 0, (vetoChamferHeight + tyvekOutThickWall -
                                                             tyvekOutThickTop) / 2.);
        tyvekOutTop = new G4UnionSolid("TyvekOutTop", tyvekOutTopConsOut, tyvekOutTopConsIn, nullptr,
                                       tyvekOutTopConsPos);
    } else {
        tyvekOutTop = new G4Tubs("TyvekOutTop", 0, tyvekOutSize.y(), tyvekOutThickTop / 2., 0, viewDeg);
    }
    G4VSolid* tyvekOut = new G4UnionSolid("TyvekOut", tyvekOutTop, tyvekOutWall, nullptr, -tyvekOutTopPos);

    tyvekOutLV = new G4LogicalVolume(tyvekOut, tyvekMat, "TyvekOutLV");
    tyvekOutPVP = new G4PVPlacement(nullptr, tyvekOutPos + tyvekOutTopPos, tyvekOutLV, "TyvekOutPVP", detContainerLV,
                                    false, 0, true);
    tyvekOutLV->SetVisAttributes(visTyvekOut);

    // Veto
    vetoSize = tyvekOutSize - G4ThreeVector(vetoThickWall, tyvekOutThickWall, vetoThickTop / 2.);
    G4ThreeVector vetoPos = tyvekOutPos;
    vetoPos[2] -= (vetoThickTop - vetoChamferHeight - vetoTopRoundedRadius) / 2;

    G4VSolid* vetoWall = new G4Tubs("VetoWall", vetoSize.x(), vetoSize.y(), vetoSize.z(), 0, viewDeg);
    G4ThreeVector vetoTopPos =
        G4ThreeVector(0, 0, vetoSize.z() + (vetoThickTop - vetoChamferHeight - vetoTopRoundedRadius) / 2.0);
    G4VSolid* vetoTop;
    if (vetoTopRoundedRadius > 0) {
        G4VSolid* vetoTopRoundingBase = new G4Torus("VetoTopRoundingBase", 0, vetoTopRoundedRadius,
                                                    vetoSize.y() - vetoTopRoundedRadius, 0, viewDeg);
        G4VSolid* cutCube = new G4Box("CutCube", vetoSize.y(), vetoSize.y(), vetoSize.y());
        G4VSolid* vetoTopRounding = new G4SubtractionSolid("VetoTopRounding", vetoTopRoundingBase, cutCube,
                                                           nullptr, G4ThreeVector(0, 0, -vetoSize.y()));
        G4VSolid* vetoTopTube = new G4Tubs("VetoTopTube", 0, vetoSize.y() - vetoTopRoundedRadius,
                                           vetoTopRoundedRadius / 2, 0, viewDeg);
        if (vetoTopRoundedRadius == vetoThickTop) {
            vetoTop = new G4UnionSolid("VetoTop", vetoTopRounding, vetoTopTube, nullptr,
                                       G4ThreeVector(0, 0, vetoTopRoundedRadius / 2.));
        } else {
            G4VSolid* vetoTopInc = new G4UnionSolid("VetoTop", vetoTopRounding, vetoTopTube, nullptr,
                                                    G4ThreeVector(0, 0, vetoTopRoundedRadius / 2.));
            G4VSolid* vetoTopBaseTube = new G4Tubs("VetoTopBaseTube", 0, vetoSize.y(),
                                                   (vetoThickTop - vetoTopRoundedRadius) / 2., 0, viewDeg);
            vetoTop = new G4UnionSolid("VetoTop", vetoTopInc, vetoTopBaseTube, nullptr,
                                       G4ThreeVector(0, 0, -(vetoThickTop - vetoTopRoundedRadius) / 2));
            vetoTopPos[2] += (vetoThickTop - vetoTopRoundedRadius) / 2;
        }
    } else if (vetoChamferHeight == vetoThickTop) {
        vetoTop = new G4Cons("VetoTop", 0, vetoSize.y(), 0, vetoSize.y() - vetoChamferHeight, vetoChamferHeight / 2, 0,
                             viewDeg);
        vetoTopPos = G4ThreeVector(0, 0, vetoSize.z() + vetoChamferHeight / 2.0);
    } else if (vetoChamferHeight > 0 and vetoChamferHeight < vetoThickTop) {
        G4VSolid* vetoTopCyl = new G4Tubs("VetoTopCyl", 0, vetoSize.y(), (vetoThickTop - vetoChamferHeight) / 2., 0,
                                          viewDeg);
        G4VSolid* vetoTopCons = new G4Cons("VetoTopCons", 0, vetoSize.y(), 0, vetoSize.y() - vetoChamferHeight,
                                           vetoChamferHeight / 2, 0, viewDeg);
        G4ThreeVector vetoTopConsPos = G4ThreeVector(0, 0, vetoThickTop / 2.);
        vetoTop = new G4UnionSolid("VetoTop", vetoTopCyl, vetoTopCons, nullptr, vetoTopConsPos);
    } else {
        vetoTop = new G4Tubs("VetoTop", 0, vetoSize.y(), vetoThickTop / 2., 0, viewDeg);
    }
    G4VSolid* veto = new G4UnionSolid("Veto", vetoTop, vetoWall, nullptr, -vetoTopPos);

    vetoLV = new G4LogicalVolume(veto, vetoMat, "VetoLV");
    vetoPVP = new G4PVPlacement(nullptr, vetoPos + vetoTopPos, vetoLV, "VetoPVP", detContainerLV, false, 0, true);
    vetoLV->SetVisAttributes(visVeto);

    // Tyvek Mid
    tyvekMidSize = vetoSize - G4ThreeVector(tyvekMidThickWall, vetoThickWall, 0);
    G4ThreeVector tyvekMidPos = vetoPos;

    G4VSolid* tyvekMidWall = new G4Tubs("TyvekMidWall", tyvekMidSize.x(), tyvekMidSize.y(), tyvekMidSize.z(), 0,
                                        viewDeg);
    G4VSolid* tyvekMidTop = new G4Tubs("TyvekMidTop", 0, tyvekMidSize.x(), tyvekMidThickTop / 2., 0, viewDeg);
    G4ThreeVector tyvekMidTopPos = G4ThreeVector(0, 0, tyvekMidSize.z() - tyvekMidThickTop / 2.0);
    G4VSolid* tyvekMid = new G4UnionSolid("TyvekMid", tyvekMidWall, tyvekMidTop, nullptr, tyvekMidTopPos);

    tyvekMidLV = new G4LogicalVolume(tyvekMid, tyvekMat, "TyvekMidLV");
    tyvekMidPVP = new G4PVPlacement(nullptr, tyvekMidPos, tyvekMidLV, "TyvekMidPVP", detContainerLV, false, 0, true);
    tyvekMidLV->SetVisAttributes(visTyvekMid);

    // Optic Layer for Veto
    G4VSolid* vetoOpticLayer = new G4Tubs("VetoOpticLayer", tyvekMidSize.x(), tyvekOutSize.y(),
                                          vetoOpticLayerHeight / 2., 0, viewDeg);
    vetoOpticLayerLV = new G4LogicalVolume(vetoOpticLayer, opticLayerMat, "VetoOpticLayerLV");
    G4ThreeVector vetoOpticLayerPos = tyvekOutPos + G4ThreeVector(0, 0, -tyvekOutSize.z() - (vetoOpticLayerHeight -
                                                                      vetoChamferHeight - vetoTopRoundedRadius) / 2.0);
    vetoOpticLayerPVP = new G4PVPlacement(nullptr, vetoOpticLayerPos, vetoOpticLayerLV, "VetoOpticLayerPVP",
                                          detContainerLV, false, 0,
                                          true);
    vetoOpticLayerLV->SetVisAttributes(visOpticLayer);
}

void Detector::ConstructShell() {
    // Rubber gasket
    G4ThreeVector rubberPos = G4ThreeVector(
                                            0, 0, (detContainerTopSize.z() - rubberHeight) / 2.0 - tyvekOutThickTop -
                                            vetoThickTop - tyvekMidThickTop);
    G4VSolid* rubber = new G4Tubs("Rubber", 0, tyvekMidSize.x(), rubberHeight / 2., 0, viewDeg);
    G4LogicalVolume* rubberLV = new G4LogicalVolume(rubber, rubberMat, "RubberLV");
    new G4PVPlacement(nullptr, rubberPos, rubberLV, "RubberPVP", detContainerLV, false, 0, true);
    rubberLV->SetVisAttributes(visRubber);

    // Core Logical volume for all core components
    coreTopSize = G4ThreeVector(0, tyvekMidSize.x(),
                                detContainerTopSize.z() - rubberHeight - tyvekOutThickTop - vetoThickTop -
                                tyvekMidThickTop + plateThick);
    G4VSolid* coreTop = new G4Tubs("CoreTop", 0, coreTopSize.y(), coreTopSize.z() / 2., 0, viewDeg);
    G4VSolid* coreBottom = new G4Tubs("CoreBottom", 0, bottomCapInnerRadius, (bottomCapHeight - bottomCapThick) / 2., 0,
                                      viewDeg);
    G4ThreeVector coreBottomPos = G4ThreeVector(0, 0, -(coreTopSize.z() + bottomCapHeight - bottomCapThick) / 2.);
    G4VSolid* core = new G4UnionSolid("Core", coreTop, coreBottom, nullptr, coreBottomPos);

    G4ThreeVector corePos = G4ThreeVector(
                                          0, 0, (detContainerTopSize.z() - coreTopSize.z()) / 2.0 - tyvekOutThickTop -
                                          vetoThickTop - tyvekMidThickTop -
                                          rubberHeight);
    coreLV = new G4LogicalVolume(core, galacticMat, "CoreLV");
    new G4PVPlacement(nullptr, corePos, coreLV, "CorePVP", detContainerLV, false, 0, true);
    coreLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Shell
    G4ThreeVector shellWallSize = coreTopSize;
    shellWallSize[0] = shellWallSize[1] - shellThickWall;
    G4VSolid* shellWall = new G4Tubs("ShellWall", shellWallSize.x(), shellWallSize.y(), shellWallSize.z() / 2., 0,
                                     viewDeg);
    G4LogicalVolume* shellWallLV = new G4LogicalVolume(shellWall, AlMat, "ShellWallLV");
    G4ThreeVector shellWallPos = G4ThreeVector(0, 0, 0);
    new G4PVPlacement(nullptr, shellWallPos, shellWallLV, "ShellWallPVP", coreLV, false, 0, true);
    shellWallLV->SetVisAttributes(visShell);

    G4VSolid* shellTop = new G4Tubs("ShellTop", 0, shellWallSize.x(), shellThickTop / 2., 0, viewDeg);
    G4ThreeVector shellTopPos = G4ThreeVector(0, 0, (shellWallSize.z() - shellThickTop) / 2.);
    G4LogicalVolume* shellTopLV = new G4LogicalVolume(shellTop, AlMat, "ShellTopLV");
    new G4PVPlacement(nullptr, shellTopPos, shellTopLV, "ShellTopPVP", coreLV, false, 0, true);
    shellTopLV->SetVisAttributes(visShell);

    G4ThreeVector shellBottomSize = G4ThreeVector(shellWallSize.x(), bottomCapInnerRadius,
                                                  (bottomCapHeight - bottomCapThick) / 2.);
    G4VSolid* shellBottom = new G4Tubs("ShellBottom", shellBottomSize.x(), shellBottomSize.y(), shellBottomSize.z(), 0,
                                       viewDeg);
    G4ThreeVector shellBottomPos = coreBottomPos;
    G4LogicalVolume* shellBottomLV = new G4LogicalVolume(shellBottom, AlMat, "ShellBottomLV");
    new G4PVPlacement(nullptr, shellBottomPos, shellBottomLV, "ShellBottomPVP", coreLV, false, 0, true);
    shellBottomLV->SetVisAttributes(visShell);

    G4ThreeVector shellTabSize = G4ThreeVector(shellWallSize.x() - shellTabLength, shellWallSize.x(),
                                               shellTabHeight / 2.);
    G4VSolid* shellTab = new G4Tubs("ShellTab", shellTabSize.x(), shellTabSize.y(), shellTabSize.z(), 0, viewDeg);
    G4double shellTabDepth = crystalHeight + crystalShellThickTop + tyvekInThickTop + crystalGlassHeight + shellThickTop
        + GasketHeight;
    G4ThreeVector shellTabPos = G4ThreeVector(0, 0, shellWallSize.z() / 2. - shellTabDepth - shellTabHeight / 2.);
    G4LogicalVolume* shellTabLV = new G4LogicalVolume(shellTab, AlMat, "ShellTabLV");
    new G4PVPlacement(nullptr, shellTabPos, shellTabLV, "ShellTabPVP", coreLV, false, 0, true);
    shellTabLV->SetVisAttributes(visAl);

    // Rubber Gasket between Crystal and SiPM Holder
    G4ThreeVector GasketSize = G4ThreeVector(shellWallSize.x() - GasketLength, shellWallSize.x(), GasketHeight / 2.);
    G4VSolid* Gasket = new G4Tubs("GasketTab", GasketSize.x(), GasketSize.y(), GasketSize.z(), 0, viewDeg);
    G4ThreeVector GasketPos = G4ThreeVector(0, 0, shellWallSize.z() / 2. - shellTabDepth + GasketHeight / 2.);
    G4LogicalVolume* GasketLV = new G4LogicalVolume(Gasket, rubberMat, "GasketLV");
    new G4PVPlacement(nullptr, GasketPos, GasketLV, "GasketPVP", coreLV, false, 0, true);
    GasketLV->SetVisAttributes(visAl);
}

void Detector::ConstructBottomVeto() {
    // Bottom veto aluminium shell
    G4ThreeVector bottomVetoShellSize = G4ThreeVector(bottomVetoRadius + tyvekBottomThickWall,
                                                      plateBottomHoleRadius - shellThickWall,
                                                      bottomVetoShellHeight / 2.);
    G4ThreeVector bottomVetoShellPos = G4ThreeVector(0, 0,
                                                     -((coreTopSize.z() - bottomVetoShellHeight) / 2 +
                                                         bottomCapHeight - bottomCapThick - holderHeight));

    G4VSolid* bottomVetoShell = new G4Tubs("BottomVetoShell", bottomVetoShellSize.x(), bottomVetoShellSize.y(),
                                           bottomVetoShellSize.z(), 0, viewDeg);
    G4LogicalVolume* bottomVetoShellLV = new G4LogicalVolume(bottomVetoShell, AlMat, "BottomVetoShellLV");
    new G4PVPlacement(nullptr, bottomVetoShellPos, bottomVetoShellLV, "BottomVetoShellPVP", coreLV, false, 0, true);
    bottomVetoShellLV->SetVisAttributes(visShell);

    G4ThreeVector bottomVetoShellTabSize = G4ThreeVector(bottomVetoShellSize.x() - bottomVetoShellTabLength,
                                                         bottomVetoShellSize.x(), bottomVetoShellTabHeight / 2.);
    G4VSolid* bottomVetoShellTab = new G4Tubs("BottomVetoShellSize", bottomVetoShellTabSize.x(),
                                              bottomVetoShellTabSize.y(), bottomVetoShellTabSize.z(), 0, viewDeg);
    G4ThreeVector bottomVetoShellTabPos = bottomVetoShellPos + G4ThreeVector(
                                                                             0, 0, (bottomVetoShellHeight -
                                                                                 bottomVetoShellTabHeight) / 2. -
                                                                             bottomVetoHeight - tyvekBottomThickTop);
    G4LogicalVolume* bottomVetoShellTabLV = new G4LogicalVolume(bottomVetoShellTab, AlMat, "BottomVetoShellTabLV");
    new G4PVPlacement(nullptr, bottomVetoShellTabPos, bottomVetoShellTabLV, "BottomVetoShellTabPVP", coreLV, false, 0,
                      true);
    bottomVetoShellTabLV->SetVisAttributes(visAl);

    // Bottom Tyvek
    tyvekBottomSize = G4ThreeVector(bottomVetoRadius, bottomVetoRadius + tyvekBottomThickWall,
                                    (bottomVetoHeight + tyvekBottomThickTop) / 2.);
    G4ThreeVector tyvekBottomPos = bottomVetoShellTabPos + G4ThreeVector(0, 0, tyvekBottomSize.z() +
                                                                         bottomVetoShellTabHeight / 2.);

    G4VSolid* tyvekBottomWall = new G4Tubs("TyvekBottomWall", tyvekBottomSize.x(), tyvekBottomSize.y(),
                                           tyvekBottomSize.z(), 0, viewDeg);
    G4VSolid* tyvekBottomTop = new G4Tubs("TyvekBottomTop", 0, tyvekBottomSize.x(), tyvekBottomThickTop / 2., 0,
                                          viewDeg);
    G4ThreeVector tyvekBottomTopPos = G4ThreeVector(0, 0, tyvekBottomSize.z() - tyvekBottomThickTop / 2.0);
    G4VSolid* tyvekBottom = new G4UnionSolid("TyvekBottom", tyvekBottomWall, tyvekBottomTop, nullptr,
                                             tyvekBottomTopPos);

    tyvekBottomLV = new G4LogicalVolume(tyvekBottom, tyvekMat, "TyvekBottomLV");
    tyvekBottomPVP = new G4PVPlacement(nullptr, tyvekBottomPos, tyvekBottomLV, "TyvekBottomPVP", coreLV, false, 0,
                                       true);
    tyvekBottomLV->SetVisAttributes(visTyvekBottom);

    // Bottom Veto
    G4ThreeVector bottomVetoPos = tyvekBottomPos - G4ThreeVector(0, 0, tyvekBottomThickTop / 2.0);
    G4VSolid* bottomVeto = new G4Tubs("BottomVeto", 0., bottomVetoRadius, bottomVetoHeight / 2., 0, viewDeg);
    bottomVetoLV = new G4LogicalVolume(bottomVeto, vetoMat, "BottomVetoLV");
    bottomVetoPVP = new G4PVPlacement(nullptr, bottomVetoPos, bottomVetoLV, "BottomVetoPVP", coreLV, false, 0, true);
    bottomVetoLV->SetVisAttributes(visVeto);

    // Optic Layer for Bottom Veto
    G4VSolid* bottomVetoOpticLayer = new G4Tubs("BottomVetoOpticLayer", 0, bottomVetoShellTabSize.x(),
                                                bottomVetoOpticLayerHeight / 2., 0, viewDeg);
    bottomVetoOpticLayerLV = new G4LogicalVolume(bottomVetoOpticLayer, opticLayerMat, "BottomVetoOpticLayerLV");
    G4ThreeVector bottomVetoOpticLayerPos = bottomVetoShellTabPos + G4ThreeVector(
         0, 0, (bottomVetoShellTabHeight - bottomVetoOpticLayerHeight) / 2.0);
    bottomVetoOpticLayerPVP = new G4PVPlacement(nullptr, bottomVetoOpticLayerPos, bottomVetoOpticLayerLV,
                                                "BottomVetoOpticLayerPVP", coreLV, false, 0, true);
    bottomVetoOpticLayerLV->SetVisAttributes(visOpticLayer);
}


void Detector::ConstructCrystal() {
    // Crystal Container
    crystalContSize = G4ThreeVector(0, crystalRadius + tyvekInThickWall + crystalShellThickWall,
                                    (crystalHeight + tyvekInThickTop + crystalShellThickTop + crystalGlassHeight) / 2);
    G4ThreeVector crystalContPos = G4ThreeVector(0, 0, coreTopSize.z() / 2 - crystalContSize.z() - shellThickTop);

    G4VSolid* crystalCont = new G4Tubs("CrystalContainer", crystalContSize.x(), crystalContSize.y(),
                                       crystalContSize.z(), 0, viewDeg);
    crystalContLV = new G4LogicalVolume(crystalCont, galacticMat, "CrystalContainerLV");
    new G4PVPlacement(nullptr, crystalContPos, crystalContLV, "CrystalContainerPVP", coreLV, false, 0, true);
    crystalContLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Crystal
    G4ThreeVector crystalPos = G4ThreeVector(0, 0, -crystalContSize.z() + crystalGlassHeight + crystalHeight / 2.);
    G4VSolid* crystal = new G4Tubs("Crystal", 0, crystalRadius, crystalHeight / 2., 0, viewDeg);
    crystalLV = new G4LogicalVolume(crystal, CrystalMat, "CrystalLV");
    crystalPVP = new G4PVPlacement(nullptr, crystalPos, crystalLV, "CrystalPVP", crystalContLV, false, 0, true);
    crystalLV->SetVisAttributes(visCrystal);

    // Tyvek In
    tyvekInSize = G4ThreeVector(crystalRadius, crystalRadius + tyvekInThickWall,
                                (crystalHeight + tyvekInThickTop) / 2.);
    G4ThreeVector tyvekInPos = crystalPos + G4ThreeVector(0, 0, tyvekInThickTop / 2.);

    G4VSolid* tyvekInWall = new G4Tubs("TyvekInWall", tyvekInSize.x(), tyvekInSize.y(), tyvekInSize.z(), 0, viewDeg);
    G4VSolid* tyvekInTop = new G4Tubs("TyvekInTop", 0, tyvekInSize.x(), tyvekInThickTop / 2., 0, viewDeg);
    G4ThreeVector tyvekInTopPos = G4ThreeVector(0, 0, tyvekInSize.z() - tyvekInThickTop / 2.0);
    G4VSolid* tyvekIn = new G4UnionSolid("TyvekIn", tyvekInWall, tyvekInTop, nullptr, tyvekInTopPos);

    tyvekInLV = new G4LogicalVolume(tyvekIn, tyvekMat, "TyvekInLV");
    tyvekInPVP = new G4PVPlacement(nullptr, tyvekInPos, tyvekInLV, "TyvekInPVP", crystalContLV, false, 0, true);
    tyvekInLV->SetVisAttributes(visTyvekIn);

    // Construct crystal shell
    G4ThreeVector crystalShellSize = tyvekInSize + G4ThreeVector(tyvekInThickWall, crystalShellThickWall,
                                                                 (crystalShellThickTop + crystalGlassHeight) / 2.);
    G4ThreeVector crystalShellPos = tyvekInPos + G4ThreeVector(0, 0, (crystalShellThickTop - crystalGlassHeight) / 2.);

    G4VSolid* crystalShellWall = new G4Tubs("CrystalShellWall", crystalShellSize.x(), crystalShellSize.y(),
                                            crystalShellSize.z(), 0, viewDeg);
    G4VSolid* crystalShellTop = new G4Tubs("CrystalShellTop", 0, crystalShellSize.x(), crystalShellThickTop / 2., 0,
                                           viewDeg);
    G4ThreeVector crystalShellTopPos = G4ThreeVector(0, 0, crystalShellSize.z() - crystalShellThickTop / 2.0);
    G4VSolid* crystalShell = new G4UnionSolid("CrystalShell", crystalShellWall, crystalShellTop, nullptr,
                                              crystalShellTopPos);

    G4LogicalVolume* crystalShellLV = new G4LogicalVolume(crystalShell, AlMat, "CrystalShellLV");
    new G4PVPlacement(nullptr, crystalShellPos, crystalShellLV, "CrystalShellPVP", crystalContLV, false, 0, true);
    crystalShellLV->SetVisAttributes(visShell);

    // Crystal Glass
    G4VSolid* crystalGlass = new G4Tubs("CrystalGlass", 0, crystalShellSize.x(), crystalGlassHeight / 2., 0, viewDeg);
    crystalGlassLV = new G4LogicalVolume(crystalGlass, glassMat, "CrystalGlassLV");
    G4ThreeVector crystalGlassPos = crystalShellPos + G4ThreeVector(
                                                                    0, 0, -crystalShellSize.z() + crystalGlassHeight /
                                                                    2.0);
    new G4PVPlacement(nullptr, crystalGlassPos, crystalGlassLV, "CrystalGlassPVP", crystalContLV, false, 0, true);
    crystalGlassLV->SetVisAttributes(visGlass);

    // Optic layer for crystal
    G4VSolid* crystalOpticLayer = new G4Tubs("CrystalOpticLayer", 0, crystalSize.y(), crystalOpticLayerHeight / 2., 0,
                                             viewDeg);
    crystalOpticLayerLV = new G4LogicalVolume(crystalOpticLayer, opticLayerMat, "CrystalOpticLayerLV");
    G4ThreeVector crystalOpticLayerPos = crystalContPos + G4ThreeVector(
                                                                        0, 0, -crystalContSize.z() -
                                                                        crystalOpticLayerHeight / 2.0);
    crystalOpticLayerPVP = new G4PVPlacement(nullptr, crystalOpticLayerPos, crystalOpticLayerLV, "CrystalOpticLayerPVP",
                                             coreLV, false, 0, true);
    crystalOpticLayerLV->SetVisAttributes(visOpticLayer);
}

void Detector::ConstructSiPM() {
    // SiPM Frame
    auto* SiPMFrameBase = new G4Box("SiPMFrameBase", SiPMLength / 2., SiPMWidth / 2., SiPMHeight / 2.);
    auto* SiPMHole = new G4Box("SiPMHole", SiPMLength / 2. - SiPMFrameSize, SiPMWidth / 2. - SiPMFrameSize,
                               SiPMHeight / 2. + 5 * mm);

    auto* SiPMFrame = new G4SubtractionSolid("SiPMFrame", SiPMFrameBase, SiPMHole);
    SiPMFrameLV = new G4LogicalVolume(SiPMFrame, SiPMFrameMat, "SiPMFrameLV");
    SiPMFrameLV->SetVisAttributes(visSiPMFrame);

    // SiPM Body
    auto* SiPMBody = new G4Box("SiPMBody", SiPMLength / 2. - SiPMFrameSize, SiPMWidth / 2. - SiPMFrameSize,
                               (SiPMHeight - SiPMWindowThick) / 2.);
    SiPMBodyLV = new G4LogicalVolume(SiPMBody, SiPMMat, "SiPMBodyLV");
    SiPMBodyLV->SetVisAttributes(visSiPM);

    // SiPM Optical Window
    auto* SiPMWindow = new G4Box("SiPMWindow", SiPMLength / 2. - SiPMFrameSize, SiPMWidth / 2. - SiPMFrameSize,
                                 SiPMWindowThick / 2.);
    SiPMWindowLV = new G4LogicalVolume(SiPMWindow, SiPMEncapsulantMat, "SiPMWindowLV");
    SiPMWindowLV->SetVisAttributes(visGlass);

    if (!SiPMPhotocathodeSurf) {
        SiPMPhotocathodeSurf = new G4OpticalSurface("SiPMPhotocathode");
        SiPMPhotocathodeSurf->SetModel(unified);
        SiPMPhotocathodeSurf->SetType(dielectric_metal);
        SiPMPhotocathodeSurf->SetFinish(polished);

        auto* mpt = new G4MaterialPropertiesTable();

        auto pde = Utils::ReadCSV("../OpticalParameters/SiPM_PDE.csv", 1.0, true);
        Utils::NormalizeMaxToOne(pde);
        mpt->AddProperty("EFFICIENCY", pde.E, pde.V, pde.E.size());

        std::vector refl(pde.E.size(), 0.02);
        mpt->AddProperty("REFLECTIVITY", pde.E, refl, pde.E.size());

        SiPMPhotocathodeSurf->SetMaterialPropertiesTable(mpt);
    }
}


void Detector::ConstructHolder(G4ThreeVector& refPos, const G4String& prefix) {
    // Holder body
    G4ThreeVector holderSize = G4ThreeVector(coreTopSize.y() - holderThickWall - shellThickWall,
                                             coreTopSize.y() - shellThickWall, holderHeight / 2.);
    G4VSolid* holderWall = new G4Tubs(prefix + "HolderWall", holderSize.x(), holderSize.y(), holderSize.z(), 0,
                                      viewDeg);
    G4ThreeVector holderWallPos = refPos + G4ThreeVector(0, 0, holderHeight / 2.);
    G4LogicalVolume* holderWallLV = new G4LogicalVolume(holderWall, AlMat, prefix + "HolderWallLV");
    new G4PVPlacement(nullptr, holderWallPos, holderWallLV, prefix + "HolderWallPVP", coreLV, false, 0, true);
    holderWallLV->SetVisAttributes(visHolder);

    G4VSolid* holderBottom = new G4Tubs(prefix + "HolderBottom", 0, holderSize.x(), holderThickBottom / 2., 0, viewDeg);
    G4ThreeVector holderBottomPos = refPos + G4ThreeVector(0, 0, holderThickBottom / 2.);
    G4LogicalVolume* holderBottomLV = new G4LogicalVolume(holderBottom, AlMat, prefix + "HolderBottomLV");
    new G4PVPlacement(nullptr, holderBottomPos, holderBottomLV, prefix + "HolderBottomPVP", coreLV, false, 0, true);
    holderBottomLV->SetVisAttributes(visHolder);

    // Spring holder
    G4VSolid* springHolderBase = new G4Tubs("SpringHolderBase", 0, holderSize.x(), springHolderHeight / 2., 0, viewDeg);
    G4VSolid* cutterX = new G4Box("CutterX", springHolderGapX / 2., holderSize.x() + 5 * mm,
                                  springHolderHeight / 2. + 5 * mm);
    G4VSolid* cutterY = new G4Box("CutterY", holderSize.x() + 5 * mm, springHolderGapY / 2.,
                                  springHolderHeight / 2. + 5 * mm);

    G4VSolid* springHole = new G4Tubs("SpringHole", 0, springRadius, springHolderHeight / 2. + 5 * mm, 0, 360 * deg);
    G4VSolid* tempSolid;
    G4double diff = -((springHolderGapX + springHolderGapY) / 4. + springRadius) + 0.5 * std::sqrt(
         2 * springHoleCenterRadius * springHoleCenterRadius - (springHolderGapX - springHolderGapY) * (springHolderGapX
             - springHolderGapY) / 4.);
    std::vector<std::pair<G4int, G4int>> signs = {
        {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
    };
    for (size_t i = 0; i < static_cast<size_t>(360 * deg / viewDeg); i++) {
        G4ThreeVector pos = G4ThreeVector(springHolderGapX / 2. + springRadius + diff,
                                          springHolderGapY / 2. + springRadius + diff, 0);
        pos[0] *= signs[i].first;
        pos[1] *= signs[i].second;
        tempSolid = new G4SubtractionSolid("TempSolid", springHolderBase, springHole, nullptr, pos);
        springHolderBase = tempSolid;
    }
    G4VSolid* springHolderIncomplete = new G4SubtractionSolid("SpringHolderIncomplete", tempSolid, cutterX);

    G4VSolid* springHolder = new G4SubtractionSolid("SpringHolder", springHolderIncomplete, cutterY);
    G4ThreeVector springHolderPos = refPos + G4ThreeVector(0, 0, holderThickBottom + springHolderHeight / 2.);
    G4LogicalVolume* springHolderLV = new G4LogicalVolume(springHolder, AlMat, prefix + "SpringHolderLV");
    new G4PVPlacement(nullptr, springHolderPos, springHolderLV, prefix + "SpringHolderPVP", coreLV, false, 0, true);
    springHolderLV->SetVisAttributes(visHolder);

    // Board
    G4ThreeVector boardPos = refPos + G4ThreeVector(0, 0, holderThickBottom + springLength + boardHeight / 2.);
    G4VSolid* board = new G4Tubs(prefix + "Board", 0, holderSize.x(), boardHeight / 2., 0, viewDeg);
    G4LogicalVolume* boardLV = new G4LogicalVolume(board, boardMat, prefix + "BoardLV");
    new G4PVPlacement(nullptr, boardPos, boardLV, prefix + "BoardPVP", coreLV, false, 0, true);
    boardLV->SetVisAttributes(visBoard);

    // Payload
    G4ThreeVector payloadPos = refPos + G4ThreeVector(0, 0, holderThickBottom + springLength - payloadHeight / 2.);
    G4VSolid* payload = new G4Box(prefix + "Payload", payloadLength / 2., payloadWidth / 2., payloadHeight / 2.);
    G4LogicalVolume* payloadLV = new G4LogicalVolume(payload, payloadMat, prefix + "PayloadLV");
    new G4PVPlacement(nullptr, payloadPos, payloadLV, prefix + "PayloadPVP", coreLV, false, 0, true);
    payloadLV->SetVisAttributes(visPayload);
}


void Detector::ConstructCrystalSiPM() {
    // Holder
    G4ThreeVector SiPMHolderPos = G4ThreeVector(0, 0,
                                                coreTopSize.z() / 2. - (shellThickTop + crystalHeight +
                                                    tyvekInThickTop + crystalShellThickTop +
                                                    crystalGlassHeight + crystalOpticLayerHeight + SiPMHeight +
                                                    boardHeight + springLength +
                                                    holderThickBottom));
    ConstructHolder(SiPMHolderPos, "CrystalSiPM");

    // SiPM Container
    G4ThreeVector SiPMContPos = SiPMHolderPos + G4ThreeVector(0, 0,
                                                              holderThickBottom + springLength + boardHeight +
                                                              SiPMHeight / 2.);
    G4VSolid* SiPMCont = new G4Tubs("CrystalSiPMContainer", 0, coreTopSize.y() - holderThickWall - shellThickWall,
                                    SiPMHeight / 2., 0, 360 * deg);
    G4LogicalVolume* SiPMContLV = new G4LogicalVolume(SiPMCont, galacticMat, "CrystalSiPMContainerLV");
    new G4PVPlacement(nullptr, SiPMContPos, SiPMContLV, "CrystalSiPMContainerPVP", coreLV, false, 0, true);
    SiPMContLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Crystal SiPM
    std::vector countVec = {2, 4, 4, 2}; // {2, 2};
    G4int copyN = 0;

    for (int i = 0; i < countVec.size(); i++) {
        for (int j = 0; j < countVec[i]; j++) {
            G4ThreeVector SiPMPos((countVec.size() / 2 - 0.5 - i) * /*(*/SiPMLength,               // + crystalSiPMDist)
                                  (countVec[i] / 2. - 1) * SiPMWidth - (j - 0.5) * /*(*/SiPMWidth, // + crystalSiPMDist)
                                  0);

            new G4PVPlacement(nullptr, SiPMPos, SiPMFrameLV,
                              "CrystalSiPMFramePVP", SiPMContLV, false, copyN, true);

            auto* bodyPVP = new G4PVPlacement(nullptr,
                                              SiPMPos + G4ThreeVector(0, 0, -SiPMWindowThick / 2.),
                                              SiPMBodyLV,
                                              "CrystalSiPMBodyPVP", SiPMContLV, false, copyN, true);

            auto* windowPVP = new G4PVPlacement(nullptr,
                                                SiPMPos + G4ThreeVector(0, 0, (SiPMHeight - SiPMWindowThick) / 2.),
                                                SiPMWindowLV,
                                                "CrystalSiPMWindowPVP", SiPMContLV, false, copyN, true);

            new G4LogicalBorderSurface("CrystalSiPM_Photocathode_" + std::to_string(copyN),
                                       windowPVP, bodyPVP, SiPMPhotocathodeSurf);

            copyN++;
        }
    }

    // G4int crystalEdgeSiPMCount = crystalSiPMCount - copyN;
    // G4double crystalSiPMRadius = crystalRadius -
    //     std::ceil(0.5 * std::sqrt(SiPMWidth * SiPMWidth + SiPMLength * SiPMLength)) * mm;
    // for (size_t i = 0; i < crystalEdgeSiPMCount; i++) {
    //     auto* rotMat = new G4RotationMatrix(i * 360 * deg / crystalEdgeSiPMCount, 0, 0);
    //     G4ThreeVector SiPMPos(crystalSiPMRadius * std::cos(i * 360 * deg / crystalEdgeSiPMCount),
    //                           crystalSiPMRadius * std::sin(i * 360 * deg / crystalEdgeSiPMCount),
    //                           0);
    //     new G4PVPlacement(rotMat, SiPMPos, SiPMFrameLV, "CrystalSiPMFramePVP", SiPMContLV, false, copyN + i, true);
    //
    //     auto* bodyPVP = new G4PVPlacement(rotMat,
    //                                       SiPMPos + G4ThreeVector(0, 0, -SiPMWindowThick / 2.),
    //                                       SiPMBodyLV, "CrystalSiPMBodyPVP", SiPMContLV, false, copyN + i, true);
    //
    //     auto* windowPVP = new G4PVPlacement(rotMat,
    //                                         SiPMPos + G4ThreeVector(0, 0, (SiPMHeight - SiPMWindowThick) / 2.),
    //                                         SiPMWindowLV, "CrystalSiPMWindowPVP", SiPMContLV, false, copyN + i, true);
    //
    //     new G4LogicalBorderSurface("CrystalSiPM_Photocathode_" + std::to_string(i),
    //                                windowPVP, bodyPVP, SiPMPhotocathodeSurf);
    // }
}

void Detector::ConstructVetoSiPM() {
    // Holder springs
    G4int springNumber = 4;
    G4ThreeVector refPos = G4ThreeVector(0, 0, -detContainerTopSize.z() / 2.);
    G4ThreeVector vetoSpringSize = G4ThreeVector(detContainerTopSize.y() - vetoSpringWidth, detContainerTopSize.y(),
                                                 detContainerTopSize.z() - tyvekOutThickTop - vetoHeight -
                                                 vetoOpticLayerHeight - SiPMHeight - boardHeight);
    G4double vetoSpringGap = vetoSpringSize.z() / springNumber - vetoSpringHeight;
    G4VSolid* vetoSpring = new G4Tubs("VetoSpringPart", vetoSpringSize.x(), vetoSpringSize.y(), vetoSpringHeight / 2.,
                                      0, viewDeg);
    G4LogicalVolume* vetoSpringLV = new G4LogicalVolume(vetoSpring, AlMat, "VetoSpringLV");
    for (size_t i = 0; i < springNumber; i++) {
        G4ThreeVector vetoSpringPos = refPos + G4ThreeVector(0, 0,
                                                             vetoSpringGap * (i + 1) + vetoSpringHeight * (i + 0.5));
        new G4PVPlacement(nullptr, vetoSpringPos, vetoSpringLV, "VetoSpringPVP", detContainerLV, false, i, true);
        vetoSpringLV->SetVisAttributes(visSpring);
    }

    // Board
    G4ThreeVector vetoBoardPos = refPos + G4ThreeVector(0, 0,
                                                        (vetoSpringGap + vetoSpringHeight) * springNumber +
                                                        boardHeight / 2.);
    G4VSolid* vetoBoard = new G4Tubs("VetoBoard", coreTopSize.y(), detContainerTopSize.y(), boardHeight / 2., 0,
                                     360 * deg);
    G4LogicalVolume* vetoBoardLV = new G4LogicalVolume(vetoBoard, galacticMat, "VetoBoardLV");
    new G4PVPlacement(nullptr, vetoBoardPos, vetoBoardLV, "VetoBoardPVP", detContainerLV, false, 0, true);
    vetoBoardLV->SetVisAttributes(visBoard);

    // Veto Payload
    G4ThreeVector vetoPayloadPos = vetoBoardPos + G4ThreeVector(0, 0, -(vetoPayloadHeight + boardHeight) / 2.);
    G4int vetoPayloadNumber = 4;
    G4double vetoPayloadPhase = 22.5 * deg;
    G4double vetoPayloadRadius = (detContainerTopSize.y() + coreTopSize.y()) / 2.;
    G4VSolid* vetoPayload = new G4Box("VetoPayload", vetoPayloadLength / 2., vetoPayloadWidth / 2.,
                                      vetoPayloadHeight / 2);
    G4LogicalVolume* vetoPayloadLV = new G4LogicalVolume(vetoPayload, payloadMat, "VetoPayloadLV");
    vetoPayloadLV->SetVisAttributes(visPayload);
    for (size_t i = 0; i < vetoPayloadNumber; i++) {
        G4RotationMatrix* rotMat = new G4RotationMatrix(vetoPayloadPhase + i * 360 * deg / vetoPayloadNumber, 0, 0);
        G4ThreeVector payloadPos = vetoPayloadPos + G4ThreeVector(
                                                                  vetoPayloadRadius *
                                                                  std::cos(vetoPayloadPhase + i * 360 * deg /
                                                                           vetoPayloadNumber),
                                                                  vetoPayloadRadius *
                                                                  std::sin(vetoPayloadPhase + i * 360 * deg /
                                                                           vetoPayloadNumber), 0);
        new G4PVPlacement(rotMat, payloadPos, vetoPayloadLV, "VetoPayloadPVP", detContainerLV, false, i, true);
    }

    // SiPM Container
    G4ThreeVector SiPMContPos = refPos + G4ThreeVector(0, 0,
                                                       (vetoSpringGap + vetoSpringHeight) * springNumber + boardHeight +
                                                       SiPMHeight / 2.);
    G4VSolid* SiPMCont = new G4Tubs("VetoSiPMContainer", coreTopSize.y(), detContainerTopSize.y(), SiPMHeight / 2., 0,
                                    360 * deg);
    G4LogicalVolume* SiPMContLV = new G4LogicalVolume(SiPMCont, galacticMat, "VetoSiPMContainerLV");
    new G4PVPlacement(nullptr, SiPMContPos, SiPMContLV, "VetoSiPMContainerPVP", detContainerLV, false, 0, true);
    SiPMContLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Veto SiPM
    G4double vetoSiPMRadius = (detContainerTopSize.y() + coreTopSize.y()) / 2.;
    for (size_t i = 0; i < vetoSiPMCount; i++) {
        auto* rotMat = new G4RotationMatrix(i * 360 * deg / vetoSiPMCount, 0, 0);
        G4ThreeVector SiPMPos(vetoSiPMRadius * std::cos(i * 360 * deg / vetoSiPMCount),
                              vetoSiPMRadius * std::sin(i * 360 * deg / vetoSiPMCount),
                              0);

        new G4PVPlacement(rotMat, SiPMPos, SiPMFrameLV, "VetoSiPMFramePVP", SiPMContLV, false, i, true);

        auto* bodyPVP = new G4PVPlacement(rotMat,
                                          SiPMPos + G4ThreeVector(0, 0, -SiPMWindowThick / 2.),
                                          SiPMBodyLV, "VetoSiPMBodyPVP", SiPMContLV, false, i, true);

        auto* windowPVP = new G4PVPlacement(rotMat,
                                            SiPMPos + G4ThreeVector(0, 0, (SiPMHeight - SiPMWindowThick) / 2.),
                                            SiPMWindowLV, "VetoSiPMWindowPVP", SiPMContLV, false, i, true);

        new G4LogicalBorderSurface("VetoSiPM_Photocathode_" + std::to_string(i),
                                   windowPVP, bodyPVP, SiPMPhotocathodeSurf);
    }
}

void Detector::ConstructBottomVetoSiPM() {
    // Holder
    G4ThreeVector SiPMHolderPos = G4ThreeVector(
                                                0, 0, -coreTopSize.z() / 2. - bottomCapHeight + bottomCapThick);
    ConstructHolder(SiPMHolderPos, "BottomVetoSiPM");

    // SiPM Container
    G4ThreeVector SiPMContPos = SiPMHolderPos + G4ThreeVector(0, 0,
                                                              holderThickBottom + springLength + boardHeight +
                                                              SiPMHeight / 2.);
    G4VSolid* SiPMCont = new G4Tubs("BottomVetoSiPMContainer", 0, coreTopSize.y() - holderThickWall - shellThickWall,
                                    SiPMHeight / 2., 0,
                                    360 * deg);
    G4LogicalVolume* SiPMContLV = new G4LogicalVolume(SiPMCont, galacticMat, "BottomVetoSiPMContainerLV");
    new G4PVPlacement(nullptr, SiPMContPos, SiPMContLV, "BottomVetoSiPMContainerPVP", coreLV, false, 0, true);
    SiPMContLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Bottom Veto SiPM
    G4double bottomVetoSiPMRadius = (coreTopSize.y() - holderThickWall - shellThickWall) / 2.;
    for (size_t i = 0; i < bottomVetoSiPMCount; i++) {
        auto* rotMat = new G4RotationMatrix(i * 360 * deg / bottomVetoSiPMCount, 0, 0);
        G4ThreeVector SiPMPos(bottomVetoSiPMRadius * std::cos(i * 360 * deg / bottomVetoSiPMCount),
                              bottomVetoSiPMRadius * std::sin(i * 360 * deg / bottomVetoSiPMCount),
                              0);

        new G4PVPlacement(rotMat, SiPMPos, SiPMFrameLV, "BottomVetoSiPMFramePVP", SiPMContLV, false, i, true);

        auto* bodyPVP = new G4PVPlacement(rotMat,
                                          SiPMPos + G4ThreeVector(0, 0, -SiPMWindowThick / 2.),
                                          SiPMBodyLV, "BottomVetoSiPMBodyPVP", SiPMContLV, false, i, true);

        auto* windowPVP = new G4PVPlacement(rotMat,
                                            SiPMPos + G4ThreeVector(0, 0, (SiPMHeight - SiPMWindowThick) / 2.),
                                            SiPMWindowLV, "BottomVetoSiPMWindowPVP", SiPMContLV, false, i, true);

        new G4LogicalBorderSurface("BottomVetoSiPM_Photocathode_" + std::to_string(i),
                                   windowPVP, bodyPVP, SiPMPhotocathodeSurf);
    }
}

std::vector<G4LogicalVolume*> Detector::GetSensitiveLV() const {
    return {
        crystalLV,
        vetoLV,
        bottomVetoLV,
        tyvekOutLV,
        tyvekMidLV,
        tyvekInLV,
        tyvekBottomLV,
    };
}

void Detector::AddBorderSurface(const G4String& name,
                                G4VPhysicalVolume* pvFrom,
                                G4VPhysicalVolume* pvTo,
                                G4OpticalSurface* surf) {
    if (!pvFrom || !pvTo || !surf) return;
    new G4LogicalBorderSurface(name, pvFrom, pvTo, surf);
}

void Detector::AddBidirectionalBorder(const G4String& nameAToB,
                                      const G4String& nameBToA,
                                      G4VPhysicalVolume* pvA,
                                      G4VPhysicalVolume* pvB,
                                      G4OpticalSurface* surf) {
    if (!pvA || !pvB || !surf) return;
    new G4LogicalBorderSurface(nameAToB, pvA, pvB, surf);
    new G4LogicalBorderSurface(nameBToA, pvB, pvA, surf);
}

void Detector::ConstructOpticalSurfaces() {
    auto* tyvekSurf = new G4OpticalSurface("TyvekSurface");
    tyvekSurf->SetModel(unified);
    // tyvekSurf->SetType(dielectric_metal);
    // tyvekSurf->SetFinish(polished);
    // tyvekSurf->SetSigmaAlpha(0.0);
    tyvekSurf->SetType(dielectric_dielectric);
    tyvekSurf->SetFinish(groundfrontpainted);
    tyvekSurf->SetSigmaAlpha(0.2);

    const auto refl = Utils::ReadCSV("../OpticalParameters/Tyvek_reflectivity.csv", 1.0, true);
    auto* mpt = new G4MaterialPropertiesTable();
    mpt->AddProperty("REFLECTIVITY", refl.E, refl.V, refl.E.size());

    std::vector E = {1.0 * eV, 4.0 * eV};
    std::vector spike = {0.0, 0.0};
    std::vector lobe = {0.02, 0.02};
    std::vector back = {0.0, 0.0};

    mpt->AddProperty("SPECULARSPIKECONSTANT", E.data(), spike.data(), static_cast<G4int>(E.size()), true);
    mpt->AddProperty("SPECULARLOBECONSTANT", E.data(), lobe.data(), static_cast<G4int>(E.size()), true);
    mpt->AddProperty("BACKSCATTERCONSTANT", E.data(), back.data(), static_cast<G4int>(E.size()), true);

    tyvekSurf->SetMaterialPropertiesTable(mpt);

    // 1) Crystal <-> TyvekIn
    AddBidirectionalBorder("CrystalToTyvekIn", "TyvekInToCrystal", crystalPVP, tyvekInPVP, tyvekSurf);

    // 2) Veto <-> TyvekMid
    AddBidirectionalBorder("VetoToTyvekMid", "TyvekMidToVeto", vetoPVP, tyvekMidPVP, tyvekSurf);
    // 2) Veto <-> TyvekOut
    AddBidirectionalBorder("VetoToTyvekOut", "TyvekOutToVeto", vetoPVP, tyvekOutPVP, tyvekSurf);

    // 2) Bottom Veto <-> TyvekBottom
    AddBidirectionalBorder("BottomVetoToTyvekBottom", "TyvekBottomToBottomVeto", bottomVetoPVP, tyvekBottomPVP,
                           tyvekSurf);
}
