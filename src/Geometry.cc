#include "Geometry.hh"

#include <sstream>

#include <CLHEP/Random/MixMaxRng.h>
#include <G4Point3D.hh>
#include <Randomize.hh>

using namespace Sizes;
using namespace Configuration;


Geometry::Geometry() {
    std::vector<G4String> detectorList = {"NaI", "CsI"};
    if (std::find(detectorList.begin(), detectorList.end(), detectorType) == detectorList.end()) {
        G4Exception("Geometry::ConstructDetector", "DetectorType", FatalException,
                    ("Detector not found: " + detectorType + ".\nAvailable detectors: NaI, CsI").c_str());
    }
    std::vector<G4String> crystalSiPMConfigList = {
        "2x2", "16-cross", "13-circle", "12-cross", "12-circle", "12-rhombus", "8-circle"
    };
    if (std::find(crystalSiPMConfigList.begin(), crystalSiPMConfigList.end(), crystalSiPMConfig) ==
        crystalSiPMConfigList.end()) {
        G4Exception("Geometry::ConstructDetector", "CrystalSiPMConfiguration", FatalException,
                    ("SiPM configuration not found: " + crystalSiPMConfig +
                        ".\nAvailable configurations: 2x2, 16-cross, 12-cross, 12-circle").c_str());
    }

    nist = G4NistManager::Instance();
    zeroRot = new G4RotationMatrix(0, 0, 0);

    detContainerSize = G4ThreeVector(0 * mm,
                                     modelRadius - tunaCanThickWall,
                                     (modelHeight - bottomCapThick - tunaCanThickTop) / 2.0);
    detContainerPos = G4ThreeVector(0, 0, -(plateCenterThick + tunaCanThickTop) / 2);

    worldHalfSize = std::max({modelRadius * 2, modelHeight}) * 2;
    const GenSurface genSurface = GenSurface::For(fluxDirection);
    worldHalfSize = std::max(worldHalfSize, genSurface.OuterRadius() + genSurface.Origin().mag() + GenSurface::Margin());
    worldHalfSize = std::max(worldHalfSize, GenSurface::BoundingRadius()
                                            + GenSurface::PayloadCentre().mag() + GenSurface::Margin());

    detContVisAttr = new G4VisAttributes(G4VisAttributes::GetInvisible());
    detContVisAttr->SetForceSolid(false);
}

void Geometry::ConstructDetector() {
    G4VSolid* detContTopTube = new G4Tubs("DetContTopTube", 0, detContainerSize.y(),
                                          (modelHeight - tunaCanThickTop + plateCenterThick) / 2, 0,
                                          360 * deg);
    G4VSolid* detContMidTube = new G4Tubs("DetContMidTube", 0, plateBottomHoleRadius, plateThick / 2., 0,
                                          360 * deg);
    G4VSolid* detContBottomTube = new G4Tubs("DetContBottomTube", 0, bottomCapInnerRadius,
                                             (bottomCapHeight - bottomCapThick) / 2., 0,
                                             360 * deg);
    G4ThreeVector detContMidPos = G4ThreeVector(
                                                0, 0, -(modelHeight - tunaCanThickTop + plateCenterThick + plateThick) /
                                                2);
    G4ThreeVector detContBottomPos = G4ThreeVector(
                                                   0, 0, -(modelHeight - tunaCanThickTop + plateCenterThick +
                                                       bottomCapHeight - bottomCapThick) / 2 - plateThick);
    G4VSolid* detContIncomplete = new G4UnionSolid("DetContIncomplete", detContTopTube, detContMidTube, zeroRot,
                                                   detContMidPos);

    detContainer = new G4UnionSolid("DetectorContainer", detContIncomplete, detContBottomTube, zeroRot,
                                    detContBottomPos);
    detContainerLV = new G4LogicalVolume(detContainer, worldMat, "DetectorContainerLV");
    detContainerPVPL = new G4PVPlacement(Layout::NataliaPlacement() * G4Translate3D(detContainerPos), detContainerLV,
                                         "DetectorContainerPVPL", worldLV, false, 0, true);
    detContainerLV->SetVisAttributes(detContVisAttr);

    detector = new Detector(detContainerLV, nist);
    detector->Construct();
    std::vector<G4LogicalVolume*> sensitiveLV = detector->GetSensitiveLV();
    crystalLV = sensitiveLV.at(0);
    vetoLV = sensitiveLV.at(1);
    bottomVetoLV = sensitiveLV.at(2);
    tyvekOutLV = sensitiveLV.at(3);
    tyvekMidLV = sensitiveLV.at(4);
    tyvekInLV = sensitiveLV.at(5);
    tyvekBottomLV = sensitiveLV.at(6);
}


G4VPhysicalVolume* Geometry::Construct() {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    worldMat = nist->FindOrBuildMaterial("G4_Galactic");

    auto* mptVac = new G4MaterialPropertiesTable();
    const G4int N = 2;
    G4double E[N] = {2.0 * eV, 3.4 * eV};
    G4double n1[N] = {1.0, 1.0};
    mptVac->AddProperty("RINDEX", E, n1, N);
    worldMat->SetMaterialPropertiesTable(mptVac);

    worldBox = new G4Box("World", worldHalfSize, worldHalfSize, worldHalfSize);
    worldLV = new G4LogicalVolume(worldBox, worldMat, "WorldLV");
    worldPVP = new G4PVPlacement(zeroRot, G4ThreeVector(0, 0, 0), worldLV, "WorldPVPL", nullptr, false, 0, false);
    worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

    if (Layout::HasSatellite()) {
        CubeSat_GeoScan_3U cubeSat(worldLV, nist);
        cubeSat.ConstructCubeSat();
    }
    if (!engineConfig.empty()) {
        Engine engine(worldLV, nist);
        engine.ConstructEngine();
        new G4PVPlacement(nullptr, Layout::EnginePosition(), engine.GetEngineLV(), "EnginePVP", worldLV, false, 0, true);
    }

    ConstructDetector();
    if (Layout::HasTunaCan()) {
        TunaCan tunaCan(worldLV, nist, Layout::NataliaPlacement(), detContainerSize.y());
        tunaCan.Construct();
    }

    CheckEnvelope();

    return worldPVP;
}

void Geometry::CheckEnvelope() const {
    // Every top-level volume has to lie inside the generation envelope (or the flat / isotropic normalisation
    // misses it) and inside BoundingRadius minus the margin (or primaries can start in material).
    // Surface points are drawn from a private engine, so the run's random sequence stays untouched.
    CLHEP::MixMaxRng localEngine;
    CLHEP::HepRandomEngine* runEngine = G4Random::getTheEngine();
    G4Random::setTheEngine(&localEngine);

    const G4ThreeVector centre = GenSurface::PayloadCentre();
    const G4double rMax = GenSurface::BoundingRadius() - GenSurface::Margin();
    constexpr G4double tolerance = 1e-3 * mm;
    constexpr G4int nPoints = 2000;

    std::ostringstream problems;
    for (std::size_t i = 0; i < worldLV->GetNoDaughters(); ++i) {
        const G4VPhysicalVolume* pv = worldLV->GetDaughter(i);
        const G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
        const G4Transform3D toWorld(pv->GetObjectRotationValue(), pv->GetObjectTranslation());
        const G4bool needEnvelope = !GenSurface::AllowedOutsideEnvelope(pv->GetName());

        auto outside = [&](const G4ThreeVector& local, G4ThreeVector& world) {
            const G4Point3D p = toWorld * G4Point3D(local);
            world = G4ThreeVector(p.x(), p.y(), p.z());
            return (world - centre).mag() > rMax + tolerance || (needEnvelope && !GenSurface::InEnvelope(world,
                tolerance));
        };

        G4ThreeVector pMin, pMax, world;
        solid->BoundingLimits(pMin, pMax);
        G4bool cornersInside = true;
        for (G4int c = 0; c < 8 && cornersInside; ++c) {
            const G4ThreeVector corner(c & 1 ? pMax.x() : pMin.x(), c & 2 ? pMax.y() : pMin.y(),
                                       c & 4 ? pMax.z() : pMin.z());
            cornersInside = !outside(corner, world);
        }
        if (cornersInside) continue;

        for (G4int n = 0; n < nPoints; ++n) {
            if (!outside(solid->GetPointOnSurface(), world)) continue;
            problems << "\n  " << pv->GetName() << ": point (" << world.x() / mm << ", " << world.y() / mm << ", "
                << world.z() / mm << ") mm";
            break;
        }
    }

    G4Random::setTheEngine(runEngine);

    if (!problems.str().empty()) {
        G4Exception("Geometry::CheckEnvelope", "OutsideEnvelope", FatalException,
                    ("Volumes outside the generation envelope of -gc " + geometryConfig +
                        " (update GenSurface::PayloadEnvelope / Layout::SatelliteExtentZ):" + problems.str()).c_str());
    }
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
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    auto* detectorSD = new SensitiveDetector("DetectorSD", 0, "Crystal");
    sdManager->AddNewDetector(detectorSD);
    crystalLV->SetSensitiveDetector(detectorSD);

    auto* vetoSD = new SensitiveDetector("VetoSD", 1, "Veto");
    sdManager->AddNewDetector(vetoSD);
    vetoLV->SetSensitiveDetector(vetoSD);

    auto* bottomVetoSD = new SensitiveDetector("BottomVetoSD", 2, "BottomVeto");
    sdManager->AddNewDetector(bottomVetoSD);
    bottomVetoLV->SetSensitiveDetector(bottomVetoSD);

    if (useOptics) {
        auto* sipmSD = new SiPMOpticalSD("SiPMOpticalSD");
        auto* sipmWindowLV = detector->GetSiPMWindowLV();
        sipmSD->SetSiPMWindowLV(sipmWindowLV);
        sdManager->AddNewDetector(sipmSD);
        sipmWindowLV->SetSensitiveDetector(sipmSD);
    }
}
