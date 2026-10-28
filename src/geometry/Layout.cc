#include "geometry/Layout.hh"

#include <algorithm>
#include <vector>

#include <G4Exception.hh>
#include <G4Point3D.hh>
#include <G4SystemOfUnits.hh>

#include "Configuration.hh"
#include "geometry/Sizes.hh"
#include "geometry/CubeSatSizes.hh"
#include "geometry/EngineSizes.hh"

using namespace Configuration;


namespace {
    // NATALIA (TunaCan centre) above the satellite centre
    // cubesat: upside down in the central unit
    const G4double centralUnitZ = 3 * mm;

    // tunacan: on top of the upper unit, upright (position taken over from the tunacan branch)
    G4double TopZ() {
        return 3 * CubeSatSizes::UnitFrame::halfZ - 12 * mm + Sizes::modelHeight - CubeSatSizes::Frame::bracketHeight;
    }

    G4bool Known(const G4String& value, const std::vector<G4String>& list) {
        return std::find(list.begin(), list.end(), value) != list.end();
    }

    G4String Join(const std::vector<G4String>& list) {
        G4String out;
        for (const auto& s : list) out += (out.empty() ? "" : ", ") + s;
        return out;
    }
}


void Layout::Configure() {
    const std::vector<G4String> geometries = {"default", "tunacan", "cubesat", "experiment"};
    const std::vector<G4String> engines = {"sputnix", "geoscan"};
    const std::vector<G4String> cameras = {"galper"};

    if (!Known(geometryConfig, geometries)) {
        G4Exception("Layout::Configure", "GeometryConfig", FatalException,
                    ("Geometry configuration not found: " + geometryConfig +
                        ".\nAvailable configurations: " + Join(geometries)).c_str());
    }
    if (!engineConfig.empty() && !Known(engineConfig, engines)) {
        G4Exception("Layout::Configure", "EngineConfig", FatalException,
                    ("Engine configuration not found: " + engineConfig +
                        ".\nAvailable configurations: " + Join(engines)).c_str());
    }
    if (!cameraConfig.empty() && !Known(cameraConfig, cameras)) {
        G4Exception("Layout::Configure", "CameraConfig", FatalException,
                    ("Camera configuration not found: " + cameraConfig +
                        ".\nAvailable configurations: " + Join(cameras)).c_str());
    }
    if (!HasSatellite() && (!engineConfig.empty() || !cameraConfig.empty())) {
        G4Exception("Layout::Configure", "NoSatellite", FatalException,
                    ("--engine and --camera need a satellite: -gc tunacan or -gc cubesat (given: " +
                        geometryConfig + ")").c_str());
    }
    if (!cameraConfig.empty()) {
        // Camera.cc builds the old model; its placement in the central unit (looking sideways) is not defined yet
        G4Exception("Layout::Configure", "CameraNotPlaced", FatalException,
                    "The camera is not placed yet: define its position in Layout and the placement in "
                    "Geometry::Construct, then extend SatelliteExtentZ if it sticks out of the satellite");
    }

    if (HasSatellite()) {
        // Mounting plate of the satellite version
        Sizes::plateSize = 78 * mm;
        Sizes::plateCornerSize = 10 * mm;
    }
}


G4bool Layout::HasSatellite() {
    return geometryConfig == "tunacan" || geometryConfig == "cubesat";
}


G4bool Layout::HasTunaCan() {
    return geometryConfig != "experiment";
}


G4bool Layout::BareCrystal() {
    return geometryConfig == "experiment";
}


G4Transform3D Layout::NataliaPlacement() {
    if (geometryConfig == "cubesat") {
        return G4Translate3D(0, 0, centralUnitZ) * G4RotateX3D(180 * deg);
    }
    if (geometryConfig == "tunacan") {
        return G4Translate3D(0, 0, TopZ());
    }
    return G4Transform3D::Identity;
}


G4ThreeVector Layout::EnginePosition() {
    // Both presets: lower part of the upper unit. Replace with the real positions when they are known.
    if (engineConfig == "sputnix") return {0, 0, EngineSizes::engineContainer::halfZ + 8 * cm};
    if (engineConfig == "geoscan") return {0, 0, EngineSizes::engineContainer::halfZ + 8 * cm};
    return {};
}


void Layout::SatelliteExtentZ(G4double& zMin, G4double& zMax) {
    zMin = -CubeSatSizes::CubeSat::halfZ;
    zMax = CubeSatSizes::CubeSat::halfZ;
    if (geometryConfig == "tunacan") zMax = std::max(zMax, TopZ() + Sizes::modelHeight / 2.0);
}


G4double Layout::SatelliteTopZ() {
    return CubeSatSizes::CubeSat::halfZ - CubeSatSizes::Frame::bracketHeight;
}


G4ThreeVector Layout::IsotropicCentre() {
    const G4Point3D local(0, 0, SatelliteTopZ() - TopZ());
    const G4Point3D world = NataliaPlacement() * local;
    return {world.x(), world.y(), world.z()};
}


G4String Layout::Tag() {
    G4String tag = geometryConfig == "default" ? "" : geometryConfig;
    if (!engineConfig.empty()) tag += "_engine-" + engineConfig;
    if (!cameraConfig.empty()) tag += "_camera-" + cameraConfig;
    return tag;
}
