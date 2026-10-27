#ifndef LAYOUT_HH
#define LAYOUT_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>

// What is built and where it stands, for the chosen -gc/--engine/--camera options.
// Geometry places the volumes and GenSurface builds the generation envelope from the same numbers,
// so a new position goes here and nowhere else.
namespace Layout {
    // Checks the option names and their combination, applies the configuration-dependent sizes.
    // Must be called once after argv parsing, before anything reads Sizes or builds a GenSurface.
    void Configure();

    G4bool HasSatellite();   // tunacan, cubesat
    G4bool HasTunaCan();     // TunaCan, Ti plate and bottom cap around NATALIA
    G4bool BareCrystal();    // experiment: crystal and its SiPMs only, no veto and no shells

    // NATALIA local frame (TunaCan centre at the origin, +z towards the TunaCan top) -> world
    G4Transform3D NataliaPlacement();
    G4ThreeVector EnginePosition();

    // Satellite configurations: world z range occupied by the whole model
    void SatelliteExtentZ(G4double& zMin, G4double& zMax);
    // Top of the satellite body without the brackets: what the TunaCan plate stands on
    G4double SatelliteTopZ();

    // Suffix for output names; empty for the default geometry, so its file names do not change
    G4String Tag();
}

#endif //LAYOUT_HH
