#ifndef CUBESAT_HH
#define CUBESAT_HH

#include <G4VisAttributes.hh>
#include <G4QuadrangularFacet.hh>
#include <G4Polyhedron.hh>
#include <G4TessellatedSolid.hh>
#include <G4TriangularFacet.hh>
#include <G4SubtractionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4UniformRandPool.hh>
#include <G4UnionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Element.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Trd.hh>
#include <G4LogicalVolume.hh>
// #include <G4GDMLParser.hh>
#include <G4VUserDetectorConstruction.hh>

#include "geometry/CubeSatSizes.hh"
#include "geometry/CameraSizes.hh"
#include "geometry/Camera.hh"

class CubeSat_GeoScan_3U {
public:
    CubeSat_GeoScan_3U(G4LogicalVolume* world, G4NistManager* nistManager);
    ~CubeSat_GeoScan_3U() = default;

    void ConstructCubeSat();

    G4LogicalVolume* GetCubeSatLV() {return cubeSatLV;}

private:
    G4LogicalVolume* worldLV{};
    G4NistManager* nist{};

    G4Material* frameMat{};
    G4Material* lintelMat{};
    G4Material* holderMat{};
    G4Material* coverMat{};
    G4Material* solarPanelMat{};
    G4Material* boardMat{};
    G4Material* mechanicsMat{};
    G4Material* vacuumMat{};

    G4VSolid* solarPanel{};
    G4VSolid* frame{};
    G4VSolid* cubeSat{};
    G4VSolid* lintel{};
    G4VSolid* fixator{};
    G4VSolid* serviceSystem{};
    G4VSolid* cover{};
    G4VSolid* mechanics{};
    G4VSolid* unitFrame{};
    G4VSolid* unit{};
    G4VSolid* bridge{};

    G4LogicalVolume* solarPanelHolderLV{};
    G4LogicalVolume* solarPanelLV{};
    G4LogicalVolume* frameLV{};
    G4LogicalVolume* cubeSatLV{};
    G4LogicalVolume* lintelLV{};
    G4LogicalVolume* serviceSystemLV{};
    G4LogicalVolume* coverLV{};
    G4LogicalVolume* mechanicsLV{};
    G4LogicalVolume* unitLV{};
    G4LogicalVolume* fixatorLV{};
    G4LogicalVolume* unitFrameLV{};
    G4LogicalVolume* bridgeLV{};
    G4LogicalVolume* cameraLV{};

    G4VisAttributes* visSolarPanel{};
    G4VisAttributes* visHolder{};
    G4VisAttributes* visFrame{};
    G4VisAttributes* visCover{};
    G4VisAttributes* visBoard{};
    G4VisAttributes* visMechanics{};
    G4VisAttributes* visServiceSystem{};

    void DefineMaterial();
    void DefineVisual();

    void ConstructBridge();
    void ConstructLintel();
    void ConstructFixator();
    void ConstructUnitFrame();
    void ConstructFrame();
    void ConstructSolarPanelHolder();
    void ConstructSolarPanel();
    void ConstructCover();
    void ConstructMechanics();
    void ConstructServiceSystem();
};


#endif // CUBESAT_HH
