#ifndef CAMERA_HH
#define CAMERA_HH

#include <G4VisAttributes.hh>
#include <G4QuadrangularFacet.hh>
#include <G4Polyhedron.hh>
#include <G4Polyhedra.hh>
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
#include <G4Cons.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <G4Trd.hh>
#include <G4LogicalVolume.hh>
// #include <G4GDMLParser.hh>
#include <G4VUserDetectorConstruction.hh>

#include "geometry/CameraSizes.hh"

class Camera {
public:
    Camera(G4LogicalVolume* world, G4NistManager* nistManager);
    ~Camera() = default;

    void ConstructCamera();

    G4LogicalVolume* GetCameraLV() {
        return cameraLV;
    }

private:
    G4LogicalVolume* worldLV{};
    G4NistManager* nist{};

    G4Material* lenseMat{};
    G4Material* lenseShellMat{};
    G4Material* retainerMat{};
    G4Material* lenseHolderMat{};
    G4Material* electronicsMat{};
    G4Material* matrixMat{};
    G4Material* electronicsHolderMat{};
    G4Material* cameraPlateMat{};
    G4Material* vacuumMat{};

    G4VSolid* lenseTop{};
    G4VSolid* lenseBottom{};
    G4VSolid* lenseShell{};
    G4VSolid* retainer{};
    G4VSolid* lenseHolder{};
    G4VSolid* electronicsHolder{};
    G4VSolid* electronics{};
    G4VSolid* matrix{};
    G4VSolid* cameraPlate{};

    G4LogicalVolume* cameraLV{};
    G4LogicalVolume* lenseTopLV{};
    G4LogicalVolume* lenseBottomLV{};
    G4LogicalVolume* lenseShellLV{};
    G4LogicalVolume* retainerLV{};
    G4LogicalVolume* lenseHolderLV{};
    G4LogicalVolume* electronicsHolderLV{};
    G4LogicalVolume* electronicsLV{};
    G4LogicalVolume* matrixLV{};
    G4LogicalVolume* cameraPlateLV{};

    G4VisAttributes* visLense{};
    G4VisAttributes* visLenseShell{};
    G4VisAttributes* visRetainer{};
    G4VisAttributes* visLenseHolder{};
    G4VisAttributes* visElectronics{};
    G4VisAttributes* visMatrix{};
    G4VisAttributes* visElectronicsHolder{};
    G4VisAttributes* visCameraPlate{};

    void DefineMaterial();
    void DefineVisual();

    void ConstructLense();
    void ConstructLenseHolder();
    void ConstructCameraPlate();
    void ConstructElectronics();
    void ConstructElectronicsHolder();
};

#endif //CAMERA_HH
