#ifndef SVETLANA2_0_ENGINE_HH
#define SVETLANA2_0_ENGINE_HH

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

#include "geometry/EngineSizes.hh"

class Engine {
public:
    Engine(G4LogicalVolume* world, G4NistManager* nistManager);
    ~Engine() = default;

    void ConstructEngine();

    G4LogicalVolume* GetEngineLV() {
        return engineContainerLV;
    }

private:
    G4LogicalVolume* worldLV{};
    G4NistManager* nist{};

    G4Material* boardMat{};
    G4Material* magneticMat{};
    G4Material* capacitorMat{};
    G4Material* vacuumMat{};
    G4Material* electronicsMat{};
    G4Material* AlMat{};

    G4VSolid* boardTop{};
    G4VSolid* boardBottom{};
    G4VSolid* capacitor{};
    G4VSolid* fieldCoilCentral{};
    G4VSolid* electronics{};
    G4VSolid* fieldCoil{};

    G4LogicalVolume* engineContainerLV{};
    G4LogicalVolume* boardTopLV{};
    G4LogicalVolume* boardBottomLV{};
    G4LogicalVolume* capacitorLV{};
    G4LogicalVolume* fieldCoilCentralLV{};
    G4LogicalVolume* fieldCoilCentralCuLV{};
    G4LogicalVolume* electronicsLV{};
    G4LogicalVolume* fieldCoilLV{};

    G4VisAttributes* visBoardTop{};
    G4VisAttributes* visBoardBottom{};
    G4VisAttributes* visMagnetic{};
    G4VisAttributes* visCapacitor{};
    G4VisAttributes* visAl{};
    G4VisAttributes* visElectronics{};

    void DefineMaterial();
    void DefineVisual();

    void ConstructBoardTop();
    void ConstructBoardBottom();
    void ConstructCapacitor();
    void ConstructFieldCoilCentral();
    void ConstructFieldCoil();
    void ConstructElectronics();
};

#endif //SVETLANA2_0_ENGINE_HH
