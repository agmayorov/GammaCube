#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <G4VisAttributes.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4UnionSolid.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4Version.hh>
#include <G4LossTableManager.hh>
#include <G4Element.hh>
#include <G4Box.hh>
#include <G4VSolid.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>

#include "Sizes.hh"


class Detector {
public:
    G4NistManager *nist;

    G4String detectorType;

    G4LogicalVolume *detContainerLV;
    G4ThreeVector detContainerSize;

    G4VisAttributes *visBlue;
    G4VisAttributes* visWhite;
    G4VisAttributes *visRed;
    G4VisAttributes *visCyan;
    G4VisAttributes *visGrey;

    Detector() = default;

    virtual ~Detector() = default;

    void init();

    virtual void DefineMaterials() = 0;
    virtual void DefineVisual();

    virtual void Construct() = 0;

    [[nodiscard]] G4String GetDetectorType() const;
    [[nodiscard]] virtual std::vector<G4LogicalVolume *>GetSensitiveLV() const = 0;
};


#endif //DETECTOR_HH
