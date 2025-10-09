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
#include <G4MultiUnion.hh>
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

    G4VisAttributes *visBlue{};
    G4VisAttributes *visWhite{};
    G4VisAttributes *visRed{};
    G4VisAttributes *visCyan{};
    G4VisAttributes *visGrey{};

    Detector(G4LogicalVolume *, const G4ThreeVector &, G4NistManager *, const Sizes &, G4double, G4int,
             const G4String &);
    ~Detector() = default;

    void DefineMaterials();
    void DefineVisual();

    void Construct();

    [[nodiscard]] G4String GetDetectorType() const;
    [[nodiscard]] std::vector<G4LogicalVolume *>GetSensitiveLV() const;

private:
    G4MultiUnion *LEDs;
    G4double gapSize;
    G4double LEDSize;

    G4double viewDeg;
    G4int nLED;

    G4double tyvekThick;
    G4double shellThick;
    G4double vetoThick;
    G4ThreeVector crystalSize;

    G4Material *vetoMat{};
    G4Material *CrystalMat{};
    G4Material *tyvekInMat{};
    G4Material *tyvekOutMat{};
    G4Material *LEDMat{};
    G4Material *shellMat{};

    G4LogicalVolume *crystalLV{};
    G4LogicalVolume *tyvekOutLV{};
    G4LogicalVolume *tyvekInLV{};
    G4LogicalVolume *vetoLV{};
    G4LogicalVolume *shellLV{};

    G4Tubs *hole{};

    void ConstructShell();
    void ConstructLED();
};


#endif //DETECTOR_HH
