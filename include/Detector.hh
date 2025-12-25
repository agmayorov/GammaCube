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
#include <G4Cons.hh>
#include <G4VSolid.hh>
#include <G4Trd.hh>
#include <G4Tubs.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>

#include "Sizes.hh"


class Detector {
public:
    G4NistManager* nist;
    G4String detectorType;

    G4ThreeVector detContainerTopSize;

    G4LogicalVolume* detContainerLV;
    G4VPhysicalVolume* detContainerPVP;
    G4double detContainerSizeZ;

    G4VisAttributes* visTyvekOut{};
    G4VisAttributes* visTyvekMid{};
    G4VisAttributes* visTyvekIn{};
    G4VisAttributes* visTyvekBottom{};
    G4VisAttributes* visCrystal{};
    G4VisAttributes* visVeto{};
    G4VisAttributes* visSpring{};
    G4VisAttributes* visSiPM{};
    G4VisAttributes* visSiPMFrame{};
    G4VisAttributes* visAl{};
    G4VisAttributes* visShell{};
    G4VisAttributes* visHolder{};
    G4VisAttributes* visGlass{};
    G4VisAttributes* visOpticLayer{};
    G4VisAttributes* visBottomVetoShell{};
    G4VisAttributes* visGasket{};
    G4VisAttributes* visBoard{};
    G4VisAttributes* visPayload{};
    G4VisAttributes* visRubber{};

    Detector(G4LogicalVolume*, G4NistManager*, G4double, G4int, const G4String&);
    ~Detector() = default;

    void DefineMaterials();
    void DefineVisual();

    void Construct();

    [[nodiscard]] G4String GetDetectorType() const;
    [[nodiscard]] std::vector<G4LogicalVolume*> GetSensitiveLV() const;

    G4LogicalVolume* GetSiPMWindowLV() const { return SiPMWindowLV; }

private:
    G4double viewDeg;
    G4int yieldScale;

    G4ThreeVector crystalSize;
    G4ThreeVector vetoSize;
    G4ThreeVector tyvekInSize;
    G4ThreeVector tyvekMidSize;
    G4ThreeVector tyvekOutSize;
    G4ThreeVector tyvekBottomSize;
    G4ThreeVector coreTopSize;
    G4ThreeVector crystalContSize;

    G4Material* vetoMat{};
    G4Material* CrystalMat{};
    G4Material* tyvekMat{};
    G4Material* SiPMMat{};
    G4Material* SiPMFrameMat{};
    G4Material* AlMat{};
    G4Material* rubberMat{};
    G4Material* boardMat{};
    G4Material* payloadMat{};
    G4Material* glassMat{};
    G4Material* SiPMGlassMat{};
    G4Material* opticLayerMat{};
    G4Material* galacticMat{};

    G4LogicalVolume* crystalLV{};
    G4LogicalVolume* coreLV{};
    G4LogicalVolume* crystalGlassLV;
    G4LogicalVolume* crystalContLV{};
    G4LogicalVolume* tyvekInLV{};
    G4LogicalVolume* tyvekOutLV{};
    G4LogicalVolume* tyvekMidLV{};
    G4LogicalVolume* tyvekBottomLV{};
    G4LogicalVolume* vetoLV{};
    G4LogicalVolume* bottomVetoLV{};
    G4LogicalVolume* SiPMBodyLV{};
    G4LogicalVolume* SiPMWindowLV{};
    G4LogicalVolume* SiPMFrameLV{};

    G4LogicalVolume* crystalOpticLayerLV{};
    G4LogicalVolume* vetoOpticLayerLV{};
    G4LogicalVolume* bottomVetoOpticLayerLV{};

    G4OpticalSurface* SiPMPhotocathodeSurf = nullptr;

    G4VPhysicalVolume* crystalPVP{};
    G4VPhysicalVolume* tyvekInPVP{};
    G4VPhysicalVolume* AlPVP{};
    G4VPhysicalVolume* rubberPVP{};
    G4VPhysicalVolume* tyvekMidPVP{};
    G4VPhysicalVolume* tyvekBottomPVP{};
    G4VPhysicalVolume* vetoPVP{};
    G4VPhysicalVolume* bottomVetoPVP{};
    G4VPhysicalVolume* tyvekOutPVP{};
    G4VPhysicalVolume* crystalOpticLayerPVP{};
    G4VPhysicalVolume* vetoOpticLayerPVP{};
    G4VPhysicalVolume* bottomVetoOpticLayerPVP{};

    void ConstructCrystal();
    void ConstructShell();
    void ConstructVeto();
    void ConstructBottomVeto();
    void ConstructHolder(G4ThreeVector&, const G4String&);

    void ConstructSiPM();
    void ConstructCrystalSiPM();
    void ConstructVetoSiPM();
    void ConstructBottomVetoSiPM();

    void AddBorderSurface(const G4String& name,
                          G4VPhysicalVolume* pvFrom,
                          G4VPhysicalVolume* pvTo,
                          G4OpticalSurface* surf);

    void AddBidirectionalBorder(const G4String& nameAToB,
                                const G4String& nameBToA,
                                G4VPhysicalVolume* pvA,
                                G4VPhysicalVolume* pvB,
                                G4OpticalSurface* surf);
    void ConstructOpticalSurfaces();
};


#endif //DETECTOR_HH
