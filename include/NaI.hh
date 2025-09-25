#ifndef NAI_HH
#define NAI_HH

#include "Detector.hh"


class NaI : public Detector {
public:
    NaI(G4LogicalVolume *, const G4ThreeVector &, G4NistManager *, const Sizes &, G4double, G4bool);

    void DefineMaterials() override;
    void Construct() override;
    [[nodiscard]] std::vector<G4LogicalVolume *> GetSensitiveLV() const override;

private:
    G4double gapSize;
    G4double LEDSize;

    G4double viewDeg;
    G4bool doubleLED;

    G4double tyvekThick;
    G4double shellThick;
    G4double vetoThick;
    G4ThreeVector NaISize;

    G4Material *vetoMat;
    G4Material *NaIMat;
    G4Material *tyvekInMat;
    G4Material *tyvekOutMat;
    G4Material *LEDMat;
    G4Material *shellMat;

    G4LogicalVolume *NaILV;
    G4LogicalVolume *tyvekOutLV;
    G4LogicalVolume *tyvekInLV;
    G4LogicalVolume *vetoLV;
    G4LogicalVolume *shellLV;

    G4Tubs *hole;

    void ConstructShell();
    void ConstructLED();
};

#endif //NAI_HH
