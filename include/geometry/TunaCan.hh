#ifndef TUNACAN_HH
#define TUNACAN_HH

#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4Transform3D.hh>
#include <G4VisAttributes.hh>

// Al TunaCan around the detector container, the Ti mounting plate and the bottom cap.
// Built in the NATALIA local frame (TunaCan centre at the origin) and placed with one transform.
class TunaCan {
public:
    TunaCan(G4LogicalVolume* mother, G4NistManager* nistManager, const G4Transform3D& placement,
            G4double innerRadius);
    ~TunaCan() = default;

    void Construct();

private:
    G4LogicalVolume* motherLV;
    G4NistManager* nist;
    G4Transform3D placement;
    G4double innerRadius;

    G4VisAttributes* tunaCanVisAttr;

    void Place(G4LogicalVolume* lv, const G4ThreeVector& pos, const G4String& name) const;
};

#endif //TUNACAN_HH
