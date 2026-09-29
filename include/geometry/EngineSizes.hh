#ifndef SVETLANA2_0_ENGINESIZES_HH
#define SVETLANA2_0_ENGINESIZES_HH

namespace EngineSizes
{
    namespace engineContainer
    {
        const G4double halfX = 48.0 * mm;
        const G4double halfY = 48.0 * mm;
        const G4double halfZ = 26.1 * mm;
    }
    namespace BoardBottom
    {
        const G4double halfX = 48.0 * mm;
        const G4double halfY = 48.0 * mm;
        const G4double halfZ = 1.0 * mm;

        const G4double cutoutHalfZ = halfZ;
        const G4double cornerCutoutHalfX = 4.0 * mm;
        const G4double cornerCutoutHalfY = 5.5 * mm;
        const G4double edgeCutoutHalfX = 9.0 * mm;
        const G4double edgeCutoutHalfY = 1.5 * mm;
    }

    namespace Capacitor
    {
        const G4double halfX = 21.05 * mm;
        const G4double halfY = 12.2 * mm;
        const G4double halfZ = 22.1 * mm;
    }

    namespace BoardTop
    {
        const G4double halfX = 0.765 * mm;
        const G4double halfY = 16.0 * mm;
        const G4double halfZ = 22.0 * mm;
    }

    namespace FieldCoilCentral
    {
        const G4double radius = 16.0 * mm;
        const G4double halfHeight = 22.0 * mm;
    }

    namespace FieldCoil
    {
        const G4double outerRadius = 11.5 * mm;
        const G4double innerRadius = 6.5 * mm;
        const G4double halfHeight = 3.5 * mm;
        const G4double wire = 0.25 * mm;
    }

    namespace Electronics
    {
        const G4double halfX = 1.65 * mm;
        const G4double halfY = 12.7 * mm;
        const G4double halfZ = 12.2 * mm;
    }


}

#endif //SVETLANA2_0_ENGINESIZES_HH
