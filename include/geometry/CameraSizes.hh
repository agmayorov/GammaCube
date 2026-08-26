#ifndef CAMERASIZES_HH
#define CAMERASIZES_HH

namespace CameraSizes
{
    namespace CameraBox
    {
        const G4double halfX = 49.4 * mm;
        const G4double halfY = 49.4 * mm;
        const G4double halfZ = 25.4 * mm;

        const G4double bottomHalfX = 30.0 * mm;
        const G4double bottomHalfY = 30.0 * mm;
        const G4double bottomHalfZ = 8.25 * mm;

        const G4double ledge = 4.25 * mm;
    }

    namespace Lenses
    {
        const G4double curvitureRadiusTop = 28.57 * mm;
        const G4double radiusTop = 14.6 * mm;
        const G4double thicknessTop = 0.25 * mm;
        const G4double offsetZTop = 0.27 * mm;

        const G4double radiusBottom = 5.2 * mm;
        const G4double thicknessBottom = 0.25 * mm;
        const G4double offsetZBottom = 0.97 * mm;

        const G4double lenseOffsetZ = 0.79 * mm;
    }

    namespace LenseShellTop
    {
        const G4double radius = 20.0 * mm;
        const G4double height = 2.275 * mm;

        const G4double hole1Deep = 0.1 * mm;

        const G4double hole2Deep = 1.55 * mm;
        const G4double hole2Radius = 17.02 * mm;

        const G4double hole3Deep = 0.075 * mm;
        const G4double hole3Radius = 17.50 * mm;

        const G4double hole4Deep = 0.15 * mm;
        const G4double hole4Radius = 17.80 * mm;

        const G4double hole5Deep = 0.1 * mm;
        const G4double hole5Radius = 17.85 * mm;

        const G4double hole6Deep = 0.3 * mm;
        const G4double hole6Radius = 18.45 * mm;

        const G4double notchRadius = 19.6 * mm;
        const G4double notchHeight = 0.4 * mm;

        const G4double delta = 0.05 * mm;
    }

    namespace LenseShellBody
    {
        const G4double radius1 = 20.5 * mm;
        const G4double height1 = 3.85 * mm;

        const G4double radius2 = 14.25 * mm;
        const G4double height2 = 7.15 * mm;

        const G4double radius3 = 10.5 * mm;
        const G4double height3 = 0.5 * mm;

        const G4double radius4 = 11.0 * mm;
        const G4double height4 = 4.505 * mm;

        const G4double radius5 = 11.5 * mm;
        const G4double height5 = 3.1 * mm;

        const G4double radius6 = 15.0 * mm;
        const G4double height6 = 0.65 * mm;

        const G4double radius7 = 12.0 * mm;
        const G4double height7 = 0.4 * mm;

        const G4double radius8 = 12.7 * mm;
        const G4double height8 = 1.6 * mm;

        const G4double thickness = 0.25 * mm;

        const G4double wholeHeight = height1 + height2 + height3 + height4 + height5 + height6 + height7 + height8;

        const G4double retainerRadius = 16.9 * mm;
        const G4double retainerHeight = 2.9 * mm;
    }

    namespace LenseShellBottom
    {
        const G4double radius = 11.02 * mm;
        const G4double height = 0.545 * mm;

        const G4double hole1Deep = 0.16 * mm;

        const G4double hole2Deep = 0.3 * mm;
        const G4double hole2Radius = 6.24 * mm;

        const G4double hole3Deep = 0.085 * mm;
        const G4double hole3Radius = 8.6 * mm;

        const G4double hole4Radius = 9.05 * mm;
    }

    namespace LenseHolder
    {
        const G4double radius = 20.25 * mm;
        const G4double height = 2.5 * mm;
        const G4double length = 22.5 * mm;
        const G4double widht = 22.5 * mm;

        const G4double legHeight = 19.9 * mm;
        const G4double legRadiusInner = 1.75 * mm;
        const G4double legRadiusOuter = 2.25 * mm;

        const G4double holeOffsetX = 4.5 * mm;
        const G4double holeOffsetY = 4.5 * mm;
    }

    namespace CameraPlate
    {
        const G4double halfX = 49.4 * mm;
        const G4double halfY = 49.4 * mm;
        const G4double halfZ = 3.0 * mm;

        const G4double ledge = 4.25 * mm;

        const G4double lenseHolderHolesDeep = 1.0 * mm;
        const G4double electronicsHolderHolesDeep = 2.0 * mm;

        const G4double holeRadius = 12.7 * mm;

        const G4double recessRadius = 28.6 * mm;
        const G4double recessDeep = 1.0 * mm;
    }

    namespace ElectronicsHolder
    {
        const G4double height = 1.5 * mm;
        const G4double length = 30.0 * mm;
        const G4double widht = 30.0 * mm;

        const G4double legHeight = 8.75 * mm;
        const G4double legRadiusInner = 1.5 * mm;
        const G4double legRadiusOuter = 2.0 * mm;

        const G4double holeOffsetX = 5.0 * mm;
        const G4double holeOffsetY = 5.0 * mm;
    }

    namespace Electronics
    {
        const G4double height = 0.51 * mm;
        const G4double outerRadius = 28.58 * mm;
        const G4double innerRadius = 12.7 * mm;

        const G4double pentagonSide = 26.75 * mm;
        const G4double pentagonRadius = 22.75 * mm;
        const G4double pentagonHeight = 0.445 * mm;

        const G4double matrixHeight = 0.125 * mm;

        const G4double offsetZ = 0.79 * mm;
    }
}

#endif //CAMERASIZES_HH
