#ifndef CUBESATSIZES_HH
#define CUBESATSIZES_HH

#include <G4SystemOfUnits.hh>
#include <algorithm>

#include "CameraSizes.hh"

namespace CubeSatSizes
{
    namespace UnitFrame
    {
        const G4double halfX = 4.25 * mm;
        const G4double halfY = 4.25 * mm;
        const G4double halfZ = 49.8 * mm;

        const G4double ledgeX = 42.2 * mm;
        const G4double ledgeY = 34.2 * mm;
        const G4double ledgeHeight = 0.75 * mm;
        const G4double ledgeThickness = 0.3 * mm;

        const G4double thickness = 0.6 * mm;
    }


    namespace Unit
    {
        const G4double halfX = 50 * mm;
        const G4double halfY = 50 * mm;
        const G4double halfZ = UnitFrame::halfZ;
    }

    namespace Bridge
    {
        const G4double thickness = 0.6 * mm;
        const G4double length = Unit::halfX;
        const G4double width = 3.8 * mm;

        const G4double bracketX = 10.5 * mm;
        const G4double bracketY = 1.9 * mm;

        const G4double fixLedgeX = 4.25 * mm;

        const G4double fixLedgeHeight = 3.0 * mm;
        const G4double fixLedgeLength = 3.05 * mm;
        const G4double fixLedgeWidth = 3.2 * mm;

        const G4double ledgeThickness = 0.6 * mm;
    }

    namespace Lintel
    {
        const G4double thickness = 0.6 * mm;
        const G4double length = 41.3 * mm;
        const G4double width = 3.8 * mm;

        const G4double bracketX = 3.35 * mm;
        const G4double bracketY = 8.55 * mm;

        const G4double ledgeThickness = 0.75 * mm;
        const G4double ledgeHeight = 1.0 * mm;

        const G4double fixLedgeHeight = 3.1 * mm;
        const G4double fixLedgeXThickness = 0.5 * mm;
        const G4double fixLedgeZThickness = 0.9 * mm;

        const G4double gap = 0.1 * mm;
    }

    namespace Frame
    {
        const G4double gap = 5.425 * mm;
        const G4double bracketHeight = 5.0 * mm;
        const G4double thickness = 1.0 * mm;
    }

    namespace CubeSat
    {
        const G4double halfX = 55.4 * mm;
        const G4double halfY = 55.4 * mm;
        const G4double halfZ = 170.25 * mm;
    }

    namespace Fixator
    {
        const G4double length = 11.5 * mm;
        const G4double width = 12.05 * mm;
        const G4double thickness = 0.85 * mm;

        const G4double holeX = 5.5 * mm;
        const G4double holeY = 5.75 * mm;
        const G4double cutX = 1.6 * mm;
        const G4double cutY = 5.55 * mm;

        const G4double recessGapX = 0.6 * mm;
        const G4double recessGapY = 0.7 * mm;
        const G4double recessDeep = 0.5 * mm;
        const G4double recessX = 2.4 * mm;
        const G4double recessY = 5.0 * mm;
    }

    namespace SolarPanel
    {
        const G4double thickness = 0.295 * mm;
        const G4double deep = 0.505 * mm;
        const G4double gap = 0.1 * mm;
    }

    namespace Holder
    {
        const G4double halfX = 41.0 * mm;
        const G4double halfY = 3.0 * mm;
        const G4double halfZ = 48.5 * mm;

        const G4double ledgeX = 28.8 * mm;
        const G4double ledgeY = 1.5 * mm;
        const G4double ledgeZ = 3.25 * mm;

        const G4double ledgeThickness = 0.3 * mm;
        const G4double ledgeBackThickness = 0.6 * mm;

        const G4double recess1Deep = 0.9 * mm;
        const G4double recess1Ledge = 2.5 * mm;
        const G4double recess1Radius = 15.0 * mm;

        const G4double recess2X = 39.4 * mm;
        const G4double recess2Y = 0.6 * mm;
        const G4double recess2Z = 40.5 * mm;

        const G4double recessBackDeep = 0.9 * mm;
        const G4double recessBackLedgeX = 6.485 * mm;
        const G4double recessBackLedgeY = 3.25 * mm;
        const G4double recessBackGap = 12.835 * mm;

        const G4double holeGapZ = 16.75 * mm;
        const G4double holeX1 = 9.0 * mm;
        const G4double holeZ1 = 5.75 * mm;
        const G4double holeX2 = 7.0 * mm;
        const G4double holeZ2 = 4.5 * mm;
        const G4double holeZ = 12.05 * mm;

        const G4double stripThickness = 3.0 * mm;
        const G4double stripX = 27.0 * mm;
        const G4double stripY = 27.0 * mm;
        const G4double stripGap = 4.0 * mm;
        const G4double stripRadius = 6.0 * mm;
    }

    namespace Mechanics
    {
        const G4double halfX = 48.0 * mm;
        const G4double halfY = 41.0 * mm;
        const G4double halfZ = 44.7 * mm;

        const G4double ledgeX = 40.0 * mm;
        const G4double ledgeY = 34.5 * mm;
        const G4double ledgeZ = 10.0 * mm;

        const G4double cut1X = 3.5 * mm;
        const G4double cut1Y = 12.25 * mm;
        const G4double cut1Z = 2.5 * mm;

        const G4double cut2X = 4.0 * mm;
        const G4double cut2Y = 12.25 * mm;
        const G4double cut2Z = 1.25 * mm;

        const G4double channalWidth = 4.0 * mm;

        const G4double serviceSystemGap = 19.795 * mm;
    }

    namespace Cover
    {
        const G4double halfX = 48.0 * mm;
        const G4double halfY = 32.0 * mm;
        const G4double halfZ = 6.0 * mm;

        const G4double ledgeX = 3.5 * mm;
        const G4double ledgeY = 4.0 * mm;
        const G4double ledgeZ = 0.75 * mm;

        const G4double holeX = 5.25 * mm;
        const G4double holeY = 18.0 * mm;

        const G4double stripThickness = 0.75 * mm;
        const G4double stripHeigth = 1.5 * mm;

        const G4double bottomRecessDeep = 2.5 * mm;
        const G4double bottomRecessLength = 39.5 * mm;
        const G4double bottomRecessWidth = 30.425 * mm;

        const G4double cut1X = 37.5 * mm;
        const G4double cut1Y = 5.0 * mm;
        const G4double cut1Z = 1.0 * mm;

        const G4double cut2X = 28.0 * mm;
        const G4double cut2Y = 20.5 * mm;
        const G4double cut2Z = 1.5 * mm;
    }

    namespace Boards
    {
        const G4double half1X = 48.0 * mm;
        const G4double half1Y = 32.0 * mm;
        const G4double half1Z = 0.765 * mm;

        const G4double half2X = 27.5 * mm;
        const G4double half2Y = 20.0 * mm;
        const G4double half2Z = 0.54 * mm;

        const G4double gap1X = 0.2 * mm;
        const G4double gap2X = 0.3 * mm;
        const G4double gapY = 0.5 * mm;
        const G4double gapZ = 0.505 * mm;
    }

    // const G4double displacement = CameraSizes::CameraBox::halfZ - Frame::bracketHeight;
    const G4double displacement = 0.0 * mm;
}

#endif // CUBESATSIZES_HH
