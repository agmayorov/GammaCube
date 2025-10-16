#ifndef SIZES_HH
#define SIZES_HH

#include <G4Types.hh>
#include <G4SystemOfUnits.hh>

namespace  Sizes {
    inline G4double modelHeight{78 * mm};
    inline G4double modelRadius{40 * mm};

    inline G4double tunaCanThickWall{2 * mm};
    inline G4double tunaCanThickTop{6 * mm};
    inline G4double tunaCanThickBottom{0 * mm};

    inline G4double gapSizeWall{2 * mm};
    inline G4double gapSizeTop{2 * mm};
    inline G4double gapSizeBottom{2 * mm};

    inline G4double tyvekOutThickWall{1 * mm};
    inline G4double tyvekOutThickTop{1 * mm};
    inline G4double tyvekOutThickBottom{1 * mm};

    inline G4double vetoThickWall{7 * mm};
    inline G4double vetoThickTop{7 * mm};
    inline G4double vetoThickBottom{7 * mm};

    inline G4double vetoChamferHeigh{0 * mm};

    inline G4double tyvekMidThickWall{1 * mm};
    inline G4double tyvekMidThickTop{1 * mm};
    inline G4double tyvekMidThickBottom{1 * mm};

    inline G4double rubberRadius{5 * mm};
    inline G4double rubberHeight{5 * mm};

    inline G4double AlThickWall{0.5 * mm};
    inline G4double AlThickTop{0.5 * mm};

    inline G4double AlCapThickWall{0.5 * mm};
    inline G4double AlCapThickBottom{0.5 * mm};

    inline G4double tyvekInThickWall{1 * mm};
    inline G4double tyvekInThickTop{1 * mm};
    inline G4double tyvekInThickBottom{1 * mm};

    inline G4int crystalLEDCount{4};
    inline G4double crystalLEDRadius{3 * mm};
    inline G4double crystalLEDHeight{1 * mm};

    inline G4int vetoLEDCount{4};
    inline G4double vetoLEDRadius{3 * mm};
    inline G4double vetoLEDHeight{3 * mm};

    inline G4int vetoBottomLEDCount{4};
    inline G4double vetoBottomLEDRadius{3 * mm};
    inline G4double vetoBottomLEDHeight{3 * mm};

    inline G4double wireRadius{1 * mm};
    inline G4double wireInsulationThick{0.2 * mm};

    inline G4int pinCount{4};
    inline G4double pinRadius{0.5 * mm};

    inline G4double boardHeight{5 * mm};
    inline G4double boardLength{5 * mm};
    inline G4double boardWidth{5 * mm};

    inline G4double boardSpace{5 * mm};

    inline G4double tunaCanAllSize{5 * mm};
    inline G4double tyvekInAllSize{5 * mm};
    inline G4double tyvekOutAllSize{5 * mm};
    inline G4double tyvekMidAllSize{5 * mm};
    inline G4double vetoAllSize{5 * mm};
    inline G4double AlAllSize{5 * mm};
    inline G4double rubberAllSize{5 * mm};
    inline G4double crystalLEDAllSize{5 * mm};
    inline G4double vetoLEDAllSize{5 * mm};
    inline G4double vetoBottomLEDAllSize{5 * mm};

};

#endif //SIZES_HH