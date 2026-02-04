#ifndef SIZES_HH
#define SIZES_HH

#include <G4Types.hh>
#include <G4SystemOfUnits.hh>

namespace  Sizes {
    inline G4double modelHeight{63 * mm};
    inline G4double modelRadius{37 * mm};

    inline G4double tunaCanThickWall{5 * mm};
    inline G4double tunaCanThickTop{5 * mm};

    inline G4double plateSize{83 * mm};
    inline G4double plateCornerSize{8.5 * mm};
    inline G4double plateThick{3 * mm};

    inline G4double plateOuterHoleRadius{37 * mm};
    inline G4double plateInnerHoleRadius{32 * mm};
    inline G4double plateBottomHoleRadius{23 * mm};
    inline G4double plateCenterThick{4.1 * mm};

    inline G4double bottomCapInnerRadius{25 * mm};
    inline G4double bottomCapHeight{6.3 * mm};
    inline G4double bottomCapThick{1.9 * mm};

    inline G4double tyvekOutThickWall{1 * mm};
    inline G4double tyvekOutThickTop{1 * mm};

    inline G4double tyvekBottomThickWall{1 * mm};
    inline G4double tyvekBottomThickTop{1 * mm};

    inline G4double vetoHeight{53.5 * mm};
    inline G4double vetoRadius{31 * mm};
    inline G4double vetoThickWall{7 * mm};
    inline G4double vetoThickTop{7 * mm};

    inline G4double vetoChamferHeight{0 * mm};
    inline G4double vetoTopRoundedRadius{0 * mm};

    inline G4double vetoOpticLayerHeight{1 * mm};

    inline G4double bottomVetoRadius{18.5 * mm};
    inline G4double bottomVetoHeight{5 * mm};

    inline G4double supportLedgeHoleRadius{2.25 * mm};
    inline G4double supportLedgeHeight{2 * mm};

    inline G4double supportLedgeDistanceX{19.75 * mm};
    inline G4double supportLedgeDistanceY{12.75 * mm};

    inline G4double bottomVetoShellHeight{8 * mm};

    inline G4double bottomVetoShellTabLength{1 * mm};
    inline G4double bottomVetoShellTabHeight{1 * mm};

    inline G4double bottomVetoOpticLayerHeight{1 * mm};

    inline G4double tyvekMidThickWall{1 * mm};
    inline G4double tyvekMidThickTop{1 * mm};

    inline G4double rubberHeight{0.5 * mm};

    inline G4double shellThickWall{2 * mm};
    inline G4double shellThickTop{2 * mm};

    inline G4double shellTabLength{1.5 * mm};
    inline G4double shellTabHeight{1.5 * mm};

    inline G4double GasketHeight{0.5 * mm};
    inline G4double GasketLength{2 * mm};

    inline G4double crystalRadius{19 * mm};
    inline G4double crystalHeight{27 * mm};

    inline G4double tyvekInThickWall{1 * mm};
    inline G4double tyvekInThickTop{1 * mm};

    inline G4double crystalShellThickWall{1 * mm};
    inline G4double crystalShellThickTop{1 * mm};

    inline G4double crystalGlassHeight{1 * mm};
    inline G4double crystalOpticLayerHeight{1 * mm};

    inline G4double holderThickWall{1.5 * mm};
    inline G4double holderThickBottom{2 * mm};
    inline G4double holderHeight{9 * mm};

    inline G4double springHolderHeight{2 * mm};
    inline G4double springHolderGapX{12.75 * mm};
    inline G4double springHolderGapY{19.75 * mm};

    inline G4double springLength{5.4 * mm};
    inline G4double springRadius{2.25 * mm};
    inline G4double springHoleCenterRadius{16 * mm};

    inline G4double payloadHeight{2 * mm};
    inline G4double payloadLength{28 * mm};
    inline G4double payloadWidth{19.5 * mm};

    inline G4double vetoPayloadHeight{2 * mm};
    inline G4double vetoPayloadLength{7 * mm};
    inline G4double vetoPayloadWidth{7 * mm};

    inline G4double vetoSpringHeight{0.5 * mm};
    inline G4double vetoSpringWidth{0.5 * mm};

    inline G4double boardHeight{2 * mm};

    inline G4double SiPMHeight{0.6 * mm};
    inline G4double SiPMWindowThick{0.21 * mm};
    inline G4double SiPMLength{7 * mm};
    inline G4double SiPMWidth{7 * mm};
    inline G4double SiPMFrameSize{0.5 * mm};

    inline G4int crystalSiPMCount{16};
    inline G4int vetoSiPMCount{8};
    inline G4int bottomVetoSiPMCount{4};

}

#endif //SIZES_HH