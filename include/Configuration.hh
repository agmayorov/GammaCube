#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <G4String.hh>
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>


namespace Configuration
{
    inline G4String detectorType{"CsI"};

    inline G4String fluxType{"Uniform"};
    inline G4String fluxDirection{"isotropic"};

    inline G4double eCrystalThreshold{0 * MeV};
    inline G4double eVetoThreshold{0 * MeV};

    inline G4bool useOptics{false};
    inline G4int yieldScale{1};
    inline G4int oCrystalThreshold{0};
    inline G4int oVetoThreshold{0};
    inline G4int oBottomVetoThreshold{0};

    inline G4String crystalSiPMConfig{"12-cross"};
    inline G4bool polishedTyvek{false};
    inline G4double viewDeg{360 * deg};

    inline G4int nBins{1000};
    inline G4String outputFile{"GammaCube.root"};
    inline G4bool saveSecondaries{false};
    inline G4bool savePhotons{false};

    inline G4double weight4{1.};
    inline G4double weight8{1.};

    inline G4double Emin{-1};
    inline G4double Emax{-1};

    inline G4double N_MIN{10};
    inline G4double N_MIN_opt{10};

    inline G4bool isLogBin{false};
}


namespace DEG
{
    inline G4double tyvekOut{360};
    inline G4double veto{360};
    inline G4double tyvekMid{360};
    inline G4double opticLayerVeto{360};
    inline G4double rubber{360};
    inline G4double shell{360};
    inline G4double bottomVetoShell{360};
    inline G4double tyvekBottom{360};
    inline G4double bottomVeto{360};
    inline G4double opticLayerBottomVeto{360};
    inline G4double crystalContainer{360};
    inline G4double crystal{360};
    inline G4double tyvekIn{360};
    inline G4double crystalShell{360};
    inline G4double crystallGlass{360};
    inline G4double opticLayerCrystall{360};
    inline G4double holder{360};
    inline G4double crystalSiPM{360};
    inline G4double vetoSiPM{360};
    inline G4double vetoSpring{360};
    inline G4double vetoBoard{360};
    inline G4double bottomVetoSiPM{360};
    inline G4double tunaCan{360};
    inline G4double plate{360};
    inline G4double plateHole{360};
}

namespace shift
{
    inline G4double tyvekOut{0};
    inline G4double veto{0};
    inline G4double tyvekMid{0};
    inline G4double opticLayerVeto{0};
    inline G4double rubber{0};
    inline G4double shell{0};
    inline G4double bottomVetoShell{0};
    inline G4double tyvekBottom{0};
    inline G4double bottomVeto{0};
    inline G4double opticLayerBottomVeto{0};
    inline G4double crystalContainer{0};
    inline G4double crystal{0};
    inline G4double tyvekIn{0};
    inline G4double crystalShell{0};
    inline G4double crystallGlass{0};
    inline G4double opticLayerCrystall{0};
    inline G4double holder{0};
    inline G4double crystalSiPM{0};
    inline G4double vetoSiPM{0};
    inline G4double vetoSpring{0};
    inline G4double vetoBoard{0};
    inline G4double bottomVetoSiPM{0};
    inline G4double tunaCan{0};
    inline G4double plate{0};
    inline G4double plateHole{0};
}


#endif //CONFIGURATION_HH
