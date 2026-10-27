#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <G4String.hh>
#include <G4Types.hh>
#include <G4SystemOfUnits.hh>


namespace Configuration
{
    inline G4String detectorType{"CsI"};

    inline G4String geometryConfig{"default"};   // default | tunacan | cubesat | experiment
    inline G4String engineConfig{""};            // "" (no engine) | sputnix | geoscan
    inline G4String cameraConfig{""};            // "" (no camera) | galper

    inline G4String fluxType{"Uniform"};
    inline G4String fluxDirection{"isotropic"};
    inline G4double beamTheta{0.};
    inline G4double beamPhi{0.};
    inline G4double sourceZ{0.};
    inline G4double sourceR{-1.};
    inline G4double coneAngle{180. * deg};

    inline G4long seed{0};

    inline G4double eCrystalThreshold{0 * MeV};
    inline G4double eVetoThreshold{0 * MeV};

    inline G4bool useOptics{false};
    inline G4int yieldScale{1};
    inline G4int oCrystalThreshold{0};
    inline G4int oVetoThreshold{0};
    inline G4int oBottomVetoThreshold{0};
    inline G4double photonTimeLimit{-1.};   // readout gate; <= 0 means disabled

    inline G4int sipmCells{0};                   // microcells per channel; 0 disables saturation
    inline G4double sipmRecoveryTime{95 * ns};   // microcell recharge time constant
    inline G4double sipmCrosstalk{0.};           // optical crosstalk probability; 0 disables it

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


#endif //CONFIGURATION_HH
