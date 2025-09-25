#ifndef PRMIARYGENERATIONACTION_HH
#define PRMIARYGENERATIONACTION_HH

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4VVisManager.hh>
#include <G4Event.hh>
#include <G4Circle.hh>
#include <G4EventManager.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <cmath>
#include <fstream>
#include <numeric>

#include "EventAction.hh"
#include "Geometry.hh"


enum class PdfType {
    PerE,
    PerLog10E,
    IntegralAboveE
};

struct Row {
    double E_MeV;
    double flux;
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction(const Geometry *, G4bool);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event *evt) override;

private:
    G4ParticleGun *particleGun = nullptr;

    G4double radius;
    G4ThreeVector center;
    G4ThreeVector detectorHalfSize;

    G4double Emin;
    G4double Emax;

    G4bool usePLAWForGammas;
    G4double alpha;
    G4double A;

    std::vector<G4double> bandE;
    std::vector<G4double> bandCDF;

    G4bool verticalFlux;

    void GenerateOnSphere(G4ThreeVector &pos, G4ThreeVector &dir) const;
    void PickParticle(G4ParticleDefinition *&pdef, G4String &name);

    G4double BandNofE(G4double E) const;
    void BuildBandCDF();
    G4double SampleBandEnergyCDF() const;

    G4double SamplePLAWEnergy() const;

    G4bool useCSVForProtons = false;
    G4String csvPath = "../Data_Sheet_2.CSV";
    int csvYear = 1999;
    int csvOrder = 15;

    G4double csvEnergyToMeV;
    std::vector<G4double> csvE, csvCDF;
    PdfType csvPdfType;

    void BuildCSVFluxCDF();
    G4double SampleLogPolyEnergyCDF() const;
};

#endif //PRMIARYGENERATIONACTION_HH
