#include "Detector.hh"


void Detector::init() {
    DefineMaterials();
    DefineVisual();
}


G4String Detector::GetDetectorType() const {
    return detectorType;
}


void Detector::DefineVisual() {
    visWhite = new G4VisAttributes(G4Color(1.0, 1.0, 1.0));
    visWhite->SetForceSolid(true);
    visBlue = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
    visBlue->SetForceSolid(true);
    visRed = new G4VisAttributes(G4Color(1.0, 0.0, 0.0));
    visRed->SetForceSolid(true);
    visCyan = new G4VisAttributes(G4Color(0.0, 1.0, 1.0));
    visCyan->SetForceSolid(true);
    visGrey = new G4VisAttributes(G4Color(0.73, 0.746, 0.7578));
    visGrey->SetForceSolid(true);
}
