#include "PrimaryGeneratorAction.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(G4String fDir, const G4String& fluxType, const G4double cThreshold)
    : particleGun(new G4ParticleGun(1)),
      fluxDirection(std::move(fDir)),
      eCrystalThreshold(cThreshold) {
    genSurface = GenSurface::For(fluxDirection);

    std::vector<G4String> fluxTypeList = {"Uniform", "PLAW", "COMP", "SEP", "Galactic", "Table"};
    if (std::find(fluxTypeList.begin(), fluxTypeList.end(), fluxType) == fluxTypeList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxType", FatalException,
                    ("Flux type not found: " + fluxType + ".\nAvailable flux types: Uniform, PLAW, SEP, Galactic")
                    .
                    c_str());
    }

    if (fluxType == "Uniform") {
        flux = new UniformFlux();
    } else if (fluxType == "PLAW") {
        flux = new PLAWFlux();
    } else if (fluxType == "COMP") {
        flux = new COMPFlux();
    } else if (fluxType == "SEP") {
        flux = new SEPFlux();
    } else if (fluxType == "Galactic") {
        flux = new GalacticFlux();
    } else if (fluxType == "Table") {
        flux = new TableFlux();
    }
}


PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt) {
    G4ThreeVector x, v;
    genSurface.Sample(x, v);
    ParticleInfo info = flux->GenerateParticle();

    particleGun->SetParticleDefinition(info.def);
    particleGun->SetParticleEnergy(info.energy);
    particleGun->SetParticlePosition(x);
    particleGun->SetParticleMomentumDirection(v);
    particleGun->SetParticleTime(0.0 * ns);
    particleGun->GeneratePrimaryVertex(evt);

    if (auto* ea = dynamic_cast<EventAction*>(G4EventManager::GetEventManager()->GetUserEventAction())) {
        PrimaryRec rec;
        rec.index = static_cast<int>(ea->primBuf.size());
        rec.pdg = info.pdg;
        rec.name = info.name;
        rec.E_MeV = info.energy / MeV;
        rec.dir = v;
        rec.pos_mm = x / mm;
        rec.t0_ns = 0.0;
        ea->primBuf.emplace_back(std::move(rec));
    }
}
