#include <utility>

#include "PrimaryGeneratorAction.hh"


PrimaryGeneratorAction::PrimaryGeneratorAction(const Geometry *geometry, const G4bool vFlux, const G4String &fluxType)
    : particleGun(new G4ParticleGun(1)), verticalFlux(vFlux) {
    center = geometry->detectorPos;
    detectorHalfSize = geometry->modelSize;
    const G4ThreeVector tempVec = G4ThreeVector(0,
                                                geometry->modelSize.y(),
                                                geometry->modelSize.z());
    radius = sqrt(tempVec.y() * tempVec.y() + tempVec.z() * tempVec.z()) + 5 * mm;

    std::vector<G4String> fluxTypeList = {"Uniform", "PLAW", "SEP", "Galactic"};
    if (std::find(fluxTypeList.begin(), fluxTypeList.end(), fluxType) == fluxTypeList.end()) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "FluxType", FatalException,
                    ("Flux type not found: " + fluxType + ".\nAvailable flux types: Uniform, PLAW, SEP, Galactic").
                    c_str());
    }

    if (fluxType == "Uniform") {
        flux = new UniformFlux();
    } else if (fluxType == "PLAW") {
        flux = new PLAWFlux(0.01 * MeV, 100 * MeV, 1.411103);
    } else if (fluxType == "SEP") {
        flux = new SEPFlux(0.1 * MeV, 1000 * MeV, 1998, 15);
    } else if (fluxType == "Galactic") {
        flux = new GalacticFlux(1 * GeV, 1000 * GeV, 600);
    }
}


PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}


void PrimaryGeneratorAction::GenerateOnSphere(G4ThreeVector &pos, G4ThreeVector &dir) const {
    // const G4double u = 2.0 * G4UniformRand() - 1.0; // cos(theta) ~ U[-1,1]
    const G4double u = G4UniformRand(); // cos(theta) ~ U[0,1]
    const G4double phi = 2.0 * M_PI * G4UniformRand();
    const G4double l = std::sqrt(std::max(0.0, 1.0 - u * u));
    const G4ThreeVector rhat(l * std::cos(phi), l * std::sin(phi), u);

    pos = center + radius * rhat;

    const G4ThreeVector z = rhat.unit();
    const G4ThreeVector a = std::fabs(z.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector x = z.cross(a).unit();
    const G4ThreeVector y = z.cross(x).unit();

    const G4double ksi = G4UniformRand();
    const G4double sinTh = std::sqrt(ksi);
    const G4double cosTh = std::sqrt(1.0 - ksi);
    const G4double phi2 = 2.0 * M_PI * G4UniformRand();

    const G4ThreeVector v_local(sinTh * std::cos(phi2), sinTh * std::sin(phi2), cosTh);

    dir = -(v_local.x() * x + v_local.y() * y + v_local.z() * z);
    dir = dir.unit();
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event *evt) {
    G4ThreeVector x, v;
    if (verticalFlux) {
        v = G4ThreeVector(0., 0., -1.);
        const G4double x_ = 2 * (G4UniformRand() - 0.5) * detectorHalfSize.y();
        const G4double y_ = 2 * (G4UniformRand() - 0.5) * detectorHalfSize.y();
        const G4double z_ = radius;
        x = G4ThreeVector(x_, y_, z_);
    } else {
        GenerateOnSphere(x, v);
    }
    ParticleInfo info = flux->GenerateParticle();

    particleGun->SetParticleDefinition(info.def);
    particleGun->SetParticleEnergy(info.energy);
    particleGun->SetParticlePosition(x);
    particleGun->SetParticleMomentumDirection(v);
    particleGun->SetParticleTime(0.0 * ns);
    particleGun->GeneratePrimaryVertex(evt);

    if (auto *ea = dynamic_cast<EventAction *>(G4EventManager::GetEventManager()->GetUserEventAction())) {
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
