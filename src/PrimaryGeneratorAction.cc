#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(Geometry *geometry, G4bool vFlux) : particleGun(new G4ParticleGun(1)),
    verticalFlux(vFlux) {
    center = geometry->detectorPos;
    detectorHalfSize = geometry->detContainerSize;
    const G4ThreeVector tempVec = G4ThreeVector(0,
                                                geometry->modelSize.y() + geometry->sizes.lidThick,
                                                geometry->modelSize.z() + geometry->sizes.lidThick);
    radius = sqrt(tempVec.y() * tempVec.y() + tempVec.z() * tempVec.z()) + 5 * mm;

    alpha = -1.30;
    beta = -2.34;
    E0 = 0.18 / (2 + alpha);
    useBandForGammas = true;
    if (useBandForGammas) {
        Emin = 0.001 * MeV;
        Emax = 500 * MeV;
        BuildBandCDF();
    }
}


PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}


G4double PrimaryGeneratorAction::BandNofE(const G4double E) const {
    const G4double Eb = (alpha - beta) * E0;
    constexpr G4double E_100 = 0.1 * MeV;

    if (E < Eb) {
        return std::pow(E / E_100, alpha) * std::exp(-E / E0);
    }
    return std::pow((alpha - beta) * E0 / E_100, alpha - beta)
           * std::exp(beta - alpha) * std::pow(E / E_100, beta);
}


void PrimaryGeneratorAction::BuildBandCDF() {
    constexpr int N = 8192;
    bandE.resize(N);
    bandCDF.assign(N, 0.0);

    const G4double logR = std::log(Emax / Emin);

    for (int i = 0; i < N; ++i) {
        const G4double u = i / static_cast<G4double>(N - 1);
        bandE[i] = Emin * std::exp(u * logR);
    }

    std::vector<G4double> f(N);
    for (int i = 0; i < N; ++i) f[i] = BandNofE(bandE[i]);

    G4double cum = 0.0;
    bandCDF[0] = 0.0;
    for (int i = 0; i < N - 1; ++i) {
        const G4double dE = bandE[i + 1] - bandE[i];
        const G4double trap = 0.5 * (f[i] + f[i + 1]) * dE;
        cum += trap;
        bandCDF[i + 1] = cum;
    }

    if (cum <= 0) {
        for (int i = 0; i < N; ++i) bandCDF[i] = i / static_cast<G4double>(N - 1);
        return;
    }
    for (int i = 0; i < N; ++i) bandCDF[i] /= cum;
}


G4double PrimaryGeneratorAction::SampleBandEnergyCDF() const {
    const G4double u = G4UniformRand();

    const auto it = std::lower_bound(bandCDF.begin(), bandCDF.end(), u);
    if (it == bandCDF.begin()) return bandE.front();
    if (it == bandCDF.end()) return bandE.back();

    const size_t j = std::distance(bandCDF.begin(), it);
    const size_t i = j - 1;

    const G4double c0 = bandCDF[i];
    const G4double c1 = bandCDF[j];
    const G4double t = (u - c0) / std::max(c1 - c0, 1e-16);

    return bandE[i] + t * (bandE[j] - bandE[i]);
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


void PrimaryGeneratorAction::PickParticle(G4ParticleDefinition *&pdef, G4String &name) {
    auto *pt = G4ParticleTable::GetParticleTable();
    if (useBandForGammas) {
        pdef = pt->FindParticle("gamma");
        name = "gamma";
        return;
    }
    if (const G4double r = G4UniformRand(); r < 0.25) {
        pdef = pt->FindParticle("gamma");
        name = "gamma";
        Emin = 0.001 * MeV;
        Emax = 500 * MeV;
    } else if (r < 0.50) {
        pdef = pt->FindParticle("e-");
        name = "e-";
        Emin = 0.001 * MeV;
        Emax = 100. * MeV;
    } else if (r < 0.75) {
        pdef = pt->FindParticle("proton");
        name = "proton";
        Emin = 1. * MeV;
        Emax = 10000. * MeV;
    } else {
        pdef = pt->FindParticle("alpha");
        name = "alpha";
        Emin = 10. * MeV;
        Emax = 10000. * MeV;
    }
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event *evt) {
    G4ParticleDefinition *pdef = nullptr;
    G4String pname;
    PickParticle(pdef, pname);

    G4ThreeVector x, v;
    if (verticalFlux) {
        v = G4ThreeVector(0., 0., -1.);
        const G4double x_ = (G4UniformRand() - 0.5) * detectorHalfSize.x();
        const G4double y_ = (G4UniformRand() - 0.5) * detectorHalfSize.y();
        const G4double z_ = radius;
        x = G4ThreeVector(x_, y_, z_);
    } else {
        GenerateOnSphere(x, v);
    }
    G4double energy = 0.0;
    if (useBandForGammas) {
        energy = SampleBandEnergyCDF();
    } else {
        energy = Emin * std::pow(Emax / Emin, G4UniformRand());
    }

    particleGun->SetParticleDefinition(pdef);
    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticlePosition(x);
    particleGun->SetParticleMomentumDirection(v);
    particleGun->SetParticleTime(0.0 * ns);
    particleGun->GeneratePrimaryVertex(evt);

    if (auto *ea = dynamic_cast<EventAction *>(G4EventManager::GetEventManager()->GetUserEventAction())) {
        PrimaryRec rec;
        rec.index = static_cast<int>(ea->primBuf.size());
        rec.pdg = pdef->GetPDGEncoding();
        rec.name = pname;
        rec.E_MeV = energy;
        rec.dir = v;
        rec.pos_mm = x / mm;
        rec.t0_ns = 0.0;
        ea->primBuf.emplace_back(std::move(rec));
    }
}
