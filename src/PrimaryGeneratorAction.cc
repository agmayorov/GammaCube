#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(const Geometry *geometry, const G4bool vFlux) : particleGun(
        new G4ParticleGun(1)),
    verticalFlux(vFlux) {
    center = geometry->detectorPos;
    detectorHalfSize = geometry->modelSize;
    const G4ThreeVector tempVec = G4ThreeVector(0,
                                                geometry->modelSize.y(),
                                                geometry->modelSize.z());
    radius = sqrt(tempVec.y() * tempVec.y() + tempVec.z() * tempVec.z()) + 5 * mm;

    alpha = 1.411103;
    A = 9.360733;
    usePLAWForGammas = false;
    if (usePLAWForGammas) {
        Emin = 0.01 * MeV;
        Emax = 100 * MeV;
    }

    useCSVForProtons = true;
    csvPath = "../Data_Sheet_2.CSV";
    csvYear = 1998;
    csvOrder = 15;

    csvEnergyToMeV = 1.;
    csvPdfType = PdfType::PerE;

    if (useCSVForProtons) {
        Emin = 0.1 * MeV;
        Emax = 1000 * MeV;
        BuildCSVFluxCDF();
    }
}


PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
}


G4double PrimaryGeneratorAction::SamplePLAWEnergy() const {
    const G4double u = G4UniformRand();

    if (std::abs(alpha - 1.0) < 1e-12) {
        return Emin * std::pow(Emax / Emin, u);
    }

    const G4double EminPow = std::pow(Emin, 1.0 - alpha);
    const G4double EmaxPow = std::pow(Emax, 1.0 - alpha);
    const G4double val = EminPow + u * (EmaxPow - EminPow);
    return std::pow(val, 1.0 / (1.0 - alpha));
}


static std::vector<double> extract_numbers(const std::string &line) {
    static const std::regex re(R"(([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?))");
    std::vector<double> out;
    for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end; ++it) {
        out.push_back(std::stod((*it)[1].str()));
    }
    return out;
}

void PrimaryGeneratorAction::BuildCSVFluxCDF() {
    csvE.clear();
    csvCDF.clear();

    std::ifstream in(csvPath.c_str());
    if (!in) {
        G4Exception("PrimaryGeneratorAction::BuildCSVFluxCDF", "CSV_OPEN_FAIL",
                    JustWarning, ("Cannot open " + csvPath).c_str());
        csvE = {1. * MeV, 10. * MeV};
        csvCDF = {0.0, 1.0};
        return;
    }

    std::vector<Row> rows;
    rows.reserve(1884);

    std::string line;
    bool headerSkipped = false;
    while (std::getline(in, line)) {
        if (!headerSkipped) { headerSkipped = true; continue; }
        auto nums = extract_numbers(line);
        if (nums.size() < 4) continue;

        const int yr = static_cast<int>(std::llround(nums[0]));
        const double E_csv = nums[1];
        const double flx = nums[2];
        const int ord = static_cast<int>(std::llround(nums.back()));

        if (yr == csvYear && ord == csvOrder) {
            const double E_G4 = E_csv * csvEnergyToMeV * MeV;
            rows.push_back({E_G4, flx});
        }
    }
    in.close();

    if (rows.size() < 2) {
        G4Exception("PrimaryGeneratorAction::BuildCSVFluxCDF", "CSV_NO_ROWS",
                    JustWarning, "No matching rows for given year/order.");
        csvE = {1. * MeV, 10. * MeV};
        csvCDF = {0.0, 1.0};
        return;
    }

    std::sort(rows.begin(), rows.end(),
              [](const Row &a, const Row &b) { return a.E_MeV < b.E_MeV; });

    {
        std::vector<Row> uniq;
        uniq.reserve(rows.size());
        for (const auto &r: rows) {
            if (!uniq.empty() && std::abs(r.E_MeV - uniq.back().E_MeV) < 1e-12 * MeV) continue;
            uniq.push_back(r);
        }
        rows.swap(uniq);
    }

    const size_t Ntot = rows.size();
    std::vector<double> E(Ntot), lE(Ntot), f_perE(Ntot, 0.0);
    for (size_t i = 0; i < Ntot; ++i) {
        E[i]  = rows[i].E_MeV;
        lE[i] = std::log(E[i]);
    }

    auto clamp_pos = [](const double x) { return (x > 0.0 && std::isfinite(x)) ? x : 0.0; };

    if (csvPdfType == PdfType::PerE) {
        for (size_t i = 0; i < Ntot; ++i) f_perE[i] = clamp_pos(rows[i].flux);
    } else if (csvPdfType == PdfType::PerLog10E) {
        constexpr double LN10 = 2.302585092994046;
        for (size_t i = 0; i < Ntot; ++i) f_perE[i] = clamp_pos(rows[i].flux / (E[i] * LN10));
    } else {
        if (Ntot >= 2) {
            auto dphi_dlnE = [&](size_t i)->double {
                if (i == 0)         return (rows[1].flux      - rows[0].flux)      / (lE[1]      - lE[0]);
                if (i + 1 == Ntot)  return (rows[Ntot-1].flux - rows[Ntot-2].flux) / (lE[Ntot-1] - lE[Ntot-2]);
                return (rows[i+1].flux - rows[i-1].flux) / (lE[i+1] - lE[i-1]);
            };
            for (size_t i = 0; i < Ntot; ++i) f_perE[i] = clamp_pos( -dphi_dlnE(i) / E[i] );
        }
    }

    bool useBounds = (Emin > 0.0 && Emax > 0.0 && std::isfinite(Emin) &&
                      std::isfinite(Emax) && Emax > Emin);

    size_t iLo = 0, iHi = Ntot - 1;
    if (useBounds) {
        auto itLo = std::lower_bound(E.begin(), E.end(), Emin);
        if (itLo == E.begin()) {
            iLo = 0;
        } else if (itLo == E.end()) {
            iLo = Ntot - 1;
        } else {
            size_t j = static_cast<size_t>(itLo - E.begin());
            iLo = (Emin - E[j-1] <= E[j] - Emin) ? (j - 1) : j;
        }

        auto itHi = std::upper_bound(E.begin(), E.end(), Emax);
        if (itHi == E.begin()) {
            iHi = 0;
        } else if (itHi == E.end()) {
            iHi = Ntot - 1;
        } else {
            size_t j = static_cast<size_t>(itHi - E.begin());
            size_t jm1 = j - 1;
            if (j < Ntot) {
                iHi = (Emax - E[jm1] <= E[j] - Emax) ? jm1 : j;
            } else {
                iHi = jm1;
            }
        }

        if (iLo > iHi) std::swap(iLo, iHi);
        if (iHi == iLo && Ntot >= 2) {
            if (iHi + 1 < Ntot) ++iHi; else if (iLo > 0) --iLo;
        }
    }

    const size_t N = iHi - iLo + 1;
    if (N < 2) {
        G4Exception("PrimaryGeneratorAction::BuildCSVFluxCDF", "CSV_RANGE_TOO_NARROW",
                    JustWarning, "Energy range too narrow (fewer than 2 points).");
        csvE = {1. * MeV, 10. * MeV};
        csvCDF = {0.0, 1.0};
        return;
    }

    std::vector<double> Es(N), lEs(N), f_sub(N);
    for (size_t k = 0; k < N; ++k) {
        Es[k]   = E[iLo + k];
        lEs[k]  = std::log(Es[k]);
        f_sub[k]= f_perE[iLo + k];
    }

    csvE = Es;
    csvCDF.assign(N, 0.0);
    long double acc = 0.0L;
    for (size_t i = 0; i + 1 < N; ++i) {
        const long double gi   = static_cast<long double>(f_sub[i])   * static_cast<long double>(Es[i]);
        const long double gip1 = static_cast<long double>(f_sub[i+1]) * static_cast<long double>(Es[i+1]);
        const long double dlnE = lEs[i+1] - lEs[i];
        const long double wi   = 0.5L * (gi + gip1) * dlnE;
        acc += wi;
        csvCDF[i + 1] = static_cast<double>(acc);
    }

    if (acc <= 0.0L || !std::isfinite(static_cast<double>(acc))) {
        for (size_t i = 0; i < N; ++i) csvCDF[i] = static_cast<double>(i) / static_cast<double>(N - 1);
        return;
    }

    for (auto &v: csvCDF) v /= static_cast<double>(acc);
    csvCDF.front() = 0.0;
    csvCDF.back()  = 1.0;
}


G4double PrimaryGeneratorAction::SampleLogPolyEnergyCDF() const {
    if (csvE.size() < 2) return 1.0 * MeV;

    const G4double u = G4UniformRand();

    const auto it = std::upper_bound(csvCDF.begin(), csvCDF.end(), u);
    if (it == csvCDF.begin()) return csvE.front();
    if (it == csvCDF.end()) return csvE.back();
    const size_t j = std::distance(csvCDF.begin(), it);
    const size_t i = j - 1;

    const double C0 = csvCDF[i];
    const double C1 = csvCDF[j];
    const double t = (u - C0) / std::max(C1 - C0, 1e-12);

    const double lnE0 = std::log(csvE[i]);
    const double lnE1 = std::log(csvE[j]);
    const double lnEu = lnE0 + t * (lnE1 - lnE0);

    return std::exp(lnEu);
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
    if (usePLAWForGammas) {
        pdef = pt->FindParticle("gamma");
        name = "gamma";
        return;
    }
    if (useCSVForProtons) {
        pdef = pt->FindParticle("proton");
        name = "proton";
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
        const G4double x_ = 2 * (G4UniformRand() - 0.5) * detectorHalfSize.y();
        const G4double y_ = 2 * (G4UniformRand() - 0.5) * detectorHalfSize.y();
        const G4double z_ = radius;
        x = G4ThreeVector(x_, y_, z_);
    } else {
        GenerateOnSphere(x, v);
    }
    G4double energy = 0.0;
    if (usePLAWForGammas) {
        energy = SamplePLAWEnergy();
    } else if (useCSVForProtons && pname == "proton") {
        energy = SampleLogPolyEnergyCDF();
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
        rec.E_MeV = energy / MeV;
        rec.dir = v;
        rec.pos_mm = x / mm;
        rec.t0_ns = 0.0;
        ea->primBuf.emplace_back(std::move(rec));
    }
}
