#include "Flux/Flux.hh"

ParticleInfo Flux::GenerateParticle() {
    auto *pt = G4ParticleTable::GetParticleTable();

    ParticleInfo info;
    info.energy = SampleEnergy();
    info.def = pt->FindParticle(name);
    info.name = name;
    info.pdg = info.def->GetPDGEncoding();
    return info;
}