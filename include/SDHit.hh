#ifndef SDHIT_HH
#define SDHIT_HH

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>


class SDHit : public G4VHit {
public:
    explicit SDHit(G4int volID) : volumeID(volID) {}

    SDHit() = default;
    ~SDHit() = default;

    inline void *operator new(size_t);
    inline void operator delete(void *);

    G4int volumeID = -1;
    G4double edep = 0.0;
    G4double tmin = DBL_MAX;

    void AddEdep(G4double);
    void UpdateTmin(G4double t);
    G4double GetEdep() const;

    void Draw() override {}
    void Print() override {}
};


using SDHitCollection = G4THitsCollection<SDHit>;

extern G4ThreadLocal G4Allocator<SDHit> *SDHitAllocator;

inline void *SDHit::operator new(size_t) {
    if (!SDHitAllocator) SDHitAllocator = new G4Allocator<SDHit>;
    return (void *) SDHitAllocator->MallocSingle();
}

inline void SDHit::operator delete(void *hit) {
    SDHitAllocator->FreeSingle((SDHit *) hit);
}

#endif //SDHIT_HH