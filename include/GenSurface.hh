#ifndef GENSURFACE_HH
#define GENSURFACE_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "Sizes.hh"


class GenSurface {
public:
    enum class Shape { Sphere, Flat, PointBeam, PointIso };

    static GenSurface For(const G4String& fluxDirection);

    static void PayloadExtent(G4double& radius, G4double& zMin, G4double& zMax);
    static G4ThreeVector PayloadCentre();
    static G4double BoundingRadius();
    static G4ThreeVector Arrival(G4double theta, G4double phi);
    static G4String DirectionTag();
    static G4double Margin() { return 5. * mm; }

    void Sample(G4ThreeVector& pos, G4ThreeVector& dir) const;

    [[nodiscard]] G4bool IsIsotropic() const { return shape == Shape::Sphere; }
    [[nodiscard]] G4bool IsPoint() const { return shape == Shape::PointBeam || shape == Shape::PointIso; }
    [[nodiscard]] G4bool IsCone() const { return shape == Shape::PointIso; }
    [[nodiscard]] G4int Hemisphere() const { return hemisphere; }
    [[nodiscard]] G4double OuterRadius() const;

    [[nodiscard]] G4double Radius() const { return radius; }
    [[nodiscard]] G4double CylinderRadius() const { return cylinderRadius; }
    [[nodiscard]] G4double CylinderHalfHeight() const { return cylinderHalfHeight; }
    [[nodiscard]] G4double HalfLength() const { return halfLength; }
    [[nodiscard]] G4double CapHalfU() const { return capHalfU; }
    [[nodiscard]] G4double HalfU() const { return halfU; }
    [[nodiscard]] G4double HalfV() const { return halfV; }
    [[nodiscard]] G4double Standoff() const { return standoff; }
    [[nodiscard]] G4double Theta() const { return theta; }
    [[nodiscard]] G4double Phi() const { return phi; }
    [[nodiscard]] const G4ThreeVector& Origin() const { return origin; }
    [[nodiscard]] const G4ThreeVector& SourcePosition() const { return sourcePosition; }
    [[nodiscard]] G4double SourceRadius() const { return sourceRadius; }
    [[nodiscard]] G4double ConeAngle() const { return coneAngle; }
    [[nodiscard]] G4double ConeSolidAngleFraction() const;

    [[nodiscard]] const G4ThreeVector& Axis() const { return axis; }
    [[nodiscard]] const G4ThreeVector& U() const { return u; }
    [[nodiscard]] const G4ThreeVector& V() const { return v; }

    [[nodiscard]] G4bool InSilhouette(G4double a, G4double b) const;

    [[nodiscard]] G4double SPerp_cm2() const;
    [[nodiscard]] G4double GeomFactor_cm2sr() const;
    [[nodiscard]] G4double Norm_cm2() const;

    [[nodiscard]] G4String ShapeName() const;

private:
    Shape shape{Shape::Sphere};
    G4int hemisphere{0};

    G4double radius{0.};
    G4double cylinderRadius{0.};
    G4double cylinderHalfHeight{0.};
    G4double halfLength{0.};
    G4double capHalfU{0.};
    G4double halfU{0.};
    G4double halfV{0.};
    G4double standoff{0.};
    G4double theta{0.};
    G4double phi{0.};
    G4double sourceRadius{0.};
    G4double coneAngle{0.};

    G4ThreeVector origin{0., 0., 0.};
    G4ThreeVector sourcePosition{0., 0., 0.};
    G4ThreeVector axis{0., 0., -1.};
    G4ThreeVector u{1., 0., 0.};
    G4ThreeVector v{0., 1., 0.};

    void SampleSphere(G4ThreeVector& pos, G4ThreeVector& dir) const;
    void SampleFlat(G4ThreeVector& pos, G4ThreeVector& dir) const;
    void SamplePointIso(G4ThreeVector& pos, G4ThreeVector& dir) const;
};

#endif //GENSURFACE_HH
