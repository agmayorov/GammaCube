#include "GenSurface.hh"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <G4Exception.hh>
#include <Randomize.hh>

#include "Configuration.hh"

using namespace Sizes;
using namespace Configuration;


void GenSurface::PayloadExtent(G4double& radius, G4double& zMin, G4double& zMax) {
    radius = modelRadius;
    zMax = modelHeight / 2.0;
    zMin = -modelHeight / 2.0 - plateThick - plateCenterThick - bottomCapHeight;
}


G4ThreeVector GenSurface::PayloadCentre() {
    G4double r, zMin, zMax;
    PayloadExtent(r, zMin, zMax);

    return {0., 0., 0.5 * (zMin + zMax)};
}


G4double GenSurface::BoundingRadius() {
    G4double r, zMin, zMax;
    PayloadExtent(r, zMin, zMax);
    const G4double zc = 0.5 * (zMin + zMax);

    const G4double cylinder = std::hypot(r, 0.5 * (zMax - zMin));

    const G4double plateHalf = plateSize / 2.0 + plateCornerSize;
    const G4double plateZ = std::max(std::fabs(-modelHeight / 2.0 - plateThick - zc),
                                     std::fabs(-modelHeight / 2.0 - zc));
    const G4double plate = std::sqrt(2.0 * plateHalf * plateHalf + plateZ * plateZ);

    return std::max(cylinder, plate) + Margin();
}


G4ThreeVector GenSurface::Arrival(const G4double theta, const G4double phi) {
    return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
}


G4String GenSurface::DirectionTag() {
    if (fluxDirection.find("isotropic") != std::string::npos) return fluxDirection;

    std::ostringstream ss;
    ss << fluxDirection << "_t" << beamTheta / deg << "_p" << beamPhi / deg;
    if (fluxDirection == "point_beam" || fluxDirection == "point_iso") ss << "_z" << sourceZ / mm;
    if (fluxDirection == "point_iso") ss << "_c" << Configuration::coneAngle / deg;
    return ss.str();
}


GenSurface GenSurface::For(const G4String& fluxDirection) {
    GenSurface s;
    s.origin = PayloadCentre();
    s.standoff = BoundingRadius();
    s.radius = s.standoff;

    G4double zMin, zMax;
    PayloadExtent(s.cylinderRadius, zMin, zMax);
    s.cylinderHalfHeight = 0.5 * (zMax - zMin);

    if (fluxDirection == "isotropic" || fluxDirection == "isotropic_up" || fluxDirection == "isotropic_down") {
        s.shape = Shape::Sphere;
        if (fluxDirection == "isotropic_up") s.hemisphere = 1;
        else if (fluxDirection == "isotropic_down") s.hemisphere = -1;
        return s;
    }

    if (fluxDirection != "flat" && fluxDirection != "point_beam" && fluxDirection != "point_iso") {
        G4Exception("GenSurface::For", "FluxDirection", FatalException,
                    ("Flux direction is not implemented: " + fluxDirection +
                        ".\nAvailable flux directions: isotropic, isotropic_up, isotropic_down, flat,"
                        " point_beam, point_iso (aliases: vertical_up, vertical_down, horizontal)").c_str());
    }

    s.theta = beamTheta;
    s.phi = beamPhi;
    s.axis = (-Arrival(s.theta, s.phi)).unit();

    if (fluxDirection == "point_beam" || fluxDirection == "point_iso") {
        if (Configuration::coneAngle <= 0. || Configuration::coneAngle > 180. * deg) {
            G4Exception("GenSurface::For", "ConeAngle", FatalException, "Cone angle must be in (0, 180] deg");
        }
        s.shape = fluxDirection == "point_beam" ? Shape::PointBeam : Shape::PointIso;
        s.sourceRadius = BoundingRadius();
        s.sourcePosition = G4ThreeVector(s.sourceRadius * std::cos(s.phi), s.sourceRadius * std::sin(s.phi), sourceZ);
        s.coneAngle = s.shape == Shape::PointIso ? Configuration::coneAngle : 0.;
        return s;
    }

    s.shape = Shape::Flat;

    const G4ThreeVector up(0., 0., 1.);
    const G4ThreeVector upPerp = up - s.axis * up.dot(s.axis);
    s.u = upPerp.mag() > 1e-9 ? upPerp.unit() : G4ThreeVector(1., 0., 0.);
    s.v = s.axis.cross(s.u).unit();

    const G4double sinTh = upPerp.mag() > 1e-12 ? upPerp.mag() : 0.;
    const G4double cosTh = std::fabs(s.axis.z()) > 1e-12 ? std::fabs(s.axis.z()) : 0.;

    s.halfLength = s.cylinderHalfHeight * sinTh;
    s.capHalfU = s.cylinderRadius * cosTh;
    s.halfU = s.halfLength + s.capHalfU;
    s.halfV = s.cylinderRadius;

    return s;
}


G4bool GenSurface::InSilhouette(const G4double a, const G4double b) const {
    if (std::fabs(b) > cylinderRadius) return false;

    const G4double over = std::fabs(a) - halfLength;
    if (over <= 0.) return true;
    if (capHalfU <= 0.) return false;

    const G4double x = over / capHalfU;
    const G4double y = b / cylinderRadius;
    return x * x + y * y <= 1.;
}


void GenSurface::Sample(G4ThreeVector& pos, G4ThreeVector& dir) const {
    if (shape == Shape::Sphere) {
        SampleSphere(pos, dir);
    } else if (shape == Shape::Flat) {
        SampleFlat(pos, dir);
    } else if (shape == Shape::PointIso) {
        SamplePointIso(pos, dir);
    } else {
        pos = sourcePosition;
        dir = axis;
    }
}


void GenSurface::SampleSphere(G4ThreeVector& pos, G4ThreeVector& dir) const {
    G4double cosT = 2.0 * G4UniformRand() - 1.0;
    if (hemisphere > 0) cosT = G4UniformRand();
    else if (hemisphere < 0) cosT = -G4UniformRand();

    const G4double psi = twopi * G4UniformRand();
    const G4double l = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const G4ThreeVector rhat(l * std::cos(psi), l * std::sin(psi), cosT);

    pos = origin + radius * rhat;

    const G4ThreeVector z = rhat.unit();
    const G4ThreeVector a = std::fabs(z.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector x = z.cross(a).unit();
    const G4ThreeVector y = z.cross(x).unit();

    const G4double ksi = G4UniformRand();
    const G4double sinTh = std::sqrt(ksi);
    const G4double cosTh = std::sqrt(1.0 - ksi);
    const G4double psi2 = twopi * G4UniformRand();

    dir = -(sinTh * std::cos(psi2) * x + sinTh * std::sin(psi2) * y + cosTh * z);
    dir = dir.unit();
}


void GenSurface::SampleFlat(G4ThreeVector& pos, G4ThreeVector& dir) const {
    G4double a, b;
    do {
        a = (2.0 * G4UniformRand() - 1.0) * halfU;
        b = (2.0 * G4UniformRand() - 1.0) * halfV;
    } while (!InSilhouette(a, b));

    dir = axis;
    pos = origin + u * a + v * b - axis * standoff;
}


void GenSurface::SamplePointIso(G4ThreeVector& pos, G4ThreeVector& dir) const {
    const G4double cosB = 1.0 - G4UniformRand() * (1.0 - std::cos(coneAngle));
    const G4double sinB = std::sqrt(std::max(0.0, 1.0 - cosB * cosB));
    const G4double psi = twopi * G4UniformRand();

    const G4ThreeVector a = std::fabs(axis.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector x = axis.cross(a).unit();
    const G4ThreeVector y = axis.cross(x).unit();

    pos = sourcePosition;
    dir = (sinB * std::cos(psi) * x + sinB * std::sin(psi) * y + cosB * axis).unit();
}


G4double GenSurface::ConeSolidAngleFraction() const {
    if (shape == Shape::PointBeam) return 0.;
    return 0.5 * (1.0 - std::cos(coneAngle));
}


G4double GenSurface::OuterRadius() const {
    if (shape == Shape::Sphere) return radius;
    if (IsPoint()) return (sourcePosition - origin).mag();
    return std::sqrt(standoff * standoff + halfU * halfU + halfV * halfV);
}


G4double GenSurface::SPerp_cm2() const {
    if (IsPoint()) return 0.0;
    if (shape == Shape::Sphere) {
        return pi * (radius / cm) * (radius / cm);
    }
    return pi * (cylinderRadius / cm) * (capHalfU / cm) + 4.0 * (cylinderRadius / cm) * (halfLength / cm);
}


G4double GenSurface::GeomFactor_cm2sr() const {
    if (shape != Shape::Sphere) return 0.0;
    const G4double r_cm = radius / cm;
    return (hemisphere != 0 ? 2.0 : 4.0) * pi * pi * r_cm * r_cm;
}


G4double GenSurface::Norm_cm2() const {
    if (IsPoint()) return 1.0;
    return shape == Shape::Sphere ? GeomFactor_cm2sr() : SPerp_cm2();
}


G4String GenSurface::ShapeName() const {
    if (shape == Shape::Sphere) return hemisphere != 0 ? "hemisphere" : "sphere";
    if (shape == Shape::PointBeam) return "point_beam";
    if (shape == Shape::PointIso) return "point_cone";
    return capHalfU > 0. && halfLength > 0. ? "cylinder_silhouette" : (halfLength > 0. ? "rectangle" : "disk");
}
