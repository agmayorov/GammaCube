#include "GenSurface.hh"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <G4Exception.hh>
#include <Randomize.hh>

#include "Configuration.hh"
#include "geometry/CubeSatSizes.hh"

using namespace Sizes;
using namespace Configuration;


// The Ti plate of the bare NATALIA sticks out of the cylinder and is left out of the generation area on purpose;
// in the satellite configurations the box contains it.
static G4bool PlateOutsideEnvelope() {
    return Layout::HasTunaCan() && !Layout::HasSatellite();
}


GenSurface::Envelope GenSurface::PayloadEnvelope() {
    Envelope e;
    if (Layout::HasSatellite()) {
        e.box = true;
        e.halfX = CubeSatSizes::CubeSat::halfX;
        e.halfY = CubeSatSizes::CubeSat::halfY;
        Layout::SatelliteExtentZ(e.zMin, e.zMax);
    } else if (Layout::BareCrystal()) {
        // detector container without the TunaCan
        e.halfX = e.halfY = modelRadius - tunaCanThickWall;
        e.zMax = modelHeight / 2.0 - tunaCanThickTop;
        e.zMin = -modelHeight / 2.0 - plateThick - plateCenterThick - bottomCapHeight + bottomCapThick;
    } else {
        e.halfX = e.halfY = modelRadius;
        e.zMax = modelHeight / 2.0;
        e.zMin = -modelHeight / 2.0 - plateThick - plateCenterThick - bottomCapHeight;
    }
    return e;
}


G4ThreeVector GenSurface::PayloadCentre() {
    const Envelope e = PayloadEnvelope();

    return {0., 0., 0.5 * (e.zMin + e.zMax)};
}


G4double GenSurface::BoundingRadius() {
    const Envelope e = PayloadEnvelope();
    const G4double zc = 0.5 * (e.zMin + e.zMax);

    if (e.box) return std::sqrt(e.halfX * e.halfX + e.halfY * e.halfY + 0.25 * (e.zMax - e.zMin) * (e.zMax - e.zMin))
                      + Margin();

    const G4double cylinder = std::hypot(e.halfX, 0.5 * (e.zMax - e.zMin));
    if (!PlateOutsideEnvelope()) return cylinder + Margin();

    const G4double plateHalf = plateSize / 2.0 + plateCornerSize;
    const G4double plateZ = std::max(std::fabs(-modelHeight / 2.0 - plateThick - zc),
                                     std::fabs(-modelHeight / 2.0 - zc));
    const G4double plate = std::sqrt(2.0 * plateHalf * plateHalf + plateZ * plateZ);

    return std::max(cylinder, plate) + Margin();
}


G4double GenSurface::MinSourceRadius(const G4double z) {
    const Envelope e = PayloadEnvelope();
    const G4double r = e.box ? std::hypot(e.halfX, e.halfY) : e.halfX;

    const G4double plateTop = -modelHeight / 2.0;
    const G4double plateBottom = plateTop - plateThick;
    const G4double plateHalf = plateSize / 2.0 + plateCornerSize;

    G4double rMin = 0.;
    if (z > e.zMin - Margin() && z < e.zMax + Margin()) rMin = r + Margin();
    if (PlateOutsideEnvelope() && z > plateBottom - Margin() && z < plateTop + Margin())
        rMin = std::sqrt(2.0) * plateHalf + Margin();
    return rMin;
}


G4bool GenSurface::InEnvelope(const G4ThreeVector& p, const G4double tolerance) {
    const Envelope e = PayloadEnvelope();
    if (p.z() < e.zMin - tolerance || p.z() > e.zMax + tolerance) return false;
    if (e.box) return std::fabs(p.x()) <= e.halfX + tolerance && std::fabs(p.y()) <= e.halfY + tolerance;
    return std::hypot(p.x(), p.y()) <= e.halfX + tolerance;
}


G4bool GenSurface::AllowedOutsideEnvelope(const G4String& physicalVolumeName) {
    if (!PlateOutsideEnvelope()) return false;
    return physicalVolumeName == "PlatePVPL" || physicalVolumeName == "PlateStripLeftPVPL" ||
           physicalVolumeName == "PlateStripRightPVPL";
}


G4ThreeVector GenSurface::Arrival(const G4double theta, const G4double phi) {
    return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
}


G4String GenSurface::DirectionTag() {
    if (fluxDirection.find("isotropic") != std::string::npos) return fluxDirection;

    std::ostringstream ss;
    ss << fluxDirection << "_t" << beamTheta / deg << "_p" << beamPhi / deg;
    if (fluxDirection == "point_beam" || fluxDirection == "point_iso") {
        ss << "_r" << (sourceR < 0. ? BoundingRadius() : sourceR) / mm << "_z" << sourceZ / mm;
    }
    if (fluxDirection == "point_iso") ss << "_c" << Configuration::coneAngle / deg;
    return ss.str();
}


GenSurface GenSurface::For(const G4String& fluxDirection) {
    GenSurface s;
    s.origin = PayloadCentre();
    s.standoff = BoundingRadius();
    s.radius = s.standoff;

    const Envelope e = PayloadEnvelope();
    s.box = e.box;
    s.cylinderRadius = e.box ? 0. : e.halfX;
    s.cylinderHalfHeight = 0.5 * (e.zMax - e.zMin);
    if (e.box) s.boxHalfSize = G4ThreeVector(e.halfX, e.halfY, s.cylinderHalfHeight);

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
        s.sourceRadius = sourceR < 0. ? BoundingRadius() : sourceR;
        if (s.sourceRadius < MinSourceRadius(sourceZ)) {
            std::ostringstream msg;
            msg << "Point source at r = " << s.sourceRadius / mm << " mm, z = " << sourceZ / mm
                << " mm is inside or too close to the model. Minimal radius at this height: "
                << MinSourceRadius(sourceZ) / mm << " mm";
            G4Exception("GenSurface::For", "SourceRadius", FatalException, msg.str().c_str());
        }
        s.sourcePosition = G4ThreeVector(s.sourceRadius * std::cos(s.phi), s.sourceRadius * std::sin(s.phi), sourceZ);
        s.coneAngle = s.shape == Shape::PointIso ? Configuration::coneAngle : 0.;
        return s;
    }

    s.shape = Shape::Flat;

    const G4ThreeVector up(0., 0., 1.);
    const G4ThreeVector upPerp = up - s.axis * up.dot(s.axis);
    s.u = upPerp.mag() > 1e-9 ? upPerp.unit() : G4ThreeVector(1., 0., 0.);
    s.v = s.axis.cross(s.u).unit();

    if (s.box) {
        const G4ThreeVector& h = s.boxHalfSize;
        s.halfU = h.x() * std::fabs(s.u.x()) + h.y() * std::fabs(s.u.y()) + h.z() * std::fabs(s.u.z());
        s.halfV = h.x() * std::fabs(s.v.x()) + h.y() * std::fabs(s.v.y()) + h.z() * std::fabs(s.v.z());
        return s;
    }

    const G4double sinTh = upPerp.mag() > 1e-12 ? upPerp.mag() : 0.;
    const G4double cosTh = std::fabs(s.axis.z()) > 1e-12 ? std::fabs(s.axis.z()) : 0.;

    s.halfLength = s.cylinderHalfHeight * sinTh;
    s.capHalfU = s.cylinderRadius * cosTh;
    s.halfU = s.halfLength + s.capHalfU;
    s.halfV = s.cylinderRadius;

    return s;
}


G4bool GenSurface::InSilhouette(const G4double a, const G4double b) const {
    if (box) {
        // the ray through origin + a*u + b*v along the axis crosses the box (slab test)
        const G4ThreeVector p = u * a + v * b;
        G4double tMin = -DBL_MAX;
        G4double tMax = DBL_MAX;
        for (G4int i = 0; i < 3; ++i) {
            if (std::fabs(axis[i]) < 1e-12) {
                if (std::fabs(p[i]) > boxHalfSize[i]) return false;
                continue;
            }
            const G4double t1 = (-boxHalfSize[i] - p[i]) / axis[i];
            const G4double t2 = (boxHalfSize[i] - p[i]) / axis[i];
            tMin = std::max(tMin, std::min(t1, t2));
            tMax = std::min(tMax, std::max(t1, t2));
        }
        return tMin <= tMax;
    }

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
    if (box) {
        const G4ThreeVector h = boxHalfSize / cm;
        return 4.0 * (h.y() * h.z() * std::fabs(axis.x()) + h.x() * h.z() * std::fabs(axis.y()) +
                      h.x() * h.y() * std::fabs(axis.z()));
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
    if (box) {
        const G4int faces = (std::fabs(axis.x()) > 1e-12) + (std::fabs(axis.y()) > 1e-12) +
                            (std::fabs(axis.z()) > 1e-12);
        return faces > 1 ? "box_silhouette" : "rectangle";
    }
    return capHalfU > 0. && halfLength > 0. ? "cylinder_silhouette" : (halfLength > 0. ? "rectangle" : "disk");
}
