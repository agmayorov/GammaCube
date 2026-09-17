#include "GenSurfaceVis.hh"

#include <algorithm>
#include <cmath>
#include <vector>

#include <G4Polyline.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>
#include <G4ExtrudedSolid.hh>
#include <G4Sphere.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>
#include <G4Circle.hh>

#include "Configuration.hh"

using namespace Configuration;


static void Segment(G4VVisManager* vis, const G4ThreeVector& a, const G4ThreeVector& b, const G4Colour& colour) {
    G4Polyline line;
    line.push_back(a);
    line.push_back(b);
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


static void Loop(G4VVisManager* vis, const std::vector<G4ThreeVector>& pts, const G4Colour& colour) {
    G4Polyline line;
    for (const auto& p : pts) line.push_back(p);
    line.push_back(pts.front());
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


static void Arc(G4VVisManager* vis, const G4ThreeVector& centre, G4double radius,
                const G4ThreeVector& e1, const G4ThreeVector& e2,
                G4double from, G4double to, const G4Colour& colour) {
    constexpr G4int nSteps = 72;
    G4Polyline line;
    for (G4int i = 0; i <= nSteps; ++i) {
        const G4double a = from + (to - from) * i / nSteps;
        line.push_back(centre + radius * (std::cos(a) * e1 + std::sin(a) * e2));
    }
    line.SetVisAttributes(G4VisAttributes(colour));
    vis->Draw(line);
}


static void Arrow(G4VVisManager* vis, const G4ThreeVector& tail, const G4ThreeVector& dir,
                  G4double length, const G4ThreeVector& e1, const G4ThreeVector& e2, const G4Colour& colour) {
    const G4ThreeVector tip = tail + dir * length;
    const G4double barb = 0.2 * length;

    Segment(vis, tail, tip, colour);
    Segment(vis, tip, tip - dir * barb + e1 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb - e1 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb + e2 * barb * 0.35, colour);
    Segment(vis, tip, tip - dir * barb - e2 * barb * 0.35, colour);
}


static std::vector<G4TwoVector> SilhouetteOutline(const GenSurface& s) {
    std::vector<G4TwoVector> pts;
    if (s.CapHalfU() <= 0.) {
        pts.emplace_back(-s.HalfLength(), -s.HalfV());
        pts.emplace_back(s.HalfLength(), -s.HalfV());
        pts.emplace_back(s.HalfLength(), s.HalfV());
        pts.emplace_back(-s.HalfLength(), s.HalfV());
        return pts;
    }

    constexpr G4int nSteps = 72;
    for (G4int i = 0; i < nSteps; ++i) {
        const G4double t = twopi * (i + 0.5) / nSteps;
        const G4double c = std::cos(t);
        pts.emplace_back((c >= 0. ? s.HalfLength() : -s.HalfLength()) + s.CapHalfU() * c,
                         s.CylinderRadius() * std::sin(t));
    }
    return pts;
}


void GenSurfaceVis::DrawPayloadCylinder(G4VVisManager* vis) {
    G4double r, zMin, zMax;
    GenSurface::PayloadExtent(r, zMin, zMax);
    const G4ThreeVector ex(1., 0., 0.);
    const G4ThreeVector ey(0., 1., 0.);
    const G4Colour colour = G4Colour::Grey();

    Arc(vis, {0., 0., zMin}, r, ex, ey, 0., twopi, colour);
    Arc(vis, {0., 0., zMax}, r, ex, ey, 0., twopi, colour);

    constexpr G4int nLines = 8;
    for (G4int i = 0; i < nLines; ++i) {
        const G4double psi = twopi * i / nLines;
        const G4ThreeVector p = r * (std::cos(psi) * ex + std::sin(psi) * ey);
        Segment(vis, p + G4ThreeVector(0., 0., zMin), p + G4ThreeVector(0., 0., zMax), colour);
    }
}


void GenSurfaceVis::DrawFlat(G4VVisManager* vis, const GenSurface& s) {
    const G4ThreeVector centre = s.Origin() - s.Axis() * s.Standoff();
    const std::vector<G4TwoVector> outline = SilhouetteOutline(s);

    const G4RotationMatrix rotation(s.U(), s.V(), s.Axis());
    const G4double halfThick = 0.002 * s.Standoff();
    const G4ExtrudedSolid sheet("GenSurfaceSheet", outline, halfThick, G4TwoVector(), 1., G4TwoVector(), 1.);
    G4VisAttributes sheetAttributes(G4Colour(1., 1., 0., 0.15));
    sheetAttributes.SetForceSolid(true);
    vis->Draw(sheet, sheetAttributes, G4Transform3D(rotation, centre));

    std::vector<G4ThreeVector> loop;
    for (const auto& p : outline) loop.push_back(centre + s.U() * p.x() + s.V() * p.y());
    Loop(vis, loop, G4Colour::Yellow());

    constexpr G4int nShort = 6;
    const G4double spacing = 2.0 * std::max(std::min(s.HalfU(), s.HalfV()), 1e-3 * mm) / nShort;
    const G4int nU = std::clamp(static_cast<G4int>(std::lround(2.0 * s.HalfU() / spacing)), 2, 16);
    const G4int nV = std::clamp(static_cast<G4int>(std::lround(2.0 * s.HalfV() / spacing)), 2, 16);
    const G4double length = 0.45 * s.Standoff();

    for (G4int i = 0; i < nU; ++i) {
        for (G4int j = 0; j < nV; ++j) {
            const G4double a = (2.0 * (i + 0.5) / nU - 1.0) * s.HalfU();
            const G4double b = (2.0 * (j + 0.5) / nV - 1.0) * s.HalfV();
            if (!s.InSilhouette(a, b)) continue;
            Arrow(vis, centre + s.U() * a + s.V() * b, s.Axis(), length, s.U(), s.V(), G4Colour::Cyan());
        }
    }
}


void GenSurfaceVis::DrawSphere(G4VVisManager* vis, const GenSurface& s) {
    const G4ThreeVector centre = s.Origin();
    const G4double r = s.Radius();
    const G4ThreeVector up(0., 0., 1.);
    const G4ThreeVector ex(1., 0., 0.);
    const G4ThreeVector ey(0., 1., 0.);
    const G4Colour colour = G4Colour::Yellow();

    G4double from = 0.;
    G4double to = pi;
    if (s.Hemisphere() > 0) to = halfpi;
    else if (s.Hemisphere() < 0) from = halfpi;

    const G4Sphere shell("GenSurfaceShell", 0., r, 0., twopi, from, to - from);
    G4VisAttributes shellAttributes(G4Colour(1., 1., 0., 0.15));
    shellAttributes.SetForceSolid(true);
    vis->Draw(shell, shellAttributes, G4Transform3D(G4RotationMatrix(), centre));

    Arc(vis, centre, r, ex, ey, 0., twopi, colour);

    constexpr G4int nMeridians = 12;
    for (G4int i = 0; i < nMeridians; ++i) {
        const G4double psi = twopi * i / nMeridians;
        const G4ThreeVector radial = std::cos(psi) * ex + std::sin(psi) * ey;
        Arc(vis, centre, r, up, radial, from, to, colour);
    }
}


void GenSurfaceVis::DrawPoint(G4VVisManager* vis, const GenSurface& s) {
    const G4ThreeVector p = s.SourcePosition();
    const G4ThreeVector axis = s.Axis();
    const G4double length = 0.5 * GenSurface::BoundingRadius();

    G4Circle marker(p);
    marker.SetScreenSize(8.);
    marker.SetFillStyle(G4VMarker::filled);
    marker.SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
    vis->Draw(marker);

    const G4ThreeVector a = std::fabs(axis.z()) < 0.999 ? G4ThreeVector(0, 0, 1) : G4ThreeVector(1, 0, 0);
    const G4ThreeVector e1 = axis.cross(a).unit();
    const G4ThreeVector e2 = axis.cross(e1).unit();

    Arrow(vis, p, axis, length, e1, e2, G4Colour::Cyan());
    if (!s.IsCone()) return;

    const G4double alpha = s.ConeAngle();
    const G4ThreeVector rimCentre = p + axis * length * std::cos(alpha);
    const G4double rimRadius = length * std::sin(alpha);
    const G4Colour colour = G4Colour::Yellow();

    Arc(vis, rimCentre, rimRadius, e1, e2, 0., twopi, colour);

    constexpr G4int nLines = 12;
    for (G4int i = 0; i < nLines; ++i) {
        const G4double psi = twopi * i / nLines;
        Segment(vis, p, rimCentre + rimRadius * (std::cos(psi) * e1 + std::sin(psi) * e2), colour);
    }
}


void GenSurfaceVis::Draw() {
    G4VVisManager* vis = G4VVisManager::GetConcreteInstance();
    if (!vis) return;

    const GenSurface s = GenSurface::For(fluxDirection);

    DrawPayloadCylinder(vis);
    if (s.IsIsotropic()) DrawSphere(vis, s);
    else if (s.IsPoint()) DrawPoint(vis, s);
    else DrawFlat(vis, s);
}


G4VisExtent GenSurfaceVis::Extent() {
    const GenSurface s = GenSurface::For(fluxDirection);
    return {G4Point3D(s.Origin()), s.OuterRadius() + GenSurface::Margin()};
}
