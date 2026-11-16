#ifndef GEOMAGNETIC_HH
#define GEOMAGNETIC_HH

#include <string>
#include <vector>

#include <G4ThreeVector.hh>

namespace GeoMagnetic {
    constexpr double kEarthRadius = 6371.2;
    constexpr double kEquatorialRadius = 6378.137;
    constexpr double kFlattening = 1. / 298.257223563;
    constexpr double kGM = 398600.4418;
    constexpr double kEarthRotation = 7.2921159e-5;
    constexpr double kDipoleMomentAP8 = 31165.3;

    double FractionalYear(const std::string& date);

    G4ThreeVector GeodeticToECEF(double latitudeDeg, double longitudeDeg, double altitudeKm);

    class IGRF {
    public:
        IGRF(const std::string& shcPath, double year);

        [[nodiscard]] G4ThreeVector Field(const G4ThreeVector& ecefKm) const;
        [[nodiscard]] double DipoleMoment() const;
        [[nodiscard]] double Year() const { return year; }

    private:
        int nmax{0};
        double year{0.};
        std::vector<double> g, h;

        [[nodiscard]] int Index(const int n, const int m) const { return n * (nmax + 1) + m; }
    };

    struct LShell {
        double L{0.};
        double B{0.};
        double Bmin{0.};
        double I{0.};
    };

    LShell McIlwainL(const IGRF& field, const G4ThreeVector& ecefKm, double dipoleMoment, double maxL = 15.);

    double McIlwainFunction(double I, double Bm, double dipoleMoment);
}

#endif //GEOMAGNETIC_HH
