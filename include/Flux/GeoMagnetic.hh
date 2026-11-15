#ifndef GEOMAGNETIC_HH
#define GEOMAGNETIC_HH

#include <string>
#include <vector>

#include <G4ThreeVector.hh>

// Геомагнитное поле IGRF и L-оболочка Мак-Илвейна для потока радиационных поясов (Belt).
// Координаты — геоцентрические декартовы, связанные с Землёй (ECEF), в км; поле в нТ.
namespace GeoMagnetic {
    constexpr double kEarthRadius = 6371.2;          // км, опорный радиус IGRF и единица L
    constexpr double kEquatorialRadius = 6378.137;   // км, WGS-84
    constexpr double kFlattening = 1. / 298.257223563;
    constexpr double kGM = 398600.4418;              // км³/с²
    constexpr double kEarthRotation = 7.2921159e-5;  // рад/с
    constexpr double kDipoleMomentAP8 = 31165.3;     // нТ·R_E³ — постоянный момент, в котором заданы L и B0 моделей AP8/AE8

    // "2026-01-01" или "2026.5" -> дробный год
    double FractionalYear(const std::string& date);

    // Геодезические широта и долгота (градусы), высота над эллипсоидом WGS-84 (км) -> ECEF, км
    G4ThreeVector GeodeticToECEF(double latitudeDeg, double longitudeDeg, double altitudeKm);

    class IGRF {
    public:
        // Коэффициенты из файла .shc (формат IGRF-14), линейная интерполяция по времени
        IGRF(const std::string& shcPath, double year);

        [[nodiscard]] G4ThreeVector Field(const G4ThreeVector& ecefKm) const;   // вектор B, нТ
        [[nodiscard]] double DipoleMoment() const;                             // √(g10²+g11²+h11²), нТ·R_E³
        [[nodiscard]] double Year() const { return year; }

    private:
        int nmax{0};
        double year{0.};
        std::vector<double> g, h;   // [n·(nmax+1) + m], нТ

        [[nodiscard]] int Index(const int n, const int m) const { return n * (nmax + 1) + m; }
    };

    struct LShell {
        double L{0.};      // L Мак-Илвейна для частицы, отражающейся в этой точке; +inf — линия уходит за maxL
        double B{0.};      // |B| в точке, нТ
        double Bmin{0.};   // минимум |B| на силовой линии, нТ
        double I{0.};      // второй инвариант I = ∫√(1 − B/B_m) ds между точками отражения, R_E
    };

    // Трассировка силовой линии (RK4) от точки до сопряжённой точки отражения, интеграл I и L по
    // аппроксимации функции Мак-Илвейна (та же, что в GTsimulation/IRBEM). L считается с моментом
    // dipoleMoment: для моделей AP8/AE8 это kDipoleMomentAP8. Трассировка обрывается, если линия
    // уходит дальше maxL радиусов Земли (тогда L = +inf).
    LShell McIlwainL(const IGRF& field, const G4ThreeVector& ecefKm, double dipoleMoment, double maxL = 15.);

    // L по I и B_m (нТ, R_E): L³·B_m/M = 1 + exp(Y(ln(I³·B_m/M)))
    double McIlwainFunction(double I, double Bm, double dipoleMoment);
}

#endif //GEOMAGNETIC_HH
