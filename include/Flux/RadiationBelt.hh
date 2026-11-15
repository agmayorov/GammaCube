#ifndef RADIATIONBELT_HH
#define RADIATIONBELT_HH

#include <map>
#include <string>
#include <vector>

// Поток Belt: захваченные частицы радиационных поясов по трёхмерным картам J(α_eq, E, L)
// (../RadiationBelts/<модель>.txt, их пишет Processing/convert_belt_mat.py).
//
// В точке с полем B на оболочке L частица с локальным питч-углом α имеет на экваторе угол
// sin²α_eq = sin²α · B0/B, B0 = M/L³, и по теореме Лиувилля J_loc(E, α) = J(E, α_eq, L).
// Изотропная генерация получает усреднённую по направлениям интенсивность
//     <J>(E) = (1/4π) ∫ J_loc dΩ = ∫₀^{π/2} J(E, α_eq(α), L) sin α dα,   см⁻²·с⁻¹·ср⁻¹·МэВ⁻¹,
// а для орбиты — её среднее по времени. Спектр пишется в CSV, дальше поток работает как Table.
namespace RadiationBelt {
    class Model {
    public:
        explicit Model(const std::string& path);

        // <J>(E) на сетке Energies() для оболочки L и отношения B/B0 (≥ 1 у захваченных частиц)
        [[nodiscard]] std::vector<double> MeanIntensity(double L, double BoverB0) const;

        [[nodiscard]] const std::string& Name() const { return name; }
        [[nodiscard]] const std::string& Particle() const { return particle; }
        [[nodiscard]] const std::vector<double>& Energies() const { return energy; }
        [[nodiscard]] double Lmin() const { return L.front(); }
        [[nodiscard]] double Lmax() const { return L.back(); }

    private:
        std::string name, particle;
        std::vector<double> alpha, energy, L;   // α_eq в градусах, E в МэВ
        std::vector<double> J;                  // [(iα·nE + iE)·nL + iL]

        [[nodiscard]] double At(const size_t ia, const size_t ie, const size_t il) const {
            return J[(ia * energy.size() + ie) * L.size() + il];
        }
    };

    struct Result {
        std::string model, particle;
        std::string position;          // point | geo | orbit
        std::string description;       // параметры положения одной строкой для info_*.txt
        std::string tag;               // для имён файлов: PSB97_orbit_h500_i51.6
        std::string spectrumPath;      // CSV «E, <J>» — вход для TableFlux и CountRates
        std::string profilePath;       // CSV по точкам орбиты (только для orbit)
        double year{0.};
        double dipoleMoment{0.};
        double L{0.}, B{0.}, BoverB0{0.};   // для point и geo
        int points{1};
        double trappedFraction{1.};          // доля точек орбиты с ненулевым потоком
        double omniFlux{0.};                 // 4π∫<J>dE в диапазоне E_min…E_max, см⁻²·с⁻¹
        double modelEmin{0.}, modelEmax{0.}; // диапазон энергий модели, МэВ
    };

    // Читает ../Flux_config/Belt_params.txt, считает спектр и пишет CSV в текущий каталог.
    // EminMeV/EmaxMeV — уже ограниченный порогом диапазон запуска (для интеграла в Result).
    Result Prepare(const std::string& configPath, double EminMeV, double EmaxMeV);

    // Параметры «ключ: значение» (# — комментарий); общий разбор для Loader и Prepare
    std::map<std::string, std::string> ReadParams(const std::string& path);
}

#endif //RADIATIONBELT_HH
