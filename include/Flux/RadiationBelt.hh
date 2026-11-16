#ifndef RADIATIONBELT_HH
#define RADIATIONBELT_HH

#include <map>
#include <string>
#include <vector>

namespace RadiationBelt {
    class Model {
    public:
        explicit Model(const std::string& path);

        [[nodiscard]] std::vector<double> MeanIntensity(double L, double BoverB0) const;

        [[nodiscard]] const std::string& Name() const { return name; }
        [[nodiscard]] const std::string& Particle() const { return particle; }
        [[nodiscard]] const std::vector<double>& Energies() const { return energy; }
        [[nodiscard]] double Lmin() const { return L.front(); }
        [[nodiscard]] double Lmax() const { return L.back(); }

    private:
        std::string name, particle;
        std::vector<double> alpha, energy, L;
        std::vector<double> J;

        [[nodiscard]] double At(const size_t ia, const size_t ie, const size_t il) const {
            return J[(ia * energy.size() + ie) * L.size() + il];
        }
    };

    struct Result {
        std::string model, particle;
        std::string position;
        std::string description;
        std::string tag;
        std::string spectrumPath;
        std::string profilePath;
        double year{0.};
        double dipoleMoment{0.};
        double L{0.}, B{0.}, BoverB0{0.};
        int points{1};
        double trappedFraction{1.};
        double omniFlux{0.};
        double modelEmin{0.}, modelEmax{0.};
    };

    struct ManualPoint {
        double L{-1.}, B{-1.}, BoverB0{-1.};
    };

    Result Prepare(const std::string& configPath, double EminMeV, double EmaxMeV, const ManualPoint& manual = {});

    std::map<std::string, std::string> ReadParams(const std::string& path);
}

#endif //RADIATIONBELT_HH
