#include "Flux/GeoMagnetic.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace GeoMagnetic {
    double FractionalYear(const std::string& date) {
        int y = 0, mo = 1, d = 1;
        char dash1 = 0, dash2 = 0;
        std::istringstream in(date);
        if (date.find('-') == std::string::npos) {
            double year = 0.;
            if (!(in >> year)) throw std::runtime_error("Неверная дата: " + date);
            return year;
        }
        if (!(in >> y >> dash1 >> mo >> dash2 >> d) || dash1 != '-' || dash2 != '-' || mo < 1 || mo > 12 || d < 1 ||
            d > 31)
            throw std::runtime_error("Неверная дата (нужно ГГГГ-ММ-ДД): " + date);

        const bool leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
        int days[] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        int dayOfYear = d - 1;
        for (int i = 0; i < mo - 1; ++i) dayOfYear += days[i];
        return y + dayOfYear / (leap ? 366. : 365.);
    }

    G4ThreeVector GeodeticToECEF(const double latitudeDeg, const double longitudeDeg, const double altitudeKm) {
        const double lat = latitudeDeg * M_PI / 180.;
        const double lon = longitudeDeg * M_PI / 180.;
        const double e2 = kFlattening * (2. - kFlattening);
        const double N = kEquatorialRadius / std::sqrt(1. - e2 * std::sin(lat) * std::sin(lat));
        return {(N + altitudeKm) * std::cos(lat) * std::cos(lon), (N + altitudeKm) * std::cos(lat) * std::sin(lon),
                (N * (1. - e2) + altitudeKm) * std::sin(lat)};
    }

    // ---------------------------------------------------------------- IGRF

    IGRF::IGRF(const std::string& shcPath, const double requestedYear) : year(requestedYear) {
        std::ifstream in(shcPath);
        if (!in) throw std::runtime_error("Не удалось открыть " + shcPath);

        std::string line;
        auto nextDataLine = [&]() {
            while (std::getline(in, line))
                if (line.find_first_not_of(" \t\r") != std::string::npos && line[line.find_first_not_of(" \t")] != '#')
                    return true;
            return false;
        };

        // заголовок: nmin nmax ntimes spline_order steps [start end], затем строка эпох
        int nmin = 0, ntimes = 0;
        if (!nextDataLine()) throw std::runtime_error(shcPath + ": нет заголовка");
        std::istringstream(line) >> nmin >> nmax >> ntimes;
        if (nmax < 1 || ntimes < 1) throw std::runtime_error(shcPath + ": неверный заголовок");

        std::vector<double> epochs(ntimes);
        if (!nextDataLine()) throw std::runtime_error(shcPath + ": нет строки эпох");
        std::istringstream times(line);
        for (auto& t : epochs) times >> t;

        if (year < epochs.front() || year > epochs.back()) {
            std::ostringstream msg;
            msg << "Дата " << year << " вне диапазона IGRF " << epochs.front() << "–" << epochs.back();
            throw std::runtime_error(msg.str());
        }
        size_t k = 0;
        while (k + 2 < epochs.size() && year > epochs[k + 1]) ++k;
        const double w = epochs.size() > 1 ? (year - epochs[k]) / (epochs[k + 1] - epochs[k]) : 0.;

        g.assign((nmax + 1) * (nmax + 1), 0.);
        h.assign((nmax + 1) * (nmax + 1), 0.);
        int rows = 0;
        while (nextDataLine()) {
            std::istringstream row(line);
            int n = 0, m = 0;
            row >> n >> m;
            std::vector<double> v(ntimes);
            for (auto& x : v) row >> x;
            if (!row || n < 1 || n > nmax || std::abs(m) > n) throw std::runtime_error(shcPath + ": неверная строка: " + line);
            const double value = epochs.size() > 1 ? v[k] + w * (v[k + 1] - v[k]) : v[0];
            (m >= 0 ? g : h)[Index(n, std::abs(m))] = value;
            ++rows;
        }
        if (rows != nmax * (nmax + 2)) throw std::runtime_error(shcPath + ": неполный набор коэффициентов");
    }

    double IGRF::DipoleMoment() const {
        return std::sqrt(g[Index(1, 0)] * g[Index(1, 0)] + g[Index(1, 1)] * g[Index(1, 1)] +
                         h[Index(1, 1)] * h[Index(1, 1)]);
    }

    // B = −∇V, V = a Σ (a/r)^{n+1} Σ (g cos mφ + h sin mφ) P_n^m(cos θ), P — квазинормированные по Шмидту
    G4ThreeVector IGRF::Field(const G4ThreeVector& x) const {
        const double r = x.mag();
        const double cosT = x.z() / r;
        const double sinT = std::max(std::sqrt(std::max(0., 1. - cosT * cosT)), 1e-12);
        const double phi = std::atan2(x.y(), x.x());

        const int size = (nmax + 1) * (nmax + 1);
        thread_local std::vector<double> P, dP, cosM, sinM;
        P.assign(size, 0.);
        dP.assign(size, 0.);
        cosM.resize(nmax + 1);
        sinM.resize(nmax + 1);
        for (int m = 0; m <= nmax; ++m) {
            cosM[m] = std::cos(m * phi);
            sinM[m] = std::sin(m * phi);
        }

        P[Index(0, 0)] = 1.;
        for (int n = 1; n <= nmax; ++n) {
            // диагональ P_n^n
            if (n == 1) {
                P[Index(1, 1)] = sinT;
                dP[Index(1, 1)] = cosT;
            } else {
                const double k = std::sqrt((2. * n - 1.) / (2. * n));
                P[Index(n, n)] = k * sinT * P[Index(n - 1, n - 1)];
                dP[Index(n, n)] = k * (cosT * P[Index(n - 1, n - 1)] + sinT * dP[Index(n - 1, n - 1)]);
            }
            // P_n^m, m < n
            for (int m = 0; m < n; ++m) {
                const double denom = std::sqrt(static_cast<double>(n * n - m * m));
                const double a = (2. * n - 1.) / denom;
                const double b = n >= 2 ? std::sqrt(static_cast<double>((n - 1) * (n - 1) - m * m)) / denom : 0.;
                const double p2 = n >= 2 && m <= n - 2 ? P[Index(n - 2, m)] : 0.;
                const double dp2 = n >= 2 && m <= n - 2 ? dP[Index(n - 2, m)] : 0.;
                P[Index(n, m)] = a * cosT * P[Index(n - 1, m)] - b * p2;
                dP[Index(n, m)] = a * (cosT * dP[Index(n - 1, m)] - sinT * P[Index(n - 1, m)]) - b * dp2;
            }
        }

        double Br = 0., Bt = 0., Bp = 0.;
        const double ratio = kEarthRadius / r;
        double ar = ratio * ratio;   // (a/r)^{n+2}, n = 0
        for (int n = 1; n <= nmax; ++n) {
            ar *= ratio;
            for (int m = 0; m <= n; ++m) {
                const int i = Index(n, m);
                const double gh = g[i] * cosM[m] + h[i] * sinM[m];
                Br += (n + 1) * ar * gh * P[i];
                Bt -= ar * gh * dP[i];
                Bp += ar * m * (g[i] * sinM[m] - h[i] * cosM[m]) * P[i] / sinT;
            }
        }

        const double cosP = std::cos(phi), sinP = std::sin(phi);
        const double rho = Br * sinT + Bt * cosT;   // компонента в экваториальной плоскости
        return {rho * cosP - Bp * sinP, rho * sinP + Bp * cosP, Br * cosT - Bt * sinT};
    }

    // ---------------------------------------------------------------- L Мак-Илвейна

    double McIlwainFunction(const double I, const double Bm, const double dipoleMoment) {
        // Коэффициенты аппроксимации Y(X) функции Мак-Илвейна по шести диапазонам X (как в GTsimulation/IRBEM)
        static const double a[10][6] = {
            {3.0062102e-1, 6.2337691e-1, 6.228644e-1, 6.222355e-1, 2.0007187, -3.0460681},
            {3.33338e-1, 4.3432642e-1, 4.3352788e-1, 4.3510529e-1, -1.8461796e-1, 1},
            {0, 1.5017245e-2, 1.4492441e-2, 1.2817956e-2, 1.2038224e-1, 0},
            {0, 1.3714667e-3, 1.1784234e-3, 2.1680398e-3, -6.7310339e-3, 0},
            {0, 8.2711096e-5, 3.8379917e-5, -3.2077032e-4, 2.170224e-4, 0},
            {0, 3.2916354e-6, -3.3408822e-6, 7.9451313e-5, -3.8049276e-6, 0},
            {0, 8.1048663e-8, -5.3977642e-7, -1.2531932e-5, 2.8212095e-8, 0},
            {0, 1.0066362e-9, -2.1997983e-8, 9.9766148e-7, 0, 0},
            {0, 8.3232531e-13, 2.3028767e-9, -3.958306e-8, 0, 0},
            {0, -8.1537735e-14, 2.6047023e-10, 6.3271665e-10, 0, 0}};

        const double X = I > 0. ? std::log(I * I * I * Bm / dipoleMoment) : -50.;
        const int col = X < -22. ? 0 : X < -3. ? 1 : X < 3. ? 2 : X < 11.7 ? 3 : X < 23. ? 4 : 5;
        double Y = 0., power = 1.;
        for (const auto& row : a) {
            Y += row[col] * power;
            power *= X;
        }
        return std::cbrt(dipoleMoment / Bm * (1. + std::exp(Y)));
    }

    LShell McIlwainL(const IGRF& field, const G4ThreeVector& start, const double dipoleMoment, const double maxL) {
        LShell out;
        const G4ThreeVector b0 = field.Field(start);
        const double Bm = b0.mag();
        out.B = Bm;
        out.Bmin = Bm;

        // идём вдоль линии туда, где |B| убывает; если растёт в обе стороны — точка в минимуме, I = 0
        const double probe = 1.;   // км
        double sign = 1.;
        if (field.Field(start + probe * b0.unit()).mag() >= Bm) {
            if (field.Field(start - probe * b0.unit()).mag() >= Bm) {
                out.L = McIlwainFunction(0., Bm, dipoleMoment);
                return out;
            }
            sign = -1.;
        }

        G4ThreeVector x = start;
        G4ThreeVector bx = b0;
        double Bprev = Bm, qprev = 0., I = 0.;
        auto direction = [&](const G4ThreeVector& p) { return sign * field.Field(p).unit(); };

        for (int step = 0; step < 100000; ++step) {
            const double r = x.mag();
            if (r > maxL * kEarthRadius) {
                out.L = std::numeric_limits<double>::infinity();
                return out;
            }
            const double ds = std::clamp(0.005 * r, 1., 100.);   // км: ~32 км у Земли
            const G4ThreeVector k1 = sign * bx.unit();
            const G4ThreeVector k2 = direction(x + 0.5 * ds * k1);
            const G4ThreeVector k3 = direction(x + 0.5 * ds * k2);
            const G4ThreeVector k4 = direction(x + ds * k3);
            const G4ThreeVector xn = x + ds / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
            const G4ThreeVector bn = field.Field(xn);
            const double Bn = bn.mag();

            if (Bn >= Bm) {
                // сопряжённая точка отражения внутри шага: подынтегральное √(1 − B/B_m) спадает до нуля как √s
                const double t = Bn > Bprev ? std::clamp((Bm - Bprev) / (Bn - Bprev), 0., 1.) : 1.;
                I += 2. / 3. * qprev * t * ds;
                out.I = I / kEarthRadius;
                out.L = McIlwainFunction(out.I, Bm, dipoleMoment);
                return out;
            }
            const double q = std::sqrt(1. - Bn / Bm);
            // первый шаг начинается в точке отражения (q = 0, рост как √s), дальше — трапеции
            I += step == 0 ? 2. / 3. * q * ds : 0.5 * (qprev + q) * ds;
            out.Bmin = std::min(out.Bmin, Bn);
            x = xn;
            bx = bn;
            Bprev = Bn;
            qprev = q;
        }
        throw std::runtime_error("McIlwainL: силовая линия не замкнулась за 100000 шагов");
    }
}
