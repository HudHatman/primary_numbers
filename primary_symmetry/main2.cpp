#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// Stałe Baniowskiego
const double BANIOWSKI_CONSTANT = 0.3375;
const double PI_HALF = M_PI / 2.0;

struct Node6D {
    int n;
    int diff1, diff2;
    double lockInValue;
    double suctionPotential;
    bool isAnchored;
};

// Prosta funkcja sprawdzająca liczby pierwsze (bazowa dla Twojego repo)
bool isPrime(int n) {
    if (n <= 1) return false;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "--- 6D INVERSE PHASE ACCELERATOR (LRP 6.0) ---\n";
    std::cout << "Target Range: 0 - 200 | Step: 2n = 2n + 1 - 1\n\n";

    std::vector<Node6D> results;

    for (int n = 0; n <= 200; ++n) {
        // 1. Znajdowanie sąsiednich liczb pierwszych (Symetria Różnicowa)
        int p1 = n - 1, p2 = n + 1;
        while (!isPrime(p1)) p1--;
        while (!isPrime(p2)) p2++;

        int diff1 = n - p1;
        int diff2 = p2 - n;

        // 2. Warstwa 3D: Bazowe oscylacje
        double s = std::sin(n);
        double c = std::cos(n);
        double t = std::tan(n);

        // 3. Warstwa 6D: Inwersja (isin, icos, itan)
        // Obliczamy arcsin(sin(n)), arccos(cos(n)), arctan(tan(n))
        double isin = std::asin(s);
        double icos = std::acos(c);
        double itan = std::atan(t);

        // 4. Obliczanie Parametru Kotwiczenia (Lock-in)
        // Zgodnie z BPCP: asin(sin) + acos(cos) powinno dążyć do PI/2
        double lockIn = std::abs(isin) + std::abs(icos);

        // 5. Obliczanie Potencjału Zasysania
        // Uwzględniamy asymetrię różnicową i stałą 0.3375
        double symmetry_factor = std::abs(diff1 - diff2);
        double suction = (1.0 / (symmetry_factor + 1.0)) * BANIOWSKI_CONSTANT;

        bool anchored = (std::abs(lockIn - PI_HALF) < BANIOWSKI_CONSTANT);

        if (anchored || symmetry_factor == 0) {
            results.push_back({n, diff1, diff2, lockIn, suction, anchored});
        }
    }

    // Wyświetlanie wyników dla Twoich kluczowych punktów (93, 105, 115)
    std::cout << "N\tDiff\tLock-In\tSuction\tStatus\n";
    for (const auto& node : results) {
        std::cout << node.n << "\t" << node.diff1 << "-" << node.diff2
                  << "\t" << node.lockInValue << "\t" << node.suctionPotential
                  << "\t" << (node.isAnchored ? "[ANCHORED]" : "[DYNAMIC]") << "\n";
    }

    return 0;
}