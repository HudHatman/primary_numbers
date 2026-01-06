#include <iostream>



int main() {
    std::cout << "--- BANIOWSKI WAVE-FUNCTION PRIME SEARCH ---\n";
    std::cout << "Metoda: Detekcja miejsc zerowych (Vortex Nodes)\n\n";
    std::cout << std::left << std::setw(10) << "NODE (n)"
              << std::setw(20) << "Psi(n) Voltage"
              << std::setw(15) << "Density (Sn)"
              << "Status\n";
    std::cout << "------------------------------------------------------------\n";

    // Zakres badawczy (np. wokół Twoich wysp stabilności)
    for (long long n = 2; n <= 200; ++n) {
        long double psi = calculate_wave_function(n);
        int sn = get_vortex_density(psi);

        // KRYTERIUM BANIOWSKIEGO:
        // Liczba pierwsza manifestuje się, gdy napięcie fazowe Psi(n)
        // spada poniżej progu harmonicznego skorelowanego z B_CONST.
        long double threshold = B_CONST / (long double)n;

        if (psi < threshold * 2.5) { // Filtr rezonansowy
            std::cout << std::left << std::setw(10) << n
                      << std::setw(20) << std::fixed << std::setprecision(10) << psi
                      << std::setw(15) << sn;

            if (sn > 25) { // Wysoka gęstość informacyjna
                std::cout << "[VORTEX PRIME]";
            } else {
                std::cout << "[CANDIDATE]";
            }
            std::cout << "\n";
        }
    }

    return 0;
}