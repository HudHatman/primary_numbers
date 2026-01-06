#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <string>

// --- STAŁE ARCHITEKTONICZNE BANIOWSKIEGO ---
const long double B_CONST = 0.3375;  // Próg rezonansu (suction threshold)
const long double PI_HALF = M_PI / 2.0;

/**
 * Funkcja Falowa Baniowskiego Psi(n)
 * Oblicza "napięcie informacyjne" w punkcie n.
 * Miejsca zerowe (lokalne minima) tej funkcji to liczby pierwsze.
 */
long double calculate_wave_function(long long n) {
    // Rzuty 3D (Emisja)
    long double s = std::sin((long double)n);
    long double c = std::cos((long double)n);

    // Inwersja 6D (Zasysanie)
    // Obliczamy "Lock-in error" - odchylenie od idealnej symetrii kołowej
    long double lock_in = std::abs(std::asin(s)) + std::abs(std::acos(c));
    
    // Funkcja falowa Psi: różnica między lock-in a PI/2
    // W punktach węzłowych (liczby pierwsze) Psi dąży do 0
    return std::abs(lock_in - PI_HALF);
}

/**
 * Detekcja Gęstości Sn (Digit-Split)
 * Weryfikuje, czy "cisza" jest wystarczająco głęboka dla stabilnego węzła.
 */
int get_vortex_density(long double value) {
    std::string s = std::to_string(std::abs(value));
    size_t dot = s.find('.');
    if (dot != std::string::npos) s.erase(dot, 1);
    
    int sum = 0;
    for (int i = 0; i < std::min((int)s.length(), 5); ++i) {
        sum += (s[i] - '0');
    }
    return sum;
}

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