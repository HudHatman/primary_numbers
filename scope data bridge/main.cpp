#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

// Stałe rezonansu
const long double B_CONST = 0.3375;
const long double PI_HALF = M_PI / 2.0;

struct ScopeFrame {
    long long n;
    double ch1_emission;  // Sin(n) - surowa fala 3D
    double ch2_inversion; // Asin(Sin(n)) - echo 6D
    double ch3_voltage;   // Psi(n) - napięcie Baniowskiego
    double ch4_density;   // Sn - gęstość informacyjna
};

int main() {
    std::string filename = "vortex_scope_data.csv";
    std::ofstream file(filename);

    // Nagłówek zgodny ze standardem importu oscyloskopów (CSV)
    file << "Time,CH1,CH2,CH3,CH4\n";

    std::cout << "Generowanie danych dla oscyloskopu (N=2 do 400)..." << std::endl;

    for (long long n = 2; n <= 400; ++n) {
        long double s = std::sin((long double)n);
        long double c = std::cos((long double)n);
        
        // Obliczanie kanałów
        double ch1 = (double)s;
        double ch2 = (double)(std::abs(std::asin(s)) + std::abs(std::acos(c)));
        double ch3 = (double)std::abs(ch2 - (double)PI_HALF);
        
        // Gęstość Sn skalowana dla czytelności na oscyloskopie (0-1V)
        double ch4 = (double)(std::abs(s) * B_CONST + std::abs(c) * (1.0 - B_CONST));

        // Zapis do pliku: n jako pseudo-czas (oś X)
        file << n << "," 
             << std::fixed << std::setprecision(6) 
             << ch1 << "," 
             << ch2 << "," 
             << ch3 << "," 
             << ch4 << "\n";
    }

    file.close();
    std::cout << "Eksport zakończony: " << filename << std::endl;
    return 0;
}