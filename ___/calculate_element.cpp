#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>

using namespace std;

// --- ELEMENTARNE FUNKCJE LOGICZNE ---

// "Wybór odpowiedniej liczby" z bloku (Suma cyfr)
// Przekształca kod wizualny (np. "2222") na wartość obliczeniową (8)
int getBlockValue(string block) {
    int sum = 0;
    for (char c : block) {
        if (isdigit(c)) {
            sum += (c - '0');
        }
    }
    return sum;
}

// Sprawdzanie czy liczba jest pierwsza
bool isPrime(int n) {
    if (n <= 1) return false;
    for (int i = 2; i <= sqrt(n); i++) {
        if (n % i == 0) return false;
    }
    return true;
}

// --- STRUKTURA DANYCH PIERWIASTKA ---


struct ElementData {
    int angle;
    std::vector<std::vector<std::string>> table;
    int table_sum_provided;
    bool is_sum_prime;

    int getProtons() const { return abs(angle); }
    int getElectrons() const { return abs(angle); }

    int getNeutrons() const {
        int n_sum = 0;
        for (int i = 0; i < table.size() && i < 3; ++i) {
            for (int j = 0; j < table[i].size(); ++j) {
                for (int k = 0; k < table[i][j].size(); ++k) {
                    n_sum += table[i][j][k];
                }
            }
        }
        return n_sum;
    }

    double getAtomicMass() const {
        double protonMass = 1.007276;
        double neutronMass = 1.008665;
        double electronMass = 0.000548;

        int P = getProtons();
        int N = getNeutrons();
        int E = getElectrons();

        // Suma mas składników (bez uwzględnienia deficytu masy / energii wiązania)
        return (P * protonMass) + (N * neutronMass) + (E * electronMass);
    }
};


// --- FUNKCJA RAPORTUJĄCA ---

void generateElementReport(const ElementData& e) {
    cout << ">>> RAPORT DLA PIERWIASTKA (Angle: " << e.angle << ") <<<" << endl;
    cout << "Liczba Protonów (P):  " << e.getProtons() << endl;
    cout << "Liczba Elektronów (E): " << e.getElectrons() << endl;
    cout << "Masa atomowa: " << e.getAtomicMass() << endl;
    cout << "Liczba Neutronów (N):  " << e.getNeutrons() << " [Wyliczone z bloków jądrowych]" << endl;
    cout << "Suma jest Pierwsza:    " << (isPrime(e.table_sum_provided) ? "TAK" : "NIE") << endl;
    cout << "-----------------------------------------------" << endl << endl;
}

int main() {
    // PRZYKŁAD 1: Dane z Angle -200 (Pierwszy rekord z data.txt)
    ElementData eNeg200;
    eNeg200.angle = -200;
    eNeg200.table_sum_provided = 59;
    eNeg200.table = {
        {"3", "22", "2", "0"},       // Rząd 0 (Jądro)
        {"333", "3", "333", "33"},   // Rząd 1 (Jądro)
        {"0", "3", "33", "3"},       // Rząd 2 (Jądro)
        {"2", "0", "0", "0"},       // Rząd 3 (Warstwa e-)
        {"1", "33", "0", "0"},       // Rząd 4 (Warstwa e-)
        {"2", "0", "0", "0"}        // Rząd 5 (Warstwa e-)
    };

    // PRZYKŁAD 2: Dane z Angle 10 (Neon)
    ElementData eNeon;
    eNeon.angle = 14;
    eNeon.table_sum_provided = 29;
    eNeon.table = {
        {"0", "2"},
        {"333", "7"},
        {"2", "22"},
        {"1", "22"},
        {"1", "22"},
        {"1", "22"}
    };

    generateElementReport(eNeg200);
    generateElementReport(eNeon);

    return 0;
}