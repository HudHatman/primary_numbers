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
    int angle;                      // Numer porządkowy (Z)
    vector<vector<string>> table;   // Tabela 2D bloków
    int table_sum_provided;         // Dostarczona suma kontrolna
    bool is_sum_prime;

    // Protony i Elektrony = Wartość bezwzględna numeru angle
    int getProtons() const { return abs(angle); }
    int getElectrons() const { return abs(angle); }

    // Obliczanie Neutronów (N)
    // Według algorytmu, jądro (neutrony) jest zakodowane w pierwszych 3 rzędach tabeli
    int getNeutrons() const {
        int n_sum = 0;
        // Sumujemy wartości (sumy cyfr) bloków z rzędów 0, 1 i 2
        for (int i = 0; i < 3 && i < table.size(); ++i) {
            for (const string& block : table[i]) {
                n_sum += getBlockValue(block);
            }
        }
        return n_sum;
    }

    // Weryfikacja spójności danych: czy suma cyfr wszystkich bloków = table_sum
    bool verifyTableSum() const {
        int total = 0;
        for (const auto& row : table) {
            for (const string& block : row) {
                total += getBlockValue(block);
            }
        }
        return total == table_sum_provided;
    }
};

// --- FUNKCJA RAPORTUJĄCA ---

void generateElementReport(const ElementData& e) {
    cout << ">>> RAPORT DLA PIERWIASTKA (Angle: " << e.angle << ") <<<" << endl;
    cout << "Liczba Protonów (P):  " << e.getProtons() << endl;
    cout << "Liczba Elektronów (E): " << e.getElectrons() << endl;
    cout << "Liczba Neutronów (N):  " << e.getNeutrons() << " [Wyliczone z bloków jądrowych]" << endl;
    cout << "Suma Kontrolna:        " << (e.verifyTableSum() ? "ZGODNA" : "BŁĄD")
         << " (" << e.table_sum_provided << ")" << endl;
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
    eNeon.angle = 10;
    eNeon.table_sum_provided = 29;
    eNeon.table = {
        {"1", "7"},
        {"333", "2222"},
        {"0", "1"},
        {"1", "0"},
        {"1", "0"},
        {"1", "0"}
    };

    generateElementReport(eNeg200);
    generateElementReport(eNeon);

    return 0;
}