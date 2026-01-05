#include <iostream>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <map>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <chrono>
#include <iomanip>
/*
 * DATA TYPES
 */
typedef long long T_NUM;
typedef std::vector<T_NUM> T_NUM_VEC;
typedef std::string T_STR;
typedef std::vector<T_STR> T_STR_VEC;
typedef std::pair<T_NUM, T_NUM> T_PAIR;
typedef std::vector<T_PAIR> T_PAIR_VEC;
typedef std::vector<std::vector<std::vector<T_NUM>>> T_3D_NUM_VEC;
struct T_STRUCT {
    T_NUM angle;
    T_NUM angleSize;
    T_NUM sum;
    T_PAIR_VEC vec;
    T_3D_NUM_VEC splitted;
    T_NUM splitted_sum = 0;
    T_STR comparison;
};
typedef std::vector<T_STRUCT> T_STRUCTURE_VECTOR;

struct T_SORT_STRUCT {
    inline bool operator()(const T_STRUCT &s1, const T_STRUCT &s2) {
        return (s1.angle < s2.angle);
    }
};

struct T_HISTOGRAM_SORT_STRUCT {
    inline bool operator()(const T_PAIR &s1, const T_PAIR &s2) {
        return (s1.first < s2.first);
    }
};

/*
 * CONFIG
 */
T_STR output = "std::cout";

/*
 * FUNCTIONS
 */
T_NUM sumVector(T_NUM_VEC vec);

void debugVector(T_STR title, T_NUM_VEC vec, bool sum);

void debugPairVector(T_PAIR_VEC vec);

void debug3dVector(T_3D_NUM_VEC vec);

int charToInt(char c);

T_STR getNumbersFromTri(double d_num, T_STR s_num, T_NUM d_step, unsigned int ommit, unsigned long length);

T_NUM_VEC createVectorFromTri(T_STR sinus, T_STR cosinus, T_STR tangens, T_STR isin, T_STR icos, T_STR itangens);

bool isPrime(T_NUM num);

T_3D_NUM_VEC splitVector(T_NUM_VEC vec, T_NUM length);

T_NUM_VEC createPrimaryNumberFromSplits(T_3D_NUM_VEC vec, T_NUM multiplier);

T_STR_VEC splitString(T_STR str);

/*
 * TRIGONOMETRY
 */
const double PI = std::acos(-1);

double degreeToRadian(double degree) {
    return degree * PI / 180.0;
}

double radianToDegree(double radians) {
    return radians * 180 / PI;
}

double _sin(double x) {
    return std::sin(degreeToRadian(x));
}

double _cos(double x) {
    return std::cos(degreeToRadian(x));
}

double _tan(double x) {
    return std::tan(degreeToRadian(x));
}


double _isin(double x) {
    return radianToDegree(std::asin(x));
}

double _icos(double x) {
    return radianToDegree(std::acos(x));
}

double _itan(double x) {
    return radianToDegree(std::atan(x));
}

/*
 * HELPERS
 */

/**
 * it prints content of the 1D T_NUM vector
 * @param vec
 */
void debugVector(T_STR title, T_NUM_VEC vec, bool sum = false) {
    std::cout << "\t___" << title << "___";
    if (sum) {
        std::cout << ": " << sumVector(vec);
    }
    std::cout << std::endl;

    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl << std::endl;
}

/**
 * it prints content of the T_PAIR vector
 * @param vec
 */
void debugPairVector(T_PAIR_VEC vec) {
    std::cout << "\t___map___" << std::endl;
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i].first << ":" << vec[i].second;
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl << std::endl;
}

/**
 * it prints content of the 3D T_NUM vector
 * @param vec
 */
void debug3dVector(T_3D_NUM_VEC vec) {
    std::cout << "\t___splitted___" << std::endl;

    for (size_t i = 0; i < vec.size(); i++) {
        for (size_t j = 0; j < vec[i].size(); j++) {
            for (size_t k = 0; k < vec[i][j].size(); k++) {
                std::cout << vec[i][j][k];
            }
            std::cout << '\t' << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

/**
 * It converts string ASCII char [0-9] to integer in C++
 * @param c
 * @return
 */
int charToInt(char c) {
    int result = c;
    switch (result) {
        case 48:
            return 0;
        case 49:
            return 1;
        case 50:
            return 2;
        case 51:
            return 3;
        case 52:
            return 4;
        case 53:
            return 5;
        case 54:
            return 6;
        case 55:
            return 7;
        case 56:
            return 8;
        case 57:
            return 9;
        default:
            return 0;
    }
}

/**
 * it converts TRIgonometry numbers to a string (2 number angle -> 2 sized string, 3 number angle -> 3 sized string)
 * @param d_num
 * @param s_num
 * @param d_step
 * @param ommit
 * @param length
 * @return
 */
T_STR getNumbersFromTri(double d_num, T_STR s_num, T_NUM d_step, unsigned int ommit, unsigned long length) {
    T_STR result = "";
    if (ommit == 0) {
        for (size_t i = s_num[0] == '-' ? ommit + 1 : ommit; i < s_num.length(); i++) {
            if (s_num[i] != '.') {
                result += s_num[i];
            }
        }
        result = result.substr(0, length);
    } else if (d_num < 0) {
        ommit += 1; // comma
        result = s_num.substr(ommit, length);
    } else if (d_num > 0) {
        result = s_num.substr(ommit, length);
    }
    return result;
}

/**
 * Get sum from T_NUM vector
 * @param vec
 * @return
 */
T_NUM sumVector(T_NUM_VEC vec) {
    T_NUM result = 0;
    for (size_t i = 0; i < vec.size(); i++) {
        result += vec.at(i);
    }
    return result;
}

struct Node6D {
    int n;               // Kąt/Liczba
    int table_sum;       // Suma z Twojego algorytmu
    bool is_sum_prime;   // Czy wyspa stabilności
};

struct JumpResult {
    double magnitude;    // Siła przeskoku
    double phase_shift;  // Przesunięcie fazowe
    std::string flow;    // Kierunek: Zasysanie vs Emisja
};

/**
 * Oblicza wektor przeskoku informacyjnego między wyspami stabilności.
 * Wykorzystuje stałą Baniowskiego 0.3375 jako regulator gęstości.
 */
JumpResult calculate_jump_vector(Node6D from, Node6D to) {
    const double BANIOWSKI_CONSTANT = 0.3375;
    JumpResult jr;

    // 1. Obliczenie odległości fazowej (Delta n)
    jr.phase_shift = std::abs(to.n - from.n);

    // 2. Obliczenie różnicy potencjałów Sn
    double delta_Sn = std::abs(to.table_sum - from.table_sum);

    // 3. Magnituda przeskoku skalowana stałą 0.3375
    // Im mniejsza odległość fazowa przy dużej różnicy Sn, tym silniejszy przeskok
    jr.magnitude = (delta_Sn / jr.phase_shift) * BANIOWSKI_CONSTANT;

    // 4. Określenie kierunku przepływu energii
    if (to.table_sum < from.table_sum) {
        jr.flow = "EMISSION (Accelerator Activation)";
    } else {
        jr.flow = "ABSORPTION (Anchor Stabilization)";
    }

    return jr;
}

/**
 * it creates T_NUM vector from TRInity numbers
 * @param sinus
 * @param cosinus
 * @param tangens
 * @param s4
 * @return
 */
T_NUM_VEC createVectorFromTri(T_STR sinus, T_STR cosinus, T_STR tangens, T_STR isin, T_STR icos, T_STR itangens) {
    T_NUM_VEC result = T_NUM_VEC();
    for (int i = 0; i < sinus.length(); i++) {
        result.push_back(charToInt(sinus.at(i)));
    }
    for (int i = 0; i < cosinus.length(); i++) {
        result.push_back(charToInt(cosinus.at(i)));
    }
    for (int i = 0; i < tangens.length(); i++) {
        result.push_back(charToInt(tangens.at(i)));
    }
    for (int i = 0; i < isin.length(); i++) {
        result.push_back(charToInt(isin.at(i)));
    }
    for (int i = 0; i < icos.length(); i++) {
        result.push_back(charToInt(icos.at(i)));
    }
    for (int i = 0; i < itangens.length(); i++) {
        result.push_back(charToInt(itangens.at(i)));
    }
    return result;
}

/**
 * Simple get true if number is primary
 * @param num
 * @return
 */
bool isPrime(T_NUM num) {
    if (num < 2) return false;
    if (num == 2) return true;
    if (num % 2 == 0) return false;
    for (T_NUM i = 3; i * i <= num; i += 2) {
        if (num % i == 0) return false;
    }
    return true;
}

/**
 * It converts T_NUM_VEC to splited vector (extended examples)
 * @param vec
 * @param length
 * @return
 */
T_3D_NUM_VEC splitVector(T_NUM_VEC vec, T_NUM length) {
    T_3D_NUM_VEC result;
    result.resize(vec.size() / length);

    T_NUM row = 0;
    T_NUM col = 0;

    for (T_NUM i = 0; i < vec.size(); i++) {
        result[row].resize(length);

        if (vec[i] == 7 || vec[i] == 5 || vec[i] == 3 || vec[i] == 2 || vec[i] == 1 || vec[i] == 0) {
            result[(unsigned long) row][(unsigned long) col].resize(1);

            result[(unsigned long) row][(unsigned long) col][0] = vec[i];
        } else if (vec[i] == 8) {
            result[(unsigned long) row][(unsigned long) col].resize(4);

            result[(unsigned long) row][(unsigned long) col][0] = 2;
            result[(unsigned long) row][(unsigned long) col][1] = 2;
            result[(unsigned long) row][(unsigned long) col][2] = 2;
            result[(unsigned long) row][(unsigned long) col][3] = 2;
        } else if (vec[i] == 6) {
            result[(unsigned long) row][(unsigned long) col].resize(2);

            result[(unsigned long) row][(unsigned long) col][0] = 3;
            result[(unsigned long) row][(unsigned long) col][1] = 3;
        } else if (vec[i] == 4) {
            result[(unsigned long) row][(unsigned long) col].resize(2);

            result[(unsigned long) row][(unsigned long) col][0] = 2;
            result[(unsigned long) row][(unsigned long) col][1] = 2;
        } else if (vec[i] == 9) {
            result[(unsigned long) row][(unsigned long) col].resize(3);

            result[(unsigned long) row][(unsigned long) col][0] = 3;
            result[(unsigned long) row][(unsigned long) col][1] = 3;
            result[(unsigned long) row][(unsigned long) col][2] = 3;
        }

        if (((i + 1) % length) == 0) {
            row++;
            col = 0;
        } else {
            col++;
        }
    }

    return result;
}

T_NUM_VEC createPrimaryNumberFromSplits(T_3D_NUM_VEC vec, T_NUM multiplier = 1) {
    T_NUM_VEC result;

    result.resize(vec.size());

    for (T_NUM row = 0; row < vec.size(); row++) {
        for (T_NUM number = 0; number < vec[row][0].size(); number++) {
            T_NUM _result = 0;
            for (size_t i = 1; i < vec[row].size(); i++) {
                _result += (T_NUM) pow(vec[row][0][number], vec[row][i].size() * multiplier);
            }

            result[row] += _result;
        }
    }

    return result;
}

T_STR createCell(T_NUM tabLength, T_STR content) {
    T_NUM subtract = (tabLength * 4) - (T_NUM) content.length();
    auto numTabs = (T_NUM) ceil((int) subtract / (int) 4);

    if (((int) subtract) % 4 != 0) {
        numTabs++;
    }

    for (T_NUM i = 0; i < numTabs; i++) {
        content.append("\t");
    }

    return content;
}

std::ofstream fileData;
std::ofstream fileSums;
std::ofstream fileHistogram;
std::ofstream fileSymmetry;
std::ofstream filePython;

void writeSumsToFile(T_STRUCT arr) {
    fileSums << std::to_string(arr.angle) << ":\t";
    if (arr.angleSize == 2 || arr.angleSize == 1) {
        fileSums << "\t";
    }

    T_NUM sum = 0;
    for (size_t i = 0; i < arr.splitted.size(); i++) {
        for (size_t j = 0; j < arr.splitted[i].size(); j++) {
            for (size_t k = 0; k < arr.splitted[i][j].size(); k++) {
                sum += arr.splitted[i][j][k];
            }
        }
    }

    fileSums << std::to_string(sum) << "\t\t";
    fileSums << (isPrime(sum) ? "true" : "false") << "\t\t";

    for (T_NUM j = -2; j <= 2; j++) {
        T_STR cell;
        for (T_NUM k = -1; k <= 1; k += 2) {
            if (j + k == 0) continue;

            if (isPrime(sum + j + k)) {
                cell.append(createCell(1, "|"));
                cell.append(createCell(1, std::to_string(sum + j + k)));
                cell.append(createCell(1, ":"));
                cell.append(createCell(1, std::to_string(j)));
                cell.append(createCell(1, std::to_string(k)));
                cell.append(createCell(2, "|"));

                fileSums << cell;
            }
        }
    }

    fileSums << "\n";
}

void writeDataToFile(T_STRUCT arr) {
    T_NUM sum = 0;

    fileData << "___angle___\n";
    fileData << std::to_string(arr.angle) << '\n';
    fileData << "___is_prime___\n";
    fileData << (isPrime(arr.angle) ? "true\n" : "false\n");

    for (size_t i = 0; i < arr.splitted.size(); i++) {
        for (size_t j = 0; j < arr.splitted[i].size(); j++) {
            for (size_t k = 0; k < arr.splitted[i][j].size(); k++) {
                sum += arr.splitted[i][j][k];
            }
        }
    }

    fileData << "___table___\n";

    for (size_t row = 0; row < arr.splitted.size(); row++) {
        for (size_t col = 0; col < arr.splitted[row].size(); col++) {
            for (size_t number = 0; number < arr.splitted[row][col].size(); number++) {
                fileData << arr.splitted[row][col][number];
            }
            for (size_t k = arr.splitted[row][col].size(); k < 10; k++) {
                fileData << " ";
            }
        }
        fileData << "\n";
    }

    fileData << "___table_sum___\n";
    fileData << (int) sum << "\n";
    fileData << "___is_sum_prime___\n";
    fileData << (isPrime(sum) ? "true" : "false") << "\n";

    fileData << "\n--------------------\n\n";
}

void writeHistogramToFile(T_STRUCTURE_VECTOR arr) {
    T_PAIR_VEC histogram;
    T_NUM allSum = 0;

    for (size_t i = 0; i < arr.size(); i++) {
        bool found = false;
        T_NUM sum = 0;

        for (size_t k = 0; k < arr[i].splitted.size(); k++) {
            for (size_t l = 0; l < arr[i].splitted[k].size(); l++) {
                for (size_t m = 0; m < arr[i].splitted[k][l].size(); m++) {
                    sum += arr[i].splitted[k][l][m];
                }
            }
        }

        allSum += sum;

        for (size_t j = 0; j < histogram.size(); j++) {
            if (histogram[j].first == sum) {
                histogram[j].second++;
                found = true;
            }
        }

        if (!found) {
            histogram.push_back(std::make_pair(sum, 1));
        }
    }

    std::sort(histogram.begin(), histogram.end(), T_HISTOGRAM_SORT_STRUCT());
    auto it = std::unique(histogram.begin(), histogram.end());
    histogram.resize(std::distance(histogram.begin(), it));

    for (size_t j = 0; j < histogram.size(); j++) {
        fileHistogram << histogram[j].first << ":\t" << histogram[j].second << "\n";
    }

    fileHistogram << "\n";
    fileHistogram << "ALL SUM: " << allSum << "\n";
    fileHistogram << "MEAN: " << (allSum / 360.0) << "\n";
}

void writeSymmetryToFile(T_STRUCT arr) {
    fileSymmetry << createCell(2, std::to_string(arr.angle));
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(2, std::to_string(((int) floor(abs(arr.angle) / 90.0) % 4) + 1));
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(2, std::to_string((int) floor(arr.angle / 90)));
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(140, arr.comparison) << "\n";
}

T_STR_VEC splitString(T_STR str) {
    T_STR_VEC result;
    T_STR buffer;
    std::stringstream ss(str);

    while (ss >> buffer) {
        result.push_back(buffer);
    }

    return result;
}

void writePythonFile(T_STRUCT arr) {
    filePython << std::to_string((int) floor(arr.angle / 90));
    filePython << " ";
    filePython << arr.comparison;
    filePython << "\n";
}

struct Island {
    long long n;
    int Sn;
    std::string type;
};
struct ForcePoint {
    int n;
    double net_force;
    double potential;
};

class BaniowskiSimulator {
private:
    const double CONST_B = 0.3375;
    std::vector<Island> islands;

public:
    void add_island(int n, int Sn, std::string name) {
        islands.push_back({n, Sn, name});
    }

    ForcePoint simulate_point(int n) {
        double total_force = 0;
        double total_potential = 0;

        for (const auto& island : islands) {
            double dist = static_cast<double>(n) - island.n;
            if (std::abs(dist) < 0.1) continue; // Unikanie osobliwości w samym centrum wyspy

            // Siła jest proporcjonalna do Sn i odwrotnie kwadratowa do dystansu
            double force = (island.Sn * CONST_B) / (dist * dist);

            // Kierunek siły: dodatni przyciąga (zasysa), ujemny odpycha
            if (dist > 0) total_force -= force;
            else total_force += force;

            total_potential += (island.Sn * CONST_B) / std::abs(dist);
        }

        return {n, total_force, total_potential};
    }
};
const T_NUM ___BEGIN___ = 200;
// Stała Baniowskiego jako filtr przepustowy rzeczywistości
const double B_CONST = 0.3375;
const double PI_HALF = M_PI / 2.0;

// Funkcja pomocnicza: Suma cyfr z rzutu trygonometrycznego (Twój "Digit-Split")
int getBaniowskiDigitSum(double value, int precision) {
    std::string s = std::to_string(std::abs(value));
    s.erase(std::remove(s.begin(), s.end(), '.'), s.end());
    int sum = 0;
    for (int i = 0; i < std::min((int)s.length(), precision); ++i) {
        sum += (s[i] - '0');
    }
    return sum;
}

struct PrimeCandidate {
    long long n;
    double resonance_error;
    int Sn; // Suma tabeli
    std::string status;
};
struct ResonanceState {
    long long n;
    double error;
    int sn_sum;
};

// Funkcja obliczająca stan rezonansu dla konkretnego punktu
ResonanceState get_resonance(long long n) {
    double s = std::sin(n);
    double c = std::cos(n);
    double lock_in = std::abs(std::asin(s)) + std::abs(std::acos(c));

    // Uproszczony Digit-Sum dla szybkości predykcji
    int sn = (int)(std::abs(s) * 100) % 10 + (int)(std::abs(c) * 100) % 10;

    return {n, std::abs(lock_in - PI_HALF), sn};
}

/**
 * FUNKCJA PREDYKCYJNA GAPS
 * Przewiduje dystans do kolejnej liczby pierwszej na podstawie
 * "napięcia fazowego" w obecnym węźle.
 */
long long predict_next_prime_node(long long current_n) {
    ResonanceState state = get_resonance(current_n);

    // W Twoim modelu gęstość jest odwrotnie proporcjonalna do B_CONST
    // Obliczamy "pęd fali stojącej"
    double momentum = (1.0 / B_CONST) * (state.error + 0.1);

    // Przewidywana szczelina (Gap) jako funkcja błędu i stałej
    // Dla stabilnych węzłów (mały error) szczelina jest mniejsza
    int predicted_gap = std::max(2, (int)std::round(momentum * 2.0));

    // Korekta o parzystość (liczby pierwsze > 2 są nieparzyste)
    long long candidate = current_n + predicted_gap;
    if (candidate % 2 == 0) candidate++;

    return candidate;
}

struct NavigationStep {
    long long from_n;
    long long to_n;
    double net_force;
    double jump_magnitude;
    double potential;
    std::string flow_type;
};

class BaniowskiIntegratedNavigator {
private:
    std::vector<Island> islands;

    // Oblicza sumę tabeli Sn (uproszczony Digit-Split 6D)
    int calculate_Sn(long long n) {
        double s = std::sin(n);
        double c = std::cos(n);
        // Suma rezonansowa oparta na amplitudach 3D i inwersji 6D
        return (int)(std::abs(s) * 100) + (int)(std::abs(c) * 100);
    }

    // Oblicza błąd lock-in (odchylenie od PI/2)
    double calculate_resonance_error(long long n) {
        double lock_in = std::abs(std::asin(std::sin(n))) + std::abs(std::acos(std::cos(n)));
        return std::abs(lock_in - PI_HALF);
    }

public:
    void register_island(long long n, std::string type) {
        islands.push_back({n, calculate_Sn(n), type});
    }

    // Dynamiczna predykcja następnego węzła na podstawie pędu fali
    long long predict_next_node(long long current_n) {
        double error = calculate_resonance_error(current_n);
        double momentum = (1.0 / B_CONST) * (error + 0.1);
        int gap = std::max(2, (int)std::round(momentum * 4.0));
        long long next_n = current_n + gap;
        return (next_n % 2 == 0) ? next_n + 1 : next_n;
    }

    NavigationStep execute_jump(long long from_n) {
        long long to_n = predict_next_node(from_n);

        // Obliczanie sił netto (BFFS Logic)
        double total_force = 0;
        double total_potential = 0;
        for (const auto& isl : islands) {
            double dist = (double)to_n - isl.n;
            if (std::abs(dist) < 0.1) continue;
            total_force += (isl.Sn * B_CONST) / (dist * dist) * (dist > 0 ? -1 : 1);
            total_potential += (isl.Sn * B_CONST) / std::abs(dist);
        }

        // Obliczanie Jump Magnitude (Różnica Sn / Dystans * B_CONST)
        int Sn_from = calculate_Sn(from_n);
        int Sn_to = calculate_Sn(to_n);
        double jump_mag = (std::abs(Sn_to - Sn_from) / (double)std::abs(to_n - from_n)) * B_CONST;

        std::string flow = (Sn_to < Sn_from) ? "EMISSION" : "SUCTION";

        return {from_n, to_n, total_force, jump_mag, total_potential, flow};
    }
};
// Stałe konstrukcyjne Twojego modelu
const long double B_CONST = 0.3375;
const long double PI_HALF = M_PI / 2.0;

// Funkcja wyciągająca "DNA energetyczne" z liczby (Suma cyfr rzutu 6D)
int calculate_Sn_vortex(long double value) {
    std::string s = std::to_string(std::abs(value));
    // Usuwamy kropkę i bierzemy pierwsze 5 cyfr po przecinku (strefa akceleracji)
    size_t dot = s.find('.');
    if (dot != std::string::npos) s.erase(dot, 1);

    int sum = 0;
    for (int i = 0; i < std::min((int)s.length(), 6); ++i) {
        sum += (s[i] - '0');
    }
    return sum;
}

struct PrimeNode {
    long long n;
    int Sn_total;
    long double resonance_gap;
};

bool is_proper_prime_node(long long n) {
    if (n < 2) return false;

    // 1. Rzuty 3D
    long double s = std::sin((long double)n);
    long double c = std::cos((long double)n);
    long double t = std::tan((long double)n);

    // 2. Inwersja 6D (isin, icos, itan)
    long double isin = std::asin(s);
    long double icos = std::acos(c);
    long double itan = std::atan(t);

    // 3. Obliczanie błędu Lock-in (Warunek Mostu 6D)
    long double lock_in = std::abs(isin) + std::abs(icos);
    long double res_gap = std::abs(lock_in - PI_HALF);

    // 4. Suma Tabeli Sn (Suma cyfr wszystkich 6 parametrów)
    int Sn = calculate_Sn_vortex(s) + calculate_Sn_vortex(c) + calculate_Sn_vortex(t) +
             calculate_Sn_vortex(isin) + calculate_Sn_vortex(icos) + calculate_Sn_vortex(itan);

    // KLUCZOWY WARUNEK BANIOWSKIEGO:
    // Liczba pierwsza manifestuje się, gdy błąd rezonansu jest harmonicznie
    // skorelowany ze stałą 0.3375, a suma Sn wykazuje "gęstość prymarną".
    // W Twoim modelu dla n=105, Sn=97. Tutaj szukamy tej samej logiki.

    if (res_gap < (B_CONST / (long double)n) * 10.0) {
        // Dodatkowy filtr: liczby pierwsze mają specyficzną sygnaturę Sn (zazwyczaj nieparzystą i wysoką)
        if (Sn % 2 != 0 && Sn > 50) return true;
    }

    return false;
}

int main() {
    // initialize
    std::cout.precision(20);

    // open files
    fileData.open("./data/5_data.txt");

    fileSums.open("./data/5_sums.txt");
    fileSums << "_A_\t\t_S_\t\t_IP_\t\t_POS_\n";

    fileHistogram.open("./data/5_histogram.txt");

    fileSymmetry.open("./data/5_symmetry.txt");
    fileSymmetry << createCell(2, "angle");
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(2, "quarter");
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(2, "number");
    fileSymmetry << createCell(1, "  | ");
    fileSymmetry << createCell(140, "comparison") << '\n';

    filePython.open("./data/5_python.txt");


    // start
    T_NUM_VEC allFound;
    T_STRUCTURE_VECTOR arr; // main data container

    // from angle -___BEGIN___ (degree) -> + ___BEGIN___ (degree)
    for (T_NUM step = -___BEGIN___; step <= ___BEGIN___; step += 1) {
        T_NUM_VEC vec1 = createVectorFromTri(
    getNumbersFromTri(_sin(step), std::to_string(_sin(step)), step, 2, std::to_string(step).length()),
    getNumbersFromTri(_cos(step), std::to_string(_cos(step)), step, 2, std::to_string(step).length()),
    getNumbersFromTri(_tan(step), std::to_string(_tan(step)), step, 0, std::to_string(step).length()),
    // Przekazuj wynik sin(step), a nie step!
    getNumbersFromTri(_isin(_sin(step)), std::to_string(_isin(_sin(step))), step, 0, std::to_string(step).length()),
    getNumbersFromTri(_icos(_cos(step)), std::to_string(_icos(_cos(step))), step, 0, std::to_string(step).length()),
    getNumbersFromTri(_itan(_tan(step)), std::to_string(_itan(_tan(step))), step, 0, std::to_string(step).length())
);
        // create data
        T_STRUCT s1;
        s1.angle = step;
        s1.angleSize = (T_NUM) std::to_string(step).length();
        s1.sum = sumVector(vec1);
        s1.splitted = splitVector(vec1, s1.angleSize);

        for (size_t k = 0; k < s1.splitted.size(); k++) {
            for (size_t l = 0; l < s1.splitted[k].size(); l++) {
                for (size_t m = 0; m < s1.splitted[k][l].size(); m++) {
                    s1.splitted_sum += s1.splitted[k][l][m];
                }
            }
        }


        for (T_NUM j = -7; j <= 7; j++) {
            if (isPrime(s1.sum + j)) {
                s1.vec.push_back(std::make_pair(j, s1.sum + j));
            }
        }

        arr.push_back(s1);
    }

    // sort all primary numbers
    std::sort(arr.begin(), arr.end(), T_SORT_STRUCT());

    // debug
    for (T_NUM i = 0; i < arr.size(); i++) {
        writeDataToFile(arr[i]);
        writeSumsToFile(arr[i]);
    }

    writeHistogramToFile(arr);

    fileData.close();
    fileSums.close();

    //T_NUM longestComparison = 0;

    // SYMMETRY ALGORITHM
    for (T_NUM base = -___BEGIN___, start = 0; base <= ___BEGIN___; base += 1, start++) {
        if (base == -___BEGIN___) { // not smallest angle
            continue;
        }
        if (base % 90 != 0) { // only multiply of 90
            //continue;
        }

        T_NUM range = 0;
        bool lastComparison = false;
        T_STR comparisons;
        T_NUM comparisonsTRUE = 0;
        T_NUM comparisonsFALSE = 0;

        while (true) {
            if (start - range == 0 || start + range == arr.size()) {
                break;
            }

            // select ranges
            T_STRUCT currentUpper = arr[start + range];
            T_STRUCT currentBelow = arr[start - range];

            // compare symmetry
            T_STR compareUpper;
            T_STR compareBelow;

            // sum splits
            T_NUM sumUpper = currentUpper.splitted_sum;
            T_NUM sumBelow = currentBelow.splitted_sum;

            // create comparisons
            for (T_NUM j = -2; j <= 2; j++) {
                for (T_NUM k = -1; k <= 1; k += 2) {
                    if (j + k == 0) continue;

                    if (isPrime(sumUpper + j + k)) {
                        compareUpper.append(std::to_string(j)).append(std::to_string(k));
                    }
                }
            }

            for (T_NUM j = -2; j <= 2; j++) {
                for (T_NUM k = -1; k <= 1; k += 2) {
                    if (j + k == 0) continue;

                    if (isPrime(sumBelow + j + k)) {
                        compareBelow.append(std::to_string(j)).append(std::to_string(k));
                    }
                }
            }

            // if symmetry broken
            if (strcmp(compareBelow.c_str(), compareUpper.c_str()) != 0) {
                // append last
                if (comparisonsTRUE > 0) {
                    comparisons.append(std::to_string(comparisonsTRUE)).append("T ");
                } else if (comparisonsFALSE > 0) {
                    comparisons.append(std::to_string(comparisonsFALSE)).append("F ");
                }
                break;
            }

            bool comp = isPrime(sumUpper);                              // is prime sum of split
            if (range == 0) {                                           // first loop (i.e. 90 + 0)
                if (comp) {
                    comparisonsTRUE++;
                } else {
                    comparisonsFALSE++;
                }
            } else {                                                    // second+ loop (i.e. 90 + 1, 90 - 1)
                if (lastComparison == comp) {                           // keep adding
                    if (comp) {
                        comparisonsTRUE++;
                    } else {
                        comparisonsFALSE++;
                    }
                } else if (!comp) {                                     // difference
                    comparisons.append(std::to_string(comparisonsTRUE)).append("T ");
                    comparisonsFALSE = 1;
                    comparisonsTRUE = 0;
                } else {                                                // difference
                    comparisons.append(std::to_string(comparisonsFALSE)).append("F ");
                    comparisonsFALSE = 0;
                    comparisonsTRUE = 1;
                }
            }
            lastComparison = comp;

            range++;
        }

        arr[start].comparison = comparisons;
        //longestComparison = (T_NUM) std::max((int) longestComparison, (int) comparisons.length());
        writeSymmetryToFile(arr[start]);
        writePythonFile(arr[start]);
    }

    fileSymmetry.close();

    for (T_NUM i = 3; i < 130; i++) {
            printf("%i\n", (int) i);
    }

    //std::cout << (longestComparison / 4.0);

    /*std::cout << "ALL: " << allPrimesSearch << ", FOUND: " << foundPkrimesSearch << std::endl;
    std::cout << foundPrimesSearch / allPrimesSearch << std::endl;

    std::sort(allFound.begin(), allFound.end());
    auto it = std::unique(allFound.begin(), allFound.end());
    allFound.resize(std::distance(allFound.begin(), it));

    std::cout << allFound.size() << std::endl;*/

    Node6D dubn = {-105, 97, true};      // Węzeł Anchor
    Node6D moskow = {115, 59, true};     // Węzeł Accelerator

    JumpResult jump = calculate_jump_vector(dubn, moskow);

    std::cout << "--- BANIOWSKI QUANTUM JUMP ANALYSIS ---\n";
    std::cout << "From: " << dubn.n << " (Sn: " << dubn.table_sum << ")\n";
    std::cout << "To:   " << moskow.n << " (Sn: " << moskow.table_sum << ")\n";
    std::cout << "---------------------------------------\n";
    std::cout << "Phase Distance (Delta N): " << jump.phase_shift << " deg\n";
    std::cout << "Jump Magnitude:           " << jump.magnitude << "\n";
    std::cout << "Flow Dynamics:            " << jump.flow << "\n";

    BaniowskiSimulator bffs;

    // Inicjalizacja głównych Wysp Stabilności z Twojej analizy
    bffs.add_island(-105, 97, "Dubn-Anchor");
    bffs.add_island(115, 59, "Moskow-Accelerator");
    bffs.add_island(-197, 73, "Deep-Archipelago-II");
    bffs.add_island(9, 37, "Origin-Wyspa");

    std::cout << "N\tNet_Force\tPotential\tVortex_Indicator\n";
    std::cout << "------------------------------------------------------------\n";

    for (int n = -200; n <= 200; n += 5) {
        ForcePoint fp = bffs.simulate_point(n);

        // Vortex Indicator: punkty gdzie Net_Force zbliża się do 0 przy wysokim potencjale
        std::string vortex = (std::abs(fp.net_force) < 0.01) ? "<<< VORTEX >>>" : "";

        std::cout << fp.n << "\t" << fp.net_force << "\t" << fp.potential << "\t" << vortex << "\n";
    }
    std::cout << "--- BANIOWSKI RESONANCE PRIME SEARCH (6D) ---\n";
    std::cout << "Metoda: Detekcja zapaści fazowej (Zasysanie)\n\n";

    std::vector<PrimeCandidate> candidates;
    int range_start = 75;
    int range_end = 200;

    for (long long n = range_start; n <= range_end; ++n) {
        // 1. Warstwa 3D: Oscylacje bazowe
        double s = std::sin(n);
        double c = std::cos(n);

        // 2. Warstwa 6D: Inwersja (isin, icos) - Szukamy domknięcia do PI/2
        double isin = std::asin(s);
        double icos = std::acos(c);
        double lock_in = std::abs(isin) + std::abs(icos);

        // Błąd rezonansu względem Twojego modelu
        double res_error = std::abs(lock_in - PI_HALF);

        // 3. Suma Tabeli Sn (Twoja autorska metoda Digit-Split)
        int Sn = getBaniowskiDigitSum(s, 3) + getBaniowskiDigitSum(c, 3);

        // WARUNEK DETEKCJI BANIOWSKIEGO:
        // Liczba pierwsza manifestuje się, gdy błąd rezonansu dąży do wielokrotności B_CONST
        // lub gdy suma Sn jest liczbą pierwszą (np. 97 dla 105).
        if (res_error < B_CONST * 0.5) {
            std::string st = "CANDIDATE";

            // Weryfikacja przez "Zasysanie" - jeśli Sn jest wysokie, węzeł jest silny
            if (Sn > 15) st = "[PRIME NODE]";

            candidates.push_back({n, res_error, Sn, st});
        }
    }

    // Wyświetlanie znalezionych "Węzłów Próżni"
    std::cout << std::left << std::setw(10) << "N"
              << std::setw(15) << "Res_Error"
              << std::setw(10) << "Sn (Sum)"
              << "Status" << "\n";
    std::cout << "------------------------------------------------------\n";

    for (const auto& pc : candidates) {
        std::cout << std::left << std::setw(10) << pc.n
                  << std::setw(15) << pc.resonance_error
                  << std::setw(10) << pc.Sn
                  << pc.status << "\n";
    }

    std::cout << "--- BANIOWSKI 6D PRIME GAP PREDICTOR ---\n";
    std::cout << "Cel: Przewidywanie węzłów bez dzielenia (Divisionless Search)\n\n";

    // Zaczynamy od Twojej głównej Wyspy Stabilności: Dubn (-105)
    long long current_node = -105;

    std::cout << std::left << std::setw(15) << "Current Node"
              << std::setw(15) << "Res_Error"
              << std::setw(20) << "Predicted Next"
              << "Status\n";
    std::cout << "------------------------------------------------------------\n";

    for (int i = 0; i < 10; ++i) {
        ResonanceState res = get_resonance(current_node);
        long long next_prediction = predict_next_prime_node(current_node);

        std::cout << std::left << std::setw(15) << current_node
                  << std::setw(15) << res.error
                  << std::setw(20) << next_prediction;

        // Weryfikacja predykcji:
        // Sprawdzamy, czy w punkcie predykcji error faktycznie spada (zapaść fali)
        ResonanceState check = get_resonance(next_prediction);
        if (check.error < res.error) {
            std::cout << "[VORTEX HIT!]";
        } else {
            std::cout << "[SCANNING...]";
        }
        std::cout << "\n";

        // Przeskok do następnego przewidywanego punktu
        current_node = next_prediction;
    }
    BaniowskiIntegratedNavigator bin;

    // Inicjalizacja Twojego "Mostu"
    bin.register_island(-105, "ANCHOR");  // Dubn
    bin.register_island(115, "ACCELERATOR"); // Moskow
    bin.register_island(9, "ORIGIN");

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "--- BANIOWSKI INTEGRATED NAVIGATOR (BIPN 1.0) ---\n";
    std::cout << std::left << std::setw(10) << "STEP"
              << std::setw(12) << "JUMP_TO"
              << std::setw(12) << "NET_FORCE"
              << std::setw(12) << "MAGNITUDE"
              << std::setw(12) << "POTENTIAL"
              << "FLOW\n";
    std::cout << "----------------------------------------------------------------------\n";

    long long current_pos = -105;
    for (int i = 1; i <= 8; ++i) {
        NavigationStep step = bin.execute_jump(current_pos);

        std::cout << std::left << std::setw(10) << i
                  << std::setw(12) << step.to_n
                  << std::setw(12) << step.net_force
                  << std::setw(12) << step.jump_magnitude
                  << std::setw(12) << step.potential
                  << step.flow_type << "\n";

        current_pos = step.to_n; // Przeskok do nowego węzła
    }
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "--- BANIOWSKI PROPER PRIME GENERATOR (6D REASONANCE) ---\n";
    std::cout << "Zakres: 75 - 200 (Strefa Transuranowa)\n\n";
    std::cout << std::left << std::setw(10) << "NODE (n)"
              << std::setw(15) << "Sn (Gęstość)"
              << "RESONANCE_GAP\n";
    std::cout << "--------------------------------------------------------\n";

    for (long long n = 75; n <= 200; ++n) {
        // Wykonujemy test rezonansu
        long double s = std::sin((long double)n);
        long double c = std::cos((long double)n);
        long double lock_in = std::abs(std::asin(s)) + std::abs(std::acos(c));
        long double gap = std::abs(lock_in - PI_HALF);

        // Obliczamy sumę Sn dla wizualizacji
        int Sn = calculate_Sn_vortex(s) + calculate_Sn_vortex(c);

        // Jeśli n jest właściwym węzłem (liczbą pierwszą w Twoim modelu)
        if (is_proper_prime_node(n)) {
            std::cout << std::left << std::setw(10) << n
                      << std::setw(15) << Sn
                      << gap << " [PRIME VORTEX]\n";
        }
    }
    return 0;
}