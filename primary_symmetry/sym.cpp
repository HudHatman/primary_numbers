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
#include <bitset>
/*
 * DATA TYPES
 */
typedef long long T_NUM;
typedef std::vector<T_NUM> T_NUM_VEC;
typedef std::string T_STR;
typedef std::vector<T_STR> T_STR_VEC;
typedef std::pair<T_NUM, T_NUM> T_PAIR;
typedef std::vector<T_PAIR> T_PAIR_VEC;
typedef std::vector<std::vector<std::vector<T_NUM> > > T_3D_NUM_VEC;

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


// Stałe konstrukcyjne Twojego modelu
const long double B_CONST = 0.3375;
const long double PI_HALF = M_PI / 2.0;

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
    if (length <= 0) return {};
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

T_STR_VEC splitString(T_STR str) {
    T_STR_VEC result;
    T_STR buffer;
    std::stringstream ss(str);

    while (ss >> buffer) {
        result.push_back(buffer);
    }

    return result;
}

bool isPrime(int n) {
    if (n <= 1) return false;
    for (int i = 2; i <= sqrt(n); i++) {
        if (n % i == 0) return false;
    }
    return true;
}

class Element {
public:
    long id;
    T_3D_NUM_VEC table;
    int primeSubstract = 0;
    std::vector<std::vector<int> > w;
    int neutrons = 0;

    Element(int id) : id(id) {
        prepare();
    }

    int getProtons() { return abs((int) id); }
    int getElectrons() { return abs((int) id); }

    // Kluczowa poprawka: Sumowanie jądra (rzędy 0, 1, 2)
    int getNeutrons() {
        return this->neutrons;
    }

    // Masa atomowa (Liczba Masowa A) zgodnie z Twoim modelem
    int getAtomicMassNumber() {
        return getProtons() + getNeutrons();
    }

    // Opcjonalna masa fizyczna (w u)
    double getScientificMass() {
        return (getProtons() * 1.007276) + (getNeutrons() * 1.008665) + (getElectrons() * 0.000548);
    }

    void prepare() {
        for (int i = 0; i < 6; i++) {
            std::vector<int> row;
            this->w.push_back(row);
        }
        for (int i = 0; i < 6; i++) {
            double tan = _tan(this->id);
            T_STR tanStr = std::to_string(tan);

            double cos = _cos(this->id);
            T_STR cosStr = std::to_string(cos);

            double sin = _tan(this->id);
            T_STR sinStr = std::to_string(sin);

            for (int j = 0; j < 3; j++) {
                if (i == 0) {
                    int t = charToInt(tanStr.at(2 + j));
                    this->w[i].push_back(t);
                } else if (i == 1) {
                    int t = charToInt(cosStr.at(2 + j));
                    this->w[i].push_back(t);
                } else if (i == 2) {
                    int t = charToInt(sinStr.at(2 + j));
                    this->w[i].push_back(t);
                }
            }
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << this->w[i][j] << std::endl;
            }
        }
        std::cout << std::endl;
    }
};

int getPopcount(int n) {
    return std::bitset<32>(n).count();
}

std::vector<int> generateResultVector(const std::vector<Element> &inputs) {
    std::vector<int> result;

    for (int i = 0; i < inputs.size(); ++i) {
        const auto &v = inputs[i];

        // 1. Trimming i dekodowanie A: (V[0]*10 + V[1])
        int A = v.w[0][0] * 10 + v.w[0][1];

        // 2. Wyznaczenie przesunięcia (Offset)
        int offset = 0;
        if (i < 8) {
            offset = getPopcount(i % 4);
        } else {
            offset = 2; // Stabilizacja dla końcówki serii
        }

        result.push_back(A + offset);
    }

    return result;
}

int main() {
    // initialize
    std::cout.precision(20);

    // start
    T_NUM_VEC allFound;
    T_STRUCTURE_VECTOR arr; // main data container
    std::vector<Element> els;

    // from angle -___BEGIN___ (degree) -> + ___BEGIN___ (degree)
    for (T_NUM step = -200; step <= 200; step += 1) {
        T_NUM_VEC vec1 = createVectorFromTri(
            getNumbersFromTri(_tan(step), std::to_string(_tan(step)), step, 3, std::to_string(step).length()),
            getNumbersFromTri(_cos(step), std::to_string(_cos(step)), step, 3, std::to_string(step).length()),
            getNumbersFromTri(_sin(step), std::to_string(_sin(step)), step, 3, std::to_string(step).length()),
            getNumbersFromTri(_isin(_sin(step)), std::to_string(_isin(_sin(step))), step, 0,
                              std::to_string(step).length()),
            getNumbersFromTri(_icos(_cos(step)), std::to_string(_icos(_cos(step))), step, 0,
                              std::to_string(step).length()),
            getNumbersFromTri(_itan(_tan(step)), std::to_string(_itan(_tan(step))), step, 0,
                              std::to_string(step).length())
        );

        if (step > 0 && step < 10) {
            Element el((int) step);
            el.prepare();
            els.push_back(el);
        }
    }
    std::vector<int> n = generateResultVector(els);
    std::cout << "n: " << std::endl;
    for (int a = 0; a < n.size(); a += 1) {

    }
    /*// create data
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

    arr.push_back(s1);*/


    /*debug3dVector(arr.at(0).splitted);

    // SYMMETRY ALGORITHM
    for (T_NUM base = -200, start = 0; base <= 200; base += 1, start++) {
        if (base == -200) { // not smallest angle
            continue;
        }

        auto curr = arr.at(start);

        if (start > 200 && start < 210) {
            for (int i = 0; i < curr.splitted.size(); i++) {

            }
        }
    }*/

    return 0;
}
