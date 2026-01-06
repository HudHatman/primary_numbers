cmake_minimum_required(VERSION 3.13)
project(primary_numbers)

set(CMAKE_CXX_STANDARD 14)

include_directories(__SIEVE)
include_directories(not_primary)
include_directories(not_primary_old)
include_directories(not_primary_working)
include_directories(sequence_vectorized)
include_directories(sieve)

add_executable(primary_numbers
        __SIEVE/bigInt.cpp
        __SIEVE/bigInt.h
        __SIEVE/main.cpp
        not_primary/bigInt.cpp
        not_primary/bigInt.h
        not_primary/main.cpp
        not_primary_old/bigInt.cpp
        not_primary_old/bigInt.h
        not_primary_old/main.cpp
        not_primary_working/BigInteger.cpp
        not_primary_working/BigInteger.h
        not_primary_working/BigIntegerAlgorithms.cpp
        not_primary_working/BigIntegerAlgorithms.h
        not_primary_working/BigIntegerLibrary.h
        not_primary_working/BigIntegerUtils.cpp
        not_primary_working/BigIntegerUtils.h
        not_primary_working/BigUnsigned.cpp
        not_primary_working/BigUnsigned.h
        not_primary_working/BigUnsignedInABase.cpp
        not_primary_working/BigUnsignedInABase.h
        not_primary_working/main.cpp
        not_primary_working/NumberlikeArray.h
        primary_1/main.cpp
        primary_symmetry/main.cpp
        sequence/main.cpp
        sequence_vectorized/bigInt.h
        sequence_vectorized/main.cpp
        sieve/bigInt.cpp
        sieve/bigInt.h
        sieve/main.cpp
        test/main.cpp
        untitled2/main.cpp
        Riemann_Hypothesis/main.cpp
        ___/calculate_element.cpp)
