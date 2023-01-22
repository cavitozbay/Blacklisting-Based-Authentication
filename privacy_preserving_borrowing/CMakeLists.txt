cmake_minimum_required(VERSION 3.0.0)
project(ppb VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(ppb 
    ppb/main.cpp 
    MIRACL/source/curve/pairing/zzn12a.cpp 
    MIRACL/source/curve/pairing/zzn4.cpp 
    MIRACL/source/curve/pairing/ecn2.cpp 
    MIRACL/source/curve/pairing/zzn2.cpp 
    MIRACL/source/curve/poly.cpp 
    MIRACL/source/big.cpp 
    MIRACL/source/zzn.cpp 
    MIRACL/source/ecn.cpp
    MIRACL/source/crt.cpp
    ppb/helpers.cpp
    ppb/Acc.cpp
    ppb/SPSEQ.cpp
    ppb/Com.cpp
    ppb/TLP.cpp
    ppb/PPBIssuer.cpp
    ppb/PPBHolder.cpp
    MIRACL/source/curve/pairing/pairing_bn.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)

include_directories("MIRACL/include" "MIRACL/source/curve" "MIRACL/source/curve/pairing")

find_library(
    mymiracl miracl
    NAMES miracl
    HINTS 
        "MIRACL"
    NO_DEFAULT_PATH
)

target_link_libraries(ppb ${mymiracl})