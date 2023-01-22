cmake_minimum_required(VERSION 3.0.0)
project(ppb VERSION 0.1.0)

include(CTest)
enable_testing()

set(mrsource .)


include_directories(
    # /home/cavit/another_flutter_exc/ffifunc
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/include
    ${mrsource}/MIRACL/include
    ${mrsource}/MIRACL/source/curve
    ${mrsource}/MIRACL/source/curve/pairing
    ./ntl-11.5.1/include
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/source/curve
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/source/curve/pairing
    # /home/cavit/dartffi/dartffiprac/ppb
)

add_library(
    MIRACL_CMAKE SHARED
    ${mrsource}/MIRACL/source/mrcore.c
    ${mrsource}/MIRACL/source/mrarth0.c
    ${mrsource}/MIRACL/source/mrarth1.c
    ${mrsource}/MIRACL/source/mrarth2.c
    ${mrsource}/MIRACL/source/mralloc.c
    ${mrsource}/MIRACL/source/mrsmall.c
    ${mrsource}/MIRACL/source/mrio1.c
    ${mrsource}/MIRACL/source/mrio2.c
    ${mrsource}/MIRACL/source/mrgcd.c
    ${mrsource}/MIRACL/source/mrjack.c
    ${mrsource}/MIRACL/source/mrxgcd.c
    ${mrsource}/MIRACL/source/mrarth3.c
    ${mrsource}/MIRACL/source/mrbits.c
    ${mrsource}/MIRACL/source/mrrand.c
    ${mrsource}/MIRACL/source/mrprime.c
    ${mrsource}/MIRACL/source/mrcrt.c
    ${mrsource}/MIRACL/source/mrscrt.c
    ${mrsource}/MIRACL/source/mrmonty.c
    ${mrsource}/MIRACL/source/mrpower.c
    ${mrsource}/MIRACL/source/mrsroot.c
    ${mrsource}/MIRACL/source/mrcurve.c
    ${mrsource}/MIRACL/source/mrfast.c
    ${mrsource}/MIRACL/source/mrshs.c
    ${mrsource}/MIRACL/source/mrshs256.c
    ${mrsource}/MIRACL/source/mrshs512.c
    ${mrsource}/MIRACL/source/mrsha3.c
    ${mrsource}/MIRACL/source/mrfpe.c
    ${mrsource}/MIRACL/source/mraes.c
    ${mrsource}/MIRACL/source/mrgcm.c
    ${mrsource}/MIRACL/source/mrlucas.c
    ${mrsource}/MIRACL/source/mrzzn2.c
    ${mrsource}/MIRACL/source/mrzzn2b.c
    ${mrsource}/MIRACL/source/mrzzn3.c
    ${mrsource}/MIRACL/source/mrzzn4.c
    ${mrsource}/MIRACL/source/mrecn2.c
    ${mrsource}/MIRACL/source/mrstrong.c
    ${mrsource}/MIRACL/source/mrbrick.c
    ${mrsource}/MIRACL/source/mrebrick.c
    ${mrsource}/MIRACL/source/mrec2m.c
    ${mrsource}/MIRACL/source/mrgf2m.c
    ${mrsource}/MIRACL/source/mrflash.c
    ${mrsource}/MIRACL/source/mrfrnd.c
    ${mrsource}/MIRACL/source/mrdouble.c
    ${mrsource}/MIRACL/source/mrround.c
    ${mrsource}/MIRACL/source/mrbuild.c
    ${mrsource}/MIRACL/source/mrflsh1.c
    ${mrsource}/MIRACL/source/mrpi.c
    ${mrsource}/MIRACL/source/mrflsh2.c
    ${mrsource}/MIRACL/source/mrflsh3.c
    ${mrsource}/MIRACL/source/mrflsh4.c
)

target_include_directories(
    MIRACL_CMAKE
    PUBLIC /home/cavit/another_flutter_exc/MIRACL_SHAREDcomba/include
    PUBLIC /home/cavit/another_flutter_exc/MIRACL_SHAREDcomba/source/curve
    PUBLIC /home/cavit/another_flutter_exc/MIRACL_SHAREDcomba/source/curve/pairing
)

# find_library(
#     SHARED
#     mymiracl MIRACL_SHARED
#     NAMES MIRACL_SHARED
#     HINTS 
#         ${mrsource}/MIRACL/build
#     NO_CACHE
#     NO_DEFAULT_PATH
# )

find_library(
    myntl ntl
    NAMES ntl
    HINTS "./ntl-11.5.1/src"
    NO_CACHE
    NO_DEFAULT_PATH
)

find_library(
    mygmp gmp
    NAMES gmp
)

find_library(
    mym m
    NAMES m
)

add_executable(ppb 
    blacw/main.cpp 
    ${mrsource}/MIRACL/source/curve/pairing/zzn12a.cpp 
    ${mrsource}/MIRACL/source/curve/pairing/zzn4.cpp 
    ${mrsource}/MIRACL/source/curve/pairing/ecn2.cpp 
    ${mrsource}/MIRACL/source/curve/pairing/zzn2.cpp 
    ${mrsource}/MIRACL/source/curve/poly.cpp 
    ${mrsource}/MIRACL/source/curve/pairing/pairing_bn.cpp
    ${mrsource}/MIRACL/source/big.cpp 
    ${mrsource}/MIRACL/source/zzn.cpp 
    ${mrsource}/MIRACL/source/ecn.cpp
    ${mrsource}/MIRACL/source/crt.cpp
    blacw/helpers.cpp
    blacw/Acc.cpp
    blacw/SPSEQ.cpp
    blacw/Com.cpp
    blacw/TLP.cpp
    blacw/BLACWIssuer.cpp
    blacw/BLACWHolder.cpp
    )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)

target_link_libraries(ppb MIRACL_CMAKE)
target_link_libraries(ppb ${myntl})
target_link_libraries(ppb ${mygmp})
target_link_libraries(ppb ${mym})

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread")
