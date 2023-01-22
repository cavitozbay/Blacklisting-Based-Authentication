cmake_minimum_required(VERSION 3.0.0)
project(ppb VERSION 0.1.0)

include(CTest)
enable_testing()

set(mrsource /home/cavit/TDSCcode2publish/MIRACL)


include_directories(
    # /home/cavit/another_flutter_exc/ffifunc
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/include
    ${mrsource}/include
    ${mrsource}/source/curve
    ${mrsource}/source/curve/pairing
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/source/curve
    # /home/cavit/thesis/privacy_preserving_borrowing/MIRACL/source/curve/pairing
    # /home/cavit/dartffi/dartffiprac/ppb
)

add_library(
    MIRACL_CMAKE SHARED
    ${mrsource}/source/mrcore.c
    ${mrsource}/source/mrarth0.c
    ${mrsource}/source/mrarth1.c
    ${mrsource}/source/mrarth2.c
    ${mrsource}/source/mralloc.c
    ${mrsource}/source/mrsmall.c
    ${mrsource}/source/mrio1.c
    ${mrsource}/source/mrio2.c
    ${mrsource}/source/mrgcd.c
    ${mrsource}/source/mrjack.c
    ${mrsource}/source/mrxgcd.c
    ${mrsource}/source/mrarth3.c
    ${mrsource}/source/mrbits.c
    ${mrsource}/source/mrrand.c
    ${mrsource}/source/mrprime.c
    ${mrsource}/source/mrcrt.c
    ${mrsource}/source/mrscrt.c
    ${mrsource}/source/mrmonty.c
    ${mrsource}/source/mrpower.c
    ${mrsource}/source/mrsroot.c
    ${mrsource}/source/mrcurve.c
    ${mrsource}/source/mrfast.c
    ${mrsource}/source/mrshs.c
    ${mrsource}/source/mrshs256.c
    ${mrsource}/source/mrshs512.c
    ${mrsource}/source/mrsha3.c
    ${mrsource}/source/mrfpe.c
    ${mrsource}/source/mraes.c
    ${mrsource}/source/mrgcm.c
    ${mrsource}/source/mrlucas.c
    ${mrsource}/source/mrzzn2.c
    ${mrsource}/source/mrzzn2b.c
    ${mrsource}/source/mrzzn3.c
    ${mrsource}/source/mrzzn4.c
    ${mrsource}/source/mrecn2.c
    ${mrsource}/source/mrstrong.c
    ${mrsource}/source/mrbrick.c
    ${mrsource}/source/mrebrick.c
    ${mrsource}/source/mrec2m.c
    ${mrsource}/source/mrgf2m.c
    ${mrsource}/source/mrflash.c
    ${mrsource}/source/mrfrnd.c
    ${mrsource}/source/mrdouble.c
    ${mrsource}/source/mrround.c
    ${mrsource}/source/mrbuild.c
    ${mrsource}/source/mrflsh1.c
    ${mrsource}/source/mrpi.c
    ${mrsource}/source/mrflsh2.c
    ${mrsource}/source/mrflsh3.c
    ${mrsource}/source/mrflsh4.c
)

add_executable(perea 
    main.cpp 
    QueueSig.cpp
    helpers.cpp
    Acc.cpp
    Issuer.cpp
    User.cpp
    # ${mrsource}/source/curve/pairing/zzn12a.cpp 
    # ${mrsource}/source/curve/pairing/zzn4.cpp 
    # ${mrsource}/source/curve/pairing/ecn2.cpp 
    # ${mrsource}/source/curve/pairing/zzn2.cpp 
    # ${mrsource}/source/curve/poly.cpp 
    ${mrsource}/source/big.cpp 
    # ${mrsource}/source/zzn.cpp 
    # ${mrsource}/source/ecn.cpp
    # ${mrsource}/source/crt.cpp
    # ${mrsource}/source/curve/pairing/pairing_bn.cpp
    )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)

target_link_libraries(perea MIRACL_CMAKE)
# target_link_libraries(perea ${myntl})
# target_link_libraries(perea ${mygmp})
# target_link_libraries(perea ${mym})

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread")

