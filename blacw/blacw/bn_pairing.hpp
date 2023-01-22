#ifndef BN_PAIRING_H
#define BN_PAIRING_H

#include <iostream>
#include <fstream>
#include <string.h>
#include "ecn.h"
#include <ctime>
#include "ecn2.h"
#include "zzn12a.h"

namespace BN_Pairing
{

    #ifdef MR_AFFINE_ONLY
        #define AFFINE
    #else
        #define PROJECTIVE
    #endif

    // Using SHA-256 as basic hash algorithm

    #define HASH_LEN 32

    //
    // R-ate Pairing Code
    //

    class BNPairing
    {
        public:
        miracl* mip;
        ZZn2 X;
        ZZn Beta;
        Big p,q,x,cf,t,BB[4][4],WB[4],SB[2][2],W[2];
        int B;

        void set_frobenius_constant(ZZn2 &X);
        
        BNPairing(miracl* mip);
        ECn G1_mult(ECn &P,Big &e);
        ECn2 G2_mult(ECn2 &P,Big &e);
        ZZn12 GT_pow(ZZn12 &res,Big &e);
        BOOL fast_pairing(ECn2& P, ECn& Q,ZZn12& res);
        void cofactor(ECn2& S,ZZn2 &F);
        BOOL ecap(ECn2& P,ECn& Q,ZZn12& r);
        BOOL member(ZZn12 r);


    };

}

#endif