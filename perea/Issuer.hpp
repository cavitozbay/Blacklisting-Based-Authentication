#ifndef ISSUER_H
#define ISSUER_H

#include <iostream>
#include "big.h"

#include "Acc.hpp"
#include "QueueSig.hpp"

class Issuer
{
private:
    /* data */
public:
    int lN = 2048;
    int lm = 256;
    int lp = 128; 
    int lpp = 128;
    int le = 128;
    int ls = 128;
    Big two;
    Big mrange;
    Big mtrange;
    Big ntrange;
    Big crange;
    Big g,h,N;

    Acc acc;
    QueueSig qsig;
    Issuer(size_t K);
    RegResp _register(ReqSig req);
    AuthResp auth(InitAuth req);
    BOOL verify_kos(KoSProof prf);
    BOOL verify_range(RangeProof prf, Big E, Big a, Big b);

};




#endif