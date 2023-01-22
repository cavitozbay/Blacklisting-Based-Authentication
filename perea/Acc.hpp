#ifndef ACC_H
#define ACC_H

#include <iostream>
#include <fstream>
#include <big.h>
#include <zzn.h>
#include "helpers.hpp"


class Acc
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


    Big N;
    Big g;
    Big V;
    Big sk_acc;
    Big up;
    Big g1, h1, h;
    int lx;
    Acc(/* args */);
    void add(Big t);
    Witness create_nonmem(Big t);
    BOOL verify_nonmem(Big t, Witness wit);
    NonMemProof prove_nonmem(Big t, Witness wit);
    BOOL verify_nonmem_proof(NonMemProof prf);
    RangeProof prove_range(Big E, Big x, Big r, Big a, Big b);
    BOOL verify_range(RangeProof prf, Big E, Big a, Big b);
};

#endif