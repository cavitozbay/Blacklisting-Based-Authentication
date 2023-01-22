#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <map>
#include "pairing_bn.h"

#define LOG

struct Sigma
{
    map<size_t, G2> tetas_hat;
    map<size_t, G1> tetas;
};

struct InitRegResp 
{
    G2 Cx0, CT0;
    G2 Cx0_t, CT0_t;

    Big x_z, rx0_z, T0_z, rT0_z;
};

struct RegResp
{
    Sigma sigma;
    vector<G2> M;
};

struct InitAuthResp 
{
    G2 Cx0p;
    G2 Cx0p_t;
    Big x_z, rx0pp_z;

    G2 CT0p;
    G2 CT0p_t;
    Big T0_z, rT0pp_z;

    G2 CR0p;
    G2 CR0p_t;
    Big R0_z, rR0pp_z;

    G2 CR1;
    G2 CR1_t;
    Big R0p_z, rR1_z;

    G2 CT1;
    G2 CT1_t;
    Big T1_z, rT1_z;

    Big rtetaP1_z, rtetaP2_z, rtetaP5_z, rx0p_z, rT0p_z, rR0p_z, rW_z;

    G2 CtetaP1, CtetaP2, CtetaP5, tetaP7, tetaP4, CP0p;
    G1 tetaP3, tetaP6;
    G1 CW;
    GT A_t, B_t, NM_t;

};

struct AuthResp
{
    Sigma sigma;
    size_t t;
    vector<G2> M;
};
#endif