#ifndef AHO_H
#define AHO_H

#include <vector>
#include <map>
#include "pairing_bn.h"
#include "helpers.hpp"

class AHO {
private:
    void rndmz(G1& x, G2& y);

public:
    size_t L;
    PFC& bn;
    AHO(PFC& _bn);
    GT A, B;
    G1 Gr, Hr, Gz, Hz;
    G2 g_hat;
    vector<G1> Gs;
    vector<G1> Hs;

    Big alpha_a, alpha_b, mu_z, nu_z;
    vector<Big> mus;
    vector<Big> nus;
    void key_gen(size_t _L);
    Sigma sign(vector<G2>& M);
    BOOL verify(Sigma sigma, vector<G2> M);
    Sigma randomize(Sigma sigma);
};

#endif