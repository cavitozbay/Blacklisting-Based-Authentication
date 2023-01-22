#ifndef ACC_H
#define ACC_H

#include <vector>
#include <NTL/ZZ_pX.h>
#include "pairing_bn.h"
#include "poly.h"
#include "helpers.hpp"

struct Witness 
{
    G1 w1_hat;
    G2 w2;
};


class Acc
{
private:
    PFC& bn;
    
    Big alpha;
    size_t a_size;
    
public:
    GT gT;
    G2& g;
    G1& g_hat;
    G1 V_hat;
    vector<Big> elements;
    vector<ZZn> elementsZZn;
    vector<G2> rs;
    vector<G1> rs_hat;
    Poly f;
    NTL::ZZ_pX fNTL;
    Acc(PFC& _bn, G2& _g, G1& _g_hat);
    ~Acc();
    void setup(size_t _a_size);
    void add(Big& y);
    Witness create_nonmembership(vector<ZZn>& nonmembers);
    BOOL verify_nonmembership(G2& Y, Witness& wit);
    void update_nonmembership(Big& y, Witness& wit);
    void create_NM_proof_param(vector<ZZn>& nonmembers, G2& g_h, G1& g_hat_h0);
};


#endif