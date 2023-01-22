#ifndef ACC_H
#define ACC_H
#include "pairing_bn.h"
#include <map>
#include <vector>

class Acc
{
private:
    /* data */
public:
    size_t n;
    G1 g;
    G2 g_hat;
    PFC& bn;
    map<size_t, G1> rs;
    map<size_t, G2> rs_hat;
    Acc(PFC& _bn);
    void setup(size_t n);
    G1 gen(vector<size_t>& V);
    G1 wit_gen(vector<size_t>& U, vector<size_t>& V);
    BOOL verify(vector<size_t>& U, G1 acc_v, G1 wit);
    ~Acc() {};
};
#endif