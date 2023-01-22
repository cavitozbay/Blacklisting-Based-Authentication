#include "pairing_bn.h"
#include "helpers.hpp"
#include "acc.hpp"
#include "aho.hpp"
#include "serviceprovider.hpp"

class User
{
public:
    PFC& bn;
    Acc& acc;
    AHO& aho;
    Big x, T0, rx0, rT0, rR0, R0;
    Big R1, rR1, T1, rT1;

    G2 g_hat, h_hat;
    G1 h;
    G2 P0, Cx0, CT0, CP0, CR0;
    Sigma sigma;
    vector<size_t> LU;

    User(ServiceProvider& sp);
    InitRegResp init_register();
    void finish_register(RegResp resp);
    InitAuthResp init_auth(vector<size_t> V);
    void finish_auth(AuthResp resp);
};