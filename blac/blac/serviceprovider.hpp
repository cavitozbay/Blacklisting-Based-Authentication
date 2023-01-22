#ifndef SERVICEPROVIDER_H
#define SERVICEPROVIDER_H

#include "pairing_bn.h"
#include "acc.hpp"
#include "aho.hpp"
#include "helpers.hpp"
#define ACC_SIZE 50

class ServiceProvider 
{
public:
    PFC& bn;
    Acc acc;
    AHO aho;
    G2 g_hat, h_hat;
    G1 h;
    ServiceProvider(PFC& _bn);
    size_t t = 1;
    RegResp _register(InitRegResp resp);
    AuthResp auth(InitAuthResp resp, vector<size_t> V);
};

#endif