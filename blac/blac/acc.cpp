#include "acc.hpp"

Acc::Acc(PFC& _bn): bn(_bn)
{
    bn.random(g);
    bn.random(g_hat);
}

void Acc::setup(size_t _n)
{
    n = _n;
    Big gamma = rand(bn.order());
    G1 tmp = bn.mult(g, gamma);
    G2 tmp_hat = bn.mult(g_hat, gamma);
    for (size_t i = 1; i <= n; i++)
    {
        rs[i] = tmp;
        tmp = bn.mult(tmp, gamma);   
        rs_hat[i] = tmp_hat;
        tmp_hat = bn.mult(tmp_hat, gamma);
    }

    tmp = bn.mult(tmp, gamma);   
    tmp_hat = bn.mult(tmp_hat, gamma);

    for (size_t i = n+2; i <= 2*n; i++)
    {
        rs[i] = tmp;
        tmp = bn.mult(tmp, gamma);   
        rs_hat[i] = tmp_hat;
        tmp_hat = bn.mult(tmp_hat, gamma);
    }
    
}

G1 Acc::gen(vector<size_t>& V)
{
    G1 acc_v = G1();

    for(auto i: V)
    {
        acc_v = acc_v + rs[n+1-i];
    }
    return acc_v;
}

G1 Acc::wit_gen(vector<size_t>& U, vector<size_t>& V)
{
    G1 wit = G1();

    for(auto j: U)
    {
        for(auto i: V)
        {
            if(i != j) wit = wit + rs[n+1+j-i];
        }
    }
    return wit;
}

BOOL Acc::verify(vector<size_t>& U, G1 acc_v, G1 wit)
{
    G2 gs_hat = G2();

    for(auto i: U)
    {
        gs_hat = gs_hat + rs_hat[i];
    }

    return (bn.pairing(gs_hat, acc_v) == bn.pairing(g_hat, wit));
}