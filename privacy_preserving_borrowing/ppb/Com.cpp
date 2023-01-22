#include "Com.hpp"

Com::Com(PFC& _bn, size_t _k): bn(_bn)
{
    k = _k;
    G2 tmp;
    for (size_t i = 0; i < k; i++)
    {
        bn.random(tmp);
        gens.push_back(tmp);
    }
}

Com::Com(PFC& _bn): bn(_bn)
{
}

Com::~Com()
{
}

CmOp Com::commit(vector<Big>& m)
{
    CmOp res;
    res.m = m;

    res.C = bn.mult(gens[0], m[0]);

    for (size_t i = 1; i < k; i++)
    {
        res.C = res.C + bn.mult(gens[i], m[i]);
    }
    
    return res;
}

void Com::init_gens(vector<G2>& _gens)
{
    gens = _gens;
    k = gens.size();
    for (size_t i = 0; i < k; i++)
    {
        bn.precomp_for_mult(gens[i]);
    }
    
}

void Com::pok_init(vector<Big>& alpha_ts, G2& u_t)
{
    alpha_ts.push_back(rand(bn.order()));
    u_t = bn.mult(gens[0], alpha_ts[0]);

    for (size_t i = 1; i < k; i++)
    {
        alpha_ts.push_back(rand(bn.order()));
        u_t = u_t + bn.mult(gens[i], alpha_ts[i]);
    }
}

void Com::pok_response(vector<Big>& alphas, vector<Big>& alpha_ts, Big& c, vector<Big>& z)
{

    for (size_t i = 0; i < k; i++)
    {
        z.push_back((alpha_ts[i] - modmult(alphas[i], c, bn.order()))%bn.order());
    }
}

BOOL Com::pok_verify(G2& u, G2& u_t, vector<Big>& z, Big& c)
{
    G2 l = bn.mult(gens[0], z[0]);

    for (size_t i = 1; i < k; i++)
    {
        l = l + bn.mult(gens[i], z[i]);
    }
    l = l + bn.mult(u, c);
    G2 r = u_t;

    return l == r;
}