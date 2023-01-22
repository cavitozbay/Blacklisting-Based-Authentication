#include "SPSEQ.hpp"

SPSEQ::SPSEQ(size_t _l, PFC& _bn, G2& _g, G1& _g_hat) : bn(_bn)
{
    l = _l;
    g = _g;
    g_hat = _g_hat;
    pk_hat = new G1*[l];
    for (size_t i = 0; i < l; i++)
    {
        pk_hat[i] = new G1;
    }
    
}

SPSEQ::~SPSEQ()
{
    for (size_t i = 0; i < l; i++)
    {
        delete pk_hat[i];
    }
    
    delete[] pk_hat;
}

vector<Big> SPSEQ::keygen()
{
    vector<Big> sk;
    Big tmp;
    for (size_t i = 0; i < l; i++)
    {
        tmp = rand(bn.order());
        sk.push_back(tmp);
        *(pk_hat[i]) = bn.mult(g_hat, tmp);
    }
    return sk;
}

Sigma SPSEQ::sign(vector<Big>& sk, G2* M)
{
    Sigma sigma;
    Big y = rand(bn.order());
    sigma.Z = bn.mult(M[0], sk[0]);

    for (size_t i = 1; i < l; i++)
    {
        sigma.Z = sigma.Z + bn.mult(M[i], sk[i]);
    }

    sigma.Z = bn.mult(sigma.Z, y);
    
    Big y_inv = inverse(y,bn.order());

    sigma.Y = bn.mult(g, y_inv);
    sigma.Y_hat = bn.mult(g_hat, y_inv);

    return sigma;
}

BOOL SPSEQ::verify(G2* M, Sigma& sigma)
{
    GT v1, _v1, v2, _v2;
    G2* Mp[3];
    Mp[0] = &(M[0]);
    Mp[1] = &(M[1]);
    Mp[2] = &(M[2]);
    v1 = bn.multi_pairing(l, Mp, pk_hat);
    /*v1 = bn.pairing(pk_hat[0], M[0]);
    for (size_t i = 1; i < l; i++)
    {
        v1 = v1 * bn.pairing(pk_hat[i], M[i]);
    }*/
    
    
    _v1 = bn.pairing(sigma.Z, sigma.Y_hat);
    //_v1 = bn.final_exp(_v1);

    v2 = bn.pairing(g, sigma.Y_hat);
    _v2 = bn.pairing(sigma.Y, g_hat);
    if(v1 == _v1 && v2 == _v2)
        return TRUE;
    return FALSE;
}

Sigma SPSEQ::chgrep(Sigma& sigma, G2* M, Big& mu)
{
    Big psi = rand(bn.order());
    Big psi_inv = inverse(psi, bn.order());
    Sigma adapted;

    adapted.Z = bn.mult(sigma.Z, modmult(mu, psi, bn.order()));
    adapted.Y = bn.mult(sigma.Y, psi_inv);
    adapted.Y_hat = bn.mult(sigma.Y_hat, psi_inv);

    for (size_t i = 0; i < l; i++)
    {
        M[i] = bn.mult(M[i], mu);
    }
    
    return adapted;
}