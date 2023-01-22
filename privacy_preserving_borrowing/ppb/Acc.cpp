#include "Acc.hpp"

Acc::Acc(PFC& _bn, G2& _g, G1& _g_hat) : bn(_bn), g(_g), g_hat(_g_hat)
{
    //bn = _bn;
    modulo(*bn.mod);
    /*g = _g;
    g_hat = _g_hat;*/
    gT = bn.pairing(g, g_hat);
    bn.precomp_for_power(gT);
    //gT = bn.final_exp(gT);
}

Acc::~Acc()
{

}

void Acc::setup(size_t _a_size)
{
    modulo(*bn.mod);
    a_size = _a_size;
    
    G2 t = g;
    G1 t_hat = g_hat;
    V_hat = g_hat;

    alpha = rand(bn.order());

    for (size_t i = 0; i < a_size; i++)
    {
        rs.push_back(t);
        rs_hat.push_back(t_hat);
        bn.precomp_for_mult(rs[i]);
        bn.precomp_for_mult(rs_hat[i]);
        t = bn.mult(t, alpha);
        t_hat = bn.mult(t_hat, alpha);
    }

    modulo(bn.order());
    for (size_t i = 0; i < BLACKLISTED_COUNT; i++)
    {
        Big y = rand(bn.order());
        elements.push_back(y);
        ZZn yZ(y);
        elementsZZn.push_back(yZ);
    }

    f = build_from_roots(elementsZZn, bn.order());
    modulo(bn.order());
    V_hat = mult_poly_hat(bn, rs_hat, f);
    modulo(*bn.mod);
}

void Acc::add(Big& y)
{   
    modulo(bn.order());
    elements.push_back(y);
    ZZn yZ(y);
    //cout << "==Acc Add: " << yZ << endl;
    elementsZZn.push_back(yZ);
    f = build_from_roots(elementsZZn, bn.order());
    modulo(bn.order());
    //cout << "  f: " << f << endl;
    V_hat = mult_poly_hat(bn, rs_hat, f);
    //V_hat *= alpha - y;
    modulo(*bn.mod);

}

Witness Acc::create_nonmembership(vector<ZZn>& nonmembers)
{
    modulo(*bn.mod);

    if(nonmembers.empty())
        return Witness({g_hat, bn.mult(g, Big(0))});
    //cout << "==Acc Nonmem: " << endl;
    //Poly f = build_from_roots(elementsZZn, bn.order());
    Poly f0 = build_from_roots(nonmembers, bn.order());
    //cout << f0 << endl;
    modulo(bn.order());
    //cout << "  f: " << f << endl;
    //cout << "  f0: " << f0 << endl;
    
    Poly result[3];
    egcd(result, f, f0);
    
    
    ZZn divisor = result[0].coeff(0);
    result[0] /= divisor;
    result[1] /= divisor;  
    result[2] /= divisor;

    //result[2] *= ZZn(-1);
    //cout << "f: " << f << endl;
    /*cout << "f0: " << f0 << endl;
    cout << "gcd: " << result[0] << endl;
    cout << "h: " << result[1] << endl;
    cout << "h0: " << result[2] << endl;*/
    G1 w1_hat = mult_poly_hat(bn, rs_hat,result[2]);
    G2 w2 = mult_poly(bn, rs, result[1]);
    modulo(bn.order());
    //cout << "  GCD check: " << f*result[1] + f0*result[2] << endl;

    modulo(*bn.mod);
    return Witness({w1_hat,w2});
} 

BOOL Acc::verify_nonmembership(G2& Y, Witness& wit)
{
    
    modulo(*bn.mod);
    BOOL Ok;
    GT res1, res2;

    res1 = bn.pairing(Y, wit.w1_hat); 
    res2 = bn.pairing(wit.w2, V_hat);

    BOOL res = res2 * res1 == gT;

#ifdef LOG
    cout << "==Acc Verif: " << endl;
    cout << "  " << res << endl;
#endif

    modulo(*bn.mod);
    return res;
}

void Acc::update_nonmembership(Big& y, Witness& wit)
{
}
