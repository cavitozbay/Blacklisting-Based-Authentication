#ifndef COM_H
#define COM_H

#include <vector>
#include "pairing_bn.h"

using namespace std;

struct CmOp
{
    G2 C;
    vector<Big> m;
};

class Com
{
private:
    PFC& bn;
    
public:
    vector<G2> gens;
    size_t k;
    Com(PFC& _bn);
    Com(PFC& _bn, size_t _k);
    ~Com();
    CmOp commit(vector<Big>& m);
    void init_gens(vector<G2>& _gens);
    void pok_init(vector<Big>& alpha_ts, G2& u_t);
    void pok_response(vector<Big>& alphas, vector<Big>& alpha_ts, Big& c, vector<Big>& z);
    BOOL pok_verify(G2& u, G2& u_t, vector<Big>& z, Big& c);
};


#endif