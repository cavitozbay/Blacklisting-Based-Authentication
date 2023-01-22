#ifndef SPSEQ_H
#define SPSEQ_H

#include <vector>
#include "pairing_bn.h"

using namespace std;

struct Sigma 
{
    G2 Z;
    G2 Y;
    G1 Y_hat;
};

class SPSEQ
{
private:
    G2 g;
    G1 g_hat;
    size_t l;
    PFC& bn;
    //vector<Big> sk;

public:
    G1** pk_hat;
    SPSEQ(size_t _l, PFC& _bn, G2& _g, G1& _g_hat);
    vector<Big> keygen();
    Sigma sign(vector<Big>& sk, G2* M);
    BOOL verify(G2* M, Sigma& sigma);
    Sigma chgrep(Sigma& sigma, G2* M, Big& mu);
    ~SPSEQ();
};

#endif