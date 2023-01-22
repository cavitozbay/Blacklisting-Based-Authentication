#include <iostream>
#include <deque>
#include "big.h"

#include "helpers.hpp"
#include "Issuer.hpp"

class User
{
private:
    /* data */
public:
    int lN = 2048;
    int lm = 256;
    int lp = 128; 
    int lpp = 128;
    int le = 128;
    int ls = 128;
    Big two;
    Big mrange;
    Big mtrange;
    Big ntrange;
    Big crange;
    Big g,h,N;
    deque<Big> Q;
    Big r;
    Sig sig;
    deque<Witness> wits;
    Acc& acc;
    QueueSig& qsig;
    User(Issuer& issuer);
    ReqSig init_register();
    void finish_register(RegResp resp);
    InitAuth init_auth();
    void finish_auth(AuthResp resp);
    KoSProof prove_kos(Sig sig, deque<Big> Q);
    BOOL verify_kos(KoSProof prf);
    RangeProof prove_range(Big E, Big x, Big r, Big a, Big b);
    BOOL verify_range(RangeProof prf, Big E, Big a, Big b);

};



