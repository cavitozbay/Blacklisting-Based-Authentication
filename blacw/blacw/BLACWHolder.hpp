#ifndef BLACWHOLDER_H
#define BLACWHOLDER_H

#include <iterator>
#include <algorithm>
#include "pairing_bn.h"
#include "big.h"
#include "BLACWIssuer.hpp"
#include "Com.hpp"
#include "Acc.hpp"
#include "SPSEQ.hpp"

using namespace std;

struct Credential
{
    G2 dsid, C[3];
    Big esk, dsrnd0, dsrnd1, z, t, y, v, u;
    Sigma sigma;
    vector<Big> BIDS_u;
    vector<ZZn> BIDS_uZZn;
    CmOp cmop;
};

class BLACWHolder
{
private:
    PFC& bn;
    Com& com;
    Acc& acc;
    SPSEQ& spseq;
    TLP& tlp;
    Credential cred;
    Big usk, uS;
    G2& g;
    G2& g_tilde;
    G2& h_tilde;
    G1& g_hat;
    G1& g_hat_p;
    G1& h_hat_p;
    Poly tmpf0;

    vector<G1> pcomp_h0s;
    Poly pcomp_h;
    G2 pcomp_w2;
    size_t pcomp_d, pcomp_rem;
    vector<Big> pcomp_bids;
    vector<ZZn> pcomp_bids_ZZn;

    void prepare_ZKP1(ZKP1& zkp, CmOp& cmop, G2& G);
    void prepare_ZKP2(ZKP2& zkp, Credential& credS, vector<G2> C_pre, Big& bidS, Witness& wit);
    void prepare_ZKP3(ZKP3& zkp, Big& yS, Big& s, Big& bid, G2& C_pr1, G2& C_1S);
    void new_bid(vector<Big>& vecbig, vector<ZZn>& veczzn);
    void remove_bid(Big& bid, vector<Big>& vecbig, vector<ZZn>& veczzn);
    void add_bid(Big& bid, vector<Big>& vecbig, vector<ZZn>& veczzn);


public:
    G2 upk;
    BLACWHolder(PFC& _bn, BLACWIssuer& issuer);
    ~BLACWHolder();
    vector<G2> receive_cred1(ZKP1& zkp);
    BOOL receive_cred2(vector<G2>& C_pre, Big& v, Sigma& _sigma, Big& esk_issuerS, G2& w);
    void borrow1(Big& vS, Big& gamma, Big& c0, Big& c1, Witness& wit, vector<G2>& C_pre, G2* _C, Sigma& _sigma, Big& bidS, ZKP2& zkp);
    BOOL borrow2(vector<G2>& C_pre, Sigma& _sigma_pr, Big& esk_issuerS, G2& w);
    void repay1(Big& bid, G2& C_1S, G2* C_pr, Sigma& sigma_pr, Big& s, ZKP3& zkp);
    BOOL repay2(Sigma& sigma_prD, G2& C_1S, Big& s, Big& v_p);

    void precompute_h0s(size_t d);
    void pcomp_borrow1(Big& vS, Big& gamma, Big& c0, Big& c1, Witness& wit, vector<G2>& C_pre, G2* _C, Sigma& _sigma, Big& bidS, ZKP2& zkp);
 
};


#endif