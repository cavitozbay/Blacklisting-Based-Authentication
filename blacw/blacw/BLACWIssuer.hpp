#ifndef BLACW_H
#define BLACW_H

#include <vector>
#include "pairing_bn.h"
#include "SPSEQ.hpp"
#include "Acc.hpp"
#include "Com.hpp"
#include "TLP.hpp"

using namespace std;

class BLACWIssuer
{
private:
    vector<Big> com_trapdoors;
    vector<Big> sk_spseq;
public:
    G2& g;
    G2 w, g_tilde, h_tilde;
    G1& g_hat;
    G1 g_hat_p, h_hat_p;
    //PublicParameters pp;

    PFC& bn;
    SPSEQ spseq;
    Acc acc;
    Com com;
    TLP& tlp;
    BLACWIssuer(PFC& _bn, G2& _g, G1& _g_hat, TLP& _tlp);
    ~BLACWIssuer();
    BOOL issue_cred(vector<G2>& C_pre, Big& esk_issuerS, Sigma& sigma, Big& v, ZKP1& zkp);
    BOOL lend(vector<G2>& C_pre, G2* C, Big& c0, Big& c1, Sigma& sigma, Big& esk_issuerS, Sigma& sigma_pr, ZKP2& zkp);
    BOOL credit(Big& bid, G2& C_1S, G2* C_pr, Sigma& sigma_pr, Sigma& sigma_prD, Big v_p, ZKP3& zkp);
    BOOL verify_ZKP1(ZKP1& zkp, vector<G2>& C_pre);
    BOOL verify_ZKP2(ZKP2& zkp, vector<G2>& C_pre, G2* C);
    BOOL verify_ZKP3(ZKP3& zkp, G2& C_pre1, G2& C_1S, Big& bid);
};



#endif