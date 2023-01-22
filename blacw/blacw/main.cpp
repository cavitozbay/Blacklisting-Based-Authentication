#include <iostream>
#include <iomanip>
#include <fstream>
#include <experimental/filesystem>
#include <string.h>
#include <ctime>
#include <cstring>
#include <map>
#include "ecn.h"
#include "ecn2.h"
#include "zzn12a.h"
#include "pairing_bn.h"
#include "Acc.hpp"
#include "SPSEQ.hpp"
#include "TLP.hpp"
#include "BLACWIssuer.hpp"
#include "BLACWHolder.hpp"
#include "helpers.hpp"

using namespace std;
/*#ifdef MR_COUNT_OPS
extern "C"
{
    int fpc=0;
    int fpa=0;
    int fpx=0;
	int fpm2=0;
	int fpi2=0;
	int fpmq=0;
	int fpsq=0;
	int fpaq=0;
}
#endif

#if MIRACL==64
Miracl precision(100,0); 
#else
Miracl precision(100,0);
//Miracl precision=100;
#endif

#ifdef MR_AFFINE_ONLY
    #define AFFINE
#else
    #define PROJECTIVE
#endif*/

/* define minimum duration of each timing, and min. number of iterations */

//Miracl precision(256,0);

void test_acc()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    Acc acc(bn, g, g_hat);

    //modulo(bn.order());

    acc.setup(ACC_SIZE);
    cout << "Setup is done..." << endl;

    for (size_t i = 0; i < 5; i++)
    {
        Big tmp = rand(bn.order());
        acc.add(tmp);
    }
    cout << "Addition is done..." << endl;

    //modulo(bn.order());
    vector<ZZn> nonmembers;

    nonmembers.push_back(ZZn(7));
    // nonmembers.push_back(ZZn(8));

    Witness wit = acc.create_nonmembership(nonmembers);
    modulo(bn.order());
    Big zero = Big(0);
    Big one = Big(1);
    modulo(*bn.mod);
    // cout << "wit1: " << (wit.w1_hat == g_hat) << endl;
    // cout << "wit2: " << (wit.w2.g.iszero()) << endl;
    // cout << "rs: " << (acc.rs[0] == g) << endl;
    Poly f0 = build_from_roots(nonmembers, bn.order());
    modulo(*bn.mod);


    G2 Y = mult_poly(bn, acc.rs, f0);
    //modulo(*bn.mod);

    acc.verify_nonmembership(Y, wit);

}

void test_build_from_roots(G2& g, G1& g_hat, PFC& bn)
{
    GT l = bn.pairing(g,bn.mult(g_hat,Big(5)));
    GT r = bn.pairing(bn.mult(g,Big(5)),g_hat);
    BOOL res = l == r;
    cout << "Build From Roots Test: " << res << endl;

    Acc acc(bn, g, g_hat);

    l = bn.pairing(g,bn.mult(g_hat,Big(5)));
    r = bn.pairing(bn.mult(g,Big(5)),g_hat);
    res = l == r;
    cout << "Build From Roots Test: " << res << endl;
    //modulo(bn.order());

    acc.setup(100);
    cout << "Setup is done..." << endl;

    /*for (size_t i = 0; i < 100; i++)
    {
        Big tmp = rand(bn.order());
        acc.add(tmp);
    }
    cout << "Addition is done..." << endl;

    //modulo(bn.order());*/
    vector<ZZn> nonmembersZZn;
    vector<Big> nonmembersBig;
    //modulo(bn.order());
    modulo(bn.order());
    Big tmp(7);
    cout << "hey";
    nonmembersZZn.push_back(ZZn(tmp));
    nonmembersBig.push_back(tmp);
    modulo(*bn.mod);

    //nonmembers.push_back(ZZn(8));

    //Witness wit = acc.create_nonmembership(nonmembers);
    //modulo(bn.order());
    Poly f0 = build_from_roots(nonmembersZZn, bn.order());
    

    cout << "poly: " << f0 << endl;
    modulo(*bn.mod);
    //modulo(bn.order());

    G2 X =  bn.mult(g, nonmembersBig[0]);
    X = X + acc.rs[1];
    G1 Y = mult_poly_hat(bn, acc.rs_hat, f0);
    
    l = bn.pairing(g,Y);
    r = bn.pairing(X,g_hat);
    //cout << l.g << endl;
    //cout << r.g << endl;
    res = l == r;
    cout << "Build From Roots Test: " << res << endl;
    //modulo(bn.order());
    //acc.verify_nonmembership(Y, wit);
}

void test_recursive_bfr()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    GT l = bn.pairing(g,bn.mult(g_hat,Big(5)));
    GT r = bn.pairing(bn.mult(g,Big(5)),g_hat);
    BOOL res = l == r;
    cout << "Build From Roots Test: " << res << endl;

    Acc acc(bn, g, g_hat);

    l = bn.pairing(g,bn.mult(g_hat,Big(5)));
    r = bn.pairing(bn.mult(g,Big(5)),g_hat);
    res = l == r;
    cout << "Build From Roots Test: " << res << endl;
    //modulo(bn.order());

    acc.setup(100);
    cout << "Setup is done..." << endl;

    /*for (size_t i = 0; i < 100; i++)
    {
        Big tmp = rand(bn.order());
        acc.add(tmp);
    }
    cout << "Addition is done..." << endl;

    //modulo(bn.order());*/
    vector<ZZn> nonmembersZZn;
    vector<Big> nonmembersBig;
    //modulo(bn.order());
    modulo(bn.order());
    Big tmp(7);
    cout << "hey";
    nonmembersZZn.push_back(ZZn(tmp));
    nonmembersBig.push_back(tmp);
    modulo(*bn.mod);

    //nonmembers.push_back(ZZn(8));

    //Witness wit = acc.create_nonmembership(nonmembers);
    //modulo(bn.order());
    Poly f0 = recursive_bfr(nonmembersZZn.size(), 0, nonmembersZZn, bn.order());
    //build_from_roots(nonmembersZZn, bn.order());
    

    cout << "poly: " << f0 << endl;
    modulo(*bn.mod);
    //modulo(bn.order());

    G2 X =  bn.mult(g, nonmembersBig[0]);
    X = X + acc.rs[1];
    G1 Y = mult_poly_hat(bn, acc.rs_hat, f0);
    
    l = bn.pairing(g,Y);
    r = bn.pairing(X,g_hat);
    //cout << l.g << endl;
    //cout << r.g << endl;
    res = l == r;
    cout << "Build From Roots Test: " << res << endl;
    //modulo(bn.order());
    //acc.verify_nonmembership(Y, wit);
}

void test_spseq(G2& g, G1& g_hat, PFC& bn)
{
    size_t l = 3;
    G2* M = new G2[l];

    for (size_t i = 0; i < l; i++)
    {
        bn.random(M[i]);
    }
    
    SPSEQ spseq(l, bn, g, g_hat);

    vector<Big> sk_spseq = spseq.keygen();

    Sigma sigma = spseq.sign(sk_spseq, M);
    //bn.random(M[0]);
    cout << "Test SPSEQ verification: " << spseq.verify(M, sigma) << endl;

    Big mu = rand(bn.order());

    sigma = spseq.chgrep(sigma, M, mu);

    cout << "Test SPSEQ chgrep verification: " << spseq.verify(M, sigma) << endl;

    delete[] M;
}

void test_pcomp_protocol()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);
    
    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
    Witness wit;
    
    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
    
    vS = rand(bn.order());
    gamma = rand(bn.order());
    vector<G2> C_pre_borrow;
    G2 C[3], C_pr[3], C_1S;

    //Borrow1
    holder.precompute_h0s(1);
    ZKP2 zkp2;
    holder.pcomp_borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
    issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
    holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);

    Big bid_;
    vS = rand(bn.order());
    gamma = rand(bn.order());
    vector<G2> C_pre_borrow_;
    G2 C_[3];
    
    //Borrow2
    ZKP2 zkp2_;
    holder.pcomp_borrow1(vS, gamma, c0, c1, wit, C_pre_borrow_, C_, sigma, bid_, zkp2_);
    issuer.lend(C_pre_borrow_, C_, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2_);
    holder.borrow2(C_pre_borrow_, sigma_pr, esk_issuerS, issuer.w);

    //Repay1
    ZKP3 zkp3;
    v_p = rand(bn.order());
    holder.repay1(bid, C_1S, C_pr, sigma_pr, s, zkp3);
    issuer.credit(bid, C_1S, C_pr, sigma_pr, sigma_prD, v_p, zkp3);
    holder.repay2(sigma_prD, C_1S, s, v_p);

    //Repay2
    ZKP3 zkp3_;
    v_p = rand(bn.order());
    G2 C_pr_[3], C_1S_;
    Sigma sigma_pr_, sigma_prD_;
    holder.repay1(bid_, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
    issuer.credit(bid_, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
    holder.repay2(sigma_prD_, C_1S_, s, v_p);

    return;

}


void test_protocol()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);
    
    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
    Witness wit;
    
    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);

    ofstream f;
    f.open("issue_cred1.bin");
    f << zkp1.B_t.g;
    f << zkp1.D_t.g;
    f << zkp1.E.g;
    f << zkp1.E_t.g;
    f << zkp1.F.g;
    f << zkp1.F_t.g;
    f << zkp1.G_t.g;
    f << zkp1.M_t.g;
    f << zkp1.z_delta;
    f << zkp1.z_gamma;
    f << zkp1.z_mult;
    f << zkp1.z_tmp;
    f << zkp1.z_u;
    f << zkp1.z_y;
    for (size_t i = 0; i < zkp1.com_alpha_zs.size(); i++)
    {
        f << zkp1.com_alpha_zs[i];
    }
    for (size_t i = 0; i < C_pre_issue.size(); i++)
    {
        f << C_pre_issue[i].g;
    }
    f << v;
    f << sigma.Y.g;
    f << sigma.Y_hat.g;
    f << sigma.Z.g;
    f << esk_issuerS;
    
    vS = rand(bn.order());
    gamma = rand(bn.order());
    vector<G2> C_pre_borrow;
    G2 C[3], C_pr[3], C_1S;

    //Borrow1
    ZKP2 zkp2;
    holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
    issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
    holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);

    Big bid_;
    vS = rand(bn.order());
    gamma = rand(bn.order());
    vector<G2> C_pre_borrow_;
    G2 C_[3];

    ofstream f1;
    f1.open("borrow.bin");
    
    f1 << zkp2.D_t.g;
    f1 << zkp2.E.g;
    f1 << zkp2.E_t.g;
    f1 << zkp2.F.g;
    f1 << zkp2.F_old_t.g;
    f1 << zkp2.F_t.g;
    f1 << zkp2.G_t.g;
    f1 << zkp2.M2_t.g;
    f1 << zkp2.M_bid_t.g;
    f1 << zkp2.M_t.g;
    f1 << zkp2.N_t.g;
    f1 << zkp2.z_bidS;
    f1 << zkp2.z_delta;
    f1 << zkp2.z_gamma;
    f1 << zkp2.z_m_bidS;
    f1 << zkp2.z_mS;
    f1 << zkp2.z_mult;
    f1 << zkp2.z_r_bidS;
    f1 << zkp2.z_r_y;
    f1 << zkp2.z_r_yS;
    f1 << zkp2.z_tmp;
    f1 << zkp2.z_tmp_bidS;
    f1 << zkp2.z_tmpS;
    f1 << zkp2.z_u;
    f1 << zkp2.z_y;
    f1 << zkp2.z_yS;
    f1 << zkp2.nonmember.C_h.g;
    f1 << zkp2.nonmember.C_hat_h0.g;
    f1 << zkp2.nonmember.D_a.g;
    f1 << zkp2.nonmember.D_a_t.g;
    f1 << zkp2.nonmember.D_b.g;
    f1 << zkp2.nonmember.D_b_t.g;
    f1 << zkp2.nonmember.E_t.g;
    f1 << zkp2.nonmember.M_t.g;
    f1 << zkp2.nonmember.z_a;
    f1 << zkp2.nonmember.z_b;
    f1 << zkp2.nonmember.z_m_acc;
    f1 << zkp2.nonmember.z_r_a;
    f1 << zkp2.nonmember.z_r_b;
    f1 << zkp2.nonmember.z_tmp_acc;
    
    for (size_t i = 0; i < zkp2.com_alpha_zs.size(); i++)
    {
        f1 << zkp2.com_alpha_zs[i];
    }
    for (size_t i = 0; i < zkp2.com_old_alpha_zs.size(); i++)
    {
        f1 << zkp2.com_old_alpha_zs[i];
    }
    for (size_t i = 0; i < C_pre_borrow.size(); i++)
    {
        f1 << C_pre_borrow[i].g;
    }
    for (size_t i = 0; i < 3; i++)
    {
        f1 << C[i].g;
    }
    
    f1 << sigma_pr.Y.g;
    f1 << sigma_pr.Y_hat.g;
    f1 << sigma_pr.Z.g;
    f1 << sigma.Y.g;
    f1 << sigma.Y_hat.g;
    f1 << sigma.Z.g;
    f1 << esk_issuerS;
    

    //Borrow2
    ZKP2 zkp2_;
    holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow_, C_, sigma, bid_, zkp2_);
    issuer.lend(C_pre_borrow_, C_, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2_);
    holder.borrow2(C_pre_borrow_, sigma_pr, esk_issuerS, issuer.w);


    //Repay1
    ZKP3 zkp3;
    v_p = rand(bn.order());
    holder.repay1(bid, C_1S, C_pr, sigma_pr, s, zkp3);
    issuer.credit(bid, C_1S, C_pr, sigma_pr, sigma_prD, v_p, zkp3);
    holder.repay2(sigma_prD, C_1S, s, v_p);

    ofstream f2;
    f2.open("repay.bin");
    
    f2 << bid;
    f2 << zkp3.C_y.g;
    f2 << zkp3.C_y_t.g;
    f2 << zkp3.C_yS.g;
    f2 << zkp3.C_yS_t.g;
    f2 << zkp3.R_t.g;
    f2 << zkp3.z_r_y;
    f2 << zkp3.z_r_yS;
    f2 << zkp3.z_y;
    f2 << zkp3.z_yS;
    f2 << C_1S.g;
    for (size_t i = 0; i < 3; i++)
    {
        f2 << C_pr[i].g;
    }
    f2 << sigma_pr.Y.g;
    f2 << sigma_pr.Y_hat.g;
    f2 << sigma_pr.Z.g;
    f2 << sigma.Y.g;
    f2 << sigma.Y_hat.g;
    f2 << sigma.Z.g;

    //Repay2
    ZKP3 zkp3_;
    v_p = rand(bn.order());
    G2 C_pr_[3], C_1S_;
    Sigma sigma_pr_, sigma_prD_;
    holder.repay1(bid_, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
    issuer.credit(bid_, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
    holder.repay2(sigma_prD_, C_1S_, s, v_p);

    f.close();
    f1.close();
    f2.close();
}

void test_com(G2& g, G1& g_hat, PFC& bn)
{
    Com com(bn, 2);
    vector<Big> alphas, alpha_ts, z;
    G2 u_t;

    for (size_t i = 0; i < 2; i++)
    {
        alphas.push_back(rand(bn.order()));
    }
    
    CmOp cmop = com.commit(alphas);

    com.pok_init(alpha_ts, u_t);
    Big c = rand(bn.order());
    com.pok_response(alphas, alpha_ts, c, z);
    cout << "Com trial: " << com.pok_verify(cmop.C, u_t, z, c) << endl;

}

void bmark_multG1(G2& g, G1& g_hat, PFC& bn)
{
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    G1 tmp;
    Big m = rand(bn.order());
    do {
        tmp = bn.mult(g_hat, m);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Mult G1" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_multG2(G2& g, G1& g_hat, PFC& bn)
{
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    G2 tmp;
    Big m = rand(bn.order());
    do {
        tmp = bn.mult(g, m);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Mult G2" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_issuance(G2& g, G1& g_hat, PFC& bn)
{
    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
    Witness wit;
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    do {
        vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
        issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
        holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Issuance" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_borrowing(G2& g, G1& g_hat, PFC& bn)
{
    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
    Witness wit;
    
    int iterations = 0;
    double elapsed = 0;
    clock_t start = clock();
    clock_t total = 0;
    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
    clock_t start_ = clock();
    G2 C[3], C_pr[3], C_1S;

    ZKP2 zkp2;
    do {
        start = clock();
        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
        iterations++;
        total += clock() - start;
        elapsed=(total-start_)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Borrowing" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
}

void bmark_borrowing_advanced(G2& g, G1& g_hat, PFC& bn)
{
    TLP tlp;

    
    
    int iterations = 0;
    double elapsed = 0;
    clock_t start;
    clock_t total = 0;

    do {
        ZKP1 zkp1;
        Sigma sigma, sigma_pr, sigma_prD;
        Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
        Witness wit;
        ZKP2 zkp2;

        BLACWIssuer issuer(bn, g, g_hat, tlp);
        //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
        BLACWHolder holder(bn, issuer);

        vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
        issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
        holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
        
        G2 C[3], C_pr[3], C_1S;

        start = clock();
        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
        iterations++;
        total += clock() - start;
        elapsed = (total)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Borrowing" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_borrowing_advanced2(size_t init_debt_count)
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    TLP tlp;
    cout << "DEBT: " << init_debt_count << ", " << "BLACKLISTED: " << BLACKLISTED_COUNT << endl;
    
    int iterations = 0;
    double elapsed = 0;
    clock_t start;
    clock_t total = 0;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v;
    Witness wit;

    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);


    for (size_t i = 0; i < init_debt_count; i++)
    {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    }
    

    do {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        start = clock();
        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
        iterations++;
        total += clock() - start;
        elapsed = (total)/(double)CLOCKS_PER_SEC;

        ZKP3 zkp3_;
        v_p = rand(bn.order());
        G2 C_pr_[3], C_1S_;
        Sigma sigma_pr_, sigma_prD_;
        holder.repay1(bid, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
        issuer.credit(bid, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
        holder.repay2(sigma_prD_, C_1S_, s, v_p);

    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Borrowing" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_borrowing_advanced3(size_t init_debt_count)
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    TLP tlp;
    cout << "DEBT: " << init_debt_count << ", " << "BLACKLISTED: " << BLACKLISTED_COUNT << endl;
    
    int iterations = 0;
    double elapsed_b1 = 0;
    double elapsed_b2 = 0;
    double elapsed_l = 0;
    clock_t start_b1;
    clock_t start_b2;
    clock_t start_l;
    clock_t total_b1 = 0;
    clock_t total_b2 = 0;
    clock_t total_l = 0;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v;
    Witness wit;

    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);


    for (size_t i = 0; i < init_debt_count; i++)
    {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    }
    

    do {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        start_b1 = clock();
        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        total_b1 += clock() - start_b1;
        elapsed_b1 = (total_b1)/(double)CLOCKS_PER_SEC;
        
        start_l = clock();
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        total_l += clock() - start_l;
        elapsed_l = (total_l)/(double)CLOCKS_PER_SEC;

        start_b2 = clock();
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
        total_b2 += clock() - start_b2;
        elapsed_b2 = (total_b2)/(double)CLOCKS_PER_SEC;

        iterations++;
        ZKP3 zkp3_;
        v_p = rand(bn.order());
        G2 C_pr_[3], C_1S_;
        Sigma sigma_pr_, sigma_prD_;
        holder.repay1(bid, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
        issuer.credit(bid, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
        holder.repay2(sigma_prD_, C_1S_, s, v_p);

    } while (elapsed_b1<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed_b1=1000.0*elapsed_b1/iterations;
    elapsed_l=1000.0*elapsed_l/iterations;
    elapsed_b2=1000.0*elapsed_b2/iterations;
    cout << setw(20) << "Borrowing 1" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed_b1 << " ms per iteration" << endl;
    cout << setw(20) << "Lending" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed_l << " ms per iteration" << endl;
    cout << setw(20) << "Borrowing 2" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed_b2 << " ms per iteration" << endl;
    
}

void bmark_all()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    TLP tlp;
    //cout << "DEBT: " << INITIAL_DEBT_COUNT << ", " << "BLACKLISTED: " << BLACKLISTED_COUNT << endl;
    int iterations = 0;
    double elapsed_b1 = 0;
    double elapsed_b2 = 0;
    double elapsed_l = 0;
    double elapsed_r1 = 0;
    double elapsed_r2 = 0;
    double elapsed_c = 0;
    clock_t start_b1;
    clock_t start_b2;
    clock_t start_l;
    clock_t start_r1;
    clock_t start_r2;
    clock_t start_c;
    clock_t total_b1 = 0;
    clock_t total_b2 = 0;
    clock_t total_l = 0;
    clock_t total_r1 = 0;
    clock_t total_r2 = 0;
    clock_t total_c = 0;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v;
    Witness wit;

    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
    BOOL res = TRUE;

    for (size_t i = 0; i < INITIAL_DEBT_COUNT; i++)
    {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        res = issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    }
    

    do {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        start_b1 = clock();
        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        total_b1 += clock() - start_b1;
        elapsed_b1 = (total_b1)/(double)CLOCKS_PER_SEC;
        
        start_l = clock();
        res = res && issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        total_l += clock() - start_l;
        elapsed_l = (total_l)/(double)CLOCKS_PER_SEC;

        start_b2 = clock();
        res = res && holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
        total_b2 += clock() - start_b2;
        elapsed_b2 = (total_b2)/(double)CLOCKS_PER_SEC;

        ZKP3 zkp3_;
        G2 C_pr_[3], C_1S_;
        Sigma sigma_pr_, sigma_prD_;

        start_r1 = clock();
        v_p = rand(bn.order());
        holder.repay1(bid, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
        total_r1 += clock() - start_r1;
        elapsed_r1 = (total_r1)/(double)CLOCKS_PER_SEC;

        start_c = clock();
        res = res && issuer.credit(bid, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
        total_c += clock() - start_c;
        elapsed_c = (total_c)/(double)CLOCKS_PER_SEC;

        start_r2 = clock();
        res = res && holder.repay2(sigma_prD_, C_1S_, s, v_p);
        total_r2 += clock() - start_r2;
        elapsed_r2 = (total_r2)/(double)CLOCKS_PER_SEC;

        iterations++;
        if(res != TRUE)
        {
            cout << "ERROR IN " << iterations << endl;
            break;
        }
    } while (elapsed_b1<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed_b1=1000.0*elapsed_b1/iterations;
    elapsed_l=1000.0*elapsed_l/iterations;
    elapsed_b2=1000.0*elapsed_b2/iterations;
    elapsed_r1=1000.0*elapsed_r1/iterations;
    elapsed_c=1000.0*elapsed_c/iterations;
    elapsed_r2=1000.0*elapsed_r2/iterations;
    cout << setw(20) << "B1 - L - B2" << " -;" << setw(5) << INITIAL_DEBT_COUNT << ";" << setw(5) << BLACKLISTED_COUNT << ";" << setw(8) << setprecision(2) << fixed << elapsed_b1 << "; " << setw(8) << setprecision(2) << fixed << elapsed_l << "; " << setw(8) << setprecision(2) << fixed << elapsed_b2 << "; " << endl;
    
    //cout << setw(20) << "Borrowing 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Lending" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_l << "; ms per iteration" << endl;
    //cout << setw(20) << "Borrowing 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b2 << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Credit" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_c << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r2 << "; ms per iteration" << endl;
    
}

void bmark_all_ext()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    TLP tlp;
    //cout << "DEBT: " << INITIAL_DEBT_COUNT << ", " << "BLACKLISTED: " << BLACKLISTED_COUNT << endl;
    int iterations = 0;
    double elapsed_b1 = 0;
    double elapsed_b2 = 0;
    double elapsed_l = 0;
    double elapsed_r1 = 0;
    double elapsed_r2 = 0;
    double elapsed_c = 0;
    clock_t start_b1;
    clock_t start_b2;
    clock_t start_l;
    clock_t start_r1;
    clock_t start_r2;
    clock_t start_c;
    map<size_t,clock_t> total_b1;
    map<size_t,clock_t> total_b2;
    map<size_t,clock_t> total_l;

    for (size_t i = 0; i <= 200; i+=10)
    {
        total_b1[i] = 0;
        total_b2[i] = 0;
        total_l[i] = 0;
    }
    
    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v;
    Witness wit;

    BOOL res = TRUE;

    // for (size_t i = 0; i < INITIAL_DEBT_COUNT; i++)
    // {
    //     ZKP2 zkp2;    
    //     G2 C[3], C_pr[3], C_1S;
    //     Big bid, vS, v_p, c0, c1, s, gamma;

    //     vS = rand(bn.order());
    //     gamma = rand(bn.order());
    //     vector<G2> C_pre_borrow;
    //     holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
    //     res = issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
    //     holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    // }
    

    do {
        BLACWHolder holder(bn, issuer);
        vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
        issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
        holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);

        for (size_t i = 0; i <= 200; i+=10)
        {
            ZKP1 zkp1;
            Sigma sigma, sigma_pr, sigma_prD;
            Big esk_issuerS, v;
            Witness wit;

            ZKP2 zkp2;    
            G2 C[3], C_pr[3], C_1S;
            Big bid, vS, v_p, c0, c1, s, gamma;

            vS = rand(bn.order());
            gamma = rand(bn.order());
            vector<G2> C_pre_borrow;
            C_pre_borrow.clear();
            start_b1 = clock();
            holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
            total_b1[i] += clock() - start_b1;
            
            start_l = clock();
            res = res && issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
            total_l[i] += clock() - start_l;

            start_b2 = clock();
            res = res && holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
            total_b2[i] += clock() - start_b2;

            for (size_t j = 0; j < 9; j++)
            {
                ZKP1 zkp1;
                Sigma sigma, sigma_pr, sigma_prD;
                Big esk_issuerS, v;
                Witness wit;

                ZKP2 zkp2;    
                G2 C[3], C_pr[3], C_1S;
                Big bid, vS, v_p, c0, c1, s, gamma;
                vS = rand(bn.order());
                gamma = rand(bn.order());
                C_pre_borrow.clear();
                holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
                res = res && issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
                res = res && holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
                // cout << "j:" << j << endl;
            }
                // cout << "i:" << i << endl;

            
        }
        // cout << "itr:" << iterations << endl;
        
        iterations++;

        if(res != TRUE)
        {
            cout << "ERROR IN " << iterations << endl;
            break;
        }
    } while (iterations<MIN_ITERS);
    
    for (size_t i = 0; i <= 200; i+=10)
    {
        elapsed_b1=1000.0*((total_b1[i])/(double)CLOCKS_PER_SEC)/iterations;
        elapsed_l=1000.0*((total_l[i])/(double)CLOCKS_PER_SEC)/iterations;
        elapsed_b2=1000.0*((total_b2[i])/(double)CLOCKS_PER_SEC)/iterations;
        cout << setw(20) << "B1 - L - B2" << " -;" << setw(5) << i << ";" << setw(5) << BLACKLISTED_COUNT << ";" << setw(8) << setprecision(2) << fixed << elapsed_b1 << "; " << setw(8) << setprecision(2) << fixed << elapsed_l << "; " << setw(8) << setprecision(2) << fixed << elapsed_b2 << "; " << endl;
    }
    
    
    //cout << setw(20) << "Borrowing 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Lending" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_l << "; ms per iteration" << endl;
    //cout << setw(20) << "Borrowing 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b2 << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Credit" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_c << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r2 << "; ms per iteration" << endl;
    
}



void bmark_pcomp()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);


    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);

    //cout << "DEBT: " << INITIAL_DEBT_COUNT << ", " << "BLACKLISTED: " << BLACKLISTED_COUNT << endl;
    int iterations = 0;
    double elapsed_pc_b1[PCOMP_SIZE];
    double elapsed_pc_l[PCOMP_SIZE];
    double elapsed_pc_b2[PCOMP_SIZE];
    clock_t start_pc;
    clock_t start_pc_b1[PCOMP_SIZE];
    clock_t start_pc_l[PCOMP_SIZE];
    clock_t start_pc_b2[PCOMP_SIZE];
    clock_t total_pc_b1[PCOMP_SIZE];
    clock_t total_pc_b2[PCOMP_SIZE];
    clock_t total_pc_l[PCOMP_SIZE];

    for (size_t i = 0; i < PCOMP_SIZE; i++)
    {
        elapsed_pc_b1[i] = 0;
        elapsed_pc_l[i] = 0;
        elapsed_pc_b2[i] = 0;
        total_pc_b1[i] = 0;
        total_pc_b2[i] = 0;
        total_pc_l[i] = 0;
    }
    
    double elapsed_pc = 0;
    clock_t total_pc = 0;

    do {
        //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
        BLACWHolder holder(bn, issuer);

        ZKP1 zkp1;
        Sigma sigma, sigma_pr, sigma_prD;
        Big esk_issuerS, v;
        Witness wit;

        vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
        issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
        holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
        BOOL res = TRUE;
        
        start_pc = clock();
        holder.precompute_h0s(PCOMP_SIZE);
        total_pc += clock() - start_pc;
        elapsed_pc = (total_pc)/(double)CLOCKS_PER_SEC;

        for (size_t i = 0; i < PCOMP_SIZE; i++)
        {
            ZKP2 zkp2;    
            G2 C[3], C_pr[3], C_1S;
            Big bid, vS, v_p, c0, c1, s, gamma;

            vS = rand(bn.order());
            gamma = rand(bn.order());
            vector<G2> C_pre_borrow;
            start_pc_b1[i] = clock();
            holder.pcomp_borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
            total_pc_b1[i] += clock() - start_pc_b1[i];
            elapsed_pc_b1[i] = (total_pc_b1[i])/(double)CLOCKS_PER_SEC;
            
            start_pc_l[i] = clock();
            res = res && issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
            total_pc_l[i] += clock() - start_pc_l[i];
            elapsed_pc_l[i] = (total_pc_l[i])/(double)CLOCKS_PER_SEC;

            start_pc_b2[i] = clock();
            res = res && holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
            total_pc_b2[i] += clock() - start_pc_b2[i];
            elapsed_pc_b2[i] = (total_pc_b2[i])/(double)CLOCKS_PER_SEC;
        }
        
        

        iterations++;
        if(res != TRUE)
        {
            cout << "ERROR IN " << iterations << endl;
            break;
        }
    } while (elapsed_pc<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed_pc=1000.0*elapsed_pc/iterations;
    for (size_t i = 0; i < PCOMP_SIZE; i++)
    {
        elapsed_pc_b1[i]=1000.0*elapsed_pc_b1[i]/iterations;
        elapsed_pc_l[i]=1000.0*elapsed_pc_l[i]/iterations;
        elapsed_pc_b2[i]=1000.0*elapsed_pc_b2[i]/iterations;   
        cout << setw(20) << "PC - B1 - L - B2" << " -;"  << setw(5) << PCOMP_SIZE << ";" << setw(5) << BLACKLISTED_COUNT << ";" << setw(8) << setprecision(2) << fixed << elapsed_pc << "; " << setw(8) << setprecision(2) << fixed << elapsed_pc_b1[i] << "; " << setw(8) << setprecision(2) << fixed << elapsed_pc_l[i] << "; " << setw(8) << setprecision(2) << fixed << elapsed_pc_b2[i] << "; " << endl;
    }
    
    
    
    //cout << setw(20) << "Borrowing 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Lending" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_l << "; ms per iteration" << endl;
    //cout << setw(20) << "Borrowing 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_b2 << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 1" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r1 << "; ms per iteration" << endl;
    //cout << setw(20) << "Credit" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_c << "; ms per iteration" << endl;
    //cout << setw(20) << "Repayment 2" << " - " << setw(5) << iterations << " iterations;" << setw(8) << setprecision(2) << fixed << elapsed_r2 << "; ms per iteration" << endl;
    
}

void bmark_repayment_advanced2(G2& g, G1& g_hat, PFC& bn)
{
    TLP tlp;
    
    int iterations = 0;
    double elapsed = 0;
    clock_t start;
    clock_t total = 0;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v;
    Witness wit;

    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);

    for (size_t i = 0; i < INITIAL_DEBT_COUNT; i++)
    {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    }

    do {
        ZKP2 zkp2;    
        G2 C[3], C_pr[3], C_1S;
        Big bid, vS, v_p, c0, c1, s, gamma;

        vS = rand(bn.order());
        gamma = rand(bn.order());
        vector<G2> C_pre_borrow;
        holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
        issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
        holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);

        ZKP3 zkp3_;

        start = clock();
        v_p = rand(bn.order());
        G2 C_pr_[3], C_1S_;
        Sigma sigma_pr_, sigma_prD_;
        holder.repay1(bid, C_1S_, C_pr_, sigma_pr_, s, zkp3_);
        issuer.credit(bid, C_1S_, C_pr_, sigma_pr_, sigma_prD_, v_p, zkp3_);
        holder.repay2(sigma_prD_, C_1S_, s, v_p);
        
        iterations++;
        total += clock() - start;
        elapsed = (total)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);
    
    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Repayment" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}


void bmark_repayment(G2& g, G1& g_hat, PFC& bn)
{
    TLP tlp;

    BLACWIssuer issuer(bn, g, g_hat, tlp);
    //BLACWHolder holder(bn, g, g_hat, issuer.com, issuer.acc, issuer.spseq, issuer.w, issuer.g_tilde, issuer.h_tilde);
    BLACWHolder holder(bn, issuer);

    ZKP1 zkp1;
    Sigma sigma, sigma_pr, sigma_prD;
    Big esk_issuerS, v, vS, gamma, c0, c1, bid, s, v_p;
    Witness wit;
    
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    vector<G2> C_pre_issue = holder.receive_cred1(zkp1);
    issuer.issue_cred(C_pre_issue, esk_issuerS, sigma, v, zkp1);
    holder.receive_cred2(C_pre_issue, v, sigma, esk_issuerS, issuer.w);
    
    G2 C[3], C_pr[3], C_1S;

    ZKP2 zkp2;

    vS = rand(bn.order());
    gamma = rand(bn.order());
    vector<G2> C_pre_borrow;
    holder.borrow1(vS, gamma, c0, c1, wit, C_pre_borrow, C, sigma, bid, zkp2);
    issuer.lend(C_pre_borrow, C, c0, c1, sigma, esk_issuerS, sigma_pr, zkp2);
    holder.borrow2(C_pre_borrow, sigma_pr, esk_issuerS, issuer.w);
    ZKP3 zkp3;
    do {
        holder.repay1(bid, C_1S, C_pr, sigma_pr, s, zkp3);
        issuer.credit(bid, C_1S, C_pr, sigma_pr, sigma_prD, v_p, zkp3);
        holder.repay2(sigma_prD, C_1S, s, v_p);

        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Repayment" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}



void test_tlp()
{
    miracl* mip = mirsys(256,0);
    mip->IOBASE = 16;

    Big T(1000);
    Big s = rand(Big(100000));
    cout << "number: " << s << endl;
    TLP tlp(T);
    
    PuzzleAux paux = tlp.pgen(s);
    //PFC bn(128);
    ofstream d("rsa.bin");
    d << tlp.n;
    d.close();
    cout << "Evaluated puzzle: " << (s == tlp.evaluate(paux.puzzle)) << endl;

    ZKPTLP zkp = tlp.prove_wellformedness(paux);
    tlp.verify_wellformedness(paux.puzzle, zkp);

    return;
    ofstream f("BLACW_tlp.bin");

    f << zkp.l / tlp.n << zkp.l % tlp.n;
    f << zkp.t_l / tlp.n << zkp.t_l % tlp.n;
    f << zkp.t_u;
    f << zkp.t_v / tlp.n << zkp.t_v % tlp.n;
    f << zkp.z_d / tlp.n << zkp.z_d % tlp.n;
    f << zkp.z_r / tlp.n << zkp.z_d % tlp.n;
    f << zkp.z_s / tlp.n << zkp.z_d % tlp.n;
    f << paux.puzzle.u;
    f << paux.puzzle.v / tlp.n << paux.puzzle.v % tlp.n;
    f.close();
    ifstream file("BLACW_tlp.bin", ios::binary);
    const auto begin = file.tellg();
    file.seekg (0, ios::end);
    const auto end = file.tellg();
    const auto fsize = (end-begin);
    cout << fsize << endl;
    tlp.verify_wellformedness(paux.puzzle, zkp);
}

void bmark_tlp_pgen()
{
    miracl* mip = mirsys(256,0);
    mip->IOBASE = 16;

    Big T(1000);
    Big s = rand(2,RHO_SIZE);

    TLP tlp(T);
    ZKPTLP zkp;
    PuzzleAux paux;
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    do {
        paux = tlp.pgen(s);
        zkp = tlp.prove_wellformedness(paux);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "TLP Generate" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;

    
}

void bmark_tlp_verify()
{
    miracl* mip = mirsys(256,0);
    mip->IOBASE = 16;

    Big T(1000);
    Big s = rand(2,RHO_SIZE);

    TLP tlp(T);
    ZKPTLP zkp;
    PuzzleAux paux;

    paux = tlp.pgen(s);
    zkp = tlp.prove_wellformedness(paux);
    
    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    do {
        tlp.verify_wellformedness(paux.puzzle, zkp);    
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "TLP Verify" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    //printf("TLP Verify - %8d iterations             ",iterations);
    //printf(" %8.2lf ms per iteration\n",elapsed);
    
}

void bmark_build_from_roots()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    vector<ZZn> roots;    
    for (size_t i = 0; i < BLACKLISTED_COUNT; i++)
    {
        roots.push_back(ZZn(rand(bn.order())));
    }

    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    do {
        Poly f = build_from_roots(roots, bn.order());
        modulo(*bn.mod);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Build From Roots" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;
    
}

void bmark_mult_poly()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    Acc acc(bn, g, g_hat);

    acc.setup(ACC_SIZE);

    vector<ZZn> roots;    
    for (size_t i = 0; i < BLACKLISTED_COUNT; i++)
    {
        roots.push_back(ZZn(rand(bn.order())));
    }
    
    Poly f = build_from_roots(roots, bn.order());
    modulo(*bn.mod);

    int iterations = 0;
    double elapsed;
    clock_t start = clock();
    
    do {
        G2 tmp = mult_poly(bn, acc.rs, f);
        iterations++;
        elapsed=(clock()-start)/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << setw(20) << "Mult Poly" << " - " << setw(5) << iterations << " iterations" << setw(8) << setprecision(2) << fixed << elapsed << " ms per iteration" << endl;    
}

void test_fast_addition()
{
    PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);

    Acc acc(bn, g, g_hat);

    acc.setup(ACC_SIZE);
    modulo(bn.order());

    vector<ZZn> roots;    
    for (size_t i = 0; i < 2; i++)
    {
        roots.push_back(ZZn(rand(bn.order())));
    }
    roots.push_back(ZZn(0));
    Poly f0 = build_from_roots(roots, bn.order());
    modulo(*bn.mod);

    Poly f = build_from_roots(acc.elementsZZn, bn.order());
    modulo(bn.order());
    //cout << "  f: " << f << endl;
    //cout << "  f0: " << f0 << endl;
    
    Poly result[3];
    egcd(result, f, f0);
    
    
    ZZn divisor = result[0].coeff(0);
    result[0] /= divisor;
    result[1] /= divisor;  
    result[2] /= divisor;
    Variable x;
    Poly _h0p = result[2]*x;

    G1 h0 = mult_poly_hat(bn, acc.rs_hat,result[2]);
    G1 h0p = mult_poly_hat(bn, acc.rs_hat,_h0p);
    G2 Y = mult_poly(bn, acc.rs, f0);
    G2 h = mult_poly(bn, acc.rs, result[1]);
    modulo(bn.order());

    ZZn yp = roots.back();
    roots.pop_back();
    Poly f0p = build_from_roots(roots, bn.order());
    vector<ZZn> rootsp;
    rootsp.push_back(yp);

    cout << "egcd: " << result[1]*f + result[2]*f0 << endl;
    cout << "egcd: " << result[1]*f + result[2]*f0p*build_from_roots(rootsp, bn.order()) << endl;

    modulo(*bn.mod);
    GT res1 = bn.pairing(Y, h0); 
    GT res2 = bn.pairing(h, acc.V_hat);

    BOOL res = res2 * res1 == acc.gT;
    cout << res << endl;

    G2 Yp = mult_poly(bn, acc.rs, f0p);
    G1 ht0 = h0p + bn.mult(h0, yp);
    res1 = bn.pairing(Yp, ht0);

    res = res2*res1 == acc.gT;
    cout << res << endl;
}

int main()
{
    /*miracl* mip = &precision;
    mip->IOBASE = 16;*/
    
    long seed = rand()%1000000000;
    //test_fast_addition();
    //BN_Pairing::BNPairing bn(mip);
    
    /*cout << "primes: " << endl;
    cout << "    " << bn.order() << endl;
    cout << "    " << *bn.mod << endl;*/
    
    //modulo(*bn.mod);
    //irand(seed);
    
    /*PFC bn(128);
    G2 g;
    G1 g_hat;

    bn.random(g);
    bn.random(g_hat);

    bn.precomp_for_mult(g);
    bn.precomp_for_mult(g_hat);
    cout << "order: " << bn.order() << endl;
    cout << "order size: " << bits(bn.order()) << endl;*/

    //bn.precomp_for_pairing(g);

    //modulo(bn.order());
    //test_acc();
    //test_spseq(g, g_hat, bn);
    //test_pcomp_protocol();
    test_protocol();
    //test_build_from_roots(g, g_hat, bn)
    //test_recursive_bfr();
    //test_com(g, g_hat, bn);
    //test_tlp();


    //bmark_multG1(g, g_hat, bn);
    //bmark_multG2(g, g_hat, bn);
    // bmark_tlp_pgen();
    // bmark_tlp_verify();
    //bmark_issuance(g, g_hat, bn);
    //bmark_borrowing(g, g_hat, bn);
    //bmark_borrowing_advanced(g, g_hat, bn);
    //bmark_borrowing_advanced2(g, g_hat, bn);
    //bmark_borrowing_advanced3(g, g_hat, bn);
    // PCOMP_SIZE = 30;
    // BLACKLISTED_COUNT = 4000;
    // ACC_SIZE = 4031;
    //bmark_pcomp();
    // for (BLACKLISTED_COUNT = 1000, ACC_SIZE = 1150; BLACKLISTED_COUNT <= 4000; BLACKLISTED_COUNT += 1500, ACC_SIZE += 1500)
    // {
    //     bmark_all_ext();
    // }
    
    
    //bmark_pcomp();
        
    
    //bmark_repayment_advanced2(g, g_hat, bn);
    //bmark_repayment(g, g_hat, bn);
    /*for (BLACKLISTED_COUNT = 100, ACC_SIZE = 120; BLACKLISTED_COUNT < 3000; BLACKLISTED_COUNT += 300, ACC_SIZE += 300)
    {
        bmark_mult_poly();
    }
    /*for (BLACKLISTED_COUNT = 100, ACC_SIZE = 120; BLACKLISTED_COUNT < 3000; BLACKLISTED_COUNT += 300, ACC_SIZE += 300)
    {
        bmark_build_from_roots();
    }*/
    //bmark_mult_poly();
    //bmark_build_from_roots();
    return 0;
}