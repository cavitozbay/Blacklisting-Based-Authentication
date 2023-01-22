#include "PPBIssuer.hpp"

PPBIssuer::PPBIssuer(PFC& _bn, G2& _g, G1& _g_hat, TLP& _tlp) : 
bn(_bn), spseq(3, bn, _g, _g_hat), acc(bn, _g, _g_hat), com(_bn), g(_g), g_hat(_g_hat), tlp(_tlp)
{
    /*g = _g;
    g_hat = _g_hat;*/
    //cout << "g I: " << &g_tilde << endl;
    vector<G2> gens;

    acc.setup(ACC_SIZE);
    sk_spseq = spseq.keygen();
    bn.random(w);

    for (size_t i = 0; i < 4; i++)
    {
        com_trapdoors.push_back(rand(bn.order()));
        gens.push_back(bn.mult(g, com_trapdoors[i]));
    }
    G2 h7;
    bn.random(h7);
    bn.random(g_tilde);
    bn.random(h_tilde);
    bn.random(g_hat_p);
    bn.random(h_hat_p);
    bn.precomp_for_mult(g_tilde);
    bn.precomp_for_mult(h_tilde);
    bn.precomp_for_mult(g_hat_p);
    bn.precomp_for_mult(h_hat_p);
    //bn.precomp_for_mult(h7);

    gens.push_back(h7);
    com.init_gens(gens);
}

PPBIssuer::~PPBIssuer()
{
}

BOOL PPBIssuer::issue_cred(vector<G2>& C_pre, Big& esk_issuerS, Sigma& sigma, Big& v, ZKP1& zkp)
{
    
    //ZKP1
    BOOL res = verify_ZKP1(zkp, C_pre);

    
    G2 C[3];
    
    esk_issuerS = rand(bn.order());
    v = rand(bn.order());
    Big exp = modmult(com_trapdoors[1], esk_issuerS, bn.order()) + 
              modmult(com_trapdoors[2], v, bn.order())%bn.order();
    C[0] = C_pre[0] + bn.mult(C_pre[2], exp);
    C[1] = C_pre[1];
    C[2] = C_pre[2];

    sigma = spseq.sign(sk_spseq, C);

    return res;
}

BOOL PPBIssuer::lend(vector<G2>& C_pre, G2* C, Big& c0, Big& c1, Sigma& sigma, Big& esk_issuerS, Sigma& sigma_pr, ZKP2& zkp)
{

    BOOL res = spseq.verify(C, sigma);
#ifdef LOG
    cout << "Lend SPSEQ verification: " << res << endl;
#endif
    // ZKP
    res = res && verify_ZKP2(zkp, C_pre, C);

    esk_issuerS = rand(bn.order());
    vector<G2> M;
    G2 _C[3];

    Big exp = modmult(com_trapdoors[1], esk_issuerS, bn.order());
    _C[0] = C_pre[0] + bn.mult(C_pre[2], exp);
    _C[1] = C_pre[1];
    _C[2] = C_pre[2];

    sigma_pr = spseq.sign(sk_spseq, _C);

    return res;
}

BOOL PPBIssuer::credit(Big& bid, G2& C_1S, G2* C_pr, Sigma& sigma_pr, Sigma& sigma_prD, Big v_p, ZKP3& zkp)
{
    BOOL res = spseq.verify(C_pr, sigma_pr);
#ifdef LOG
    cout << "Credit SPSEQ verification: " << res << endl;
#endif

    res = res && verify_ZKP3(zkp, C_pr[1], C_1S, bid);
    G2 C[3];
    C[0] = C_pr[0] + bn.mult(C_pr[2], modmult(com_trapdoors[2], v_p, bn.order()));
    C[1] = C_1S;
    C[2] = C_pr[2];

    sigma_prD = spseq.sign(sk_spseq, C);

    return res;
}

BOOL PPBIssuer::verify_ZKP1(ZKP1& zkp, vector<G2>& C_pre)
{
    //ZKP Verification

    Big c = H2(zkp.F_t);

    //F_t
    G2 l = bn.mult(h_tilde, zkp.z_gamma);

    l = l + bn.mult(zkp.F, c);

    for (size_t i = 0; i < com.k; i++)
    {
        l = l + bn.mult(com.gens[i], zkp.com_alpha_zs[i]);
    }

    BOOL res1 = (l == zkp.F_t);

    //B_t
    G2 res2L = bn.mult(C_pre[2], c) + bn.mult(g, zkp.z_u);
    BOOL res2 = (res2L == zkp.B_t);

    //D_t
    G2 res3L = bn.mult(C_pre[0], c) + bn.mult(zkp.F, zkp.z_u) + bn.mult(h_tilde, -zkp.z_mult);
    BOOL res3 = (res3L == zkp.D_t);

    //E_t
    BOOL res4 = (bn.mult(zkp.E, c) + bn.mult(g_hat_p, zkp.z_gamma) + bn.mult(h_hat_p, zkp.z_delta) == zkp.E_t);

    //G_t
    G2 res5L = bn.mult(C_pre[1], c) + bn.mult(C_pre[2], zkp.z_y);
    BOOL res5 = (res5L == zkp.G_t);
    
    //M_t
    BOOL res6 = (bn.mult(zkp.E, zkp.z_u) + bn.mult(g_hat_p, (-zkp.z_mult)%bn.order()) + bn.mult(h_hat_p, (-zkp.z_tmp)%bn.order()) == zkp.M_t);

#ifdef LOG    
    cout << "ZKP1=========" << endl;
    cout << "  F_t: " << res1 << endl;
    cout << "  B_t: " << res2 << endl;
    cout << "  D_t: " << res3 << endl;
    cout << "  E_t: " << res4 << endl;
    cout << "  G_t: " << res5 << endl;
    cout << "  M_t: " << res6 << endl;
#endif
    return res1 && res2 && res3 && res4 && res5 && res6;
}

BOOL PPBIssuer::verify_ZKP2(ZKP2& zkp, vector<G2>& C_pre, G2* C)
{

    //ZKP Verification
    Big c = H2(zkp.F_t);

    //F_t
    G2 l = bn.mult(h_tilde, zkp.z_gamma);

    l = l + bn.mult(zkp.F, c);

    for (size_t i = 0; i < com.k; i++)
    {
        l = l + bn.mult(com.gens[i], zkp.com_alpha_zs[i]);
    }

    BOOL res1 = (l == zkp.F_t);

    //B_t
    G2 res2L = bn.mult(C_pre[2], c) + bn.mult(g, zkp.z_u);
    BOOL res2 = (res2L == zkp.B_t);

    //D_t
    G2 res3L = bn.mult(C_pre[0], c) + bn.mult(zkp.F, zkp.z_u) + bn.mult(h_tilde, -zkp.z_mult);
    BOOL res3 = (res3L == zkp.D_t);

    //E_t
    BOOL res4 = (bn.mult(zkp.E, c) + bn.mult(g_hat_p, zkp.z_gamma) + bn.mult(h_hat_p, zkp.z_delta) == zkp.E_t);

    //G_t
    //res = (bn.mult(C_pre[1], c) + bn.mult(C_pre[2], zkp.z_y) == zkp.G_t);
    //cout << "G_t: " << res << endl;

    //M_t
    BOOL res5 = (bn.mult(zkp.E, zkp.z_u) + bn.mult(g_hat_p, (-zkp.z_mult)%bn.order()) + bn.mult(h_hat_p, (-zkp.z_tmp)%bn.order()) == zkp.M_t);



    //C_bidS
    G1 TMP = bn.mult(zkp.C_hat_bidS, c) + bn.mult(g_hat, zkp.z_bidS) + bn.mult(h_hat_p, zkp.z_r_bidS);
    BOOL res6 = (TMP == zkp.C_hat_bidS_t);

    //C_y
    TMP = bn.mult(zkp.C_hat_y, c) + bn.mult(g_hat, zkp.z_y) + bn.mult(h_hat_p, zkp.z_r_y);
    BOOL res7 = TMP == zkp.C_hat_y_t;

    //C_yS
    BOOL res8 = (bn.mult(zkp.C_yS, c) + bn.mult(g_hat_p, zkp.z_yS) + bn.mult(h_hat_p, zkp.z_r_yS) == zkp.C_yS_t);

    //M2
    BOOL res9 = (bn.mult(zkp.C_yS, (-zkp.z_u)%bn.order()) + bn.mult(g_hat_p, zkp.z_mS) + bn.mult(h_hat_p, zkp.z_tmpS) == zkp.M2_t);
    //M_bidS
    TMP = bn.mult(zkp.C_hat_bidS, (-zkp.z_mS)%bn.order()) + bn.mult(g_hat, zkp.z_tmp_bidS) + bn.mult(h_hat_p, zkp.z_m_bidS);
    BOOL res10 = (TMP == zkp.M_bid_t);
    //N
    BOOL res11 = (bn.power(bn.pairing(C[1], h_hat_p), -zkp.z_m_bidS)*bn.power(bn.pairing(C[1], acc.rs_hat[1]+zkp.C_hat_bidS), zkp.z_mS)*bn.power(bn.pairing(C_pre[1], h_hat_p), zkp.z_r_y)*bn.power(bn.pairing(C_pre[1], zkp.C_hat_y), c) == zkp.N_t);
    //D_a
    BOOL res12 = (bn.mult(g_hat_p, zkp.nonmember.z_a) + bn.mult(h_hat_p, zkp.nonmember.z_r_a) + bn.mult(zkp.nonmember.D_a, c) == zkp.nonmember.D_a_t);
    //D_b
    BOOL res13 = (bn.mult(g_hat_p, zkp.nonmember.z_b) + bn.mult(h_hat_p, zkp.nonmember.z_r_b) + bn.mult(zkp.nonmember.D_b, c) == zkp.nonmember.D_b_t);
    //acc_M_t
    BOOL res14 = (bn.mult(g_hat_p, zkp.nonmember.z_m_acc) + bn.mult(h_hat_p, zkp.nonmember.z_tmp_acc) + bn.mult(zkp.nonmember.D_b, -zkp.z_y) == zkp.nonmember.M_t);
    //E
    BOOL res15 = (zkp.nonmember.E_t == bn.power(bn.pairing(C[1], zkp.nonmember.C_hat_h0), c) * bn.power(bn.pairing(zkp.nonmember.C_h, acc.V_hat), -zkp.z_y) * bn.power(bn.pairing(h_tilde, acc.V_hat), zkp.nonmember.z_m_acc) * bn.power(bn.pairing(C[1],h_hat_p), zkp.nonmember.z_a) * bn.power(acc.gT, zkp.z_y));
    //F_old_t
    BOOL res16 = com.pok_verify(C[0], zkp.F_old_t, zkp.com_old_alpha_zs, c);

#ifdef LOG
    cout << "ZKP2=========" << endl;
    cout << "  F_t: " << res1 << endl;
    cout << "  B_t: " << res2 << endl;
    cout << "  D_t: " << res3 << endl;
    cout << "  E_t: " << res4 << endl;
    cout << "  M_t: " << res5 << endl;
    cout << "  C_bidS_t: " << res6 << endl;
    cout << "  C_y_t: " << res7 << endl;
    cout << "  C_yS_t: " << res8 << endl;
    cout << "  M2_t: " << res9 << endl;
    cout << "  M_bidS_t: " << res10 << endl;
    cout << "  N_t: " << res11 << endl;
    cout << "  D_a_t: " << res12 << endl;
    cout << "  D_b_t: " << res13 << endl;
    cout << "  acc_M_t: " << res14 << endl;
    cout << "  E_t: " << res15 << endl;
    cout << "  F_old_t: " << res16 << endl;
#endif
    return res1 && res2 && res3 && res4 && res5 && res6 && res7 && res8 && res9 && res10 && res11 && res12 && res13 && res14 && res15 && res16;
}

BOOL PPBIssuer::verify_ZKP3(ZKP3& zkp, G2& C_pre1, G2& C_1S, Big& bid)
{
    Big c = H1(zkp.C_y_t);

    //C_y_t
    BOOL res1 = (zkp.C_y_t == bn.mult(zkp.C_y, c) + bn.mult(g_hat_p, zkp.z_y) + bn.mult(h_hat_p, zkp.z_r_y));

    //C_yS_t
    BOOL res2 = (zkp.C_yS_t == bn.mult(zkp.C_yS, c) + bn.mult(g_hat_p, zkp.z_yS) + bn.mult(h_hat_p, zkp.z_r_yS));

    //R_t
    BOOL res3 = (zkp.R_t == (bn.power(bn.pairing(C_1S, acc.rs_hat[1] + bn.mult(g_hat, bid)), zkp.z_y) * bn.power(bn.pairing(C_pre1, g_hat), -zkp.z_yS)));


    //ZKP Verification
#ifdef LOG
    cout << "ZKP3=========" << endl;
    cout << "  C_y_t: " << res1 << endl;
    cout << "  C_yS_t: " << res2 << endl;
    cout << "  R_t: " << res3 << endl;
#endif
    return res1 && res2 && res3;
}
