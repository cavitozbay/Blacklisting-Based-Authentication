#include "PPBHolder.hpp"

using namespace std;

PPBHolder::PPBHolder(PFC& _bn, PPBIssuer& issuer) : 
bn(_bn), com(issuer.com), acc(issuer.acc), spseq(issuer.spseq), g(issuer.g), 
g_hat(issuer.g_hat), g_tilde(issuer.g_tilde), h_tilde(issuer.h_tilde),
g_hat_p(issuer.g_hat_p), h_hat_p(issuer.h_hat_p), tlp(issuer.tlp)
{
    //cout << "g I: " << &g_tilde << endl;

    //com.init_gens(gens);
    usk = rand(bn.order());
    upk = bn.mult(issuer.w,usk);
    /*g = _g;
    g_hat = _g_hat;
    g_tilde = _g_tilde;
    h_tilde = _h_tilde;*/
    tmpf0.addterm(ZZn(1),0);
}

PPBHolder::~PPBHolder()
{
}

vector<G2> PPBHolder::receive_cred1(ZKP1& zkp)
{
    //Big esk_userS, dsrnd0S, dsrnd1S, zS, tS, yS;
    vector<Big> M;
    vector<G2> C_pre;

    cred.esk = rand(bn.order());
    cred.dsrnd0 = rand(bn.order());
    cred.dsrnd1 = rand(bn.order());
    cred.z = rand(bn.order());
    cred.t = rand(bn.order());
    cred.u = rand(bn.order());
    cred.y = rand(bn.order());

    M.push_back(usk);
    M.push_back(cred.esk);
    M.push_back(Big(0));
    M.push_back(cred.z);
    M.push_back(cred.t);
    
    cred.cmop = com.commit(M);

    C_pre.push_back(bn.mult(cred.cmop.C, cred.u));
    G2 tmpG = bn.mult(g, cred.u);
    C_pre.push_back(bn.mult(tmpG,cred.y));
    C_pre.push_back(tmpG);

    //ZKP

    prepare_ZKP1(zkp, cred.cmop, C_pre[2]);
    
    return C_pre;
}

BOOL PPBHolder::receive_cred2(vector<G2>& C_pre, Big& v, Sigma& _sigma, Big& esk_issuerS, G2& w)
{
    //G2 CS[3];
    G2 tmp1 = bn.mult(com.gens[1], modmult(cred.u,esk_issuerS,bn.order()));
    G2 tmp2 = bn.mult(com.gens[2], modmult(cred.u,v,bn.order()));
    cred.C[0] = C_pre[0] + tmp1 + tmp2;
    cred.C[1] = C_pre[1];
    cred.C[2] = C_pre[2];
    cred.cmop.C = cred.cmop.C + bn.mult(com.gens[2], v) + bn.mult(com.gens[1], esk_issuerS);
    cred.cmop.m[2] = v;
    cred.cmop.m[1] = (cred.cmop.m[1] + esk_issuerS)%bn.order();

    cred.esk = (cred.esk + esk_issuerS)%bn.order();
    cred.dsid = bn.mult(w, cred.esk);

    BOOL res = spseq.verify(cred.C, _sigma);
    //cout << "Credential Issuance: " << res << endl;
    
    Big u_inv = inverse(cred.u, bn.order());
    cred.sigma = spseq.chgrep(_sigma, cred.C, u_inv);

    return res;
}

void PPBHolder::borrow1(Big& vS, Big& gamma, Big& c0, Big& c1, Witness& wit, vector<G2>& C_pre, G2* _C, Sigma& _sigma, Big& bidS, ZKP2& zkp)
{

    //Big esk_userS, dsrnd0S, dsrnd1S, zS, tS, yS;
    vector<Big> M;
    
    Credential credS;

    credS.esk = rand(bn.order());
    credS.dsrnd0 = rand(bn.order());
    credS.dsrnd1 = rand(bn.order());
    credS.z = rand(bn.order());
    credS.t = rand(bn.order());
    credS.u = rand(bn.order());
    credS.y = rand(bn.order());
    credS.v = vS;
    //bidS = rand(bn.order());
    uS = credS.u;

    gamma = rand(bn.order());

    wit = acc.create_nonmembership(cred.BIDS_uZZn);
    modulo(*bn.mod);

    //cred.BIDS_u.push_back(bidS);
    //cred.BIDS_uZZn.push_back(ZZn(bidS));
    new_bid(cred.BIDS_u, cred.BIDS_uZZn);

    bidS = cred.BIDS_u.back();

    Poly f0 = build_from_roots(cred.BIDS_uZZn, bn.order());

    modulo(*bn.mod);

    G2 g_f0 = mult_poly(bn, acc.rs, f0);

    M.push_back(usk);
    M.push_back(credS.esk);
    M.push_back(credS.v);
    M.push_back(credS.z);
    M.push_back(credS.t);

    credS.cmop = com.commit(M);

    C_pre.push_back(bn.mult(credS.cmop.C, uS));

    C_pre.push_back(bn.mult(g_f0, modmult(credS.y, uS, bn.order())));

    C_pre.push_back(bn.mult(g, uS));

    c0 = (modmult(usk, gamma, bn.order()) + cred.dsrnd0)%bn.order(); 
    c1 = (modmult(cred.esk, gamma, bn.order()) + cred.dsrnd1)%bn.order(); 


    /////////////////
    //TLP
    /////////////////
    

    for (size_t i = 0; i < 3; i++)
    {
        _C[i] = cred.C[i];
    }
    _sigma = cred.sigma;
    //cout << "in borrow: " << spseq.verify(_C, cred.sigma) << endl;
    /////////////////
    //TODO: Put ZKP here!!!
    /////////////////
    
    prepare_ZKP2(zkp, credS, C_pre, bidS, wit);

    credS.BIDS_u = cred.BIDS_u;
    credS.BIDS_uZZn = cred.BIDS_uZZn;

    cred = credS;
    /*cred.esk = esk_userS;
    cred.dsrnd0 = dsrnd0S;
    cred.dsrnd1 = dsrnd1S;
    cred.z = zS;
    cred.t = tS;
    cred.y = yS;
    cred.v = vS;
    cred.u = uS;*/

    
    //cout << "in borrow: " << spseq.verify(_C, _sigma) << endl;

}

void PPBHolder::precompute_h0s(size_t d)
{
    pcomp_d = d;
    pcomp_rem = d;
    pcomp_h0s = vector<G1>();
    pcomp_bids = vector<Big>();
    pcomp_bids_ZZn = vector<ZZn>();
    cred.BIDS_u = vector<Big>();
    cred.BIDS_uZZn = vector<ZZn>();

    for (size_t i = 0; i < d; i++)
    {
        new_bid(pcomp_bids, pcomp_bids_ZZn);
    }

    modulo(bn.order());
    Poly f0 = build_from_roots(pcomp_bids_ZZn, bn.order());
    
    Poly result[3];
    egcd(result, acc.f, f0);
    
    
    ZZn divisor = result[0].coeff(0);
    result[0] /= divisor;
    result[1] /= divisor;  
    result[2] /= divisor;
    pcomp_h = result[1];
    //result[2] *= ZZn(-1);
    //cout << "f: " << f << endl;
    /*cout << "f0: " << f0 << endl;
    cout << "gcd: " << result[0] << endl;
    cout << "h: " << result[1] << endl;
    cout << "h0: " << result[2] << endl;*/
    Variable x;

    for (size_t i = 0; i < d+1; i++, result[2] = result[2] * x)
    {
        pcomp_h0s.push_back(mult_poly_hat(bn, acc.rs_hat,result[2]));
    }
    
    //G1 w1_hat = mult_poly_hat(bn, acc.rs_hat,result[2]);
    modulo(bn.order());
    //cout << "  GCD check: " << f*result[1] + f0*result[2] << endl;

    modulo(*bn.mod);
    pcomp_w2 = mult_poly(bn, acc.rs, pcomp_h);
    
}

void PPBHolder::pcomp_borrow1(Big& vS, Big& gamma, Big& c0, Big& c1, Witness& wit, vector<G2>& C_pre, G2* _C, Sigma& _sigma, Big& bidS, ZKP2& zkp)
{   
    if(pcomp_rem == 0)
    {
        borrow1(vS, gamma, c0, c1, wit, C_pre, _C, _sigma, bidS, zkp);
        return;
    }

    //Big esk_userS, dsrnd0S, dsrnd1S, zS, tS, yS;
    vector<Big> M;
    
    Credential credS;

    credS.esk = rand(bn.order());
    credS.dsrnd0 = rand(bn.order());
    credS.dsrnd1 = rand(bn.order());
    credS.z = rand(bn.order());
    credS.t = rand(bn.order());
    credS.u = rand(bn.order());
    credS.y = rand(bn.order());
    credS.v = vS;
    //bidS = rand(bn.order());
    uS = credS.u;

    gamma = rand(bn.order());

    //cred.BIDS_u.push_back(bidS);
    //cred.BIDS_uZZn.push_back(ZZn(bidS));

    //wit = acc.create_nonmembership(cred.BIDS_uZZn);
    modulo(bn.order());
    Poly exp_h0 = build_from_roots(pcomp_bids_ZZn, bn.order());
    modulo(*bn.mod);

    Poly ypoly = build_from_roots(cred.BIDS_uZZn, bn.order());
    //wit = acc.create_nonmembership(cred.BIDS_uZZn);
    modulo(*bn.mod);

    wit.w2 = pcomp_w2;
    wit.w1_hat = mult_poly_hat(bn, pcomp_h0s, exp_h0);
    G2 Y = mult_poly(bn, acc.rs, ypoly);

    acc.verify_nonmembership(Y, wit);

    bidS = pcomp_bids[pcomp_rem-1];
    pcomp_rem--;

    remove_bid(bidS, pcomp_bids, pcomp_bids_ZZn);
    add_bid(bidS, cred.BIDS_u, cred.BIDS_uZZn);

    modulo(bn.order());
    Poly f0 = build_from_roots(cred.BIDS_uZZn, bn.order());
    modulo(*bn.mod);
    G2 g_f0 = mult_poly(bn, acc.rs, f0);

    M.push_back(usk);
    M.push_back(credS.esk);
    M.push_back(credS.v);
    M.push_back(credS.z);
    M.push_back(credS.t);

    credS.cmop = com.commit(M);

    C_pre.push_back(bn.mult(credS.cmop.C, uS));

    C_pre.push_back(bn.mult(g_f0, modmult(credS.y, uS, bn.order())));

    C_pre.push_back(bn.mult(g, uS));

    c0 = (modmult(usk, gamma, bn.order()) + cred.dsrnd0)%bn.order(); 
    c1 = (modmult(cred.esk, gamma, bn.order()) + cred.dsrnd1)%bn.order(); 


    /////////////////
    //TLP
    /////////////////
    

    for (size_t i = 0; i < 3; i++)
    {
        _C[i] = cred.C[i];
    }
    _sigma = cred.sigma;
    //cout << "in borrow: " << spseq.verify(_C, cred.sigma) << endl;
    /////////////////
    //TODO: Put ZKP here!!!
    /////////////////
    
    prepare_ZKP2(zkp, credS, C_pre, bidS, wit);

    credS.BIDS_u = cred.BIDS_u;
    credS.BIDS_uZZn = cred.BIDS_uZZn;

    cred = credS;
    /*cred.esk = esk_userS;
    cred.dsrnd0 = dsrnd0S;
    cred.dsrnd1 = dsrnd1S;
    cred.z = zS;
    cred.t = tS;
    cred.y = yS;
    cred.v = vS;
    cred.u = uS;*/

    
    //cout << "in borrow: " << spseq.verify(_C, _sigma) << endl;

}


BOOL PPBHolder::borrow2(vector<G2>& C_pre, Sigma& _sigma_pr, Big& esk_issuerS, G2& w)
{
    G2 tmp1 = bn.mult(com.gens[1], modmult(cred.u,esk_issuerS,bn.order()));
    cred.C[0] = C_pre[0] + tmp1;
    cred.C[1] = C_pre[1];
    cred.C[2] = C_pre[2];

    cred.esk = (cred.esk + esk_issuerS)%bn.order();
    cred.dsid = bn.mult(w, cred.esk);
    
    cred.cmop.C = cred.cmop.C + bn.mult(com.gens[1], esk_issuerS);
    cred.cmop.m[1] = cred.esk;

    BOOL res = spseq.verify(cred.C, _sigma_pr);
#ifdef LOG
    cout << "BorrowinG1 SPSEQ Verification: " << res << endl;
#endif

    Big u_inv = inverse(cred.u, bn.order());
    cred.sigma = spseq.chgrep(_sigma_pr, cred.C, u_inv);
#ifdef LOG
    cout << "Borrowing chgrep: " << spseq.verify(cred.C, cred.sigma) << endl;
#endif
    return res;
}

void PPBHolder::repay1(Big& bid, G2& C_1S, G2* C_pr, Sigma& sigma_pr, Big& s, ZKP3& zkp)
{
    Big yS;

    s = rand(bn.order());
    yS = rand(bn.order());

    remove_bid(bid, cred.BIDS_u, cred.BIDS_uZZn);

    Poly f0 = build_from_roots(cred.BIDS_uZZn, bn.order());
    //cout << "qwer" <<degree(f0) << endl;
    //cout << ACC_SIZE << endl;
    C_1S = bn.mult(mult_poly(bn, acc.rs, f0), modmult(yS, s, bn.order()));

    for (size_t i = 0; i < 3; i++)
    {
        C_pr[i] = cred.C[i];
    }
    sigma_pr = spseq.chgrep(cred.sigma, C_pr, s);

    prepare_ZKP3(zkp, yS, s, bid, C_pr[1], C_1S);

    //cout << "Repayment: " << spseq.verify(C_pr, sigma_pr) << endl;

    cred.y = yS;

}

BOOL PPBHolder::repay2(Sigma& sigma_prD, G2& C_1S, Big& s, Big& v_p)
{

    cred.C[0] = bn.mult(cred.C[0], s) + bn.mult(com.gens[2], modmult(s, v_p, bn.order()));
    cred.C[1] = C_1S;
    cred.C[2] = bn.mult(g, s);
    BOOL res = spseq.verify(cred.C, sigma_prD);
#ifdef LOG
    cout << "Repayment2 SPSEQ verification: " << res << endl;
#endif
    Big s_inv = inverse(s, bn.order());
    cred.sigma = spseq.chgrep(sigma_prD, cred.C, s_inv);

    cred.v = (cred.v + v_p)%bn.order();
    cred.cmop.C = cred.cmop.C + bn.mult(com.gens[2], v_p);
    cred.cmop.m[2] = cred.v;

    return res;
}

void PPBHolder::prepare_ZKP1(ZKP1& zkp, CmOp& cmop, G2& G)
{
    vector<Big> alpha_ts;
    //vector<G2> u_ts;
    Big t_u, t_y, t_mult, t_tmp, t_gamma, t_delta, mult, tmp, gamma, delta;
    
    gamma = rand(bn.order());
    delta = rand(bn.order());
    mult = modmult(cred.u, gamma, bn.order());
    tmp = modmult(cred.u, delta, bn.order());

    zkp.E = bn.mult(g_hat_p, gamma) + bn.mult(h_hat_p, delta);
    zkp.F = cmop.C + bn.mult(h_tilde, gamma);

    t_u = rand(bn.order());
    t_y = rand(bn.order());
    t_mult = rand(bn.order());
    t_tmp = rand(bn.order());
    t_gamma = rand(bn.order());
    t_delta = rand(bn.order());

    //u_ts.push_back(G2());
    com.pok_init(alpha_ts, zkp.F_t);    
    //alpha_ts.push_back(alpha_tu);
    //alpha_ts.push_back(alpha_ty);
    //F_t
    zkp.F_t = zkp.F_t + bn.mult(h_tilde, t_gamma);
    //B_t
    zkp.B_t = bn.mult(g, t_u);
    //D_t
    zkp.D_t = bn.mult(zkp.F, t_u) + bn.mult(h_tilde, -t_mult);
    //E_t
    zkp.E_t = bn.mult(g_hat_p, t_gamma) + bn.mult(h_hat_p, t_delta);
    //G_t
    zkp.G_t = bn.mult(G, t_y);
    //M_t
    zkp.M_t = bn.mult(zkp.E, t_u) + bn.mult(g_hat_p, -t_mult) + bn.mult(h_hat_p, -t_tmp);

    //Challenge
    Big c = H2(zkp.F_t);

    //Response
    com.pok_response(cmop.m, alpha_ts, c, zkp.com_alpha_zs);

    

    zkp.z_u = (t_u - modmult(cred.u, c, bn.order()))%bn.order();
    zkp.z_y = (t_y - modmult(cred.y, c, bn.order()))%bn.order();
    zkp.z_mult = (t_mult - modmult(mult, c, bn.order()))%bn.order();
    zkp.z_tmp = (t_tmp - modmult(tmp, c, bn.order()))%bn.order();
    zkp.z_gamma = (t_gamma - modmult(gamma, c, bn.order()))%bn.order();
    zkp.z_delta = (t_delta - modmult(delta, c, bn.order()))%bn.order();
}

void PPBHolder::prepare_ZKP2(ZKP2& zkp, Credential& credS, vector<G2> C_pre, Big& bidS, Witness& wit)
{
    vector<Big> alpha_ts, old_alpha_ts;
    //vector<G2> u_ts;
    Big t_u, t_mult, t_tmp, t_gamma, t_delta, mult, tmp, gamma, delta;
    
    gamma = rand(bn.order());
    delta = rand(bn.order());
    mult = modmult(credS.u, gamma, bn.order());
    tmp = modmult(credS.u, delta, bn.order());

    zkp.E = bn.mult(g_hat_p, gamma) + bn.mult(h_hat_p, delta);
    zkp.F = credS.cmop.C + bn.mult(h_tilde, gamma);

    t_u = rand(bn.order());
    t_mult = rand(bn.order());
    t_tmp = rand(bn.order());
    t_gamma = rand(bn.order());
    t_delta = rand(bn.order());

    //u_ts.push_back(G2());
    com.pok_init(alpha_ts, zkp.F_t);    
    com.pok_init(old_alpha_ts, zkp.F_old_t);    
    //alpha_ts.push_back(alpha_tu);
    //alpha_ts.push_back(alpha_ty);
    //F_t
    zkp.F_t = zkp.F_t + bn.mult(h_tilde, t_gamma);
    //B_t
    zkp.B_t = bn.mult(g, t_u);
    //D_t
    zkp.D_t = bn.mult(zkp.F, t_u) + bn.mult(h_tilde, -t_mult);
    //E_t
    zkp.E_t = bn.mult(g_hat_p, t_gamma) + bn.mult(h_hat_p, t_delta);
    //G_t
    //zkp.G_t = bn.mult(G, t_y);
    //M_t
    zkp.M_t = bn.mult(zkp.E, t_u) + bn.mult(g_hat_p, -t_mult) + bn.mult(h_hat_p, -t_tmp);


    Big r_bidS, r_yS, r_y, mS, tmpS, t_bidS, t_r_bidS, t_r_yS, t_r_y, t_yS, t_y, t_mS, t_tmpS, t_m_bidS, m_bidS, tmp_bidS, t_tmp_bidS;
    r_bidS = rand(bn.order());
    r_yS = rand(bn.order());
    r_y = rand(bn.order());
    t_bidS = rand(bn.order());
    t_r_bidS = rand(bn.order());
    t_r_yS = rand(bn.order());
    t_r_y = rand(bn.order());
    t_yS = rand(bn.order());
    t_y = rand(bn.order());
    t_mS = rand(bn.order());
    t_m_bidS = rand(bn.order());
    t_tmpS = rand(bn.order());
    t_tmp_bidS = rand(bn.order());
    
    mS = modmult(credS.u, credS.y, bn.order());
    m_bidS = modmult(mS, r_bidS, bn.order());
    tmpS = modmult(credS.u, r_yS, bn.order());
    tmp_bidS = modmult(mS, bidS, bn.order());

    zkp.C_hat_bidS = bn.mult(g_hat, bidS) + bn.mult(h_hat_p, r_bidS);
    zkp.C_yS = bn.mult(g_hat_p, credS.y) + bn.mult(h_hat_p, r_yS);
    zkp.C_hat_y = bn.mult(g_hat, cred.y) + bn.mult(h_hat_p, r_y);

    zkp.C_hat_bidS_t = bn.mult(g_hat, t_bidS) + bn.mult(h_hat_p, t_r_bidS);
    zkp.C_yS_t = bn.mult(g_hat_p, t_yS) + bn.mult(h_hat_p, t_r_yS);
    zkp.C_hat_y_t = bn.mult(g_hat, t_y) + bn.mult(h_hat_p, t_r_y);
    zkp.M2_t = bn.mult(zkp.C_yS, -t_u) + bn.mult(g_hat_p, t_mS) + bn.mult(h_hat_p, t_tmpS);
    zkp.M_bid_t = bn.mult(zkp.C_hat_bidS, -t_mS) + bn.mult(g_hat, t_tmp_bidS) + bn.mult(h_hat_p, t_m_bidS);
    zkp.N_t = bn.power(bn.pairing(cred.C[1], h_hat_p), -t_m_bidS)*bn.power(bn.pairing(cred.C[1], acc.rs_hat[1]+zkp.C_hat_bidS), t_mS)*bn.power(bn.pairing(C_pre[1], h_hat_p), t_r_y);


    Big r_a, t_r_a, z_r_a, a, t_a, z_a, r_b, t_r_b, z_r_b, b, t_b, z_b, m_acc, t_m_acc, z_m_acc, tmp_acc, t_tmp_acc, z_tmp_acc;
    
    r_a = rand(bn.order());
    t_r_a = rand(bn.order());
    a = rand(bn.order());
    t_a = rand(bn.order());
    r_b = rand(bn.order());
    t_r_b = rand(bn.order());
    b = rand(bn.order());
    t_b = rand(bn.order());
    m_acc = modmult(b, cred.y, bn.order());
    tmp_acc = modmult(r_b, cred.y, bn.order());
    t_m_acc = rand(bn.order());
    t_tmp_acc = rand(bn.order());

    zkp.nonmember.C_hat_h0 = wit.w1_hat + bn.mult(h_hat_p, a);
    zkp.nonmember.C_h = wit.w2 + bn.mult(h_tilde, b);
    zkp.nonmember.D_a = bn.mult(g_hat_p, a) + bn.mult(h_hat_p, r_a);
    zkp.nonmember.D_b = bn.mult(g_hat_p, b) + bn.mult(h_hat_p, r_b);
    
    zkp.nonmember.M_t = bn.mult(zkp.nonmember.D_b, -t_y) + bn.mult(g_hat_p, t_m_acc) + bn.mult(h_hat_p, t_tmp_acc);
    zkp.nonmember.D_a_t = bn.mult(g_hat_p, t_a) + bn.mult(h_hat_p, t_r_a);
    zkp.nonmember.D_b_t = bn.mult(g_hat_p, t_b) + bn.mult(h_hat_p, t_r_b);
    zkp.nonmember.E_t = bn.power(bn.pairing(zkp.nonmember.C_h, acc.V_hat), -t_y) * bn.power(bn.pairing(h_tilde, acc.V_hat), t_m_acc) * bn.power(bn.pairing(cred.C[1], h_hat_p), t_a) * bn.power(acc.gT, t_y);
    G2 tmpG2 = bn.mult(cred.C[1], inverse(cred.y, bn.order()));

    //Challenge
    Big c = H2(zkp.F_t);

    //Response
    com.pok_response(credS.cmop.m, alpha_ts, c, zkp.com_alpha_zs);
    com.pok_response(cred.cmop.m, old_alpha_ts, c, zkp.com_old_alpha_zs);

    

    zkp.z_u = (t_u - modmult(credS.u, c, bn.order()))%bn.order();
    zkp.z_yS = (t_yS - modmult(credS.y, c, bn.order()))%bn.order();
    zkp.z_mult = (t_mult - modmult(mult, c, bn.order()))%bn.order();
    zkp.z_tmp = (t_tmp - modmult(tmp, c, bn.order()))%bn.order();
    zkp.z_gamma = (t_gamma - modmult(gamma, c, bn.order()))%bn.order();
    zkp.z_delta = (t_delta - modmult(delta, c, bn.order()))%bn.order(); 


    zkp.z_bidS = (t_bidS - modmult(bidS, c, bn.order()))%bn.order();
    zkp.z_r_bidS = (t_r_bidS - modmult(r_bidS, c, bn.order()))%bn.order();
    zkp.z_yS = (t_yS - modmult(credS.y, c, bn.order()))%bn.order();
    zkp.z_r_yS = (t_r_yS - modmult(r_yS, c, bn.order()))%bn.order();
    zkp.z_y = (t_y - modmult(cred.y, c, bn.order()))%bn.order();
    zkp.z_r_y = (t_r_y - modmult(r_y, c, bn.order()))%bn.order();
    zkp.z_mS = (t_mS - modmult(mS, c, bn.order()))%bn.order();
    zkp.z_m_bidS = (t_m_bidS - modmult(m_bidS, c, bn.order()))%bn.order();
    zkp.z_tmpS = (t_tmpS - modmult(tmpS, c, bn.order()))%bn.order();
    zkp.z_tmp_bidS = (t_tmp_bidS - modmult(tmp_bidS, c, bn.order()))%bn.order();

    zkp.nonmember.z_a = (t_a - modmult(a, c, bn.order()))%bn.order();
    zkp.nonmember.z_r_a = (t_r_a - modmult(r_a, c, bn.order()))%bn.order();
    zkp.nonmember.z_b = (t_b - modmult(b, c, bn.order()))%bn.order();
    zkp.nonmember.z_r_b = (t_r_b - modmult(r_b, c, bn.order()))%bn.order();
    zkp.nonmember.z_m_acc = (t_m_acc - modmult(m_acc, c, bn.order()))%bn.order();
    zkp.nonmember.z_tmp_acc = (t_tmp_acc - modmult(tmp_acc, c, bn.order()))%bn.order();

}

void PPBHolder::prepare_ZKP3(ZKP3& zkp, Big& yS, Big& s, Big& bid, G2& C_pr1, G2& C_1S)
{
    Big r_yS, r_y, r_yS_t, r_y_t, yS_t, y_t, m, m_t;

    
    r_y = rand(bn.order());
    r_yS = rand(bn.order());
    r_yS_t = rand(bn.order());
    r_y_t = rand(bn.order());
    y_t = rand(bn.order());
    yS_t = rand(bn.order());

    zkp.C_yS = bn.mult(g_hat_p, yS) + bn.mult(h_hat_p, r_yS);
    zkp.C_y = bn.mult(g_hat_p, cred.y) + bn.mult(h_hat_p, r_y);
    
    zkp.C_yS_t = bn.mult(g_hat_p, yS_t) + bn.mult(h_hat_p, r_yS_t);
    zkp.C_y_t = bn.mult(g_hat_p, y_t) + bn.mult(h_hat_p, r_y_t);
    zkp.R_t = bn.power(bn.pairing(C_1S, acc.rs_hat[1] + bn.mult(g_hat,bid)), y_t) * bn.power(bn.pairing(C_pr1, g_hat), -yS_t);

    //Challenge
    Big c = H1(zkp.C_y_t);

    //Response

    zkp.z_y = (y_t - modmult(cred.y, c, bn.order()))%bn.order();
    zkp.z_yS = (yS_t - modmult(yS, c, bn.order()))%bn.order();
    zkp.z_r_yS = (r_yS_t - modmult(r_yS, c, bn.order()))%bn.order();
    zkp.z_r_y = (r_y_t - modmult(r_y, c, bn.order()))%bn.order();

    //cout << "C_1S: " << C_1S.g << endl;
    //cout << "C_pr1: " << C_pr1.g << endl;
}

void PPBHolder::new_bid(vector<Big>& vecbig, vector<ZZn>& veczzn)
{
    modulo(bn.order());
    Big tmp = rand(bn.order());
    vecbig.push_back(tmp);
    veczzn.push_back(ZZn(tmp));
    modulo(*bn.mod);
}

void PPBHolder::add_bid(Big& bid, vector<Big>& vecbig, vector<ZZn>& veczzn)
{
    modulo(bn.order());
    vecbig.push_back(bid);
    veczzn.push_back(ZZn(bid));
    modulo(*bn.mod);
}

void PPBHolder::remove_bid(Big& bid, vector<Big>& vecbig, vector<ZZn>& veczzn)
{
    modulo(bn.order());
    vector<Big>::iterator it = find(vecbig.begin(), vecbig.end(), bid);
    if(it == vecbig.end()) return;
    vecbig.erase(it);

    vector<ZZn>::iterator itZZn = find(veczzn.begin(), veczzn.end(), ZZn(bid));
    if(itZZn == veczzn.end()) return;
    veczzn.erase(itZZn);
    modulo(*bn.mod);
}
