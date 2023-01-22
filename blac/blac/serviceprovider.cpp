#include "serviceprovider.hpp"

ServiceProvider::ServiceProvider(PFC& _bn): bn(_bn), acc(bn), aho(bn)
{
    acc.setup(ACC_SIZE);
    aho.key_gen(4);
    g_hat = acc.g_hat;
    bn.random(h_hat);
    bn.random(h);
    t = 1;
}

RegResp ServiceProvider::_register(InitRegResp resp)
{
    bn.start_hash();
    bn.add_to_hash(resp.CT0);
    bn.add_to_hash(resp.CT0_t);
    bn.add_to_hash(resp.Cx0);
    bn.add_to_hash(resp.Cx0_t);
    Big mdls = bn.order();
    Big c = bn.finish_hash_to_group();

    G2 res1R = bn.mult(resp.Cx0, c) + bn.mult(g_hat, resp.x_z) + bn.mult(h_hat, resp.rx0_z);
    BOOL res1 = (resp.Cx0_t == res1R);
    G2 res2R = bn.mult(resp.CT0, c) + bn.mult(g_hat, resp.T0_z) + bn.mult(h_hat, resp.rT0_z);
    BOOL res2 = (resp.CT0_t == res2R);

#ifdef LOG
    cout << "ZKP_REG==========" << endl;
    cout << (res1 && res2) << endl;
#endif

    G2 P0, CP0, CR0;

    Big R0, rR0;
    R0 = 1;
    rR0 = 1;

    CP0 = P0 + bn.mult(h_hat, R0);
    CR0 = bn.mult(g_hat, R0) + bn.mult(h_hat, rR0);

    RegResp respP;

    respP.M.push_back(resp.Cx0);
    respP.M.push_back(resp.CT0);
    respP.M.push_back(CP0);
    respP.M.push_back(CR0);


    respP.sigma = aho.sign(respP.M);

    return respP;
}

AuthResp ServiceProvider::auth(InitAuthResp resp, vector<size_t> V)
{
    G1 accL = acc.gen(V);

    bn.start_hash();
    bn.add_to_hash(resp.Cx0p);
    bn.add_to_hash(resp.Cx0p_t);
    bn.add_to_hash(resp.CT0p);
    bn.add_to_hash(resp.CT0p_t);
    bn.add_to_hash(resp.CR0p);
    bn.add_to_hash(resp.CR0p_t);
    bn.add_to_hash(resp.CR1);
    bn.add_to_hash(resp.CR1_t);
    bn.add_to_hash(resp.CT1);
    bn.add_to_hash(resp.CT1_t);
    bn.add_to_hash(resp.A_t);
    bn.add_to_hash(resp.B_t);
    bn.add_to_hash(resp.NM_t);
    Big c = bn.finish_hash_to_group();

    G2 res1R = bn.mult(resp.Cx0p, c) + bn.mult(g_hat, resp.x_z) + bn.mult(h_hat, resp.rx0pp_z);
    BOOL res1 = (resp.Cx0p_t == res1R);
    G2 res2R = bn.mult(resp.CT0p, c) + bn.mult(g_hat, resp.T0_z) + bn.mult(h_hat, resp.rT0pp_z);
    BOOL res2 = (resp.CT0p_t == res2R);
    G2 res3R = bn.mult(resp.CR0p, c) + bn.mult(g_hat, resp.R0_z) + bn.mult(h_hat, resp.rR0pp_z);
    BOOL res3 = (resp.CR0p_t == res3R);
    G2 res4R = bn.mult(resp.CR1, c) + bn.mult(g_hat, resp.R0_z + resp.R0p_z) + bn.mult(h_hat, resp.rR1_z);
    BOOL res4 = (resp.CR1_t == res4R);
    G2 res5R = bn.mult(resp.CT1, c) + bn.mult(g_hat, resp.T1_z) + bn.mult(h_hat, resp.rT1_z);
    BOOL res5 = (resp.CT1_t == res5R);

    GT A_C = bn.power(aho.A, -1) 
             * bn.pairing(resp.CtetaP1, aho.Gz)
             * bn.pairing(resp.CtetaP2, aho.Gr)
             * bn.pairing(resp.tetaP4, resp.tetaP3)
             * bn.pairing(resp.Cx0p, aho.Gs[0])
             * bn.pairing(resp.CT0p, aho.Gs[1])
             * bn.pairing(resp.CP0p, aho.Gs[2])
             * bn.pairing(resp.CR0p, aho.Gs[3]);
    
    GT Aright = bn.power(A_C, c)
              * bn.power(bn.pairing(h_hat, aho.Gz), resp.rtetaP1_z)
              * bn.power(bn.pairing(h_hat, aho.Gr), resp.rtetaP2_z)
              * bn.power(bn.pairing(h_hat, aho.Gs[0]), resp.rx0p_z)
              * bn.power(bn.pairing(h_hat, aho.Gs[1]), resp.rT0p_z)
              * bn.power(bn.pairing(h_hat, aho.Gs[2]), resp.R0p_z)
              * bn.power(bn.pairing(h_hat, aho.Gs[3]), resp.rR0p_z);

    BOOL res6 = (resp.A_t == Aright);

    GT B_C = bn.power(aho.B, -1) 
             * bn.pairing(resp.CtetaP1, aho.Hz)
             * bn.pairing(resp.CtetaP5, aho.Hr)
             * bn.pairing(resp.tetaP7, resp.tetaP6)
             * bn.pairing(resp.Cx0p, aho.Hs[0])
             * bn.pairing(resp.CT0p, aho.Hs[1])
             * bn.pairing(resp.CP0p, aho.Hs[2])
             * bn.pairing(resp.CR0p, aho.Hs[3]);

    GT Bright = bn.power(B_C, c)
              * bn.power(bn.pairing(h_hat, aho.Hz), resp.rtetaP1_z)
              * bn.power(bn.pairing(h_hat, aho.Hr), resp.rtetaP5_z)
              * bn.power(bn.pairing(h_hat, aho.Hs[0]), resp.rx0p_z)
              * bn.power(bn.pairing(h_hat, aho.Hs[1]), resp.rT0p_z)
              * bn.power(bn.pairing(h_hat, aho.Hs[2]), resp.R0p_z)
              * bn.power(bn.pairing(h_hat, aho.Hs[3]), resp.rR0p_z);

    BOOL res7 = (resp.B_t == Bright);

    GT NM_C = bn.pairing(resp.CP0p, accL)
              * bn.power(bn.pairing(acc.g_hat, resp.CW), -1);

    GT NMright = bn.power(NM_C, c)
               * bn.power(bn.pairing(h_hat, accL), resp.R0_z + resp.R0p_z)
               * bn.power(bn.pairing(acc.g_hat, h), -resp.rW_z);
    BOOL res8 = (resp.NM_t == NMright);

#ifdef LOG
    cout << "ZKP_AUTH==========" << endl;
    // cout << res1 << endl;
    // cout << res2 << endl;
    // cout << res3 << endl;
    // cout << res4 << endl;
    // cout << res5 << endl;
    // cout << res6 << endl;
    // cout << res7 << endl;
    // cout << res8 << endl;

    cout << (res1 && res2 && res3 && res4 && res5 && res6 && res7 && res8) << endl;
#endif


    G2 CP1 = resp.CP0p + acc.rs_hat[t];
    t++;

    G2 Cx1 = resp.Cx0p;

    vector<G2> M;
    M.push_back(Cx1);
    M.push_back(resp.CT1);
    M.push_back(CP1);
    M.push_back(resp.CR1);

    AuthResp aR;
    aR.sigma = aho.sign(M);
    aR.t = t-1;
    aR.M = M;
    return aR;
}
