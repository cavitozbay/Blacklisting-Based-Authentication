#include "user.hpp"

User::User(ServiceProvider& sp): bn(sp.bn), acc(sp.acc), aho(sp.aho)
{
    g_hat = acc.g_hat;
    h_hat = sp.h_hat;
    h = sp.h;

}

InitRegResp User::init_register()
{
    x = rand(bn.order());
    T0 = rand(bn.order());

    rx0 = rand(bn.order());
    rT0 = rand(bn.order());
    
    InitRegResp resp;

    resp.Cx0 = bn.mult(g_hat, x) + bn.mult(h_hat, rx0);
    resp.CT0 = bn.mult(g_hat, T0) + bn.mult(h_hat, rT0);

    Big x_t, T0_t, rx0_t, rT0_t;

    x_t = rand(bn.order());
    T0_t = rand(bn.order());
    rx0_t = rand(bn.order());
    rT0_t = rand(bn.order());

    resp.Cx0_t = bn.mult(g_hat, x_t) + bn.mult(h_hat, rx0_t);
    resp.CT0_t = bn.mult(g_hat, T0_t) + bn.mult(h_hat, rT0_t);

    bn.start_hash();
    bn.add_to_hash(resp.CT0);
    bn.add_to_hash(resp.CT0_t);
    bn.add_to_hash(resp.Cx0);
    bn.add_to_hash(resp.Cx0_t);
    Big mdls = bn.order();
    Big c = bn.finish_hash_to_group();

    resp.x_z = (x_t - modmult(x, c, bn.order()))%bn.order();
    resp.rx0_z = (rx0_t - modmult(rx0, c, bn.order()))%bn.order();
    resp.T0_z = (T0_t - modmult(T0, c, bn.order()))%bn.order();
    resp.rT0_z = (rT0_t - modmult(rT0, c, bn.order()))%bn.order();

    return resp;
}

void User::finish_register(RegResp resp)
{
    aho.verify(resp.sigma, resp.M);

    Cx0 = resp.M[0];
    CT0 = resp.M[1];
    CP0 = resp.M[2];
    CR0 = resp.M[3];
    P0 = G2();

    rR0 = 1;
    R0 = 1;
    sigma = resp.sigma;

}

InitAuthResp User::init_auth(vector<size_t> V)
{
    InitAuthResp resp;

    Big rx0p, rT0p, R0p, rR0p;
    rx0p = rand(bn.order());
    rT0p = rand(bn.order());
    R0p = rand(bn.order());
    rR0p = rand(bn.order());
    
    Big rx0pp = rx0 + rx0p;
    Big rT0pp = rT0 + rT0p;
    Big R0pp = R0 + R0p;
    Big rR0pp = rR0 + rR0p;

    resp.Cx0p = Cx0 + bn.mult(h_hat, rx0p);
    resp.CT0p = CT0 + bn.mult(h_hat, rT0p);
    resp.CP0p = CP0 + bn.mult(h_hat, R0p);
    resp.CR0p = CR0 + bn.mult(h_hat, rR0p);

    Sigma sigmaP = aho.randomize(sigma);

    // cout << (aho.verify(sigma,resp.M)) << endl;

    Big rtetaP1, rtetaP2, rtetaP5;

    rtetaP1 = rand(bn.order());
    rtetaP2 = rand(bn.order());
    rtetaP5 = rand(bn.order());

    resp.CtetaP1 = sigmaP.tetas_hat[1] + bn.mult(h_hat, rtetaP1);
    resp.CtetaP2 = sigmaP.tetas_hat[2] + bn.mult(h_hat, rtetaP2);
    resp.CtetaP5 = sigmaP.tetas_hat[5] + bn.mult(h_hat, rtetaP5);
    resp.tetaP3 = sigmaP.tetas[3];
    resp.tetaP4 = sigmaP.tetas_hat[4];
    resp.tetaP6 = sigmaP.tetas[6];
    resp.tetaP7 = sigmaP.tetas_hat[7];

    G1 W = acc.wit_gen(LU, V);
    Big rW = rand(bn.order());
    resp.CW = W + bn.mult(h, rW);

    ///////////
    Big R1 = R0pp;
    Big rR1 = rand(bn.order());
    Big T1 = rand(bn.order());
    Big rT1 = rand(bn.order());

    resp.CR1 = bn.mult(g_hat, R1) + bn.mult(h_hat, rR1);
    resp.CT1 = bn.mult(g_hat, T1) + bn.mult(h_hat, rT1);


    //////// ZKP ////////////

    Big x_t = rand(bn.order());
    Big rx0pp_t = rand(bn.order());
    Big T0_t = rand(bn.order());
    Big rT0pp_t = rand(bn.order());
    Big R0_t = rand(bn.order());
    Big rR0pp_t = rand(bn.order());
    Big R0p_t = rand(bn.order());
    Big rR1_t = rand(bn.order());
    Big T1_t = rand(bn.order());
    Big rT1_t = rand(bn.order());
    Big rtetaP1_t = rand(bn.order());
    Big rtetaP2_t = rand(bn.order());
    Big rx0p_t = rand(bn.order());
    Big rT0p_t = rand(bn.order());
    Big rR0p_t = rand(bn.order());
    Big rtetaP5_t = rand(bn.order());
    Big rW_t = rand(bn.order());

    resp.Cx0p_t = bn.mult(g_hat, x_t) + bn.mult(h_hat, rx0pp_t);
    resp.CT0p_t = bn.mult(g_hat, T0_t) + bn.mult(h_hat, rT0pp_t);
    resp.CR0p_t = bn.mult(g_hat, R0_t) + bn.mult(h_hat, rR0pp_t);
    resp.CR1_t = bn.mult(g_hat, R0_t + R0p_t) + bn.mult(h_hat, rR1_t);
    resp.CT1_t = bn.mult(g_hat, T1_t) + bn.mult(h_hat, rT1_t);
    resp.A_t = bn.power(bn.pairing(h_hat, aho.Gz), rtetaP1_t)
             * bn.power(bn.pairing(h_hat, aho.Gr), rtetaP2_t)
             * bn.power(bn.pairing(h_hat, aho.Gs[0]), rx0p_t)
             * bn.power(bn.pairing(h_hat, aho.Gs[1]), rT0p_t)
             * bn.power(bn.pairing(h_hat, aho.Gs[2]), R0p_t)
             * bn.power(bn.pairing(h_hat, aho.Gs[3]), rR0p_t);
    
    resp.B_t = bn.power(bn.pairing(h_hat, aho.Hz), rtetaP1_t)
             * bn.power(bn.pairing(h_hat, aho.Hr), rtetaP5_t)
             * bn.power(bn.pairing(h_hat, aho.Hs[0]), rx0p_t)
             * bn.power(bn.pairing(h_hat, aho.Hs[1]), rT0p_t)
             * bn.power(bn.pairing(h_hat, aho.Hs[2]), R0p_t)
             * bn.power(bn.pairing(h_hat, aho.Hs[3]), rR0p_t);
    
    G1 accL = acc.gen(V);
    resp.NM_t = bn.power(bn.pairing(h_hat, accL), R0_t + R0p_t)
              * bn.power(bn.pairing(acc.g_hat, h), -rW_t);

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

    resp.x_z = (x_t - modmult(x, c, bn.order()));
    resp.rx0pp_z = (rx0pp_t - modmult(rx0pp, c, bn.order()));
    resp.T0_z = (T0_t - modmult(T0, c, bn.order()));
    resp.rT0pp_z = (rT0pp_t - modmult(rT0pp, c, bn.order()));
    resp.R0_z = (R0_t - modmult(R0, c, bn.order()));
    resp.rR0pp_z = (rR0pp_t - modmult(rR0pp, c, bn.order()));
    resp.R0p_z = (R0p_t - modmult(R0p, c, bn.order()));
    resp.rR1_z = (rR1_t - modmult(rR1, c, bn.order()));
    resp.T1_z = (T1_t - modmult(T1, c, bn.order()));
    resp.rT1_z = (rT1_t - modmult(rT1, c, bn.order()));
    resp.rtetaP1_z = (rtetaP1_t - modmult(rtetaP1, c, bn.order()));
    resp.rtetaP2_z = (rtetaP2_t - modmult(rtetaP2, c, bn.order()));
    resp.rtetaP5_z = (rtetaP5_t - modmult(rtetaP5, c, bn.order()));
    resp.rx0p_z = (rx0p_t - modmult(rx0p, c, bn.order()));
    resp.rT0p_z = (rT0p_t - modmult(rT0p, c, bn.order()));
    resp.rR0p_z = (rR0p_t - modmult(rR0p, c, bn.order()));
    resp.rW_z = (rW_t - modmult(rW, c, bn.order()));


    rx0 = rx0pp;
    CT0 = resp.CT1;
    CR0 = resp.CR1;
    T0 = T1;
    rT0 = rT1;
    R0 = R1;
    rR0 = rR1;
    
    return resp;
}

void User::finish_auth(AuthResp resp)
{
    aho.verify(resp.sigma, resp.M);

    sigma = resp.sigma;
    LU.push_back(resp.t);
    Cx0 = resp.M[0];
    CP0 = resp.M[2];
}
