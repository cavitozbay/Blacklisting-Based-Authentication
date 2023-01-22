#include <iostream>
#include <iomanip>
#include <fstream>
#include "poly.h"
#include "pairing_bn.h"
#include "acc.hpp"
#include "aho.hpp"
#include "serviceprovider.hpp"
#include "user.hpp"

using namespace std;

#define MIN_ITERS 50
#define MIN_TIME 10


void test_protocol()
{
        PFC bn(128);
    
    /*
    Acc acc(bn);
    acc.setup(100);

    vector<size_t> U;
    vector<size_t> V;

    for (size_t i = 1; i < 40; i++)
    {
        V.push_back(i);
    }
    U.push_back(41);
    U.push_back(42);

    G2 acc_v_hat = acc.gen(V);
    G2 wit_hat = acc.wit_gen(U, V);
    cout << "Acc Verify: " << acc.verify(U, acc_v_hat, wit_hat) << endl;


    AHO aho(bn);
    aho.key_gen(3);

    vector<G2> M;
    G2 tmp;
    for (size_t i = 0; i < 3; i++)
    {
        bn.random(tmp);
        M.push_back(tmp);
    }
    
    Sigma sigma = aho.sign(M);

    cout << "AHO Verify: " << aho.verify(sigma, M) << endl;*/

    ServiceProvider sp(bn);
    User user(sp);

    InitRegResp iRR = user.init_register();
    RegResp rR = sp._register(iRR);
    user.finish_register(rR);
    
    vector<size_t> V;

    for (size_t i = 20; i < 40; i++)
    {
        V.push_back(i);
    }



    InitAuthResp iAR = user.init_auth(V);
    AuthResp aR = sp.auth(iAR, V);
    user.finish_auth(aR);

    InitAuthResp iAR2 = user.init_auth(V);
    AuthResp aR2 = sp.auth(iAR2, V);
    user.finish_auth(aR2);


}

void bmark_register()
{
    PFC bn(128);
    
    ServiceProvider sp(bn);
    User user(sp);

    InitRegResp iRR = user.init_register();
    RegResp rR = sp._register(iRR);
    user.finish_register(rR);


    int iterations = 0;
    double elapsed;
    double total_ir = 0;
    double total_r = 0;
    double total_fr = 0;
    clock_t start = clock();
    do {
        User user(sp);
        start = clock();
        InitRegResp iRR = user.init_register();
        total_ir += clock() - start;

        start = clock();
        RegResp rR = sp._register(iRR);
        total_r += clock() - start;

        start = clock();
        user.finish_register(rR);
        total_fr += clock() - start;

        iterations++;
        elapsed=total_ir/(double)CLOCKS_PER_SEC;
    } while (elapsed<MIN_TIME || iterations<MIN_ITERS);

    elapsed=1000.0*elapsed/iterations;
    cout << 1000*(total_ir/(double)CLOCKS_PER_SEC)/iterations << endl;
    cout << 1000*(total_r/(double)CLOCKS_PER_SEC)/iterations << endl;
    cout << 1000*(total_fr/(double)CLOCKS_PER_SEC)/iterations << endl;
}

void bmark_auth()
{
    PFC bn(128);
    size_t blacklisted_cnt;
    for (blacklisted_cnt = 1000; blacklisted_cnt <= 4000; blacklisted_cnt += 1500)
    {
        ServiceProvider sp(bn);
        User user(sp);

        InitRegResp iRR = user.init_register();
        RegResp rR = sp._register(iRR);
        user.finish_register(rR);

        vector<size_t> V;
        Big smp_size = Big(ACC_SIZE);
        for (size_t i = 1; i <= blacklisted_cnt; i++)
        {
            V.push_back(i);
        }
        


        int iterations = 0;
        map<size_t, double> total_ia;
        map<size_t, double> total_a;
        map<size_t, double> total_fa;
        clock_t start;
        clock_t tmp = clock();
        for (size_t i = 0; i <= 200 ; i+=10)
        {
            total_a[i] = 0;
            total_ia[i] = 0;
            total_fa[i] = 0;
        }
        
        do {
            User userP = user;
            sp.t = blacklisted_cnt + 5;
            for (size_t i = 0; i <= 200; i+=10)
            {
                start = clock();
                InitAuthResp iAR = userP.init_auth(V);
                total_ia[i] += clock() - start;
                
                start = clock();
                AuthResp aR = sp.auth(iAR, V);
                total_a[i] += clock() - start;

                start = clock();
                userP.finish_auth(aR);
                total_fa[i] += clock() - start;
                for (size_t j = 0; j < 9; j++)
                {
                    InitAuthResp iAR = userP.init_auth(V);
                    AuthResp aR = sp.auth(iAR, V);
                    userP.finish_auth(aR);
                }
            }
            
            
            iterations++;
        } while (iterations<MIN_ITERS);
        for (size_t i = 0; i <= 200; i+=10)
        {
            cout << setw(20) << "IA - A - FA" << " - ;" << setw(5) << i << "; " << setw(6) << blacklisted_cnt << "; " << setw(8) << setprecision(2) << fixed << 1000.0*(double)(total_ia[i]/(double)CLOCKS_PER_SEC)/iterations << "; ";
            cout << setw(8) << setprecision(2) << fixed << 1000.0*(double)(total_a[i]/(double)CLOCKS_PER_SEC)/iterations << "; "; 
            cout << setw(8) << setprecision(2) << fixed << 1000.0*(double)(total_fa[i]/(double)CLOCKS_PER_SEC)/iterations << "; " << endl;
        }
    }
    
}

void bmark_comm()
{
        PFC bn(128);
    
    /*
    Acc acc(bn);
    acc.setup(100);

    vector<size_t> U;
    vector<size_t> V;

    for (size_t i = 1; i < 40; i++)
    {
        V.push_back(i);
    }
    U.push_back(41);
    U.push_back(42);

    G2 acc_v_hat = acc.gen(V);
    G2 wit_hat = acc.wit_gen(U, V);
    cout << "Acc Verify: " << acc.verify(U, acc_v_hat, wit_hat) << endl;


    AHO aho(bn);
    aho.key_gen(3);

    vector<G2> M;
    G2 tmp;
    for (size_t i = 0; i < 3; i++)
    {
        bn.random(tmp);
        M.push_back(tmp);
    }
    
    Sigma sigma = aho.sign(M);

    cout << "AHO Verify: " << aho.verify(sigma, M) << endl;*/

    ServiceProvider sp(bn);
    User user(sp);

    ofstream fr("register.bin");

    InitRegResp iRR = user.init_register();
    RegResp rR = sp._register(iRR);
    user.finish_register(rR);
    
    fr << iRR.CT0.g;
    fr << iRR.CT0_t.g;
    fr << iRR.Cx0.g;
    fr << iRR.Cx0_t.g;
    fr << iRR.rT0_z;
    fr << iRR.rx0_z;
    fr << iRR.T0_z;
    fr << iRR.x_z;
    for (size_t i = 0; i < 4; i++)
    {
        fr << rR.M[i].g;
    }
    for(auto teta: rR.sigma.tetas)
    {
        fr << teta.second.g;
    }
    for(auto teta: rR.sigma.tetas_hat)
    {
        fr << teta.second.g;
    }
    fr.close();

    vector<size_t> V;

    for (size_t i = 1; i < 40; i++)
    {
        V.push_back(i);
    }


    ofstream fa("authentication.bin");
    InitAuthResp iAR = user.init_auth(V);
    AuthResp aR = sp.auth(iAR, V);
    user.finish_auth(aR);

    fa << iAR.A_t.g;
    fa << iAR.B_t.g;
    fa << iAR.CP0p.g;
    fa << iAR.CR0p.g;
    fa << iAR.CR0p_t.g;
    fa << iAR.CR1.g;
    fa << iAR.CR1_t.g;
    fa << iAR.CT0p.g;
    fa << iAR.CT0p_t.g;
    fa << iAR.CT1.g;
    fa << iAR.CT1_t.g;
    fa << iAR.CtetaP1.g;
    fa << iAR.CtetaP2.g;
    fa << iAR.CtetaP5.g;
    fa << iAR.CW.g;
    fa << iAR.Cx0p.g;
    fa << iAR.Cx0p_t.g;
    fa << iAR.NM_t.g;
    fa << iAR.R0_z;
    fa << iAR.R0p_z;
    fa << iAR.rR0p_z;
    fa << iAR.rR0pp_z;
    fa << iAR.rR1_z;
    fa << iAR.rT0p_z;
    fa << iAR.rT1_z;
    fa << iAR.rT0pp_z;
    fa << iAR.rtetaP1_z;
    fa << iAR.rtetaP2_z;
    fa << iAR.rtetaP5_z;
    fa << iAR.rW_z;
    fa << iAR.rx0p_z;
    fa << iAR.rx0pp_z;
    fa << iAR.T0_z;
    fa << iAR.T1_z;
    fa << iAR.tetaP3.g;
    fa << iAR.tetaP4.g;
    fa << iAR.tetaP6.g;
    fa << iAR.tetaP7.g;
    fa << iAR.x_z;

    fa << aR.t;
        for (size_t i = 0; i < 4; i++)
    {
        fa << rR.M[i].g;
    }
    for(auto teta: rR.sigma.tetas)
    {
        fa << teta.second.g;
    }
    for(auto teta: rR.sigma.tetas_hat)
    {
        fa << teta.second.g;
    }
    fa.close();

}

int main(int, char**) {

    // bmark_register();
    // bmark_comm();
    // bmark_auth();
    test_protocol();
    return 0;
}

