#include <iostream>
#include <iomanip>
#include <vector>
#include "QueueSig.hpp"
#include "Acc.hpp"
#include "helpers.hpp"
#include "Issuer.hpp"
#include "User.hpp"

#define MIN_ITER 20

void bmark_issuance()
{
    size_t K = 10;
    Issuer issuer(K);
    User user(issuer);
    clock_t start;
    clock_t elapsed_r1 = 0;
    clock_t elapsed_r = 0;
    clock_t elapsed_r2 = 0;
    int iterations = 0;

    do
    {
        start = clock();
        ReqSig req = user.init_register();
        elapsed_r1 += clock() - start;

        start = clock();
        RegResp respReg = issuer._register(req);
        elapsed_r += clock() - start;

        start = clock();
        user.finish_register(respReg);
        elapsed_r2 += clock() - start;

        iterations++;
    } while (iterations<MIN_ITER);
    
    cout << setw(20) << "R1 R R2" << " - " << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r1/(double)CLOCKS_PER_SEC)/iterations << ";" << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r/(double)CLOCKS_PER_SEC)/iterations << ";" << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r2/(double)CLOCKS_PER_SEC)/iterations << ";" << endl;
}

void bmark_auth()
{
    size_t K = 10;
    Issuer issuer(K);
    User user(issuer);

    ReqSig req = user.init_register();
    RegResp respReg = issuer._register(req);
    user.finish_register(respReg);


    clock_t start;
    clock_t elapsed_r1 = 0;
    clock_t elapsed_r = 0;
    clock_t elapsed_r2 = 0;
    int iterations = 0;

    do
    {
        start = clock();
        InitAuth iA = user.init_auth();
        elapsed_r1 += clock() - start;

        start = clock();
        AuthResp respAuth = issuer.auth(iA);
        elapsed_r += clock() - start;

        start = clock();
        user.finish_auth(respAuth);
        elapsed_r2 += clock() - start;

        iterations++;
    } while (iterations<MIN_ITER);
    
    cout << setw(20) << "A1 A A2" << " - " << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r1/(double)CLOCKS_PER_SEC)/iterations << ";" << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r/(double)CLOCKS_PER_SEC)/iterations << ";" << setw(8) << setprecision(2) << fixed << 1000.0*(elapsed_r2/(double)CLOCKS_PER_SEC)/iterations << ";" << endl;
}

void bmark_comm()
{
    miracl* mip = mirsys(256,0);
    mip->IOBASE = 16;
    size_t K = 10;
    Issuer issuer(K);
    User user(issuer);

    ofstream reg("reg.bin");
    ofstream regktms("regktms.bin");
    ReqSig req = user.init_register();
    RegResp respReg = issuer._register(req);
    user.finish_register(respReg);

    reg << req.C;
    reg << req.Ct;
    reg << req.rz;
    reg << respReg.sig.e;
    reg << respReg.sig.rp;
    reg << respReg.sig.v;
    reg << respReg.wit.a;
    reg << respReg.wit.d;
    regktms << req.zs[0];

    reg.close();
    regktms.close();

    ofstream auth("auth.bin");
    ofstream authktms("authktms.bin");

    InitAuth iA = user.init_auth();
    AuthResp respAuth = issuer.auth(iA);
    user.finish_auth(respAuth);

    auth << iA.kos.C;
    auth << iA.kos.Cct;
    auth << iA.kos.Cd;
    auth << iA.kos.Ce;
    auth << iA.kos.Cet;
    auth << iA.kos.Cs;
    auth << iA.kos.Cst;
    auth << iA.kos.Ct;
    auth << iA.kos.Cv;
    auth << iA.kos.Cw;
    auth << iA.kos.Cwt;
    auth << iA.kos.Cx;
    auth << iA.kos.Cxt;
    auth << iA.kos.Cz;
    auth << iA.kos.Czmult;
    auth << iA.kos.Czt;
    auth << iA.kos.rng1.dcap1;
    auth << iA.kos.rng1.dcap2;
    auth << iA.kos.rng1.dtld2;
    auth << iA.kos.rng1.dtld1;
    auth << iA.kos.rng1.Ecap1;
    auth << iA.kos.rng1.Ecap1t;
    auth << iA.kos.rng1.Ecap2;
    auth << iA.kos.rng1.Ecap;
    auth << iA.kos.rng1.Etld1;
    auth << iA.kos.rng1.Etld1t;
    auth << iA.kos.rng1.Etld2;
    auth << iA.kos.rng1.Etld;
    auth << iA.kos.rng1.Fcap;
    auth << iA.kos.rng1.Fcapt;
    auth << iA.kos.rng1.Ftld;
    auth << iA.kos.rng1.Ftldt;
    auth << iA.kos.rng1.Wcap;
    auth << iA.kos.rng1.Wtld;
    auth << iA.kos.rng1.zr1cap;
    auth << iA.kos.rng1.zr1tld;
    auth << iA.kos.rng1.zr2cap;
    auth << iA.kos.rng1.zr2tld;
    auth << iA.kos.rng1.zr3cap/(issuer.N*issuer.N);
    auth << (iA.kos.rng1.zr3cap%(issuer.N*issuer.N))/issuer.N;
    auth << iA.kos.rng1.zr3cap%issuer.N;
    auth << iA.kos.rng1.zr3tld/(issuer.N*issuer.N);
    auth << (iA.kos.rng1.zr3tld%(issuer.N*issuer.N))/issuer.N;
    auth << iA.kos.rng1.zr3tld%issuer.N;
    auth << iA.kos.rng1.zxcap1;
    auth << iA.kos.rng1.zxtld1;
    auth << iA.kos.rng2.dcap1;
    auth << iA.kos.rng2.dcap2;
    auth << iA.kos.rng2.dtld2;
    auth << iA.kos.rng2.dtld1;
    auth << iA.kos.rng2.Ecap1;
    auth << iA.kos.rng2.Ecap1t;
    auth << iA.kos.rng2.Ecap2;
    auth << iA.kos.rng2.Ecap;
    auth << iA.kos.rng2.Etld1;
    auth << iA.kos.rng2.Etld1t;
    auth << iA.kos.rng2.Etld2;
    auth << iA.kos.rng2.Etld;
    auth << iA.kos.rng2.Fcap;
    auth << iA.kos.rng2.Fcapt;
    auth << iA.kos.rng2.Ftld;
    auth << iA.kos.rng2.Ftldt;
    auth << iA.kos.rng2.Wcap;
    auth << iA.kos.rng2.Wtld;
    auth << iA.kos.rng2.zr1cap;
    auth << iA.kos.rng2.zr1tld;
    auth << iA.kos.rng2.zr2cap;
    auth << iA.kos.rng2.zr2tld;
    auth << iA.kos.rng2.zr3cap/(issuer.N*issuer.N);
    auth << (iA.kos.rng2.zr3cap%(issuer.N*issuer.N))/issuer.N;
    auth << iA.kos.rng2.zr3cap%issuer.N;
    auth << iA.kos.rng2.zr3tld/(issuer.N*issuer.N);
    auth << (iA.kos.rng2.zr3tld%(issuer.N*issuer.N))/issuer.N;
    auth << iA.kos.rng2.zr3tld%issuer.N;
    auth << iA.kos.rng2.zxcap1;
    auth << iA.kos.rng2.zxtld1;
    auth << iA.kos.ze;
    auth << iA.kos.zr;
    auth << iA.kos.zrc;
    auth << iA.kos.zre;
    auth << iA.kos.zrs;
    auth << iA.kos.zrw;
    auth << iA.kos.zrx;
    auth << iA.kos.zrz/issuer.N;
    auth << iA.kos.zrz%issuer.N;
    auth << iA.kos.zrzmul/issuer.N;
    auth << iA.kos.zrzmul%issuer.N;
    auth << iA.kos.zw;
    auth << iA.kos.zz;
    auth << iA.req2.C;
    auth << iA.req2.Ct;
    auth << iA.req2.rz;
    auth << iA.tK;

    auth << respAuth.sig.e;
    auth << respAuth.sig.v;
    auth << respAuth.sig.rp;
    auth << respAuth.wit.a;
    auth << respAuth.wit.d;

    authktms << iA.nmps[0].C1;
    authktms << iA.nmps[0].C1t;
    authktms << iA.nmps[0].Ca;
    authktms << iA.nmps[0].Cd;
    authktms << iA.nmps[0].Ce;
    authktms << iA.nmps[0].Cegt;
    authktms << iA.nmps[0].Cet;
    authktms << iA.nmps[0].Cw;
    authktms << iA.nmps[0].Cwt;
    authktms << iA.nmps[0].Cx;
    authktms << iA.nmps[0].Cxt;
    authktms << iA.nmps[0].Cz;
    authktms << iA.nmps[0].Czt;
    authktms << iA.nmps[0].rng2.dcap1;
    authktms << iA.nmps[0].rng2.dcap2;
    authktms << iA.nmps[0].rng2.dtld2;
    authktms << iA.nmps[0].rng2.dtld1;
    authktms << iA.nmps[0].rng2.Ecap1;
    authktms << iA.nmps[0].rng2.Ecap1t;
    authktms << iA.nmps[0].rng2.Ecap2;
    authktms << iA.nmps[0].rng2.Ecap;
    authktms << iA.nmps[0].rng2.Etld1;
    authktms << iA.nmps[0].rng2.Etld1t;
    authktms << iA.nmps[0].rng2.Etld2;
    authktms << iA.nmps[0].rng2.Etld;
    authktms << iA.nmps[0].rng2.Fcap;
    authktms << iA.nmps[0].rng2.Fcapt;
    authktms << iA.nmps[0].rng2.Ftld;
    authktms << iA.nmps[0].rng2.Ftldt;
    authktms << iA.nmps[0].rng2.Wcap;
    authktms << iA.nmps[0].rng2.Wtld;
    authktms << iA.nmps[0].rng2.zr1cap;
    authktms << iA.nmps[0].rng2.zr1tld;
    authktms << iA.nmps[0].rng2.zr2cap;
    authktms << iA.nmps[0].rng2.zr2tld;
    auth << iA.nmps[0].rng2.zr3cap/(issuer.N*issuer.N);
    auth << (iA.nmps[0].rng2.zr3cap%(issuer.N*issuer.N))/issuer.N;
    auth << iA.nmps[0].rng2.zr3cap%issuer.N;
    auth << iA.nmps[0].rng2.zr3tld/(issuer.N*issuer.N);
    auth << (iA.nmps[0].rng2.zr3tld%(issuer.N*issuer.N))/issuer.N;
    auth << iA.nmps[0].rng2.zr3tld%issuer.N;
    authktms << iA.nmps[0].rng2.zxcap1;
    authktms << iA.nmps[0].rng2.zxtld1;
    authktms << iA.nmps[0].rng1.dcap1;
    authktms << iA.nmps[0].rng1.dcap2;
    authktms << iA.nmps[0].rng1.dtld2;
    authktms << iA.nmps[0].rng1.dtld1;
    authktms << iA.nmps[0].rng1.Ecap1;
    authktms << iA.nmps[0].rng1.Ecap1t;
    authktms << iA.nmps[0].rng1.Ecap2;
    authktms << iA.nmps[0].rng1.Ecap;
    authktms << iA.nmps[0].rng1.Etld1;
    authktms << iA.nmps[0].rng1.Etld1t;
    authktms << iA.nmps[0].rng1.Etld2;
    authktms << iA.nmps[0].rng1.Etld;
    authktms << iA.nmps[0].rng1.Fcap;
    authktms << iA.nmps[0].rng1.Fcapt;
    authktms << iA.nmps[0].rng1.Ftld;
    authktms << iA.nmps[0].rng1.Ftldt;
    authktms << iA.nmps[0].rng1.Wcap;
    authktms << iA.nmps[0].rng1.Wtld;
    authktms << iA.nmps[0].rng1.zr1cap;
    authktms << iA.nmps[0].rng1.zr1tld;
    authktms << iA.nmps[0].rng1.zr2cap;
    authktms << iA.nmps[0].rng1.zr2tld;
    auth << iA.nmps[0].rng1.zr3cap/(issuer.N*issuer.N);
    auth << (iA.nmps[0].rng1.zr3cap%(issuer.N*issuer.N))/issuer.N;
    auth << iA.nmps[0].rng1.zr3cap%issuer.N;
    auth << iA.nmps[0].rng1.zr3tld/(issuer.N*issuer.N);
    auth << (iA.nmps[0].rng1.zr3tld%(issuer.N*issuer.N))/issuer.N;
    auth << iA.nmps[0].rng1.zr3tld%issuer.N;
    authktms << iA.nmps[0].rng1.zxcap1;
    authktms << iA.nmps[0].rng1.zxtld1;
    authktms << iA.nmps[0].za;
    authktms << iA.nmps[0].zr;
    authktms << iA.nmps[0].zre;
    authktms << iA.nmps[0].zrw;
    authktms << iA.nmps[0].zrx;
    authktms << iA.nmps[0].zrz/issuer.N;
    authktms << iA.nmps[0].zrz%issuer.N;
    authktms << iA.nmps[0].zt;
    authktms << iA.nmps[0].zw;
    authktms << iA.nmps[0].zz;
    authktms << iA.req2.zs[0];
    auth << iA.kos.zxz[0];

    auth.close();
    authktms.close();
}

void test_protocol()
{
    size_t K = 3;
    Issuer issuer(K);
    User user(issuer);

    ReqSig req = user.init_register();
    RegResp respReg = issuer._register(req);
    user.finish_register(respReg);

    InitAuth iA = user.init_auth();
    AuthResp respAuth = issuer.auth(iA);
    user.finish_auth(respAuth);

    iA = user.init_auth();
    respAuth = issuer.auth(iA);
    user.finish_auth(respAuth);

}

int main(int, char**) {
    miracl* mip = mirsys(256,0);
    mip->IOBASE = 16;

    // QueueSig qsig(2);
    // qsig.key_gen();
    // vector<Big> Q;
    // for (size_t i = 0; i < 3; i++)
    // {
    //     Q.push_back(rand(qsig.mrange));
    // }
    
    // Big r;
    // ReqSig rs = qsig.request_sig(Q, r);
    // Sig sig = qsig.sign(rs);

    // qsig.finalize(sig, r);


    // cout << qsig.verify(sig, Q) << endl;

    // Acc acc;
    // acc.add(rand(qsig.mrange));
    // Witness wit = acc.create_nonmem(Q[0]);

    // cout << acc.verify_nonmem(Q[0], wit) << endl;

    // NonMemProof prf = acc.prove_nonmem(Q[0], wit);

    // cout << acc.verify_nonmem_proof(prf) << endl;

    test_protocol();
    // bmark_comm();
    // bmark_issuance();
    // bmark_auth();
    return 0;
}
