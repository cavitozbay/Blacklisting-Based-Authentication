#include "Issuer.hpp"

Issuer::Issuer(size_t K) : acc(), qsig(K)
{
    two = Big(2);
    mrange = pow(two, lm);
    mtrange = pow(two, lm+lp+lpp);
    ntrange = pow(two, lN+lp+lpp);
    crange = pow(two, lp);
    g = acc.g;
    h = acc.h;
    N = acc.N;

    qsig.key_gen();
    for (size_t i = 0; i < 5; i++)
    {
        long seed1 = toint(rand(9,2));
        long seed2 = toint(rand(9,2));
        Big tmp = strongp(128,seed1, seed2);
        acc.add(tmp);
    }
    
}

RegResp Issuer::_register(ReqSig req)
{
    RegResp ret;
    ret.sig = qsig.sign(req);
    ret.wit = acc.create_nonmem(Big(5));
#ifdef LOG
    cout << "REG==========" << endl;
#endif
    return ret;
}

AuthResp Issuer::auth(InitAuth req)
{
    AuthResp ret;

    BOOL res1 = verify_kos(req.kos);
#ifdef LOG
    cout << "ZKP_AUTH==========" << endl;
    cout << res1 << endl;
#endif

    ret.sig = qsig.sign(req.req2);

    for (size_t i = 0; i < qsig.K-1; i++)
    {
        acc.verify_nonmem_proof(req.nmps[i]);
    }

    ret.wit = acc.create_nonmem(req.tK);
    //cout << "nm chck" << acc.verify_nonmem(req.tK, ret.wit);
    return ret;
}

BOOL Issuer::verify_kos(KoSProof prf)
{
    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, prf.C);
    add_to_hash(&sh, prf.Ce);
    add_to_hash(&sh, prf.Cv);
    add_to_hash(&sh, prf.Cx);
    add_to_hash(&sh, prf.Cs);
    add_to_hash(&sh, prf.Cd);
    add_to_hash(&sh, prf.Cw);
    add_to_hash(&sh, prf.Cz);
    add_to_hash(&sh, prf.Czmult);
    add_to_hash(&sh, prf.Ct);
    add_to_hash(&sh, prf.Cet);
    add_to_hash(&sh, prf.Cxt);
    add_to_hash(&sh, prf.Cct);
    add_to_hash(&sh, prf.Cst);
    add_to_hash(&sh, prf.Cwt);
    add_to_hash(&sh, prf.Czt);
    Big chl = finish_hash_to_modulus(&sh, crange);

    BOOL res1 = (prf.Ct == modmult(modmult(pow(prf.C,chl,N), pow(prf.Cv, prf.ze, N), N), pow(h, prf.zr, N), N));
    BOOL res2 = (prf.Cet == modmult(modmult(pow(prf.Ce,chl,N), pow(g, prf.ze, N), N), pow(h, prf.zre, N), N));
    Big tmp3 = modmult(pow(h, prf.zrx, N), pow(prf.Cx, chl, N), N);
    Big tmp4 = modmult(pow(h, prf.zr, N), pow(moddiv(prf.C, qsig.b, N), chl, N), N);
    tmp4 = modmult(tmp4, modmult(pow(qsig.c, prf.zrc, N), pow(g, prf.zz, N),N), N);
    for (size_t i = 0; i < qsig.K+1; i++)
    {
        tmp3 = modmult(tmp3, pow(qsig.gs[i], prf.zxz[i], N), N);
        tmp4 = modmult(tmp4, pow(qsig.gs[i], prf.zxz[i], N), N);
    }

    BOOL res3 = (prf.Cxt == tmp3);
    BOOL res4 = (prf.Cct == tmp4);
    BOOL res5 = (prf.Cst == modmult(modmult(pow(g, prf.zrc, N), pow(h, prf.zrs, N), N), pow(prf.Cs, chl, N), N));
    BOOL res6 = (prf.Czt == modmult(modmult(pow(g, prf.zz, N), pow(h, prf.zrz, N), N), pow(prf.Cz, chl, N), N));
    BOOL res7 = (prf.Cwt == modmult(modmult(pow(prf.Cw,chl,N), pow(g, prf.zw, N), N), pow(h, prf.zrw, N), N));
    BOOL res8 = (prf.Czmult == modmult(modmult(pow(prf.Ce, prf.zw, N), pow(h, prf.zrzmul, N), N), pow(prf.Cz, chl, N), N));
    BOOL res9 = verify_range(prf.rng1,prf.Ce,pow(two, le-1), pow(two,le));
    // BOOL res3 = (prf.Cet == modmult(modmult(pow(prf.Ce,chl,N), pow(prf.Cd, prf.zt, N), N), pow(h, prf.zre, N), N));
    // BOOL res4 = (prf.Cegt == modmult(modmult(pow(modmult(prf.Ce, g, N),chl,N), pow(V, prf.za, N), N), modmult(pow(g, prf.zz, N),pow(h, prf.zre, N), N), N));
    // BOOL res6 = (prf.Czt == modmult(modmult(pow(prf.Cz,chl,N), pow(prf.Cx, prf.zw, N), N), pow(h,prf.zrz,N), N));
    // BOOL res7 = verify_range(prf.rng1, prf.Cx, 0, pow(two, lm));
    // BOOL res8 = verify_range(prf.rng2, prf.Ca, 0, pow(two, lm));
    //cout << res1 << res2 << res3 << res4 << res5 << res6 << res7 << res8 << res9 << endl;// << res3 << res4 << res5 << res6 << endl;
    return res1 && res2 && res3 && res4 && res5 && res6 && res7 && res8 && res9;
}

BOOL Issuer::verify_range(RangeProof prf, Big E, Big a, Big b)
{
    Big rng = sqrt(b-a);

    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, prf.Wtld);
    add_to_hash(&sh, prf.Wcap);
    Big chl = finish_hash_to_modulus(&sh, crange);

    BOOL res1 = (prf.dtld1 > chl*rng && prf.dtld1 < pow(two, lp+lpp)*rng);
    BOOL res2 = (modmult(prf.Wtld,pow(prf.Etld2, chl, N), N) == modmult(pow(g, prf.dtld1, N), pow(h, prf.dtld2, N), N));
    BOOL res3 = (prf.dcap1 > chl*rng && prf.dcap1 < pow(two, lp+lpp)*rng);
    BOOL res4 = (prf.Wcap == modmult(pow(prf.Ecap2, chl, N),modmult(pow(g, prf.dcap1, N), pow(h, prf.dcap2, N), N), N));
    BOOL res5 = (modmult(prf.Ftldt,pow(prf.Ftld, chl, N),N) == modmult(modmult(pow(g,prf.zxtld1,N), pow(h,prf.zr2tld,N), N), 1, N));
    BOOL res6 = (modmult(prf.Etld1t,pow(prf.Etld1, chl, N),N) == modmult(modmult(pow(prf.Ftld,prf.zxtld1,N), pow(h,prf.zr3tld,N), N), 1, N));
    BOOL res7 = (modmult(prf.Fcapt,pow(prf.Fcap, chl, N),N) == modmult(modmult(pow(g,prf.zxcap1,N), pow(h,prf.zr2cap,N), N), 1, N));
    BOOL res8 = (modmult(prf.Ecap1t,pow(prf.Ecap1, chl, N),N) == modmult(modmult(pow(prf.Fcap,prf.zxcap1,N), pow(h,prf.zr3cap,N), N), 1, N));

    // cout << res1 << res2 << res3 << res4 << res5 << res6 << res7 << res8 << endl;

    return res1 && res2 && res3 && res4 && res5 && res6 && res7 && res8;
}