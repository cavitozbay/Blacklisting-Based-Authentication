#include "User.hpp"

using namespace std;
User::User(Issuer& issuer): acc(issuer.acc), qsig(issuer.qsig)
{
    two = Big(2);
    mrange = pow(two, lm);
    mtrange = pow(two, lm+lp+lpp);
    ntrange = pow(two, lN+lp+lpp);
    crange = pow(two, lp);
    g = acc.g;
    h = acc.h;
    N = acc.N;

}

ReqSig User::init_register()
{   
    for (size_t i = 0; i < qsig.K; i++)
    {
        Q.push_back(Big(5));
    }
    long seed1 = toint(rand(9,2));
    long seed2 = toint(rand(9,2));
    Big tS = strongp(qsig.lm,seed1, seed2);

    Q.push_back(tS);
    ReqSig req = qsig.request_sig(Q,r);
    return req;
}

void User::finish_register(RegResp resp)
{
    qsig.finalize(resp.sig, r);
    BOOL res1 = qsig.verify(resp.sig, Q);
    BOOL res2 = acc.verify_nonmem(Big(5), resp.wit);
    sig = resp.sig;
    for (size_t i = 0; i < qsig.K; i++)
    {
        wits.push_back(resp.wit);
    }

}

InitAuth User::init_auth()
{
    InitAuth ret;
    Big tmp;
    ret.kos = prove_kos(sig, Q);
    ret.tK = Q.back();
    Q.pop_front();
    wits.pop_front();
    long seed1 = toint(rand(9,2));
    long seed2 = toint(rand(9,2));
    Big tS = strongp(qsig.lm,seed1, seed2);
    Q.push_back(tS);

    ret.req2 = qsig.request_sig(Q,r);
    
    for (size_t i = 0; i < qsig.K-1; i++)
    {
        ret.nmps.push_back(acc.prove_nonmem(Q[i], wits[i]));
    }
    // cout << ret.tK << endl;
    return ret;
}

void User::finish_auth(AuthResp resp)
{
    qsig.finalize(resp.sig, r);
    BOOL res1 = qsig.verify(resp.sig, Q);
    BOOL res2 = acc.verify_nonmem(Q[qsig.K-1], resp.wit);
    sig = resp.sig;
    wits.push_back(resp.wit);
    // cout << res1 << res2 << endl;
}

KoSProof User::prove_kos(Sig sig, deque<Big> Q)
{
    KoSProof prf;

    Big r = rand(N);
    Big w = rand(N);
    Big re = rand(N);
    Big rx = rand(lN+lp, 2);
    Big rs = rand(lN+lp, 2);
    // Big ra = rand(N);
    Big rw = rand(lN+lp, 2);
    Big rz = rand(lN+lp, 2);
    
    prf.Cv = modmult(sig.v, pow(g,w, N), N);
    prf.C = modmult(pow(prf.Cv,sig.e,N), pow(h,r,N), N);
    prf.Ce = modmult(pow(g,sig.e,N), pow(h,re,N), N);
    prf.Cx = pow(h,rx,N);
    prf.Cs = modmult(pow(g, sig.rp, N), pow(h, rs, N), N);
    for (size_t i = 0; i < qsig.K+1; i++)
    {
        prf.Cx = modmult(prf.Cx, pow(qsig.gs[i], Q[i], N), N);
    }
    
    // prf.Ca = modmult(pow(g,wit.a,N), pow(h,ra,N), N);
    // prf.Cd = modmult(wit.d, pow(g,w,N), N);
    prf.Cw = modmult(pow(g,w,N), pow(h,rw,N), N);
    Big z = sig.e*w;
    prf.Cz = modmult(pow(g,z,N), pow(h,rz,N), N);
    // // prf.Ce = modmult(pow(prf.Cd,t,N), pow(h,re,N), N);
    


    Big et = rand(ntrange);
    Big rt = rand(ntrange);
    Big ret = rand(ntrange);
    vector<Big> xts;
    for (size_t i = 0; i < qsig.K+1; i++)
    {
        xts.push_back(rand(ntrange));
    }
    Big rxt = rand(lN+lp+2*lpp,2);
    Big rct = rand(lN+lp+2*lpp,2);
    Big rst = rand(lN+lp+2*lpp,2);

    Big wt = rand(ntrange);
    // Big rxt = rand(lN+lp+2*lpp,2);
    // Big at = rand(ntrange);
    Big zt = rand(ntrange)*rand(lm, 2);
    Big rzt = rand(2*lN+lp+2*lpp,2);
    Big rwt = rand(lN+lp+2*lpp,2);
    Big rzmult = rand(2*lN+lp+2*lpp,2);
    
    prf.Ct = modmult(pow(prf.Cv,et,N), pow(h,rt,N), N);
    prf.Cet = modmult(pow(g,et,N), pow(h,ret,N), N);
    prf.Cxt = pow(h, rxt, N);
    prf.Cct = modmult(modmult(pow(h, rt, N), pow(g, zt, N), N), pow(qsig.c, rct, N), N);
    for (size_t i = 0; i < qsig.K+1; i++)
    {
        prf.Cxt = modmult(prf.Cxt, pow(qsig.gs[i], xts[i], N), N);
        prf.Cct = modmult(prf.Cct, pow(qsig.gs[i], xts[i], N), N);
    }
    prf.Cst = modmult(pow(g, rct, N), pow(h, rst, N), N);
    prf.Czt = modmult(pow(g,zt,N), pow(h,rzt,N), N);
    prf.Cwt = modmult(pow(g,wt,N), pow(h,rwt,N), N);
    prf.Czmult = modmult(pow(prf.Ce,wt,N), pow(h,rzmult,N), N);

    // prf.Cet = modmult(pow(prf.Cd,tt,N), pow(h,ret,N), N);
    // prf.Cegt = modmult(modmult(pow(V,at,N), pow(g,zt,N), N), pow(h, ret, N), N);
    // prf.Cwt = modmult(pow(g,wt,N), pow(h,rwt,N), N);
    
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

    prf.ze = et - chl*sig.e;
    prf.zr = rt - chl*r;
    prf.zre = ret - chl*re;
    for (size_t i = 0; i < qsig.K+1; i++)
    {
        prf.zxz.push_back(xts[i]-chl*Q[i]);
    }
    prf.zrx = rxt - chl*rx;
    prf.zrc = rct - chl*sig.rp;
    // prf.za = at - chl*wit.a;
    prf.zz = zt - chl*z;
    prf.zrs = rst - chl*rs;
    prf.zw = wt - chl*w;
    prf.zrz = rzt - chl*rz;
    prf.zrw = rwt - chl*rw;
    prf.zrzmul = rzmult - chl*(rz - w*re);
    prf.rng1 = prove_range(prf.Ce, sig.e, re, pow(two, le-1), pow(two,le));

    // prf.rng1 = prove_range(prf.Cx, t, rx, 0, pow(two, lm));
    // prf.rng2 = prove_range(prf.Ca, wit.a, ra, 0, pow(two, lm));

    return prf;
}




RangeProof User::prove_range(Big E, Big x, Big r, Big a, Big b)
{
    // cout << "RNGPRF" << endl;
    RangeProof prf;
    Big rng = sqrt(b-a);
    prf.Etld = moddiv(E, pow(g, a, N), N);
    prf.Ecap = moddiv(pow(g, b, N), E, N);

    Big xtld = x - a;
    Big xcap = b - x;

    Big xtld1 = sqrt(xtld);
    Big xtld2 = xtld - xtld1*xtld1;
    Big xcap1 = sqrt(xcap);
    Big xcap2 = xcap - xcap1*xcap1;

    Big rtld1 = rand(lN+lp, 2); 
    Big rtld2 = r - rtld1; 

    Big rcap1 = rand(lN+lp, 2); 
    Big rcap2 = -r - rcap1; 

    prf.Etld1 = modmult(pow(g, xtld1*xtld1, N), pow(h, rtld1, N), N);
    prf.Ecap1 = modmult(pow(g, xcap1*xcap1, N), pow(h, rcap1, N), N);
    prf.Etld2 = moddiv(prf.Etld, prf.Etld1, N);
    prf.Ecap2 = moddiv(prf.Ecap, prf.Ecap1, N);

    Big omgtld = rand(pow(two, lp+lpp)*rng);
    Big kptld = rand(lp+2*lpp+lN, 2);
    prf.Wtld = modmult(pow(g, omgtld, N), pow(h, kptld, N), N);

    Big omgcap = rand(pow(two, lp+lpp)*rng);
    Big kpcap = rand(lp+2*lpp+lN, 2);
    prf.Wcap = modmult(pow(g, omgcap, N), pow(h, kpcap, N), N);

    Big r2tld = rand(lN+lp, 2);
    prf.Ftld = modmult(pow(g, xtld1, N), pow(h, r2tld, N), N);
    Big r3tld = rtld1 - r2tld*xtld1;

    Big xtld1t = rand(lN+lp, 2);
    Big r2tldt = rand(lN+lp, 2);
    Big r3tldt = rand(3*lN+3*lp, 2);

    prf.Etld1t = modmult(pow(prf.Ftld, xtld1t, N), pow(h, r3tldt, N), N);
    prf.Ftldt = modmult(pow(g, xtld1t, N), pow(h, r2tldt, N), N);

    Big r2cap = rand(lN+lp, 2);
    prf.Fcap = modmult(pow(g, xcap1, N), pow(h, r2cap, N), N);
    Big r3cap = rcap1 - r2cap*xcap1;

    Big xcap1t = rand(lN+lp, 2);
    Big r2capt = rand(lN+lp, 2);
    Big r3capt = rand(3*lN+3*lp, 2);

    prf.Ecap1t = modmult(pow(prf.Fcap, xcap1t, N), pow(h, r3capt, N), N);
    prf.Fcapt = modmult(pow(g, xcap1t, N), pow(h, r2capt, N), N);

    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, prf.Wtld);
    add_to_hash(&sh, prf.Wcap);
    Big chl = finish_hash_to_modulus(&sh, crange);


    prf.dtld1 = omgtld + xtld2*chl;
    prf.dtld2 = kptld + rtld2*chl;
    prf.dcap1 = omgcap - xcap2*chl;
    prf.dcap2 = kpcap - rcap2*chl;

    prf.zxtld1 = xtld1t + chl*xtld1;
    prf.zr2tld = r2tldt + chl*r2tld;
    prf.zr3tld = r3tldt + chl*r3tld;

    prf.zxcap1 = xcap1t + chl*xcap1;
    prf.zr2cap = r2capt + chl*r2cap;
    prf.zr3cap = r3capt + chl*r3cap;


    return prf;
}


