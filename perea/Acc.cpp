#include "Acc.hpp"

Acc::Acc(/* args */)
{
    ifstream public_key("../public.key");
    ifstream private_key("../private.key");
    public_key >> N;
    Big p,q;
    private_key >> p >> q;
    sk_acc = (p-1)*(q-1);
    g = rand(N);
    g = modmult(g,g,N);
    V = g;
    up = 1;
    lx = 1020;
    g1 = rand(N);
    g1 = modmult(g1,g1,N);
    h1 = rand(N);
    h1 = modmult(h1,h1,N);
    h = rand(N);
    h = modmult(h,h,N);


    two = Big(2);
    mrange = pow(two, lm);
    mtrange = pow(two, lm+lp+lpp);
    ntrange = pow(two, lN+lp+lpp);
    crange = pow(two, lp);

}

void Acc::add(Big t)
{
    V = pow(V, t, N);
    up = modmult(up, t, sk_acc);
}


Witness Acc::create_nonmem(Big t)
{
    Big s, r, z;
    
    egcd(up, t, s, r, z);

    Witness wit;
    wit.a = s;
    wit.d = pow(g, -r, N);

    return wit;
}

BOOL Acc::verify_nonmem(Big t, Witness wit)
{
    BOOL res = (pow(V, wit.a, N) == modmult(pow(wit.d, t, N), g, N));
    return res;
}

NonMemProof Acc::prove_nonmem(Big t, Witness wit)
{
    NonMemProof prf;

    Big r = rand(N);
    Big w = rand(N);
    Big rx = rand(2048+lp, 2);
    Big ra = rand(N);
    Big rw = rand(N);
    Big rz = rand(N);
    Big re = rand(N);

    prf.C1 = modmult(pow(g1,t,N), pow(h1,r,N), N);
    prf.Cx = modmult(pow(g,t,N), pow(h,rx,N), N);
    prf.Ca = modmult(pow(g,wit.a,N), pow(h,ra,N), N);
    prf.Cd = modmult(wit.d, pow(g,w,N), N);
    prf.Cw = modmult(pow(g,w,N), pow(h,rw,N), N);
    Big z = t*w;
    prf.Cz = modmult(pow(g,z,N), pow(h,rz,N), N);
    prf.Ce = modmult(pow(prf.Cd,t,N), pow(h,re,N), N);
    


    Big tt = rand(ntrange);
    Big rt = rand(ntrange);
    Big wt = rand(ntrange);
    Big rxt = rand(lN+lp+2*lpp,2);
    Big ret = rand(ntrange);
    Big at = rand(ntrange);
    Big zt = rand(ntrange)*rand(lm, 2);
    Big rzt = rand(2*lN+lp+2*lpp,2);
    Big rwt = rand(lN+lp+2*lpp,2);
    
    prf.C1t = modmult(pow(g1,tt,N), pow(h1,rt,N), N);
    prf.Cxt = modmult(pow(g,tt,N), pow(h,rxt,N), N);
    prf.Cet = modmult(pow(prf.Cd,tt,N), pow(h,ret,N), N);
    prf.Cegt = modmult(modmult(pow(V,at,N), pow(g,zt,N), N), pow(h, ret, N), N);
    prf.Czt = modmult(pow(prf.Cx,wt,N), pow(h,rzt,N), N);
    prf.Cwt = modmult(pow(g,wt,N), pow(h,rwt,N), N);
    
    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, prf.C1);
    add_to_hash(&sh, prf.Ca);
    add_to_hash(&sh, prf.Cd);
    add_to_hash(&sh, prf.Ce);
    add_to_hash(&sh, prf.Cw);
    add_to_hash(&sh, prf.Cx);
    add_to_hash(&sh, prf.Cz);
    add_to_hash(&sh, prf.C1t);
    add_to_hash(&sh, prf.Cxt);
    add_to_hash(&sh, prf.Cet);
    add_to_hash(&sh, prf.Cegt);
    add_to_hash(&sh, prf.Cwt);
    add_to_hash(&sh, prf.Czt);
    Big chl = finish_hash_to_modulus(&sh, crange);

    prf.zt = tt - chl*t;
    prf.zr = rt - chl*r;
    prf.zrx = rxt - chl*rx;
    prf.zre = ret - chl*re;
    prf.za = at - chl*wit.a;
    prf.zz = zt - chl*z;
    prf.zw = wt - chl*w;
    prf.zrz = rzt - chl*(rz - rx*w);
    prf.zrw = rwt - chl*rw;

    prf.rng1 = prove_range(prf.Cx, t, rx, 0, pow(two, lm));
    prf.rng2 = prove_range(prf.Ca, wit.a, ra, 0, pow(two, lm));

    return prf;
}

BOOL Acc::verify_nonmem_proof(NonMemProof prf)
{
        sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, prf.C1);
    add_to_hash(&sh, prf.Ca);
    add_to_hash(&sh, prf.Cd);
    add_to_hash(&sh, prf.Ce);
    add_to_hash(&sh, prf.Cw);
    add_to_hash(&sh, prf.Cx);
    add_to_hash(&sh, prf.Cz);
    add_to_hash(&sh, prf.C1t);
    add_to_hash(&sh, prf.Cxt);
    add_to_hash(&sh, prf.Cet);
    add_to_hash(&sh, prf.Cegt);
    add_to_hash(&sh, prf.Cwt);
    add_to_hash(&sh, prf.Czt);
    Big chl = finish_hash_to_modulus(&sh, crange);

    BOOL res1 = (prf.C1t == modmult(modmult(pow(prf.C1,chl,N), pow(g1, prf.zt, N), N), pow(h1, prf.zr, N), N));
    BOOL res2 = (prf.Cxt == modmult(modmult(pow(prf.Cx,chl,N), pow(g, prf.zt, N), N), pow(h, prf.zrx, N), N));
    BOOL res3 = (prf.Cet == modmult(modmult(pow(prf.Ce,chl,N), pow(prf.Cd, prf.zt, N), N), pow(h, prf.zre, N), N));
    BOOL res4 = (prf.Cegt == modmult(modmult(pow(modmult(prf.Ce, g, N),chl,N), pow(V, prf.za, N), N), modmult(pow(g, prf.zz, N),pow(h, prf.zre, N), N), N));
    BOOL res5 = (prf.Cwt == modmult(modmult(pow(prf.Cw,chl,N), pow(g, prf.zw, N), N), pow(h, prf.zrw, N), N));
    BOOL res6 = (prf.Czt == modmult(modmult(pow(prf.Cz,chl,N), pow(prf.Cx, prf.zw, N), N), pow(h,prf.zrz,N), N));
    BOOL res7 = verify_range(prf.rng1, prf.Cx, 0, pow(two, lm));
    BOOL res8 = verify_range(prf.rng2, prf.Ca, 0, pow(two, lm));
    // cout << res1 << res2 << res3 << res4 << res5 << res6 << endl;
    return res1 && res2 && res3 && res4 && res5 && res6 && res7 && res8;
}

RangeProof Acc::prove_range(Big E, Big x, Big r, Big a, Big b)
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

BOOL Acc::verify_range(RangeProof prf, Big E, Big a, Big b)
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
