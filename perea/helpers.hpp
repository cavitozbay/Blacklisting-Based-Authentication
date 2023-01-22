#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <big.h>

#define HASH_LEN 32

#define LOG

struct ReqSig
{
    Big C, Ct;
    Big rz;
    vector<Big> zs;

};

struct Sig
{
    Big rp, e, v;
};


struct Witness
{
    Big a, d;
};

struct RangeProof
{
    Big Etld, Ecap;
    Big Etld1, Ecap1, Etld2, Ecap2;
    Big Wtld, Wcap;
    Big Ftld, Fcap;
    Big Ftldt, Fcapt;
    Big Etld1t, Ecap1t;
    Big dtld1, dtld2, dcap1, dcap2;
    Big zxtld1, zr1tld, zr2tld, zr3tld;
    Big zxcap1, zr1cap, zr2cap, zr3cap;
};

struct NonMemProof
{
    Big C1, Cx, Ca, Cd, Cw, Cz, Ce;
    Big C1t, Cxt, Cet, Cegt, Cwt, Czt;
    Big zt, zr, zrx, zre, za, zz, zrz, zrw, zw;
    RangeProof rng1, rng2;
};

struct RegResp
{
    Sig sig;
    Witness wit;
};

struct KoSProof
{
    Big C, Ce, Cv, Cx, Cs, Cd, Cw, Cz, Czmult ;
    Big Ct, Cet, Cxt, Cct, Cst, Cwt, Czt;
    vector<Big> zxz;
    Big ze, zr, zre, zrx, zrc, zz, zrz, zrw, zw, zrs, zrzmul;
    RangeProof rng1, rng2;
};
struct InitAuth
{
    KoSProof kos;
    ReqSig req2;
    vector<NonMemProof> nmps;
    Big tK;
};
struct AuthResp
{
    Sig sig;
    Witness wit;
};




void start_hash(sha256* sh);
void add_to_hash(sha256* sh, Big& x);
Big finish_hash_to_modulus(sha256* sh, Big p);
Big strongp(int n,long seed1,long seed2);
#endif