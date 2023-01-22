#include "QueueSig.hpp"

QueueSig::QueueSig(size_t _K)
{
    K = _K;
    two = Big(2);
    mrange = pow(two, lm);
    mtrange = pow(two, lm+lp+lpp);
    ntrange = pow(two, lN+lp+lpp);
    crange = pow(two, lp);

}

void QueueSig::key_gen()
{
    ifstream public_key("../public.key");
    ifstream private_key("../private.key");
    
    // do
    // {
    //     Big seed1, seed2;
    //     seed1 = rand(2,9);
    //     seed2 = rand(2,9);

    //     p = strongp(1024, toint(seed1), toint(seed2));
    // } while (!prime(2*p+1));
    // cout << "found p..." << endl;
    // do
    // {
    //     Big seed1, seed2;
    //     seed1 = rand(2,9);
    //     seed2 = rand(2,9);

    //     q = strongp(1024, toint(seed1), toint(seed2));
    // } while (!prime(2*q+1));
    // cout << "found q..." << endl;
    // N = p*q;
    public_key >> N;
    private_key >> p >> q;

    sk_sig = (p-1)*(q-1);
    Big tmp;
    for (size_t i = 0; i < K+1; i++)
    {
        tmp = rand(N);
        gs.push_back(modmult(tmp,tmp,N));
    }
    tmp = rand(N);
    b = modmult(tmp,tmp,N);
    tmp = rand(N);
    c = modmult(tmp,tmp,N);
    
}

ReqSig QueueSig::request_sig(deque<Big> Q, Big& r)
{
    ReqSig ret;
    Big tmp;
    r = rand(N);
    ret.C = pow(c, r, N);
    for (size_t i = 0; i < K + 1; i++)
    {
        tmp = pow(gs[i], Q[i], N);
        ret.C = modmult(ret.C, tmp, N);
    }
    Big rt = rand(ntrange);
    vector<Big> ts;
    for (size_t i = 0; i < K+1; i++)
    {
        ts.push_back(rand(mtrange));
    }
    ret.Ct = pow(c, rt, N);
    for (size_t i = 0; i < K+1; i++)
    {
        ret.Ct = modmult(ret.Ct, pow(gs[i], ts[i], N), N);
    }
    
    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, ret.C);
    add_to_hash(&sh, ret.Ct);
    Big chl = finish_hash_to_modulus(&sh, crange);

    ret.rz = rt - chl*r;
    for (size_t i = 0; i < K+1; i++)
    {
        ret.zs.push_back(ts[i]-chl*Q[i]);
    }

    return ret;
}

Sig QueueSig::sign(ReqSig& resp)
{
    sha256 sh;
    start_hash(&sh);
    add_to_hash(&sh, resp.C);
    add_to_hash(&sh, resp.Ct);
    Big chl = finish_hash_to_modulus(&sh, crange);

    Big tmp1 = pow(c, resp.rz, N);
    for (size_t i = 0; i < K+1; i++)
    {
        tmp1 = modmult(tmp1, pow(gs[i], resp.zs[i], N), N);
    }
    BOOL res1 = (resp.Ct == modmult(tmp1, pow(resp.C, chl, N), N));
    
    
    // cout << res1 << endl;

    Sig sig;
    Big seed1, seed2;
    seed1 = rand(9,2);
    seed2 = rand(9,2);

    sig.e = strongp(le, toint(seed1), toint(seed2));
    sig.rp = rand(ls,2);
    sig.v = modmult(resp.C, pow(c, sig.rp, N), N);
    sig.v = modmult(sig.v, b, N);

    Big einv = inverse(sig.e, sk_sig);
    sig.v = pow(sig.v, einv, N);
    return sig;
}

void QueueSig::finalize(Sig& sig, Big r)
{
    sig.rp += r;
}

BOOL QueueSig::verify(Sig sig, deque<Big> Q)
{
    Big l = pow(sig.v, sig.e, N);
    Big rght = modmult(b, pow(c, sig.rp, N), N);
    for (size_t i = 0; i < K+1; i++)
    {
        rght = modmult(rght, pow(gs[i], Q[i], N), N);
    }
    return (l == rght);
}

