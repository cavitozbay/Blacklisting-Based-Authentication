#include "TLP.hpp"

TLP::TLP()
{

}

TLP::TLP(Big& _T)
{
    Big tmp;
    T = _T;

    ifstream public_key("../public.key");
    ifstream private_key("../private.key");
    
    public_key >> n;
    private_key >> p >> q;

    n2 = n*n;
    
    phi_n = (p-1)*(q-1);
    phi_n2 = phi_n*n;
    tmp = rand(n);
    g = (-pow(tmp,2,n))%n;
    tmp = pow(2, T, phi_n/2);
    h = pow(g, tmp, n);

    hn = pow(h, n, n2);

    gc = pow(rand(n), Big(2),n);
    hc = pow(rand(n), Big(2),n);

    rp_range = n2*pow(Big(2),CHL_SPC+ZK_QL);
    sp_range = n*pow(Big(2),CHL_SPC+ZK_QL-2);
    mp_range = pow(Big(2),CHL_SPC+ZK_QL+RHO_SIZE);
    cout << bits(rp_range) << endl;
    cout << bits(sp_range) << endl;
    cout << bits(mp_range) << endl;
    g_inv = inverse(g, n);
    hn_inv = inverse(hn, n2);
}

TLP::~TLP()
{
}

PuzzleAux TLP::pgen(Big& s)
{
    PuzzleAux paux;

    paux.s = s;
    paux.r = rand(n2);

    paux.puzzle.u = pow(g, paux.r, n);
    Big tmp = pow(hn, paux.r, n2);
    paux.puzzle.v = modmult(tmp, (modmult(n,s,n2)+1)%n2, n2);

    return paux;
}

Big TLP::evaluate(Puzzle& puzzle)
{
    Big w = pow(puzzle.u, pow(Big(2),toint(T)), n);
    /*for (size_t i = 1; i < T; i *= 2)
    {
        w = modmult(w,w,n);
    }*/
    Big tmp = moddiv(puzzle.v,pow(w,n,n2),n2)-1;
    tmp /= n;

    return tmp;
}

//This method actually relies on proving equality to a discrete logarithm.
//We don't include discrete logarithm part here directly, we implement it in the protocol
ZKPTLP TLP::prove_wellformedness(PuzzleAux& paux, BOOL explicit_mp, Big mp)
{
    ZKPTLP zkp;
    Big rp, sp, s;
    if(!explicit_mp) mp = rand(2*mp_range) - mp_range;
    rp = rand(2*rp_range) - rp_range;
    sp = rand(2*sp_range) - sp_range;
    s = rand(n/4);
    /*
    cout << bits(rp) << endl;
    cout << bits(sp) << endl;
    cout << bits(mp) << endl;
    */
    //zkp.t_u = pow(g, 2*abs(rp), n);
    //zkp.t_v = modmult(pow(hn, 2*abs(rp), n2), (1 + modmult(n, 2*mp, n2))%n2, n2);
    zkp.l = modmult(pow(gc,paux.s,n),pow(hc,s,n),n);
    //zkp.t_l = modmult(pow(gc,mp,n),pow(hc,sp,n),n);
    if (mp < 0) zkp.t_l = pow(gc_inv, abs(mp), n);
    else zkp.t_l = pow(gc, mp, n);

    if (sp < 0) zkp.t_l = modmult(zkp.t_l, pow(hc_inv, abs(sp), n), n);
    else zkp.t_l = modmult(zkp.t_l, pow(hc, sp, n), n);
    
    if (rp < 0)
    {
        zkp.t_u = pow(g_inv, abs(2*rp), n);
        zkp.t_v = pow(hn_inv, abs(2*rp), n2);
    } 
    else
    {
        zkp.t_u = pow(g, 2*rp, n);
        zkp.t_v = pow(hn, 2*rp, n2);
    } 

    zkp.t_v = modmult(zkp.t_v, (1 + modmult(2*mp, n, n2))%n2, n2);
    
    /*if (mp < 0)
    {
        zkp.t_l = pow(gc_inv, mp, n);
    }
    else
    {
        zkp.t_l = pow(gc, mp, n);
    }*/
    /*
    zkp.l = modmult(pow(gc,paux.s,n),pow(hc,s,n),n);
    zkp.t_l = modmult(pow(gc,mp,n),pow(hc,sp,n),n);
    */
    zkp.c = rand(pow(Big(2),CHL_SPC));

    zkp.z_s = mp - paux.s*zkp.c;
    zkp.z_r = rp - paux.r*zkp.c;
    zkp.z_d = sp - s*zkp.c;

    return zkp;
}

BOOL TLP::verify_wellformedness(Puzzle& puzzle, ZKPTLP& zkp)
{
    Big t_up, t_vp, t_lp;
    
    if (zkp.z_r < 0)
    {
        t_up = pow(g_inv, abs(2*zkp.z_r), n);
        t_vp = pow(hn_inv, abs(2*zkp.z_r), n2);
    } 
    else
    {
        t_up = pow(g, 2*zkp.z_r, n);
        t_vp = pow(hn, 2*zkp.z_r, n2);
    } 

    t_vp = modmult(t_vp, (1+modmult(2*zkp.z_s,n,n2))%n2,n2);
    
    if (zkp.z_s < 0) t_lp = pow(gc_inv, abs(zkp.z_s), n);
    else t_lp = pow(gc, zkp.z_s, n);

    if (zkp.z_d < 0) t_lp = modmult(t_lp, pow(hc_inv, abs(zkp.z_d), n), n);
    else t_lp = modmult(t_lp, pow(hc, zkp.z_d, n), n);

    BOOL v1 = (zkp.t_u == modmult(t_up, pow(puzzle.u, 2*zkp.c, n), n));
    BOOL v2 = (zkp.t_v == modmult(t_vp, pow(puzzle.v, 2*zkp.c, n2), n2));
    BOOL v3 = (zkp.t_l == modmult(t_lp, pow(zkp.l, zkp.c, n), n));
    BOOL v4 = abs(zkp.z_s) < n/4;
    
    /*
    cout << "t_u check: " << v1 << endl;
    cout << "t_v check: " << v2 << endl;
    cout << "t_l check: " << v3 << endl;
    cout << "z_m check: " << v4 << endl;
    */

    return v1&&v2&&v3;
}
/*
ZKPTLP TLP::prove_wellformedness(PuzzleAux& paux)
{
    ZKPTLP zkp;
    Big sp;

    sp = rand(n2);
    PuzzleAux t_paux = pgen(sp);
    zkp.t_puzzle = t_paux.puzzle;

    char bytes[384];
    int tmp = to_binary(paux.puzzle.u, 384, bytes);
    zkp.c = rand(n);

    zkp.z_r = (t_paux.r + modmult(zkp.c, paux.r, n2))%n2;
    zkp.z_s = (t_paux.s + modmult(zkp.c, paux.s, n2))%n2;

    return zkp;
}

BOOL TLP::verify_wellformedness(Puzzle& puzzle, ZKPTLP& zkp)
{
    PuzzleAux z_paux;

    z_paux.s = zkp.z_s;
    z_paux.r = zkp.z_r;

    z_paux.puzzle.u = pow(g, z_paux.r, n);
    Big tmp = pow(hn, z_paux.r, n2);
    z_paux.puzzle.v = modmult(tmp, pow((1 + n), z_paux.s, n2), n2);

    char bytes[384];
    int bin_res = to_binary(puzzle.u, 384, bytes);
    Big c = h1(bytes);

    BOOL v1 = (modmult(pow(puzzle.u, zkp.c, n), zkp.t_puzzle.u, n) == pow(g, zkp.z_r, n));
    //BOOL v1 = (z_paux.puzzle.u == modmult(pow(puzzle.u, c, n), zkp.t_puzzle.u, n));
    BOOL v2 = (z_paux.puzzle.v == modmult(pow(puzzle.v, c, n2), zkp.t_puzzle.v,n2));

    cout << v1 << " " << v2 << endl;

    return (v1 && v2);
}*/