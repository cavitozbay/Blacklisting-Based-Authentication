#include "helpers.hpp"

using namespace std;
size_t INITIAL_DEBT_COUNT = 0;
size_t BLACKLISTED_COUNT = 100;
size_t ACC_SIZE = 120;
size_t PCOMP_SIZE = 5;

/*Poly build_from_roots(vector<ZZn>& roots, Big modulus)
{
    modulo(modulus);
    Poly f, g;
    ZZn root(0);
    g.addterm(root, 0);
    g.addterm(ZZn(1), 1);
    f.addterm(ZZn(1),0);
    //cout << "====Root Build: " << endl;
    //cout << "      f: " << f << endl;
    for (size_t i = 0; i < roots.size(); i++)
    {
        root = roots[i] - g.coeff(0);
        //cout << "       " << root << endl;
        g.addterm(root, 0);
        f = f * g;
        //cout << "      f: " << f << endl;
    }

    return f;
}*/

Poly build_from_roots(vector<ZZn>& roots, Big modulus)
{
    modulo(modulus);

    return recursive_bfr(roots.size(), 0, roots, modulus);
}

Poly recursive_bfr(size_t num, size_t pos, vector<ZZn>& roots, Big modulus)
{
    modulo(modulus);
    if(num == 0)
    {
        Poly f;
        f.addterm(ZZn(1), 0);
        return f;
    }
    if(num == 1)
    {
        Poly f;
        ZZn root(0);
        f.addterm(roots[pos], 0);
        f.addterm(ZZn(1), 1);
        return f;
    }
    return recursive_bfr(num/2, pos, roots, modulus)*recursive_bfr(num/2 + num%2, pos+num/2,roots,modulus);
}

G2 mult_poly(PFC& bn, vector<G2>& rs, Poly& f)
{
    modulo(*bn.mod);
    G2 res, t1, t1_;
    size_t deg = degree(f);
    modulo(bn.order());
    Big tmp = Big(f.coeff(0));
    modulo(*bn.mod);
    res = bn.mult(rs[0], tmp);

    if(deg == 0) return res;

    for (size_t i = 1; i <= deg; i++)
    {
        modulo(bn.order());

        tmp = Big(f.coeff(i));
        modulo(*bn.mod);
        t1 = bn.mult(rs[i], tmp);
        t1_ = res;
        res = t1 + t1_;
    }
    modulo(*bn.mod);
    
    return res;
}

G1 mult_poly_hat(PFC& bn, vector<G1>& rs, Poly& f)
{
    modulo(*bn.mod);
    
    G1 res, t2, t2_;
    size_t deg = degree(f);
    modulo(bn.order());
    //cout << "f: " << f << " !! " << deg << endl;
    Big tmp = Big(f.coeff(0));
    modulo(*bn.mod);
    res = bn.mult(rs[0], tmp);
    if(deg == 0) return res;

    for (size_t i = 1; i <= deg; i++)
    {
        modulo(bn.order());
        tmp = Big(f.coeff(i));
        //cout << tmp << endl;
        modulo(*bn.mod);
        t2 = bn.mult(rs[i], tmp);
        t2_ = res;
        res = t2 + t2_; 
    }
    modulo(*bn.mod);
    
    return res;
}

Big generate_prime()
{
    BOOL found;
    int i,spins;
    long seed;
    Big pp[NPRIMES],q,p,t;

    //cout << "Enter 9 digit seed= ";
    //cin >> seed;
    seed = rand();
    seed %= 1024;
    
    irand(seed);
    //cout << "Enter 4 digit seed= ";
    //cin >> spins;
    spins = rand();
    spins %= 32;

    for (i=0;i<spins;i++) brand();
    pp[0]=2;
    do
    {  /* find prime p = 2.pp[1].pp[2]....+1 */
        p=2;
        for (i=1;i<NPRIMES-1;i++)
        { /* generate all but last prime */
            q=rand(i+6,10);
            pp[i]=nextprime(q);
            p*=pp[i];
        }
        do
        { /* find last prime component such that p is prime */
            q=nextprime(q);
            pp[NPRIMES-1]=q;
            t=p*pp[NPRIMES-1];
            t+=1;
        } while(!prime(t));
        p=t;
        found=TRUE;
        for (i=0;i<NPRIMES;i++)
        { /* check that PROOT is a primitive root */
            if (pow(PROOT,(p-1)/pp[i],p)==1) 
            {
                found=FALSE;
                break;
            }
        }
    } while (!found);

    return p;
}

Big H2(G2& m)
{ // hash G2 point to 160-bit big number
    //Big x, y, z
    ZZn2 x, y;
    Big x1, x2, y1, y2;
    m.g.get(x,y);
    x.get(x1,x2);
    y.get(y1,y2);
    //char* ID = new char[100];
    char id[100];
    id << x1 + x2 + y1 + y2;
    char* ID = id;
    int b;
    Big h;
    char s[20];
    sha sh;
    shs_init(&sh);
    while (*ID!=0) shs_process(&sh,*ID++);
    shs_hash(&sh,s);
    h=from_binary(20,s);

    return h;
}

Big H1(G1& m)
{ // hash G2 point to 160-bit big number
    Big x, y, z;
    m.g.getxyz(x,y,z);
    //char* ID = new char[100];
    char id[100];
    id << x + y + z;
    char* ID = id;
    int b;
    Big h;
    char s[20];
    sha sh;
    shs_init(&sh);
    while (*ID!=0) shs_process(&sh,*ID++);
    shs_hash(&sh,s);
    h=from_binary(20,s);

    return h;
}

void trivial_pairing(G2& g, G1& g_hat, PFC& bn)
{
    GT l = bn.pairing(g,bn.mult(g_hat,Big(5)));
    GT r = bn.pairing(bn.mult(g,Big(5)),g_hat);
    BOOL res = l == r;
    cout << "Trivial Pairing: " << res << endl;
}



Big h1(char *string)
{ // Hash a zero-terminated string to a number < modulus
    Big h,p;
    char s[HASH_LEN];
    int i,j; 
    sha256 sh;

    shs256_init(&sh);

    for (i=0;;i++)
    {
        if (string[i]==0) break;
        shs256_process(&sh,string[i]);
    }
    shs256_hash(&sh,s);
    p=get_modulus();
    h=1; j=0; i=1;
    forever
    {
        h*=256; 
        if (j==HASH_LEN)  {h+=i++; j=0;}
        else         h+=s[j++];
        if (h>=p) break;
    }
    h%=p;
    return h;
}
