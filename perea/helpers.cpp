#include "helpers.hpp"

void start_hash(sha256* sh)
{
    shs256_init(sh);
}

void add_to_hash(sha256* sh, Big& x)
{
	int m;
	Big a=x;
    while (a>0)
    {
        m=a%256;
        shs256_process(sh,m);
        a/=256;
    }
    //h%=p;
    //return h;
}

Big finish_hash_to_modulus(sha256* sh, Big p)
{
    Big hash;
	char s[HASH_LEN];
    shs256_hash(sh,s);
    hash=from_binary(HASH_LEN,s);
	return hash%p;
}




long randise()
{ /* get a random number */
    long seed;
    cout << "Enter 9 digit random number seed  = ";
    cin >> seed;
    return seed;
}

Big strongp(int n,long seed1,long seed2)
{ /* generate strong prime number =11 mod 12 suitable for RSA encryption */
    static Big pd,pl,ph;

    Big p;
    int r,r1,r2;
    irand(seed1);
    pd=rand(2*n/3,2);
    pd=nextprime(pd);
    ph=pow((Big)2,n-1)/pd;  
    pl=pow((Big)2,n-2)/pd;
    ph-=pl;
    irand(seed2);
    ph=rand(ph);
    ph+=pl;
    r1=pd%12;
    r2=ph%12;
    r=0;
    while ((r1*(r2+r))%12!=5) r++;
    ph+=r;
    do 
    { /* find p=2*r*pd+1 = 11 mod 12 */
        p=2*ph*pd+1;
        ph+=12;
    } while (!prime(p));
    return p;
}

