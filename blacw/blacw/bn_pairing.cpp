#include "bn_pairing.hpp"

using namespace BN_Pairing;

// Use Galbraith & Scott Homomorphism idea ...

    void galscott(Big &e,Big &r,Big WB[4],Big B[4][4],Big u[4])
    {
        int i,j;
        Big v[4],w;

        for (i=0;i<4;i++)
        {
            v[i]=mad(WB[i],e,(Big)0,r,w);
            u[i]=0;
        }

        u[0]=e;
        for (i=0;i<4;i++)
            for (j=0;j<4;j++)
                u[i]-=v[j]*B[j][i];
        return;
    }

    // GLV method

    void glv(Big &e,Big &r,Big W[2],Big B[2][2],Big u[2])
    {
        int i,j;
        Big v[2],w;
        for (i=0;i<2;i++)
        {
            v[i]=mad(W[i],e,(Big)0,r,w);
            u[i]=0;
        }
        u[0]=e;
        for (i=0;i<2;i++)
            for (j=0;j<2;j++)
                u[i]-=v[j]*B[j][i];
        return;
    }

    //
    // Add A=A+B  (or A=A+A) 
    // Return line function value
    //

    //
    // Line from A to destination C. Let A=(x,y)
    // Line Y-slope.X-c=0, through A, so intercept c=y-slope.x
    // Line Y-slope.X-y+slope.x = (Y-y)-slope.(X-x) = 0
    // Now evaluate at Q -> return (Qy-y)-slope.(Qx-x)
    //

    ZZn12 line(ECn2& A,ECn2& C,ECn2& B,ZZn2& slope,ZZn2& extra,BOOL Doubling,ZZn& Qx,ZZn& Qy)
    {
        ZZn12 w;
        ZZn4 nn,dd,cc;
        ZZn2 X,Y;

    #ifdef AFFINE
        A.get(X,Y);

        if (get_mip()->TWIST==MR_SEXTIC_M)
        {
            nn.set(txx((ZZn2)-Qy),Y-slope*X);
            cc.seth(slope*Qx);
        }
        if (get_mip()->TWIST==MR_SEXTIC_D)
        {
            nn.set((ZZn2)-Qy,Y-slope*X);
            dd.set(slope*Qx);
        }
        w.set(nn,dd,cc);

    #endif
    #ifdef PROJECTIVE
        ZZn2 Z3;

        C.getZ(Z3);

    // Thanks to A. Menezes for pointing out this optimization...
        if (Doubling)
        {
            ZZn2 Z,ZZ;
            A.get(X,Y,Z);
            ZZ=Z; ZZ*=ZZ;
            if (get_mip()->TWIST==MR_SEXTIC_M)
            { // "multiplied across" by i to simplify
                nn.set((Z3*ZZ)*txx((ZZn2)Qy),slope*X-extra);
                cc.seth(-(ZZ*slope)*Qx);
            }
            if (get_mip()->TWIST==MR_SEXTIC_D)
            {
                nn.set((Z3*ZZ)*Qy,slope*X-extra);
                dd.set(-(ZZ*slope)*Qx);
            }
        }
        else
        {
            ZZn2 X2,Y2;
            B.get(X2,Y2);
            if (get_mip()->TWIST==MR_SEXTIC_M)
            {
                nn.set(Z3*txx((ZZn2)Qy),slope*X2-Y2*Z3);
                cc.seth(-slope*Qx);
            }
            if (get_mip()->TWIST==MR_SEXTIC_D)
            {
                nn.set(Z3*Qy,slope*X2-Y2*Z3);
                dd.set(-slope*Qx);
            }
        }
        w.set(nn,dd,cc);
    #endif

        return w;
    }

    void endomorph(ECn &A,ZZn &Beta)
    { // apply endomorphism (x,y) = (Beta*x,y) where Beta is cube root of unity
        ZZn x;
        x=(A.get_point())->X;
        x*=Beta;
        copy(getbig(x),(A.get_point())->X);
    }

    void q_power_frobenius(ECn2 &A,ZZn2 &F)
    { 
    // Fast multiplication of A by q (for Trace-Zero group members only)
        ZZn2 x,y,z,w,r;

    #ifdef AFFINE
        A.get(x,y);
    #else
        A.get(x,y,z);
    #endif

        w=F*F;
        r=F;

        if (get_mip()->TWIST==MR_SEXTIC_M) r=inverse(F);  // could be precalculated
        if (get_mip()->TWIST==MR_SEXTIC_D) r=F;

        w=r*r;
        x=w*conj(x);
        y=r*w*conj(y);

    #ifdef AFFINE
        A.set(x,y);
    #else
        z.conj();
        A.set(x,y,z);

    #endif
    }

    ZZn12 g(ECn2& A,ECn2& B,ZZn& Qx,ZZn& Qy)
    {
        ZZn2 lam,extra;
        ZZn12 r;
        ECn2 P=A;
        BOOL Doubling;

    // Evaluate line from A
        Doubling=A.add(B,lam,extra);

        if (A.iszero())   return (ZZn12)1; 
        r=line(P,A,B,lam,extra,Doubling,Qx,Qy);

        return r;
    }

void BNPairing::set_frobenius_constant(ZZn2 &X)
{
    Big p=get_modulus();
    switch (get_mip()->pmod8)
    {
    case 5:
         X.set((Big)0,(Big)1); // = (sqrt(-2)^(p-1)/2     
         break;
    case 3:                    // = (1+sqrt(-1))^(p-1)/2                                
         X.set((Big)1,(Big)1);      
         break;
   case 7: 
         X.set((Big)2,(Big)1); // = (2+sqrt(-1))^(p-1)/2
    default: break;
    }
    X=pow(X,(p-1)/6);
}

BNPairing::BNPairing(miracl* precision)
{
    mip = precision;
    

    mip->IOBASE=16;

// Set Curve. Note D-Type Curves are recommended. Use BN.CPP program to generate curves.

//  Curve 1.
//	x= (char *)"6000000000160007";  // found by BN.CPP 
//  B=23;
//  mip->TWIST=MR_SEXTIC_D;"

//  Curve 2.
	x= (char *)"-4080000000000001"; 
    B=2;
    mip->TWIST=MR_SEXTIC_D;

//  Curve 3.
//	x= (char *)"408000000000967A"; 
//  B=2;
//  mip->TWIST=MR_SEXTIC_D;

//  Curve 4.
//  x= (char *)"4080000000002C77"; 
//  B=3;
//  mip->TWIST=MR_SEXTIC_M;   // map Server to point on twisted curve E(Fp2)

// See ftp://ftp.computing.dcu.ie/pub/resources/crypto/twists.pdf
// D and M-type twists require a different "untwisting" operation - see paper above

    p=36*pow(x,4)+36*pow(x,3)+24*x*x+6*x+1;
    t=6*x*x+1;
    q=p+1-t;
    cf=p-1+t;
    modulo(p);

    //	Big Lambda=-(36*pow(x,3)+18*x*x+6*x+2);  // cube root of unity mod q
	Beta=-(18*pow(x,3)+18*x*x+9*x+2);    // cube root of unity mod p

    set_frobenius_constant(X);

    // Use standard Gallant-Lambert-Vanstone endomorphism method for G1

    W[0]=6*x*x+4*x+1;      // This is first column of inverse of SB (without division by determinant) 
    W[1]=-(2*x+1);
    
    SB[0][0]=6*x*x+2*x;
    SB[0][1]=-(2*x+1);
    SB[1][0]=-(2*x+1);
    SB[1][1]=-(6*x*x+4*x+1);

    // Use Galbraith & Scott Homomorphism idea for G2 & GT ... (http://eprint.iacr.org/2008/117.pdf Example 5)

	WB[0]=2*x*x+3*x+1;     // This is first column of inverse of BB (without division by determinant)
	WB[1]=12*x*x*x+8*x*x+x;
	WB[2]=6*x*x*x+4*x*x+x;
	WB[3]=-2*x*x-x;

	BB[0][0]=x+1;   BB[0][1]=x;     BB[0][2]=x;        BB[0][3]=-2*x;
	BB[1][0]=2*x+1; BB[1][1]=-x;    BB[1][2]=-(x+1);   BB[1][3]=-x;
	BB[2][0]=2*x;   BB[2][1]=2*x+1; BB[2][2]=2*x+1;    BB[2][3]=2*x+1;
	BB[3][0]=x-1;   BB[3][1]=4*x+2; BB[3][2]=-(2*x-1); BB[3][3]=x-1;

    mip->IOBASE=16;

    #ifdef AFFINE
        ecurve((Big)0,(Big)B,p,MR_AFFINE);
    #endif
    #ifdef PROJECTIVE
        ecurve((Big)0,(Big)B,p,MR_PROJECTIVE);
    #endif
}

// Use GLV endomorphism idea for multiplication in G1


ECn BNPairing::G1_mult(ECn &P,Big &e)
{
//	return e*P;
	int i;
	ECn Q;
	Big u[2];
	glv(e,q,W,SB,u);

	Q=P;
	endomorph(Q,Beta);

	Q=mul(u[0],P,u[1],Q);
	
	return Q;
}

//.. for multiplication in G2

ECn2 BNPairing::G2_mult(ECn2 &P,Big &e)
{
//	return e*P;
	int i;
	ECn2 Q[4];
	Big u[4];
	galscott(e,q,WB,BB,u);

	Q[0]=P;
	for (i=1;i<4;i++)
	{
		Q[i]=Q[i-1]; 
		q_power_frobenius(Q[i],X);
	}

// deal with -ve multipliers
	for (i=0;i<4;i++)
	{
		if (u[i]<0)
			{u[i]=-u[i];Q[i]=-Q[i];}
	}

// simple multi-addition
	return mul4(Q,u);
}

//.. and for exponentiation in GT

ZZn12 BNPairing::GT_pow(ZZn12 &res,Big &e)
{
//	return pow(res,e);
	int i,j;
	ZZn12 Y[4];
	Big u[4];

	galscott(e,q,WB,BB,u);

	Y[0]=res;
	for (i=1;i<4;i++)
		{Y[i]=Y[i-1]; Y[i].powq(X);}

// deal with -ve exponents
	for (i=0;i<4;i++)
	{
		if (u[i]<0)
			{u[i]=-u[i];Y[i].conj();}
	}

// simple multi-exponentiation
	return pow(4,Y,u);
}

//
// R-ate Pairing G2 x G1 -> GT
//
// P is a point of order q in G1. Q(x,y) is a point of order q in G2. 
// Note that P is a point on the sextic twist of the curve over Fp^2, Q(x,y) is a point on the 
// curve over the base field Fp
//

BOOL BNPairing::fast_pairing(ECn2& P, ECn& Q,ZZn12& res)
{ 
    Big xx,yy;
    ZZn Qx,Qy;

    P.norm();
    Q.get(xx,yy); Qx=xx; Qy=yy;

    ECn2 A,KA;
    ZZn2 AX,AY;
    int i,nb;
    Big n;
    ZZn12 r;
    ZZn12 t0,t1;
    ZZn12 x0,x1,x2,x3,x4,x5;

#ifdef MR_COUNT_OPS
fpc=fpa=fpx=fpmq=fpsq=fpaq=0;
#endif

	if (x<0) n=-(6*x+2);
    else n=6*x+2;
    A=P;
    nb=bits(n);
    r=1;
// Short Miller loop
	r.mark_as_miller();

    for (i=nb-2;i>=0;i--)
    {
		r*=r;
		r*=g(A,A,Qx,Qy);
		if (bit(n,i))
            r*=g(A,P,Qx,Qy);
    }
// Combining ideas due to Longa, Aranha et al. and Naehrig
	KA=P;
	q_power_frobenius(KA,X);
	if (x<0) {A=-A; r.conj();}
	r*=g(A,KA,Qx,Qy);
	q_power_frobenius(KA,X); KA=-KA;
	r*=g(A,KA,Qx,Qy);

#ifdef MR_COUNT_OPS
cout << "Miller fpc= " << fpc << endl;
cout << "Miller fpa= " << fpa << endl;
cout << "Miller fpx= " << fpx << endl;
cout << "Miller fpmq= " << fpmq << endl;
cout << "Miller fpsq= " << fpsq << endl;
cout << "Miller fpaq= " << fpaq << endl;

fpa=fpc=fpx=fpmq=fpsq=fpaq=0;
#endif
    if (r.iszero()) return FALSE;

// The final exponentiation

    t0=r;

    r.conj();

    r/=t0;    // r^(p^6-1)
	r.mark_as_regular();  // no longer "miller"

    t0=r;
    r.powq(X); r.powq(X);
    r*=t0;    // r^[(p^6-1)*(p^2+1)]

    r.mark_as_unitary();  // from now on all inverses are just conjugates !! (and squarings are faster)

	res=r;
// Newer new idea...
// See "On the final exponentiation for calculating pairings on ordinary elliptic curves" 
// Michael Scott and Naomi Benger and Manuel Charlemagne and Luis J. Dominguez Perez and Ezekiel J. Kachisa 


    t0=res;    t0.powq(X);
    x0=t0;   x0.powq(X);

    x0*=(res*t0);
    x0.powq(X);

    x1=inverse(res);  // just a conjugation!

    x4=pow(res,-x);  // x is sparse..
    x3=x4; x3.powq(X);

    x2=pow(x4,-x);
    x5=inverse(x2);
	t0=pow(x2,-x);

	x2.powq(X); 
    x4/=x2;
   
    x2.powq(X);
   
    res=t0; res.powq(X); t0*=res;

    t0*=t0;
    t0*=x4;
    t0*=x5;
    res=x3*x5;
    res*=t0;
    t0*=x2;
    res*=res;
    res*=t0;
    res*=res;
    t0=res*x1;
    res*=x0;
    t0*=t0;
    t0*=res;

#ifdef MR_COUNT_OPS
cout << "FE fpc= " << fpc << endl;
cout << "FE fpa= " << fpa << endl;
cout << "FE fpx= " << fpx << endl;
cout << "FE fpmq= " << fpmq << endl;
cout << "FE fpsq= " << fpsq << endl;
cout << "FE fpaq= " << fpaq << endl;
fpa=fpc=fpx=fpmq=fpsq=fpaq=0;
#endif

    res= t0; 
    return TRUE;
}


//
// Faster Hashing to G2 - Fuentes-Castaneda, Knapp and Rodriguez-Henriquez
//

void BNPairing::cofactor(ECn2& S,ZZn2 &F)
{
	ECn2 T,K;
	T=S;
	T*=x;
	T.norm();
	K=(T+T)+T;
	K.norm();
	q_power_frobenius(K,F);
	q_power_frobenius(S,F); q_power_frobenius(S,F); q_power_frobenius(S,F); 
	S+=T; S+=K;
	q_power_frobenius(T,F); q_power_frobenius(T,F);
	S+=T;
	S.norm();
}

//
// ecap(.) function
//

BOOL BNPairing::ecap(ECn2& P,ECn& Q,ZZn12& r)
{
    BOOL Ok;

    Ok=fast_pairing(P,Q,r);

    if (Ok) return TRUE;
    return FALSE;
}

// test if a ZZn12 element is of order q
// test r^q = r^p+1-t =1, so test r^p=r^(t-1)

BOOL BNPairing::member(ZZn12 r)
{
	ZZn12 w=r;
	w.powq(X);
	r=pow(r,x); r=pow(r,x); r=pow(r,(Big)6); // t-1=6x^2
	if (w==r) return TRUE;
	return FALSE;
}

