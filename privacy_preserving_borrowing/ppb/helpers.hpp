#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <fstream>
#include <cstring>
#include "zzn.h"
#include "ecn.h"
#include "ecn2.h"
#include "poly.h"
#include "pairing_bn.h"
#include "big.h"   /* include MIRACL system */

#define NPRIMES 10  /* =9 for > 256 bit primes */
#define PROOT 2
#define HASH_LEN 32

//#define LOG
#define MIN_TIME 10.0
#define MIN_ITERS 30
//#define INITIAL_DEBT_COUNT 1
//#define BLACKLISTED_COUNT 100
//#define ACC_SIZE 120
extern size_t INITIAL_DEBT_COUNT;
extern size_t BLACKLISTED_COUNT;
extern size_t ACC_SIZE;
extern size_t PCOMP_SIZE;



using namespace std;

struct Puzzle
{
    Big u, v;
};

struct PuzzleAux
{
    Big r, s;
    Puzzle puzzle;
};

struct ZKP1
{
    G2 F;
    G1 E;
    G2 B_t;
    G2 D_t;
    G1 E_t;
    G2 F_t;
    G2 G_t;
    G1 M_t;
    vector<Big> com_alpha_zs;
    Big z_u, z_y, z_mult, z_tmp, z_gamma, z_delta;
};

struct ZKPNonMembership
{
    G1 C_hat_h0;
    G2 C_h;
    G1 D_a, D_a_t;
    G1 D_b, D_b_t;
    G1 M_t;
    GT E_t;

    Big z_a, z_b, z_r_a, z_r_b, z_m_acc, z_tmp_acc;
};

/*
struct ZKPTLP
{
    Puzzle t_puzzle;
    Big z_r, z_s;
    Big c;
};
*/

struct ZKPTLP
{
    Big t_u, t_v, t_l;
    Big l;
    Big c;
    Big z_s, z_r, z_d;
};

struct ZKP2
{
    G2 F;
    G1 E;
    G2 B_t;
    G2 D_t;
    G1 E_t;
    G2 F_t;
    G2 G_t;
    G1 M_t;
    
    
    G1 C_hat_bidS, C_hat_bidS_t;
    G1 C_yS, C_yS_t;
    G1 C_hat_y, C_hat_y_t;
    G1 M2_t;
    G1 M_bid_t;
    GT N_t;

    ZKPNonMembership nonmember;

    vector<Big> com_alpha_zs;
    Big z_u, z_yS, z_mult, z_tmp, z_gamma, z_delta;

    Big z_y, z_r_y, z_bidS, z_r_bidS, z_r_yS, z_mS, z_tmpS, z_m_bidS, z_tmp_bidS;

    vector<Big> com_old_alpha_zs;
    G2 F_old_t;
};

struct ZKP3
{
    G1 C_yS, C_y;
    G1 C_yS_t, C_y_t;
    GT R_t;

    Big z_y, z_r_y, z_r_yS, z_yS;
};

Poly build_from_roots(vector<ZZn>& roots, Big modulus);
Poly recursive_bfr(size_t num, size_t pos, vector<ZZn>& roots, Big modulus);
G2 mult_poly(PFC& bn, vector<G2>& rs, Poly& f);
G1 mult_poly_hat(PFC& bn, vector<G1>& rs, Poly& f);
Big generate_prime();
vector<Big> generate_RSA_keys(int n);
//Big strongp(int n,long seed1,long seed2);
Big H2(G2& m);
Big H1(G1& m);
void trivial_pairing(G2& g, G1& g_hat, PFC& bn);
Big h1(char *string);

#endif