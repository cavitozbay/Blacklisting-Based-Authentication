#ifndef TLP_H
#define TLP_H

#include <cstring>
#include <fstream>
#include "big.h"
#include "helpers.hpp"

#define CHL_SPC 128
#define ZK_QL 128
#define RHO_SIZE 254

class TLP
{
private:
    Big p, q, m, phi_n, phi_n2;
public:
    Big T, g, n, n2, h, hn, gc, hc, rho;
    Big rp_range, mp_range, sp_range;
    Big g_inv, hn_inv, gc_inv, hc_inv;

    TLP(Big& _T);
    TLP();
    PuzzleAux pgen(Big& s);
    Big evaluate(Puzzle& puzzle);
    ZKPTLP prove_wellformedness(PuzzleAux& paux, BOOL explicit_mp = FALSE, Big mp = Big(0));
    BOOL verify_wellformedness(Puzzle& puzzle, ZKPTLP& zkp);
    ~TLP();
};


#endif