#include "aho.hpp"

AHO::AHO(PFC& _bn): bn(_bn)
{

}

void AHO::key_gen(size_t _L)
{
    L = _L;
    bn.random(Gr);
    bn.random(Hr);
    bn.random(g_hat);

    mu_z = rand(bn.order());
    nu_z = rand(bn.order());
    alpha_a = rand(bn.order());
    alpha_b = rand(bn.order());
    for (size_t i = 0; i < L; i++)
    {
        mus.push_back(rand(bn.order()));
        nus.push_back(rand(bn.order()));
        Gs.push_back(bn.mult(Gr,mus[i]));
        Hs.push_back(bn.mult(Hr,nus[i]));
    }
    Gz = bn.mult(Gr, mu_z);
    Hz = bn.mult(Hr, nu_z);

    A = bn.pairing(bn.mult(g_hat, alpha_a), Gr);
    B = bn.pairing(bn.mult(g_hat, alpha_b), Hr);
}

Sigma AHO::sign(vector<G2>& M)
{
    Big beta, epsilon, eta, iota, kappa;

    beta = rand(bn.order());
    epsilon = rand(bn.order());
    eta = rand(bn.order());
    iota = rand(bn.order());
    kappa = rand(bn.order());

    Sigma sigma;

    sigma.tetas_hat[1] = bn.mult(g_hat, beta);
    sigma.tetas_hat[2] = bn.mult(g_hat, (epsilon - modmult(mu_z, beta, bn.order()))%bn.order());
    for (size_t i = 0; i < L; i++)
    {
        sigma.tetas_hat[2] = sigma.tetas_hat[2] + bn.mult(M[i], -mus[i]);
    }
    
    sigma.tetas[3] = bn.mult(Gr, eta);

    sigma.tetas_hat[4] = bn.mult(g_hat, modmult((alpha_a - epsilon), inverse(eta, bn.order()), bn.order()));

    sigma.tetas_hat[5] = bn.mult(g_hat, (iota - modmult(nu_z,beta,bn.order()))%bn.order());
    for (size_t i = 0; i < L; i++)
    {
        sigma.tetas_hat[5] = sigma.tetas_hat[5] + bn.mult(M[i], -nus[i]);
    }

    sigma.tetas[6] = bn.mult(Hr, kappa);

    sigma.tetas_hat[7] = bn.mult(g_hat, modmult((alpha_b - iota), inverse(kappa, bn.order()), bn.order()));

    return sigma;
}

BOOL AHO::verify(Sigma sigma, vector<G2> M)
{
    GT Ap = bn.pairing(sigma.tetas_hat[1], Gz) * bn.pairing(sigma.tetas_hat[2], Gr) * bn.pairing(sigma.tetas_hat[4], sigma.tetas[3]);
    for (size_t i = 0; i < L; i++)
    {
        Ap = Ap * bn.pairing(M[i], Gs[i]);
    }

    GT Bp = bn.pairing(sigma.tetas_hat[1], Hz) * bn.pairing(sigma.tetas_hat[5], Hr) * bn.pairing(sigma.tetas_hat[7], sigma.tetas[6]);
    for (size_t i = 0; i < L; i++)
    {
        Bp = Bp * bn.pairing(M[i], Hs[i]);
    }

    return (A == Ap) && (B == Bp);
}


void AHO::rndmz(G1& x, G2& y)
{
    Big gamma = rand(bn.order());
    x = bn.mult(x, gamma);
    y = bn.mult(y, inverse(gamma, bn.order()));

    if(bn.pairing(y, x) == GT())
    {
        Big b = rand(Big(2));

        if(b == Big(1))
        {
            x = G1();
            y = G2();
        }
        else
        {
            b = rand(Big(2));
            if(b == Big(1))
            {
                x = G1();
                bn.random(y);
            }
            else
            {
                bn.random(x);
                y = G2();
            }
        }
    }
}

// z, r, s, t, u, v, w
// 1, 2, 3, 4, 5, 6, 7

Sigma AHO::randomize(Sigma sigma)
{
    // cout << sigma.tetas_hat[1].g << endl;
    // cout << sigma.tetas_hat[2].g << endl;
    // cout << sigma.tetas_hat[4].g << endl;
    // cout << sigma.tetas_hat[5].g << endl;
    // cout << sigma.tetas_hat[7].g << endl;
    // cout << sigma.tetas[3].g << endl;
    // cout << sigma.tetas[6].g << endl;
    // cout << "===============================" << endl;

    // randomize 2, 3, 4
    if (sigma.tetas_hat[4].g.iszero())
    {
        sigma.tetas[3] = G1();
        bn.random(sigma.tetas_hat[4]);
    }

    Big rho = rand(bn.order());
    sigma.tetas_hat[2] = sigma.tetas_hat[2] + bn.mult(sigma.tetas_hat[4], rho);

    sigma.tetas[3] = sigma.tetas[3] + bn.mult(Gr, -rho);
    rndmz(sigma.tetas[3], sigma.tetas_hat[4]);


    // randomize 5, 6, 7
    if (sigma.tetas_hat[7].g.iszero())
    {
        sigma.tetas[6] = G1();
        bn.random(sigma.tetas_hat[7]);
    }

    rho = rand(bn.order());
    sigma.tetas_hat[5] = sigma.tetas_hat[5] + bn.mult(sigma.tetas_hat[7], rho);

    sigma.tetas[6] = sigma.tetas[6] + bn.mult(Hr, -rho);
    rndmz(sigma.tetas[6], sigma.tetas_hat[7]);


    // cout << sigma.tetas_hat[1].g << endl;
    // cout << sigma.tetas_hat[2].g << endl;
    // cout << sigma.tetas_hat[4].g << endl;
    // cout << sigma.tetas_hat[5].g << endl;
    // cout << sigma.tetas_hat[7].g << endl;
    // cout << sigma.tetas[3].g << endl;
    // cout << sigma.tetas[6].g << endl;


    return sigma;
}
