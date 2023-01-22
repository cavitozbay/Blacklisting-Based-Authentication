#ifndef QUEUESIG_H
#define QUEUESIG_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <big.h>
#include "helpers.hpp"



class QueueSig
{
private:
    /* data */
public:
    int lN = 2048;
    int lm = 256;
    int lp = 128; 
    int lpp = 128;
    int le = 128;
    int ls = 128;
    Big two;
    Big mrange;
    Big mtrange;
    Big ntrange;
    Big crange;

    size_t K;
    Big N, p, q;
    Big sk_sig;
    Big b, c;
    vector<Big> gs;
    QueueSig(size_t _K);
    void key_gen();
    ReqSig request_sig(deque<Big> Q, Big& r);
    Sig sign(ReqSig& resp);
    void finalize(Sig& sig, Big r);
    BOOL verify(Sig sig, deque<Big> Q);
};


#endif