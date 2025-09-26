#pragma once
#include "expanders.h"
#include <iostream>

// 2d array of pointers to field elements
// statically allocated memory pool for FFT operations
extern prime_field::field_element *scratch[2][100];

extern bool __encode_initialized;

// main bottleneck of program
// perf report --stdio --source
// to see breakdown

// this function at a high level performs recursive expander encoding algorithms

// field elements are representations of elements in a qudratic extension field

/*
src_p:  the pointer to the list of field elements that we want to encode
dest_p: the pointer to the list where the encoded field elements will reside
n:      number of elements in the input array
dep:    the depth of recursion, determines which level of scratch to use, etc...

return int: length of encoded output
*/

namespace Optimized {

inline int encode(const prime_field::field_element *src, prime_field::field_element *dst, long long n, int dep = 0)
{
    std::cout << "n: " << n << std::endl;
    if(!__encode_initialized)
    {
        __encode_initialized = true;
        for(int i = 0; (n >> i) > 1; ++i)
        {
            // fix #1: precedence bug with memory init
            // was:
            // scratch[0][i] = new prime_field::field_element[2 * n >> i];
            // scratch[1][i] = new prime_field::field_element[2 * n >> i];

            scratch[0][i] = new prime_field::field_element[(2 * n) >> i];
            scratch[1][i] = new prime_field::field_element[(2 * n) >> i];
        }
    }
    if(n <= distance_threshold)
    {
        for(long long i = 0; i < n; ++i)
            dst[i] = src[i];
        return n;
    }
    for(long long i = 0; i < n; ++i)
    {
        scratch[0][dep][i] = src[i];
    }
    long long R = alpha * n;
    for(long long j = 0; j < R; ++j)
        scratch[1][dep][j] = prime_field::field_element(0ULL);
    //expander mult
    for(long long i = 0; i < n; ++i)
    {
        const prime_field::field_element &val = src[i];
        for(int d = 0; d < C[dep].degree; ++d)
        {
            int target = C[dep].neighbor[i][d];
            scratch[1][dep][target] = scratch[1][dep][target] + C[dep].weight[i][d] * val;
        }
    }
    long long L = encode(scratch[1][dep], &scratch[0][dep][n], R, dep + 1);
    assert(D[dep].L = L);
    R = D[dep].R;
    for(long long i = 0; i < R; ++i)
    {
        scratch[0][dep][n + L + i] = prime_field::field_element(0ULL);
    }
    for(long long i = 0; i < L; ++i)
    {
        prime_field::field_element &val = scratch[0][dep][n + i];
        for(int d = 0; d < D[dep].degree; ++d)
        {
            long long target = D[dep].neighbor[i][d];
            scratch[0][dep][n + L + target] = scratch[0][dep][n + L + target] + val * D[dep].weight[i][d];
        }
    }
    for(long long i = 0; i < n + L + R; ++i)
    {
        dst[i] = scratch[0][dep][i];
    }
    return n + L + R;
}
}