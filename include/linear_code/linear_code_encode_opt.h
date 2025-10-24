#pragma once
#include "expanders.h"

#include <iostream>
#include <array>
#include <atomic>
#include <chrono>

/*
    Some approaches I am considering:
        1. once n (recursively) falls below a certain limit, switch from heap memory to stack
            - concern here is that even though memory pool is 'statically allocated', it still holds
              pointers, which leads to pointer chasing and bad cache locality
        2. also want to use a memory pool (?) to avoid constant re-allocation
        2. vectorize/parallelize some of these loops
        3. if there are loops that are taking a long time, maybe multithread them against other loops
            - make sure no data dependency between loops first
*/

// 2d array of pointers to field elements
// statically allocated memory pool for FFT operations
extern prime_field::field_element_optimized *scratch_opt[2][100];

// scratch[2][100]: Double-buffered memory pool for recursive expander encoding
// At each recursion depth, scratch[0][dep] holds the working output buffer while
// scratch[1][dep] stores intermediate expander results, preventing read-write conflicts

extern bool __encode_initialized_opt;

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

#define TIME_LOOP(label, loop) { \
    auto start = std::chrono::high_resolution_clock::now(); \
    loop \
    auto end = std::chrono::high_resolution_clock::now(); \
    std::cerr << label << " took " \
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() \
         << " ms\n"; \
}

// most expensive loops:
// expander mult took 218 ms
// apply D edges took 154 ms
// total runtime : 531 ms (above two operatios take 372ms 70%)

namespace Optimized {

inline int encode(const prime_field::field_element_optimized *src, prime_field::field_element_optimized *dst, long long n, int dep = 0)
{
    if(!__encode_initialized_opt)
    {
        __encode_initialized_opt = true;
        for(int i = 0; (n >> i) > 1; ++i)
        {
            scratch_opt[0][i] = new prime_field::field_element_optimized[2 * n >> i];
            scratch_opt[1][i] = new prime_field::field_element_optimized[2 * n >> i];
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
        scratch_opt[0][dep][i] = src[i];
    }
    long long R = alpha * n;
    for(long long j = 0; j < R; ++j)
        scratch_opt[1][dep][j] = prime_field::field_element_optimized(0ULL);
    
    
    //expander mult
    #pragma omp parallel for
    for(long long i = 0; i < n; ++i) {
        const prime_field::field_element_optimized &val = src[i];
        for(int d = 0; d < C[dep].degree; ++d) {
            int target = C[dep].neighbor[i][d];

            scratch_opt[1][dep][target].atomic_add(C[dep].weight[i][d] * val);

            // #pragma omp critical
            // scratch_opt[1][dep][target] = scratch_opt[1][dep][target] + C[dep].weight[i][d] * val;
        }
    }

    std::cout << "CAS failures: " << cas_failures.load() << std::endl;

    long long L = encode(scratch_opt[1][dep], &scratch_opt[0][dep][n], R, dep + 1);
    assert((D[dep].L = L));
    R = D[dep].R;
    for(long long i = 0; i < R; ++i)
    {
        scratch_opt[0][dep][n + L + i] = prime_field::field_element_optimized(0ULL);
    }

    #pragma omp parallel for
    for(long long i = 0; i < L; ++i)
    {
        prime_field::field_element_optimized &val = scratch_opt[0][dep][n + i];
        for(int d = 0; d < D[dep].degree; ++d)
        {
            long long target = D[dep].neighbor[i][d];
            scratch_opt[0][dep][n + L + target].atomic_add(val * D[dep].weight[i][d]);

            // #pragma omp critical
            // scratch_opt[0][dep][n + L + target] = scratch_opt[0][dep][n + L + target] + val * D[dep].weight[i][d];
        }
    }

    std::cout << "CAS failures: " << cas_failures.load() << std::endl;

    for(long long i = 0; i < n + L + R; ++i)
    {
        dst[i] = scratch_opt[0][dep][i];
    }
    return n + L + R;
}

} // close namespace Optimized
