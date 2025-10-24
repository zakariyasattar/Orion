// file to verify correctness of my linear code encoder aginst original

#include "linear_code/linear_code_encode.h"
#include "linear_code/linear_code_encode_opt.h"

#include <chrono>
#include <iostream>

extern bool __encode_initialized;

bool compareOutputs(
    const prime_field::field_element* a, 
    const prime_field::field_element_optimized* b, 
    int                               len
) {
    for(int i = 0; i < len; ++i) {
        if(a[i] != b[i]) {
            std::cout << "Mismatch at index " << i 
                        << ": expected (" << b[i].get_real() << "," << b[i].get_img() 
                        << ") got (" << a[i].get_real() << "," << a[i].get_img() << ")" << std::endl;
            return false;
        }
    }
    return true;
}

using uint128_t = unsigned __int128;

int main() {
    bool lock_free = std::atomic<uint128_t>::is_always_lock_free;
    std::cout << std::boolalpha << "is uint128_t lock free? " << lock_free << std::endl;

    prime_field::init();
    
    int N = (1 << 13);
    expander_init(N);
    
    int buffer_size = N * 3;
    prime_field::field_element *og_coefs_p = new prime_field::field_element[N];
    prime_field::field_element *og_dest_p  = new prime_field::field_element[buffer_size];

    prime_field::field_element_optimized *opt_coefs_p = new prime_field::field_element_optimized[N];
    prime_field::field_element_optimized *opt_dest_p  = new prime_field::field_element_optimized[buffer_size];
    
    for(int i = 0; i < N; ++i) {
        og_coefs_p[i] = prime_field::random();
        
        // Copy the same values to the optimized version
        opt_coefs_p[i].real.store(og_coefs_p[i].real, std::memory_order_relaxed);
        opt_coefs_p[i].img.store(og_coefs_p[i].img, std::memory_order_relaxed);
    }

    __encode_initialized_og = false;
    
    auto og_start = std::chrono::high_resolution_clock::now();
    int original_result_length = Original::encode(og_coefs_p, og_dest_p, N);
    auto og_end = std::chrono::high_resolution_clock::now();
    
    __encode_initialized_opt = false;

    auto opt_start = std::chrono::high_resolution_clock::now();
    int optimized_result_length = Optimized::encode(opt_coefs_p, opt_dest_p, N);
    auto opt_end = std::chrono::high_resolution_clock::now();
    
    auto og_duration  = std::chrono::duration_cast<std::chrono::milliseconds>(og_end - og_start);
    auto opt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(opt_end - opt_start);

    bool isCorrect = (original_result_length == optimized_result_length) &&
                     compareOutputs(og_dest_p, opt_dest_p, original_result_length);
    
    if(isCorrect) {
        std::cout << "Success! Optimized version matches original" << std::endl;
        std::cout << "Original implementation took " << og_duration.count() << "ms \n";
        std::cout << "Optimized implementation took " << opt_duration.count() << "ms \n"; 
    } else {
        std::cout << "Error: Results differ" << std::endl;
    }
    
    return 0;
}