// file to verify correctness of my linear code encoder aginst original

#include "linear_code/linear_code_encode.h"
#include "linear_code/linear_code_encode_opt.h"

#include <iostream>

extern bool __encode_initialized;

bool compareOutputs(
    const prime_field::field_element* a, 
    const prime_field::field_element* b, 
    int                               len
) {
    for(int i = 0; i < len; ++i) {
        if(a[i] != b[i]) {
            std::cout << "Mismatch at index " << i 
                        << ": expected (" << b[i].real << "," << b[i].img 
                        << ") got (" << a[i].real << "," << a[i].img << ")" << std::endl;
            return false;
        }
    }
    return true;
}

int main() {
    prime_field::init();
    
    int N = 16384;
    expander_init(N);
    
    int buffer_size = N * 3;
    prime_field::field_element *coefs_p    = new prime_field::field_element[N];
    prime_field::field_element *og_dest_p  = new prime_field::field_element[buffer_size];
    prime_field::field_element *opt_dest_p = new prime_field::field_element[buffer_size];
    
    for(int i = 0; i < N; ++i)
        coefs_p[i] = prime_field::random();
    
    __encode_initialized = false;
    int original_result_length = Original::encode(coefs_p, og_dest_p, N);
    
    __encode_initialized = false;
    int optimized_result_length = Optimized::encode(coefs_p, opt_dest_p, N);
    
    bool isCorrect = (original_result_length == optimized_result_length) &&
                     compareOutputs(og_dest_p, opt_dest_p, original_result_length);
    
    if(isCorrect) {
        std::cout << "Success! Optimized version matches original" << std::endl;
    } else {
        std::cout << "Error: Results differ" << std::endl;
    }
    
    return 0;
}