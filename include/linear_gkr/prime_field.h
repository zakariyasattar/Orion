#pragma once
#ifndef __prime_field
#define __prime_field

//#include <boost/multiprecision/cpp_int.hpp>
//#include <boost/random.hpp>
#include "infrastructure/constants.h"
#include <cassert>
#include <immintrin.h>
#include <vector>
#include <memory>
#include <cstring>
#include <atomic>
#include <mutex>

//using namespace boost::multiprecision;
//using namespace boost::random;

namespace prime_field
{
    //extern int512_t mod;
    extern bool initialized;
    //extern independent_bits_engine<mt19937, 256, cpp_int> gen;

    extern const unsigned long long mod;
    extern __m256i packed_mod, packed_mod_minus_one;

    void init();
    void init_random();
    inline unsigned long long myMod(unsigned long long x)
    {
        return (x >> 61) + (x & mod);
    }
    inline unsigned long long mymult(const unsigned long long x, const unsigned long long y)
    {
        //return a value between [0, 2PRIME) = x * y mod PRIME
        /*
        unsigned long long lo, hi;
        lo = _mulx_u64(x, y, &hi);
        return ((hi << 3) | (lo >> 61)) + (lo & PRIME);
        */
        unsigned long long hi;
        asm (
        "mov %[x_read], %%rdx;\n"
        "mulx %[y_read], %%r9, %%r10;"
        "shld $0x3, %%r9, %%r10;\n"
        "and %[mod_read], %%r9;\n"
        "add %%r10, %%r9;\n"
        "mov %%r9, %[hi_write]"
        : [hi_write]"=r"(hi)
        : [x_read]"r"(x), [y_read]"r"(y), [mod_read]"r"(mod)
        : "rdx", "r9", "r10"
        );
        return hi;
    }
    inline __m256i packed_mymult(const __m256i x, const __m256i y)
    {
        __m256i ac, ad, bc, bd;
        __m256i x_shift, y_shift;
        x_shift = _mm256_srli_epi64(x, 32);
        y_shift = _mm256_srli_epi64(y, 32);
        bd = _mm256_mul_epu32(x, y);
        ac = _mm256_mul_epu32(x_shift, y_shift);
        ad = _mm256_mul_epu32(x_shift, y);
        bc = _mm256_mul_epu32(x, y_shift);

        __m256i ad_bc = _mm256_add_epi64(ad, bc);
        __m256i bd_srl32 = _mm256_srli_epi64(bd, 32);
        __m256i ad_bc_srl32 = _mm256_srli_epi64(_mm256_add_epi64(ad_bc, bd_srl32), 32);
        __m256i ad_bc_sll32 = _mm256_slli_epi64(ad_bc, 32);
        __m256i hi = _mm256_add_epi64(ac, ad_bc_srl32);

        __m256i lo = _mm256_add_epi64(bd, ad_bc_sll32);


        //return ((hi << 3) | (lo >> 61)) + (lo & PRIME);
        return _mm256_add_epi64(_mm256_or_si256(_mm256_slli_epi64(hi, 3), _mm256_srli_epi64(lo, 61)), _mm256_and_si256(lo, packed_mod));
    }

    inline __m256i packed_myMod(const __m256i x)
    {
        //return (x >> 61) + (x & mod);
        __m256i srl64 = _mm256_srli_epi64(x, 61);
        __m256i and64 = _mm256_and_si256(x, packed_mod);
        return _mm256_add_epi64(srl64, and64);
    }

    using ull       = unsigned long long;
    using uint128_t = unsigned __int128;

    /*
    This defines a field
    */
    class field_element
    {
    private:
    public:
        std::atomic<unsigned long long> img, real, seq_ctr = 0;

        // this function will be problematic for use with atomics... but I think its
        // used only for debugging
        std::unique_ptr<char[]> bit_stream()
        {
            char* p = new char[sizeof(field_element)];
            memcpy(p, this, sizeof(field_element));
            return std::unique_ptr<char[]>(p);
        }

        /* 
            Must implement rule of 5 in order to facilitate moves on atomics
            and to maintain original functionality for some operators
        */

        // Copy constructor
        field_element(const field_element& other) 
            : img(other.img.load()), 
            real(other.real.load()), 
            seq_ctr(other.seq_ctr.load()) {}

        // Move constructor
        field_element(field_element&& other) 
            : img(other.img.load()), 
            real(other.real.load()), 
            seq_ctr(other.seq_ctr.load()) {}

        // Copy assignment
        field_element& operator=(const field_element& other) {
            if (this != &other) {
                img.store(other.img.load());
                real.store(other.real.load());
                seq_ctr.store(other.seq_ctr.load());
            }
            return *this;
        }

        // Move assignment
        field_element& operator=(field_element&& other) {
            if (this != &other) {
                img.store(other.img.load());
                real.store(other.real.load());
                seq_ctr.store(other.seq_ctr.load());
            }
            return *this;
        }

        ull get_real() const { return real.load(); }
        void set_real(ull newReal) { real.store(newReal); }

        ull get_img() const { return img.load(); }
        void set_img(ull newImg) { img.store(newImg); }

        int size() {return sizeof(field_element);}

        inline field_element(){
            real.store(0);
            img.store(0);
        }
        inline field_element(const unsigned long long x){
            real.store(x % mod);
            img.store(0);
        }

        inline field_element operator + (const field_element &b) const
        {
            ull b_img  = b.img.load();
            ull b_real = b.real.load();
            ull a_img  = this->img.load();
            ull a_real = this->real.load();

            field_element ret;
            ret.img.store(b_img + a_img);
            ret.real.store(b_real + a_real);
            
            if(mod <= ret.img.load())
                ret.img.store(ret.img.load() - mod);
            if(mod <= ret.real.load())
                ret.real.store(ret.real.load() - mod);
            return ret;
        }

        inline field_element operator * (const field_element &b) const
        {
            field_element ret;
            ull all_prod = mymult(img.load() + real.load(), b.img.load() + b.real.load()); //at most 6 * mod
            //unsigned long long ac, bd;
            //mymult_2vec(real, b.real, img, b.img, ac, bd);
            ull ac = mymult(real.load(), b.real.load()), bd = mymult(img.load(), b.img.load()); //at most 1.x * mod
            auto nac = ac;
            if(bd >= mod)
                bd -= mod;
            if(nac >= mod)
                nac -= mod;
            nac ^= mod; //negate
            bd ^= mod; //negate

            ull t_img = all_prod + nac + bd; //at most 8 * mod
            t_img = myMod(t_img);
            
            if(t_img >= mod) { t_img -= mod; }
            
            ret.img.store(t_img);
            
            auto t_real = ac + bd;
            while(t_real >= mod) { t_real -= mod; }

            ret.real.store(t_real);

            return ret;
        }
        inline field_element operator - (const field_element &b) const
        {
            field_element ret;
            
            ull tmp_r = b.real.load() ^ mod; //tmp_r == -b.real is true in this prime field
            ull tmp_i = b.img.load() ^ mod; //same as above
            
            ret.real.store(real.load() + tmp_r);
            ret.img.store(img.load() + tmp_i);

            if(ret.real.load() >= mod)
                ret.real.store(ret.real.load() - mod);
            if(ret.img.load() >= mod)
                ret.img.store(ret.img.load() - mod);

            return ret;
        }
        inline field_element operator - () const
        {
            field_element ret;
            
            ret.real.store((mod - real.load()) % mod); // do modular in case real = 0
            ret.img.store((mod - img.load()) % mod);

            return ret;
        }

        bool operator == (const field_element &b) const;
        bool operator != (const field_element &b) const;
    };

    

    // class alignas(16) field_element_atomic {
    // public:
    //     std::atomic<uint128_t> packed_val;

    //     inline field_element_atomic(){
    //         packed_val.store(0);
    //     }
    //     inline field_element_atomic(const unsigned long long x){
    //         ull real = x % mod, img = 0;

    //         packed_val.store(pack(img, real));
    //     }

    //     static inline uint128_t pack(const ull img, const ull real) {
    //         return (static_cast<uint128_t>(img) << 64 | real);
    //     }

    //     static inline std::pair<ull, ull> unpack(const std::atomic<uint128_t>& packed_val) {
    //         ull img = static_cast<ull>(packed_val), real = static_cast<ull>(packed_val >> 64);

    //         return std::make_pair(img, real);
    //     }

    //     // so atomic_add will add another node to this node
    //     // I dont think I need to make * atomic?
    //     // * is only used with applying the weight to the current neighbor
    //     // which is not being written to anywhere in shared memory, and weights dont change
    //     inline void atomic_add(const field_element_atomic& b) {
    //         // load the value from the atomic
    //         uint128_t a_packed_val = packed_val.load(std::memory_order_acquire);
    //         uint128_t b_packed_val = b.packed_val.load(std::memory_order_acquire); // b is const

    //         while(true) {
    //             // populate new_packed_val with the addition between a and b
    //             uint128_t new_packed_val_with_addition;
    //             add(&new_packed_val_with_addition, a_packed_val, b.packed_val);

    //             // compare_exchange_weak will update the value if comp is successful
    //             if(packed_val.compare_exchange_weak(
    //                 a_packed_val,
    //                 new_packed_val_with_addition,
    //                 std::memory_order_release,
    //                 std::memory_order_acquire
    //             )) {
    //                 break;
    //             }
    //         }
    //     }

    //     inline void add(uint128_t* result_p, const uint128_t& a_packed_val, const uint128_t& b_packed_val) {
    //         ull res_img, res_real;

    //         // unpack this node and other node into ull
    //         auto [a_img, a_real] = unpack(a_packed_val);
    //         auto [b_img, b_real] = unpack(b_packed_val);

    //         res_img = b_img + a_img;
    //         res_real = b_real + a_real;

    //         if(mod <= res_img) { res_img = res_img - mod; }
    //         if(mod <= res_real) { res_real = res_real - mod; }

    //         *result_p = pack(res_img, res_real);
    //     }

    //    // Multiply atomic with regular field_element (for weights)
    //     inline field_element operator*(const field_element& b) const {
    //         auto [a_img, a_real] = unpack(packed_val.load(std::memory_order_acquire));
    //         field_element a_elem;
    //         a_elem.img = a_img;
    //         a_elem.real = a_real;
    //         return a_elem * b;  // Returns regular field_element
    //     }

    // };

    class field_element_packed
    {
    public:
        __m256i img, real;

        inline field_element_packed()
        {
            real = _mm256_set_epi64x(0, 0, 0, 0);
            img = _mm256_set_epi64x(0, 0, 0, 0);
        }

        inline field_element_packed(const field_element &x0, const field_element &x1, const field_element &x2, const field_element &x3)
        {
            real = _mm256_set_epi64x(x3.real, x2.real, x1.real, x0.real);
            img = _mm256_set_epi64x(x3.img, x2.img, x1.img, x0.img);
        }

        inline field_element_packed operator + (const field_element_packed &b) const
        {
            field_element_packed ret;
            ret.img = b.img + img;
            ret.real = b.real + real;
            __m256i msk0, msk1;
            msk0 = _mm256_cmpgt_epi64(ret.img, packed_mod_minus_one);
            msk1 = _mm256_cmpgt_epi64(ret.real, packed_mod_minus_one);
            ret.img = ret.img - _mm256_and_si256(msk0, packed_mod);
            ret.real = ret.real - _mm256_and_si256(msk1,packed_mod);
            return ret;
        }
        inline field_element_packed operator * (const field_element_packed &b) const
        {
            field_element_packed ret;
            __m256i all_prod = packed_mymult(img + real, b.img + b.real); //at most 6 * mod
            __m256i ac = packed_mymult(real, b.real), bd = packed_mymult(img, b.img); //at most 1.x * mod
            __m256i nac = ac;
            __m256i msk;
            msk = _mm256_cmpgt_epi64(bd, packed_mod_minus_one);
            bd = _mm256_sub_epi64(bd, _mm256_and_si256(packed_mod, msk));

            msk = _mm256_cmpgt_epi64(nac, packed_mod_minus_one);
            nac = _mm256_sub_epi64(nac, _mm256_and_si256(packed_mod, msk));

            nac = _mm256_xor_si256(nac, packed_mod);
            bd = _mm256_xor_si256(bd, packed_mod);

            __m256i t_img = _mm256_add_epi64(_mm256_add_epi64(all_prod, nac), bd);
            t_img = packed_myMod(t_img);

            msk = _mm256_cmpgt_epi64(t_img, packed_mod_minus_one);
            t_img = _mm256_sub_epi64(t_img, _mm256_and_si256(packed_mod, msk));

            ret.img = t_img;
            __m256i t_real = _mm256_add_epi64(ac, bd);
            while(1)
            {
                msk = _mm256_cmpgt_epi64(t_real, packed_mod_minus_one);
                int res = _mm256_testz_si256(msk, msk);
                if(res)
                    break;
                t_real = _mm256_sub_epi64(t_real, _mm256_and_si256(packed_mod, msk));
            }

            ret.real = t_real;
            return ret;
        }
        inline field_element_packed operator - (const field_element_packed &b) const
        {
            field_element_packed ret;
            __m256i tmp_r = b.real ^ packed_mod; //tmp_r == -b.real is true in this prime field
            __m256i tmp_i = b.img ^ packed_mod; //same as above
            ret.real = real + tmp_r;
            ret.img = img + tmp_i;
            __m256i msk0, msk1;
            msk0 = _mm256_cmpgt_epi64(ret.real, packed_mod_minus_one);
            msk1 = _mm256_cmpgt_epi64(ret.img, packed_mod_minus_one);

            ret.real = ret.real - _mm256_and_si256(msk0, packed_mod);
            ret.img = ret.img - _mm256_and_si256(msk1, packed_mod);

            return ret;
        }
        __mmask8 operator == (const field_element_packed &b) const;
        __mmask8 operator != (const field_element_packed &b) const;
        inline void get_field_element(field_element *dst) const
        {
            static unsigned long long real_arr[packed_size], img_arr[packed_size];
            _mm256_store_si256((__m256i*)real_arr, real);
            _mm256_store_si256((__m256i*)img_arr, img);
            for(int i = 0; i < 4; ++i)
            {
                dst[i].real = real_arr[i];
                dst[i].img = img_arr[i];
            }
        }
    };


    const int __max_order = 62;
    field_element get_root_of_unity(int order); //return a root of unity with order 2^[order]
    field_element random_real_only();
    field_element random();
    field_element fast_pow(field_element x, __uint128_t p);
    field_element inv(field_element x);
    double self_speed_test_mult(int repeat);
    double self_speed_test_add(int repeat);
}
#endif