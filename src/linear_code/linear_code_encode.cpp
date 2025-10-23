#include "linear_code/linear_code_encode.h"

prime_field::field_element *scratch[2][100];
prime_field::field_element_optimized *scratch_opt[2][100];

bool __encode_initialized_og  = false;
bool __encode_initialized_opt = false;