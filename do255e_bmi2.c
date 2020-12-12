/*
 * Curve: do255e
 * Point format: 64-bit limbs
 * Uses ADX/BMI2 opcodes with inline assembly, but no (explicit) AVX2.
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255e
#include "support.c"
#include "gf_do255e_bmi2.c"
#include "sqrt_do255e_w64.c"
#include "padd_do255e_w64.c"
#include "icore_w64.c"
#include "scalar_do255e_w64.c"
#include "pmul_base_w64.c"
#include "pmul_do255e_w64.c"
#include "pvrfy_do255e_w64.c"
#include "pmap_do255e_w64.c"
