/*
 * Curve: do255s
 * Point format: 64-bit limbs
 * Uses 64-bit intrinsics (e.g. _addcarry_u64()) but no assembly or AVX2.
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255s
#include "support.c"
#include "gf_do255s_w64.c"
#include "sqrt_do255s_w64.c"
#include "padd_do255s_w64.c"
#include "icore_w64.c"
#include "scalar_do255s_w64.c"
#include "pmul_base_w64.c"
#include "pmul_do255s_w64.c"
#include "lagrange_do255s_w64.c"
#include "pvrfy_do255s_w64.c"
#include "pmap_do255s_w64.c"
