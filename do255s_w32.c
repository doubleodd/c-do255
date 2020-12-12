/*
 * Curve: do255s
 * Point format: 32-bit limbs
 * Only portable C code. Uses 32x32->64 multiplications.
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255s
#include "support.c"
#include "gf_do255s_w32.c"
#include "sqrt_do255s_w32.c"
#include "padd_do255s_w32.c"
#include "icore_w32.c"
#include "scalar_do255s_w32.c"
#include "pmul_base_w32.c"
#include "pmul_do255s_w32.c"
#include "lagrange_do255s_w32.c"
#include "pvrfy_do255s_w32.c"
#include "pmap_do255s_w32.c"
