/*
 * Curve: do255e
 * Point format: 32-bit limbs
 * Only portable C code. Uses 32x32->64 multiplications.
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255e
#include "support.c"
#include "gf_do255e_w32.c"
#include "sqrt_do255e_w32.c"
#include "padd_do255e_w32.c"
#include "icore_w32.c"
#include "scalar_do255e_w32.c"
#include "pmul_base_w32.c"
#include "pmul_do255e_w32.c"
#include "pvrfy_do255e_w32.c"
#include "pmap_do255e_w32.c"
