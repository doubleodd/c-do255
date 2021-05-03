/*
 * Curve: do255s
 * Point format: 32-bit limbs
 * Code for ARM Cortex M0+ (with assembly)
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255s
#include "support.c"
#include "gf_do255s_cm4.c"
#include "sqrt_do255s_w32.c"
#include "padd_do255s_arm.c"
#include "icore_arm.c"
#include "scalar_do255s_arm.c"
#include "pmul_base_arm.c"
#include "pmul_do255s_arm.c"
#include "lagrange_do255s_arm.c"
#include "pvrfy_do255s_arm.c"
#include "pmap_do255s_w32.c"
