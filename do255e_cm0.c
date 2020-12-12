/*
 * Curve: do255e
 * Point format: 32-bit limbs
 * Code for ARM Cortex M0+ (with assembly)
 * Jacobian (x,w) formulas are used for all operations.
 */

#define CURVE   do255e
#include "support.c"
#include "gf_do255e_cm0.c"
#include "sqrt_do255e_w32.c"
#include "padd_do255e_arm.c"
#include "icore_arm.c"
#include "scalar_do255e_arm.c"
#include "pmul_base_arm.c"
#include "pmul_do255e_arm.c"
#include "pvrfy_do255e_arm.c"
#include "pmap_do255e_w32.c"
