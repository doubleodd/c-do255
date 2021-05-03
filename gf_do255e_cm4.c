/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined macros CURVE and CN()
 */

#define MQ   18651
static const struct do255_int256_w32 GF_INVT508 = { {
	0x2BE1CF5D, 0xD40D5B5D, 0x7282DD51, 0x3B357398,
	0x00EED6AE, 0x3ECCB228, 0x8E8FAC0B, 0x44F35C55
} };
#include "gf_arm.c"
