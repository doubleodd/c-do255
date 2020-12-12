/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined macros CURVE and CN()
 */

#define MQ   18651
static const struct do255_int256_w64 GF_INVT508 = {
	0xD40D5B5D2BE1CF5D,
	0x3B3573987282DD51,
	0x3ECCB22800EED6AE,
	0x44F35C558E8FAC0B
};
#include "gf_bmi2.c"
