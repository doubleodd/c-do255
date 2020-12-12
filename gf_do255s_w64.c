/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined macros CURVE and CN()
 */

#define MQ   3957
static const struct do255_int256_w64 GF_INVT508 = {
	0xC7E0DEC400D7BDB6,
	0xCCABD4771F6FB10F,
	0x940F23A06B74BE6E,
	0x1C45852F33548365
};
#include "gf_w64.c"
