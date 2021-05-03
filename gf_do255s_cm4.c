/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined macros CURVE and CN()
 */

#define MQ   3957
static const struct do255_int256_w32 GF_INVT508 = { {
	0x00D7BDB6, 0xC7E0DEC4, 0x1F6FB10F, 0xCCABD477,
	0x6B74BE6E, 0x940F23A0, 0x33548365, 0x1C45852F
} };
#include "gf_arm.c"
