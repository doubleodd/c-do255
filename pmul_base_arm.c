/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - included <string.h> (for memset())
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mulgen() with a 4-bit window, and support
 * functions for window creation and lookups. It is specific to ARM
 * platforms.
 */

/*
 * Lookup an affine point among 8 values (constant-time).
 * Lookup index is between 0 and 8 (inclusive). The provided array
 * is supposed to hold 1*Q, 2*Q,... 8*Q, in that order, for some
 * point Q. If the index is 0, this returns the neutral; otherwise,
 * this returns index*Q.
 * (implemented in assembly)
 */
void CN(window_lookup_8_affine)(CN(point_affine) *P,
	const CN(point_affine) *win, size_t index);
#define window_lookup_8_affine   CN(window_lookup_8_affine)

/*
 * Scalar recoding with a 4-bit window: for a 256-bit scalar s, this
 * function computes a 64-digit integer such that all digits are in
 * the -7..+8 range, and s = \sum_i d[i]*2^(4*i). The top digit
 * (index 63) is always nonnegative (i.e. in the 0..+8 range).
 *
 * Digits are encoded in sign+mantissa format: sign bit is bit 7 in the
 * byte (1 for negative, 0 for positive).
 */
static void
recode4(uint8_t *sd, const void *scalar)
{
	i256 w;
	int i;
	unsigned db;

	/*
	 * Decode and partially reduce the scalar to ensure that it
	 * is at most 2^255-1, to prevent an extra carry at the end of
	 * decoding.
	 */
	i256_decode(&w, scalar);
	modr_reduce256_partial(&w, &w, 0);
	modr_reduce256_finish(&w, &w);

	/*
	 * Convert words one by one. Each word yields exactly eight digits.
	 */
	db = 0;
	for (i = 0; i < 8; i ++) {
		uint32_t t;
		int j;

		t = w.v[i];
		for (j = 0; j < 8; j ++, t >>= 4) {
			unsigned b, mb;

			b = (t & 0x0F) + db;
			mb = (8 - b) >> 8;
			b ^= mb & (b ^ ((16 - b) | 0x80));
			db = mb & 1;
			sd[i * 8 + j] = (uint8_t)b;
		}
	}
}
