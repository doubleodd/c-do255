/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined curve multiplication core operations
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mul() with a 5-bit window; it works with
 * any finite field implementation with 64-bit limbs.
 */

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point) P;
	CN(point_affine) win[16];
	CN(point_affine) Qa;
	int i;
	uint8_t sd[52];
	uint64_t qz;

	/*
	 * Recode the scalar.
	 */
	recode5(sd, scalar);

	/*
	 * Fill the window (normalized to affine coordinates).
	 */
	window_fill_16_affine(win, P1);

	/*
	 * Top digit is nonnegative, but it can be 0; also, the source
	 * point may be the neutral.
	 */
	window_lookup_16_affine(&Qa, win, sd[51]);
	qz = gf_iszero(&Qa.X.w64);
	P.X = Qa.X;
	P.W = Qa.W;
	P.W.w64.v0 |= qz;
	P.Z.w64.v0 = 1 - qz;
	P.Z.w64.v1 = 0;
	P.Z.w64.v2 = 0;
	P.Z.w64.v3 = 0;

	/*
	 * Process other digits from top to bottom. For each digit:
	 *  - multiply current value by 32 (5 successive doublings);
	 *  - lookup point from window; negate it if the digit is
	 *    negative;
	 *  - add point to current value.
	 */
	for (i = 50; i >= 0; i --) {
		CN(double_x)(&P, &P, 5);
		window_lookup_16_affine(&Qa, win, sd[i] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}

/* see do255.h */
void
CN(mulgen)(CN(point) *P3, const void *scalar)
{
	CN(point) P;
	CN(point_affine) Qa;
	int i;
	uint8_t sd[52];
	uint64_t qz;

	/*
	 * Recode the scalar.
	 */
	recode5(sd, scalar);

	/*
	 * We split the digits into four chunks of 13, corresponding to
	 * our four precomputed windows. First batch of lookups is
	 * specialized.
	 */

	/*
	 * Top digit of the full scalar is nonnegative, but it may be zero.
	 */
	window_lookup_16_affine(&Qa, window_G195, sd[51]);
	qz = gf_iszero(&Qa.X.w64);
	P.X = Qa.X;
	P.W = Qa.W;
	P.W.w64.v0 |= qz;
	P.Z.w64.v0 = 1 - qz;
	P.Z.w64.v1 = 0;
	P.Z.w64.v2 = 0;
	P.Z.w64.v3 = 0;

	/*
	 * Lookups and additions for the top digits of the three other
	 * chunks.
	 */
	window_lookup_16_affine(&Qa, window_G, sd[12] & 31);
	gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[12] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	window_lookup_16_affine(&Qa, window_G65, sd[25] & 31);
	gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[25] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	window_lookup_16_affine(&Qa, window_G130, sd[38] & 31);
	gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[38] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	for (i = 11; i >= 0; i --) {
		CN(double_x)(&P, &P, 5);

		window_lookup_16_affine(&Qa, window_G, sd[i] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_16_affine(&Qa, window_G65, sd[i + 13] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[i + 13] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_16_affine(&Qa, window_G130, sd[i + 26] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[i + 26] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_16_affine(&Qa, window_G195, sd[i + 39] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd[i + 39] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}
