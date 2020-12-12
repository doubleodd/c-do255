/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - included <string.h> (for memset())
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined curve multiplication core operations
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mul() with a 4-bit window; it works with
 * any finite field implementation with 32-bit limbs. A 4-bit window
 * is used because 32-bit archs may be small microcontrollers, which
 * usually have very little available RAM.
 */

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point) P;
	CN(point_affine) Qa;
	CN(point_affine) win[8];
	uint32_t qz;
	int i;
	uint8_t sd[64];

	/*
	 * Recode the scalar.
	 */
	recode4(sd, scalar);

	/*
	 * Fill the window (normalized to affine coordinates).
	 */
	window_fill_8_affine(win, P1);

	/*
	 * First lookup on top digit. Top digit is in 0..+8 range, so no
	 * conditional negation is required. If the returned point is the
	 * neutral, then we must set P.Z to 0 and P.W to a non-zero value
	 * (in non-affine representation, W is not allowed to be zero);
	 * otherwise, we set Z to 1.
	 */
	window_lookup_8_affine(&Qa, win, sd[63]);
	P.X = Qa.X;
	P.W = Qa.W;
	memset(&P.Z, 0, sizeof P.Z);
	qz = gf_iszero(&Qa.X.w32);
	P.W.w32.v[0] |= qz;
	P.Z.w32.v[0] = 1 - qz;

	/*
	 * Process other digits from top to bottom. For each digit:
	 *  - multiply current value by 16 (4 successive doublings);
	 *  - lookup point from window; negate it if the digit is
	 *    negative;
	 *  - add point to current value.
	 */
	for (i = 62; i >= 0; i --) {
		CN(double_x)(&P, &P, 4);
		window_lookup_8_affine(&Qa, win, sd[i] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}
