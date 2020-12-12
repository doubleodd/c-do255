/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations
 *  - defined curve basic operations, additions, multiplications
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_verify_helper_vartime() with a 5-bit window,
 * for curve do255s. It works with any finite field implementation with
 * 64-bit limbs.
 */

/*
 * NAF5 recoding, producing 'num' digits out of the provided 64-bit word.
 * Output contains unprocessed bits, with carries added in.
 */
static uint64_t
recode_NAF5_word(int8_t *rc, uint64_t x, int num)
{
	int i;

	for (i = 0; i < num; i ++) {
		/*
		 * We use a branchless algorithm to avoid misprediction
		 * penalties. Use of NAF5 is inherently non-constant-time.
		 *
		 * If x is even, then next digit is a zero.
		 * Otherwise:
		 *  - if the five low bits are in the 1..15 range, then
		 *    this value is the next digit;
		 *  - otherwise, the five low bits are in 17..31, and
		 *    we subtract 32 to make it a negative digit in the
		 *    -15..-1 range; this implies an extra +32 to add to
		 *    the x word (carry).
		 *  Either way, the five low bits of x are then cleared.
		 *
		 * Since x is then even in all cases, we divide it by 2.
		 */
		uint64_t m, t, c;

		m = -(uint64_t)(x & 1);
		t = x & m & (uint64_t)31;
		c = (t & (uint64_t)16) << 1;
		x = (x - t) + c;
		rc[i] = (int8_t)((int)t - (int)c);
		x >>= 1;
	}
	return x;
}

/*
 * NAF5 recoding. Returned value is 1 on carry, 0 otherwise. A carry is
 * returned if the computed digit encode a value which is 2^n lower
 * (exactly) than the intended value, where n is the length of the scalar
 * (in bits).
 * This function is for a 256-bit scalar.
 */
static uint64_t
recode_NAF5_256(int8_t *rc, const i256 *c)
{
	uint64_t x;

	/*
	 * We need to leave a bit of room for carries and look-ahead, so we
	 * must call recode_NAF5_word() five times. We use four 52-bit chunks
	 * and one final 48-bit chunk.
	 */
	x = c->v0 & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc, x, 52);
	x += ((c->v0 >> 52) | (c->v1 << 12)) & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc + 52, x, 52);
	x += ((c->v1 >> 40) | (c->v2 << 24)) & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc + 104, x, 52);
	x += ((c->v2 >> 28) | (c->v3 << 36)) & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc + 156, x, 52);
	x += c->v3 >> 16;
	x = recode_NAF5_word(rc + 208, x, 48);

	return x;
}

/*
 * NAF5 recoding. Returned value is 1 on carry, 0 otherwise. A carry is
 * returned if the computed digit encode a value which is 2^n lower
 * (exactly) than the intended value, where n is the length of the scalar
 * (in bits).
 * This function is for a 128-bit scalar.
 */
static uint64_t
recode_NAF5_128(int8_t *rc, const i128 *c)
{
	uint64_t x;

	/*
	 * We need to leave a bit of room for carries and look-ahead, so we
	 * must call recode_NAF5_word() three times. We use two 52-bit chunks
	 * and one final 24-bit chunk.
	 */
	x = c->v0 & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc, x, 52);
	x += ((c->v0 >> 52) | (c->v1 << 12)) & 0x000FFFFFFFFFFFFF;
	x = recode_NAF5_word(rc + 52, x, 52);
	x += c->v1 >> 40;
	x = recode_NAF5_word(rc + 104, x, 24);

	return x;
}

/*
 * Fill win1[i] with (2*i+1)*P1 and win2[i] with (2*i+1)*P2, both in
 * affine coordinates, for i = 0..7.
 */
static void
window_fill_16odd_x2_affine(CN(point_affine) *win1, CN(point_affine) *win2,
	const CN(point) *P1, const CN(point) *P2)
{
	CN(point) T, U;
	gf ZZ[16], MZ[16];
	int i;

	/*
	 * Compute point multiples; we store the Z coordinates in a
	 * separate array.
	 */
	win1[0].X = P1->X;
	win1[0].W = P1->W;
	ZZ[0] = P1->Z.w64;
	T = *P1;
	CN(double)(&U, &T);
	for (i = 3; i <= 15; i += 2) {
		CN(add)(&T, &T, &U);
		win1[(i - 1) >> 1].X = T.X;
		win1[(i - 1) >> 1].W = T.W;
		ZZ[(i - 1) >> 1] = T.Z.w64;
	}

	win2[0].X = P2->X;
	win2[0].W = P2->W;
	ZZ[8] = P2->Z.w64;
	T = *P2;
	CN(double)(&U, &T);
	for (i = 3; i <= 15; i += 2) {
		CN(add)(&T, &T, &U);
		win2[(i - 1) >> 1].X = T.X;
		win2[(i - 1) >> 1].W = T.W;
		ZZ[((i - 1) >> 1) + 8] = T.Z.w64;
	}

	/*
	 * Invert all Z coordinates.
	 */
	MZ[0] = ZZ[0];
	for (i = 2; i <= 16; i ++) {
		gf_mul(&MZ[i - 1], &MZ[i - 2], &ZZ[i - 1]);
	}
	gf_inv(&MZ[15], &MZ[15]);
	for (i = 16; i >= 2; i --) {
		gf zi;

		gf_mul(&zi, &MZ[i - 1], &MZ[i - 2]);
		gf_mul(&MZ[i - 2], &MZ[i - 1], &ZZ[i - 1]);
		ZZ[i - 1] = zi;
	}
	ZZ[0] = MZ[0];

	/*
	 * Convert points to affine coordinates. We have computed the
	 * inverses of the Z coordinates in ZZ. If a source point is the
	 * neutral, then this sets all points in the corresponding window
	 * to 0, which is correct.
	 */
	for (i = 1; i <= 8; i ++) {
		gf zi2;

		gf_sqr(&zi2, &ZZ[i - 1]);
		gf_mul(&win1[i - 1].X.w64, &win1[i - 1].X.w64, &zi2);
		gf_mul(&win1[i - 1].W.w64, &win1[i - 1].W.w64, &ZZ[i - 1]);
		gf_sqr(&zi2, &ZZ[8 + i - 1]);
		gf_mul(&win2[i - 1].X.w64, &win2[i - 1].X.w64, &zi2);
		gf_mul(&win2[i - 1].W.w64, &win2[i - 1].W.w64, &ZZ[8 + i - 1]);
	}
}

/* see do255.h */
int
CN(verify_helper_vartime)(const void *k0,
	const CN(point) *P, const void *k1, const void *R_enc)
{
	/*
	 * We want to verify that k0*G + k1*P = R. To do so, we apply
	 * Lagrange's algorithm to get c0 and c1 such that k1 = c0/c1 mod r;
	 * c0 and c1 are at most 128 bits in absolute value (but they
	 * can be negative). The equation is then:
	 *    k0*c1*G + c0*P - c1*R = 0
	 * We compute k0*d mod r to get a 255-bit value, which we can
	 * split into low and high halves (we have precomputed windows
	 * for both G and 2^128*G).
	 */
	CN(point) U3, U4, T;
	uint8_t c0[17], c1[17];
	i256 k2;
	i128 k3, k4;
	int8_t sd2[256], sd3[128], sd4[128];
	CN(point_affine) win3[8], win4[8];
	int i;

	/*
	 * Decode point R. If decoding fails, then it's hopeless.
	 */
	if (!CN(decode)(&U4, R_enc)) {
		return 0;
	}

	/*
	 * Apply Lagrange's algorithm on k1 and decode the two obtained
	 * values.
	 */
	reduce_basis_vartime(c0, c1, k1);
	i128_decode(&k3, c0);
	i128_decode(&k4, c1);

	/*
	 * If c1 < 0, then we want to replace the equation:
	 *   k0*c1*G + c0*P - c1*R = 0
	 * with:
	 *   k0*(-c1)*G + c0*(-P) + (-c1)*(-R) = 0
	 * so that we may have a nonnegative scalar c1.
	 * Similarly, if c0 < 0, we negate it and negate P.
	 *
	 * We thus set points U3 and U4, and scalars k3 and k4, such
	 * that the verification equation is:
	 *   (k0*k4 mod r)*G + k3*U3 + k4*U4 = 0
	 * with the following rules:
	 *   - k3 = |c0|
	 *   - k4 = |c1|
	 *   - If (c0 < 0 and c1 < 0) or (c0 >= 0 and c1 >= 0), then:
	 *         U3 = P
	 *     else:
	 *         U3 = -P
	 *   - U4 = -R  (always)
	 *
	 * We will write k2 = k0*k4 mod r.
	 */
	if ((c0[16] ^ c1[16]) != 0) {
		CN(neg)(&U3, P);
	} else {
		U3 = *P;
	}
	if (c0[16] != 0) {
		unsigned char cc;

		cc = _subborrow_u64(0, 0, k3.v0,
			(unsigned long long *)&k3.v0);
		(void)_subborrow_u64(cc, 0, k3.v1,
			(unsigned long long *)&k3.v1);
	}
	if (c1[16] != 0) {
		unsigned char cc;

		cc = _subborrow_u64(0, 0, k4.v0,
			(unsigned long long *)&k4.v0);
		(void)_subborrow_u64(cc, 0, k4.v1,
			(unsigned long long *)&k4.v1);
	}
	CN(neg)(&U4, &U4);

	/*
	 * Multiply k0 by k4 modulo r.
	 * Reduction uses the fact that r = 2^254 + r0 with
	 * r0 = 56904135270672826811114353017034461895 ~= 2^125.42
	 */
	i256_decode(&k2, k0);
	modr_mul256x128(&k2, &k2, &k4);

	/*
	 * We now have the base points G (implicit), U3 and U4, and the
	 * scalars k2 (256 bits), k3 and k4 (128 bits each). All these
	 * scalars are unsigned. We proceed to compute windows, apply
	 * NAF5 recoding on scalars, and initialize the accumulator T.
	 */
	window_fill_16odd_x2_affine(win3, win4, &U3, &U4);
	if (recode_NAF5_256(sd2, &k2)) {
		T.X = window_odd_G128[0].X;
		T.W = window_odd_G128[0].W;
		T.Z.w64 = GF_ONE;
	} else {
		T.X.w64 = GF_ZERO;
		T.W.w64 = GF_ONE;
		T.Z.w64 = GF_ZERO;
	}
	if (recode_NAF5_128(sd3, &k3)) {
		CN(add_mixed)(&T, &T, &win3[0]);
	}
	if (recode_NAF5_128(sd4, &k4)) {
		CN(add_mixed)(&T, &T, &win4[0]);
	}

	/*
	 * Perform the combined point multiplications.
	 */
	for (i = 127; i >= 0; i --) {
		int j;
		CN(point_affine) Qa;

		CN(double)(&T, &T);
		if (sd2[i] != 0) {
			if (sd2[i] > 0) {
				j = sd2[i] - 1;
				CN(add_mixed)(&T, &T, &window_G[j]);
			} else {
				j = -sd2[i] - 1;
				Qa.X = window_G[j].X;
				gf_neg(&Qa.W.w64, &window_G[j].W.w64);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd2[i + 128] != 0) {
			if (sd2[i + 128] > 0) {
				j = sd2[i + 128] >> 1;
				CN(add_mixed)(&T, &T, &window_odd_G128[j]);
			} else {
				j = -sd2[i + 128] >> 1;
				Qa.X = window_odd_G128[j].X;
				gf_neg(&Qa.W.w64, &window_odd_G128[j].W.w64);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd3[i] != 0) {
			if (sd3[i] > 0) {
				j = sd3[i] >> 1;
				CN(add_mixed)(&T, &T, &win3[j]);
			} else {
				j = -sd3[i] >> 1;
				Qa.X = win3[j].X;
				gf_neg(&Qa.W.w64, &win3[j].W.w64);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd4[i] != 0) {
			if (sd4[i] > 0) {
				j = sd4[i] >> 1;
				CN(add_mixed)(&T, &T, &win4[j]);
			} else {
				j = -sd4[i] >> 1;
				Qa.X = win4[j].X;
				gf_neg(&Qa.W.w64, &win4[j].W.w64);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
	}

	/*
	 * The equation is fulfilled if and only if we get the neutral
	 * at this point.
	 */
	return (int)CN(is_neutral)(&T);
}
