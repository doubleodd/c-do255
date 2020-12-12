/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations
 *  - defined curve basic operations, additions, multiplications
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_verify_helper_vartime() with a 4-bit window,
 * for curve do255s. It works with any finite field implementation with
 * 32-bit limbs.
 */

/*
 * Fill win1[i] with (2*i+1)*P1 and win2[i] with (2*i+1)*P2, both in
 * affine coordinates, for i = 0..3.
 */
static void
window_fill_8odd_x2_affine(CN(point_affine) *win1, CN(point_affine) *win2,
	const CN(point) *P1, const CN(point) *P2)
{
	CN(point) T, U;
	gf ZZ[8], MZ[8];
	int i;

	/*
	 * Compute point multiples; we store the Z coordinates in a
	 * separate array.
	 */
	win1[0].X = P1->X;
	win1[0].W = P1->W;
	ZZ[0] = P1->Z.w32;
	T = *P1;
	CN(double)(&U, &T);
	for (i = 3; i <= 7; i += 2) {
		CN(add)(&T, &T, &U);
		win1[(i - 1) >> 1].X = T.X;
		win1[(i - 1) >> 1].W = T.W;
		ZZ[(i - 1) >> 1] = T.Z.w32;
	}

	win2[0].X = P2->X;
	win2[0].W = P2->W;
	ZZ[4] = P2->Z.w32;
	T = *P2;
	CN(double)(&U, &T);
	for (i = 3; i <= 7; i += 2) {
		CN(add)(&T, &T, &U);
		win2[(i - 1) >> 1].X = T.X;
		win2[(i - 1) >> 1].W = T.W;
		ZZ[((i - 1) >> 1) + 4] = T.Z.w32;
	}

	/*
	 * Invert all Z coordinates.
	 */
	MZ[0] = ZZ[0];
	for (i = 2; i <= 8; i ++) {
		gf_mul(&MZ[i - 1], &MZ[i - 2], &ZZ[i - 1]);
	}
	gf_inv(&MZ[7], &MZ[7]);
	for (i = 8; i >= 2; i --) {
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
	for (i = 1; i <= 4; i ++) {
		gf zi2;

		gf_sqr(&zi2, &ZZ[i - 1]);
		gf_mul(&win1[i - 1].X.w32, &win1[i - 1].X.w32, &zi2);
		gf_mul(&win1[i - 1].W.w32, &win1[i - 1].W.w32, &ZZ[i - 1]);
		gf_sqr(&zi2, &ZZ[4 + i - 1]);
		gf_mul(&win2[i - 1].X.w32, &win2[i - 1].X.w32, &zi2);
		gf_mul(&win2[i - 1].W.w32, &win2[i - 1].W.w32, &ZZ[4 + i - 1]);
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
	CN(point_affine) win3[4], win4[4];
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

		cc = _subborrow_u32(0, 0, k3.v[0], &k3.v[0]);
		cc = _subborrow_u32(cc, 0, k3.v[1], &k3.v[1]);
		cc = _subborrow_u32(cc, 0, k3.v[2], &k3.v[2]);
		(void)_subborrow_u32(cc, 0, k3.v[3], &k3.v[3]);
	}
	if (c1[16] != 0) {
		unsigned char cc;

		cc = _subborrow_u32(0, 0, k4.v[0], &k4.v[0]);
		cc = _subborrow_u32(cc, 0, k4.v[1], &k4.v[1]);
		cc = _subborrow_u32(cc, 0, k4.v[2], &k4.v[2]);
		(void)_subborrow_u32(cc, 0, k4.v[3], &k4.v[3]);
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
	window_fill_8odd_x2_affine(win3, win4, &U3, &U4);
	if (recode_NAF4_256(sd2, &k2)) {
		T.X = window_odd_G128[0].X;
		T.W = window_odd_G128[0].W;
		T.Z.w32 = GF_ONE;
	} else {
		T.X.w32 = GF_ZERO;
		T.W.w32 = GF_ONE;
		T.Z.w32 = GF_ZERO;
	}
	if (recode_NAF4_128(sd3, &k3)) {
		CN(add_mixed)(&T, &T, &win3[0]);
	}
	if (recode_NAF4_128(sd4, &k4)) {
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
				gf_neg(&Qa.W.w32, &window_G[j].W.w32);
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
				gf_neg(&Qa.W.w32, &window_odd_G128[j].W.w32);
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
				gf_neg(&Qa.W.w32, &win3[j].W.w32);
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
				gf_neg(&Qa.W.w32, &win4[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
	}

	/*
	 * The equation is fulfilled if and only if we get the neutral
	 * at this point.
	 */
	return CN(is_neutral)(&T);
}
