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
 * for curve do255e. It works with any finite field implementation with
 * 32-bit limbs.
 */

/*
 * Fill win1[i] with (2*i+1)*P1 in affine coordinates, for i = 0..3.
 */
static void
window_fill_8odd_affine(CN(point_affine) *win1, const CN(point) *P1)
{
	CN(point) T1, T2;
	gf ZZ[4], MZ[4];
	int i;

	/*
	 * Compute point multiples; we store the Z coordinates in a
	 * separate array.
	 */
	win1[0].X = P1->X;
	win1[0].W = P1->W;
	ZZ[0] = P1->Z.w32;
	T1 = *P1;
	CN(double)(&T2, &T1);
	for (i = 3; i <= 7; i += 2) {
		CN(add)(&T1, &T1, &T2);
		win1[(i - 1) >> 1].X = T1.X;
		win1[(i - 1) >> 1].W = T1.W;
		ZZ[(i - 1) >> 1] = T1.Z.w32;
	}

	/*
	 * Invert all Z coordinates.
	 */
	MZ[0] = ZZ[0];
	for (i = 2; i <= 4; i ++) {
		gf_mul(&MZ[i - 1], &MZ[i - 2], &ZZ[i - 1]);
	}
	gf_inv(&MZ[3], &MZ[3]);
	for (i = 4; i >= 2; i --) {
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
	}
}

/* see do255.h */
int
CN(verify_helper_vartime)(const void *k0,
	const CN(point) *P, const void *k1, const void *R_enc)
{
	/*
	 * We use the endomorphism to split scalar k1; for scalar k0,
	 * we simply use a normal split since windows for G and
	 * 2^128*G are hardcoded.
	 */
	i128 k1_lo, k1_hi;
	i256 k0t;
	int loneg, hineg;
	CN(point) T;
	CN(point_affine) win_lo[4], win_hi[4];
	int8_t sd0[256], sd1_lo[128], sd1_hi[128];
	int i;
	gf Rw;

	/* Decode R_enc into its w coordinate. We do not need to verify
	   that it is a valid w coordinate of a group element; we only
	   need to check that it is in the proper range. */
	if (!gf_decode(&Rw, R_enc)) {
		return 0;
	}

	/*
	 * Split scalar k1 into two signed 128-bit values.
	 */
	split_scalar(&k1_lo, &k1_hi, k1);

	/*
	 * Make sure the two scalar halves are nonnegatives, but
	 * remember the original signs.
	 */
	loneg = (int)(k1_lo.v[3] >> 31);
	hineg = (int)(k1_hi.v[3] >> 31);
	if (loneg) {
		unsigned char cc;

		cc = _subborrow_u32(0, 0, k1_lo.v[0], &k1_lo.v[0]);
		cc = _subborrow_u32(cc, 0, k1_lo.v[1], &k1_lo.v[1]);
		cc = _subborrow_u32(cc, 0, k1_lo.v[2], &k1_lo.v[2]);
		(void)_subborrow_u32(cc, 0, k1_lo.v[3], &k1_lo.v[3]);
	}
	if (hineg) {
		unsigned char cc;

		cc = _subborrow_u32(0, 0, k1_hi.v[0], &k1_hi.v[0]);
		cc = _subborrow_u32(cc, 0, k1_hi.v[1], &k1_hi.v[1]);
		cc = _subborrow_u32(cc, 0, k1_hi.v[2], &k1_hi.v[2]);
		(void)_subborrow_u32(cc, 0, k1_hi.v[3], &k1_hi.v[3]);
	}

	/*
	 * Compute window for P (low scalar half). If the scalar half was
	 * negative, then we use -P instead.
	 */
	T = *P;
	if (loneg) {
		gf_neg(&T.W.w32, &T.W.w32);
	}
	window_fill_8odd_affine(win_lo, &T);

	/*
	 * Compute window for P (high scalar); we apply the endomorphism
	 * on the points from the low window. If the low and high scalar
	 * halves had different signs, then we must also apply a negation.
	 */
	for (i = 0; i < 4; i ++) {
		gf_neg(&win_hi[i].X.w32, &win_lo[i].X.w32);
		gf_mul(&win_hi[i].W.w32, &win_lo[i].W.w32,
			loneg == hineg ? &MINUS_ETA : &ETA);
	}

	/*
	 * Recode k0 and the two halves of k1 in NAF5.
	 * Recoding k0 may yield a carry, which we use to set the
	 * initial value of the accumulator point (T).
	 * Since k1_lo and k1_hi now fit on 127 bits each, their recoding
	 * cannot yield a carry.
	 */
	for (i = 0; i < 8; i ++) {
		k0t.v[i] = dec32le((const uint8_t *)k0 + 4 * i);
	}
	if (recode_NAF4_256(sd0, &k0t)) {
		T.X = window_odd_G128[0].X;
		T.W = window_odd_G128[0].W;
		T.Z.w32 = GF_ONE;
	} else {
		T.X.w32 = GF_ZERO;
		T.W.w32 = GF_ONE;
		T.Z.w32 = GF_ZERO;
	}
	recode_NAF4_128(sd1_lo, &k1_lo);
	recode_NAF4_128(sd1_hi, &k1_hi);

	/*
	 * Perform the combined point multiplications.
	 */
	for (i = 127; i >= 0; i --) {
		int j;
		CN(point_affine) Qa;

		CN(double)(&T, &T);
		if (sd0[i] != 0) {
			if (sd0[i] > 0) {
				j = sd0[i] - 1;
				CN(add_mixed)(&T, &T, &window_G[j]);
			} else {
				j = -sd0[i] - 1;
				Qa.X = window_G[j].X;
				gf_neg(&Qa.W.w32, &window_G[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd0[i + 128] != 0) {
			if (sd0[i + 128] > 0) {
				j = sd0[i + 128] >> 1;
				CN(add_mixed)(&T, &T, &window_odd_G128[j]);
			} else {
				j = -sd0[i + 128] >> 1;
				Qa.X = window_odd_G128[j].X;
				gf_neg(&Qa.W.w32, &window_odd_G128[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd1_lo[i] != 0) {
			if (sd1_lo[i] > 0) {
				j = sd1_lo[i] >> 1;
				CN(add_mixed)(&T, &T, &win_lo[j]);
			} else {
				j = -sd1_lo[i] >> 1;
				Qa.X = win_lo[j].X;
				gf_neg(&Qa.W.w32, &win_lo[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		if (sd1_hi[i] != 0) {
			if (sd1_hi[i] > 0) {
				j = sd1_hi[i] >> 1;
				CN(add_mixed)(&T, &T, &win_hi[j]);
			} else {
				j = -sd1_hi[i] >> 1;
				Qa.X = win_hi[j].X;
				gf_neg(&Qa.W.w32, &win_hi[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
	}

	/*
	 * Verify that the resulting point matches the encoded value.
	 * We can content ourselves with:
	 *  - decoding the provided R_enc into its w coordinate;
	 *  - check equality of that w with the w coordinate of T.
	 * We do not need to perform a full decode of R_enc; neither do
	 * we need to normalize T to affine coordinate. That way, we
	 * can avoid any expensive inversion, square root or Legendre
	 * symbol computation here.
	 *
	 * The w coordinate of T is T.W / T.Z. The w coordinate of R
	 * is in Rw. We thus check that Rw*T.Z == T.W.
	 * Value zero requires some special treatment:
	 *  - If Rw == 0, then the signature is valid if and only if
	 *    we obtain the neutral, i.e. T.Z == 0.
	 *  - Otherwise, if Rw != 0, then T.Z must be non-zero, and
	 *    Rw*T.Z must be equal to T.W.
	 */
	if (gf_iszero(&Rw)) {
		return gf_iszero(&T.Z.w32);
	} else {
		gf_mul(&Rw, &Rw, &T.Z.w32);
		return !gf_iszero(&T.Z.w32) && gf_eq(&Rw, &T.W.w32);
	}
}
