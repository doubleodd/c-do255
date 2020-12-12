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
 * for curve do255e. It works with any finite field implementation with
 * 64-bit limbs.
 */

/* see do255.h */
int
CN(verify_helper_vartime)(const void *k0,
	const CN(point) *P, const void *k1, const void *R_enc)
{
	/*
	 * We use the endomorphism to split scalar k1; for scalar k0,
	 * we simply use a normal split since windows for G and
	 * 2^130*G are hardcoded.
	 */
	i128 k1_lo, k1_hi;
	CN(point_xu) M;
	CN(point_xu) win_lo[16], win_hi[16];
	uint8_t sd0[52], sd1_lo[26], sd1_hi[26];
	uint64_t sg;
	int i;
	gf Rw;

	/* Decode R_enc into its w coordinate. We do not need to verify
	   that it is a valid w coordinate of a group element; we only
	   need to check that it is in the proper range. */
	if (!gf_decode(&Rw, R_enc)) {
		return 0;
	}

	/*
	 * Corner case: P is neutral. In that case, we can ignore k1
	 * and use mulgen(), which is faster.
	 */
	if (CN(is_neutral)(P)) {
		CN(point) P2;

		CN(mulgen)(&P2, k0);
		if (gf_iszero(&P2.Z.w64)) {
			return (int)gf_iszero(&Rw);
		} else {
			gf_mul(&Rw, &Rw, &P2.Z.w64);
			return (int)gf_eq(&Rw, &P2.W.w64);
		}
	}

	/*
	 * Recode scalar k0.
	 */
	recode5(sd0, k0);

	/*
	 * Split scalar k1 into two signed 128-bit values.
	 */
	split_scalar(&k1_lo, &k1_hi, k1);

	/*
	 * Recode the low half of k1, then compute the low window
	 * with P or -P, depending on the sign of k1.
	 * We use fractional (x,u) coordinates for the window. Note
	 * that we already handled the case of P = N, so all coordinates
	 * of P are non-zero here.
	 */
	sg = recode5_small(sd1_lo, &k1_lo);
	win_lo[0].X = P->X;
	gf_sqr_inline(&win_lo[0].Z.w64, &P->Z.w64);
	win_lo[0].U = P->Z;
	gf_condneg(&win_lo[0].T.w64, &P->W.w64, sg);
	CN(double_xu)(&win_lo[1], &win_lo[0]);
	for (i = 3; i <= 15; i += 2) {
		CN(add_xu)(&win_lo[i - 1], &win_lo[i - 2], &win_lo[0]);
		CN(double_xu)(&win_lo[i], &win_lo[((i + 1) >> 1) - 1]);
	}

	/*
	 * Recode the high half of k1, then compute the high window by
	 * applying the endomorphism on the points of the low window,
	 * with a change of sign if the low and high halves of k1 did
	 * not have the same sign.
	 */
	sg ^= recode5_small(sd1_hi, &k1_hi);
	for (i = 0; i < 16; i ++) {
		gf_neg(&win_hi[i].X.w64, &win_lo[i].X.w64);
		win_hi[i].Z = win_lo[i].Z;
		win_hi[i].U = win_lo[i].U;
		gf_mul_inline(&win_hi[i].T.w64, &win_lo[i].T.w64,
			sg ? &ETA : &MINUS_ETA);
	}

	/*
	 * Perform the combined point multiplications.
	 */
	M.X.w64 = GF_ZERO;
	M.Z.w64 = GF_ONE;
	M.U.w64 = GF_ZERO;
	M.T.w64 = GF_ONE;
	for (i = 25; i >= 0; i --) {
		CN(point_xu) Q;
		CN(point_affine_xu) Qa;
		int j;

		if (i != 25) {
			CN(double_x_xu)(&M, &M, 5);
		}
		if ((sd0[i] & 31) != 0) {
			if (sd0[i] < 0x80) {
				j = sd0[i] - 1;
				CN(add_mixed_xu)(&M, &M, &window_G_xu[j]);
			} else {
				j = (sd0[i] & 31) - 1;
				Qa.X = window_G_xu[j].X;
				gf_neg(&Qa.U.w64, &window_G_xu[j].U.w64);
				CN(add_mixed_xu)(&M, &M, &Qa);
			}
		}
		if ((sd0[i + 26] & 31) != 0) {
			if (sd0[i + 26] < 0x80) {
				j = sd0[i + 26] - 1;
				CN(add_mixed_xu)(&M, &M, &window_G130_xu[j]);
			} else {
				j = (sd0[i + 26] & 31) - 1;
				Qa.X = window_G130_xu[j].X;
				gf_neg(&Qa.U.w64, &window_G130_xu[j].U.w64);
				CN(add_mixed_xu)(&M, &M, &Qa);
			}
		}
		if ((sd1_lo[i] & 31) != 0) {
			if (sd1_lo[i] < 0x80) {
				j = sd1_lo[i] - 1;
				CN(add_xu)(&M, &M, &win_lo[j]);
			} else {
				j = (sd1_lo[i] & 31) - 1;
				Q.X = win_lo[j].X;
				Q.Z = win_lo[j].Z;
				Q.U = win_lo[j].U;
				gf_neg(&Q.T.w64, &win_lo[j].T.w64);
				CN(add_xu)(&M, &M, &Q);
			}
		}
		if ((sd1_hi[i] & 31) != 0) {
			if (sd1_hi[i] < 0x80) {
				j = sd1_hi[i] - 1;
				CN(add_xu)(&M, &M, &win_hi[j]);
			} else {
				j = (sd1_hi[i] & 31) - 1;
				Q.X = win_hi[j].X;
				Q.Z = win_hi[j].Z;
				Q.U = win_hi[j].U;
				gf_neg(&Q.T.w64, &win_hi[j].T.w64);
				CN(add_xu)(&M, &M, &Q);
			}
		}
	}

	/*
	 * Verify that the resulting point matches the encoded value.
	 * We can content ourselves with:
	 *  - decoding the provided R_enc into its w coordinate;
	 *  - check equality of that w with the w coordinate of M.
	 * We do not need to perform a full decode of R_enc; neither do
	 * we need to normalize M to affine coordinate. That way, we
	 * can avoid any expensive inversion, square root or Legendre
	 * symbol computation here.
	 *
	 * The w coordinate of M is M.T / M.U. The w coordinate of R
	 * is in Rw. We thus check that Rw*M.U == M.T.
	 * Value zero requires some special treatment:
	 *  - If Rw == 0, then the signature is valid if and only if
	 *    we obtain the neutral, i.e. M.U == 0.
	 *  - Otherwise, if Rw != 0, then Rw*M.U must equal to M.T.
	 *    Note that M.T is always non-zero.
	 */
	if (gf_iszero(&Rw)) {
		return (int)gf_iszero(&M.U.w64);
	} else {
		gf_mul(&Rw, &Rw, &M.U.w64);
		return (int)gf_eq(&Rw, &M.T.w64);
	}
}
