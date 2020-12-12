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
 * for curve do255e, on ARM platforms.
 */

/* see do255.h */
int
CN(verify_helper_vartime)(const void *k0,
	const CN(point) *P, const void *k1, const void *R_enc)
{
	CN(point_affine_xu) win_lo[8], Qa;
	union {
		CN(point) j;
		CN(point_xu) f;
	} M;
	gf Rw;
	const gf *endof, *minus_endof;
	uint8_t sd0[64], sd1_lo[32], sd1_hi[32];
	uint32_t sg;
	int i, j;

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
		CN(mulgen)(&M.j, k0);
		if (gf_iszero(&M.j.Z.w32)) {
			return gf_iszero(&Rw);
		} else {
			gf_mul(&Rw, &Rw, &M.j.Z.w32);
			return gf_eq(&Rw, &M.j.W.w32);
		}
	}

	/*
	 * Recode first input scalar.
	 */
	recode4(sd0, k0);

	/*
	 * Split second input scalar into two recoded half-size scalars.
	 */
	sg = split_recode4_scalar(sd1_lo, sd1_hi, k1);

	/*
	 * Fill M with a copy of P or -P, depending on the sign of k0,
	 * then use it to initialize the window.
	 */
	M.j.X = P->X;
	gf_condneg(&M.j.W.w32, &P->W.w32, sg & 1);
	M.j.Z = P->Z;
	window_fill_8_to_xu_affine(win_lo, &M.j);

	/*
	 * Get the proper factor for the endomorphism, depending on
	 * whether k0 and k1 have the same sign or opposite signs.
	 *
	 * Same sign: factor is eta
	 * Different signs: factor is -eta
	 */
	if ((sg & 1) == (sg >> 1)) {
		endof = &ETA;
		minus_endof = &MINUS_ETA;
	} else {
		endof = &MINUS_ETA;
		minus_endof = &ETA;
	}

	/*
	 * We do NOT apply the endomorphism on all points of win0 in
	 * order to compute a window for the second point: for the ARM
	 * implementation, we suppose that RAM is a scarce resource (as
	 * is normally the case on microcontrollers) and we instead
	 * dynamically apply the endomorphism for each relevant lookup.
	 * Extra runtime cost is low (32 field multiplications, about
	 * 48000 cycles).
	 */

	/*
	 * Accumulator starts at the neutral point.
	 */
	M.f.X.w32 = GF_ZERO;
	M.f.Z.w32 = GF_ONE;
	M.f.U.w32 = GF_ZERO;
	M.f.T.w32 = GF_ONE;

	for (i = 31; i >= 0; i --) {
		if (i != 31) {
			CN(double_x_xu)(&M.f, &M.f, 4);
		}

		if ((sd0[i] & 15) != 0) {
			if (sd0[i] < 0x80) {
				CN(add_mixed_xu)(&M.f, &M.f,
					&window_G_xu[sd0[i] - 1]);
			} else {
				j = (sd0[i] & 15) - 1;
				Qa.X = window_G_xu[j].X;
				gf_neg(&Qa.U.w32, &window_G_xu[j].U.w32);
				CN(add_mixed_xu)(&M.f, &M.f, &Qa);
			}
		}
		if ((sd0[32 + i] & 15) != 0) {
			if (sd0[32 + i] < 0x80) {
				CN(add_mixed_xu)(&M.f, &M.f,
					&window_G128_xu[sd0[32 + i] - 1]);
			} else {
				j = (sd0[32 + i] & 15) - 1;
				Qa.X = window_G128_xu[j].X;
				gf_neg(&Qa.U.w32, &window_G128_xu[j].U.w32);
				CN(add_mixed_xu)(&M.f, &M.f, &Qa);
			}
		}
		if ((sd1_lo[i] & 15) != 0) {
			if (sd1_lo[i] < 0x80) {
				CN(add_mixed_xu)(&M.f, &M.f,
					&win_lo[sd1_lo[i] - 1]);
			} else {
				j = (sd1_lo[i] & 15) - 1;
				Qa.X = win_lo[j].X;
				gf_neg(&Qa.U.w32, &win_lo[j].U.w32);
				CN(add_mixed_xu)(&M.f, &M.f, &Qa);
			}
		}
		if ((sd1_hi[i] & 15) != 0) {
			j = (sd1_hi[i] & 15) - 1;
			gf_neg(&Qa.X.w32, &win_lo[j].X.w32);
			gf_mul_inline(&Qa.U.w32, &win_lo[j].U.w32,
				sd1_hi[i] < 0x80 ? endof : minus_endof);
			CN(add_mixed_xu)(&M.f, &M.f, &Qa);
		}
	}

	/*
	 * Signature is valid if and only if we get the same w
	 * coordinate as point R. We must take care to handle the
	 * neutral point properly.
	 */
	if (gf_iszero(&M.f.U.w32)) {
		return gf_iszero(&Rw);
	} else {
		gf_mul(&M.f.U.w32, &M.f.U.w32, &Rw);
		return gf_eq(&M.f.U.w32, &M.f.T.w32);
	}
}
