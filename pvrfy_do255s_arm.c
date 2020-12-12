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
 * for curve do255s. It works for ARM architectures.
 */

/*
 * Fill win1[i] with (2*i+1)*P1 and win2[i] with (2*i+1)*P2, both in
 * affine coordinates, for i = 0..3.
 * (implemented in assembly)
 */
void CN(window_fill_8odd_x2_affine)(
	CN(point_affine) *win1, CN(point_affine) *win2,
	const CN(point) *P1, const CN(point) *P2);
#define window_fill_8odd_x2_affine   CN(window_fill_8odd_x2_affine)

/*
 * Negate a 128-bit signed integer.
 */
uint32_t CN(neg_i128)(i128 *x);
#define neg_i128   CN(neg_i128)

/*
 * Custom structure for a signed integer, that can be addressed as
 * bytes (160 bits) and as a i128 (low 128 bits).
 */
typedef union {
	uint8_t b[20];
	i128 i;
} i160;

/*
 * NAF4 recoding, producing 'num' digits out of the provided 32-bit word.
 * Output contains unprocessed bits, with carries added in.
 * 'num' must be even; each digit uses 4 bits in the output buffer.
 */
static uint32_t
recode_NAF4_word(int8_t *rc, uint32_t x, int num)
{
	int i;

	for (i = 0; i < num; i += 2) {
		unsigned lo, hi;

		if ((x & 1) == 0) {
			lo = 0;
		} else {
			lo = x & 15;
			x -= lo;
			if (lo >= 8) {
				x += 16;
			}
		}
		x >>= 1;
		if ((x & 1) == 0) {
			hi = 0;
		} else {
			hi = x & 15;
			x -= hi;
			if (hi >= 8) {
				x += 16;
			}
		}
		x >>= 1;

		*(uint8_t *)(rc ++) = lo | (hi << 4);
	}
	return x;
}

/*
 * NAF4 recoding. Returned value is 1 on carry, 0 otherwise. A carry is
 * returned if the computed digit encode a value which is 2^n lower
 * (exactly) than the intended value, where n is the length of the scalar
 * (in bits).
 * This function is for a 256-bit scalar.
 */
static uint32_t
recode_NAF4_256(int8_t *rc, const i256 *c)
{
	uint32_t x;

	/*
	 * We need to leave a bit of room for carries and look-ahead, so we
	 * must call recode_NAF4_word() nine times. We use nine 28-bit chunks
	 * and one final 4-bit chunk.
	 */
	x = c->v[0] & 0x0FFFFFFF;
	x = recode_NAF4_word(rc, x, 28);
	x += ((c->v[0] >> 28) | (c->v[1] << 4)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 14, x, 28);
	x += ((c->v[1] >> 24) | (c->v[2] << 8)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 28, x, 28);
	x += ((c->v[2] >> 20) | (c->v[3] << 12)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 42, x, 28);
	x += ((c->v[3] >> 16) | (c->v[4] << 16)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 56, x, 28);
	x += ((c->v[4] >> 12) | (c->v[5] << 20)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 70, x, 28);
	x += ((c->v[5] >> 8) | (c->v[6] << 24)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 84, x, 28);
	x += c->v[6] >> 4;
	x = recode_NAF4_word(rc + 98, x, 28);
	x += c->v[7] & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 112, x, 28);
	x += c->v[7] >> 28;
	x = recode_NAF4_word(rc + 126, x, 4);

	return x;
}

/*
 * NAF4 recoding. Returned value is 1 on carry, 0 otherwise. A carry is
 * returned if the computed digit encode a value which is 2^n lower
 * (exactly) than the intended value, where n is the length of the scalar
 * (in bits).
 * This function is for a 128-bit scalar.
 */
static uint32_t
recode_NAF4_128(int8_t *rc, const i128 *c)
{
	uint32_t x;

	/*
	 * We need to leave a bit of room for carries and look-ahead, so we
	 * must call recode_NAF4_word() five times. We use four 28-bit chunks
	 * and one final 16-bit chunk.
	 */
	x = c->v[0] & 0x0FFFFFFF;
	x = recode_NAF4_word(rc, x, 28);
	x += ((c->v[0] >> 28) | (c->v[1] << 4)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 14, x, 28);
	x += ((c->v[1] >> 24) | (c->v[2] << 8)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 28, x, 28);
	x += ((c->v[2] >> 20) | (c->v[3] << 12)) & 0x0FFFFFFF;
	x = recode_NAF4_word(rc + 42, x, 28);
	x += c->v[3] >> 16;
	x = recode_NAF4_word(rc + 56, x, 16);

	return x;
}

/*
 * First internal helper for CURVE_verify_helper_vartime().
 *
 * Windows (win3 and win4) and derived scalars (k2, k3 and k4) are
 * filled with appropriate values.
 * Returned value is 0 if the source point R_enc is invalid, 1
 * otherwise.
 */
__attribute__((noinline))
static int
vrfy_helper_1(CN(point_affine) *win3, CN(point_affine) *win4,
	i256 *k2, i160 *k3, i160 *k4,
	const void *k0, const CN(point) *P, const void *k1, const void *R_enc)
{
	CN(point) U3, U4;

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
	reduce_basis_vartime(k3->b, k4->b, k1);

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
	if ((k3->b[16] ^ k4->b[16]) != 0) {
		CN(neg)(&U3, P);
	} else {
		U3 = *P;
	}
	if (k3->b[16] != 0) {
		neg_i128(&k3->i);
	}
	if (k4->b[16] != 0) {
		neg_i128(&k4->i);
	}
	CN(neg)(&U4, &U4);

	/*
	 * Multiply k0 by k4 modulo r.
	 */
	i256_decode(k2, k0);
	modr_mul256x128(k2, k2, &k4->i);

	/*
	 * We now have the base points G (implicit), U3 and U4, and the
	 * scalars k2 (256 bits), k3 and k4 (128 bits each). All these
	 * scalars are unsigned. We proceed to compute windows.
	 */
	window_fill_8odd_x2_affine(win3, win4, &U3, &U4);

	return 1;
}

/*
 * Second internal helper for CURVE_verify_helper_vartime().
 *
 * Returned value is 1 if the verification equation is fulfilled,
 * 0 otherwise.
 */
__attribute__((noinline))
static int
vrfy_helper_2(CN(point_affine) *win3, CN(point_affine) *win4,
	const i256 *k2, const i128 *k3, const i128 *k4)
{
	CN(point) T;
	int8_t sd2[128], sd3[64], sd4[64];
	int i;

	if (recode_NAF4_256(sd2, k2)) {
		T.X = window_odd_G128[0].X;
		T.W = window_odd_G128[0].W;
		T.Z.w32 = GF_ONE;
	} else {
		T.X.w32 = GF_ZERO;
		T.W.w32 = GF_ONE;
		T.Z.w32 = GF_ZERO;
	}
	if (recode_NAF4_128(sd3, k3)) {
		CN(add_mixed)(&T, &T, &win3[0]);
	}
	if (recode_NAF4_128(sd4, k4)) {
		CN(add_mixed)(&T, &T, &win4[0]);
	}

	/*
	 * Perform the combined point multiplications.
	 */
	for (i = 127; i >= 0; i --) {
		int j, nn;
		CN(point_affine) Qa;

		CN(double)(&T, &T);
		nn = (sd2[i >> 1] >> (4 * (i & 1))) & 15;
		if (nn != 0) {
			if (nn < 8) {
				j = nn - 1;
				CN(add_mixed)(&T, &T, &window_G[j]);
			} else {
				j = 15 - nn;
				Qa.X = window_G[j].X;
				gf_neg(&Qa.W.w32, &window_G[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		nn = (sd2[(i + 128) >> 1] >> (4 * (i & 1))) & 15;
		if (nn != 0) {
			if (nn < 8) {
				j = nn >> 1;
				CN(add_mixed)(&T, &T, &window_odd_G128[j]);
			} else {
				j = (16 - nn) >> 1;
				Qa.X = window_odd_G128[j].X;
				gf_neg(&Qa.W.w32, &window_odd_G128[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		nn = (sd3[i >> 1] >> (4 * (i & 1))) & 15;
		if (nn != 0) {
			if (nn < 8) {
				j = nn >> 1;
				CN(add_mixed)(&T, &T, &win3[j]);
			} else {
				j = (16 - nn) >> 1;
				Qa.X = win3[j].X;
				gf_neg(&Qa.W.w32, &win3[j].W.w32);
				CN(add_mixed)(&T, &T, &Qa);
			}
		}
		nn = (sd4[i >> 1] >> (4 * (i & 1))) & 15;
		if (nn != 0) {
			if (nn < 8) {
				j = nn >> 1;
				CN(add_mixed)(&T, &T, &win4[j]);
			} else {
				j = (16 - nn) >> 1;
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

/* see do255.h */
int
CN(verify_helper_vartime)(const void *k0,
	const CN(point) *P, const void *k1, const void *R_enc)
{
	CN(point_affine) win3[4], win4[4];
	i256 k2;
	i160 k3, k4;

	if (!vrfy_helper_1(win3, win4, &k2, &k3, &k4, k0, P, k1, R_enc)) {
		return 0;
	}
	if (!vrfy_helper_2(win3, win4, &k2, &k3.i, &k4.i)) {
		return 0;
	}
	return 1;
}
