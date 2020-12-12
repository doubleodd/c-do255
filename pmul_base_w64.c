/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined scalar integer types
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mulgen() with a 5-bit window, and support
 * functions for window creation and lookups. It works with any finite
 * field implementation with 64-bit limbs.
 */

/*
 * Fill win[i] with (i+1)*P in affine coordinates, for i = 0..7.
 */
static void
window_fill_16_affine(CN(point_affine) *win, const CN(point) *P)
{
	CN(point) T;
	gf ZZ[16], MZ[16];
	int i, j;

	/*
	 * Compute point multiples; we store the Z coordinates in a
	 * separate array.
	 */
	win[0].X = P->X;
	win[0].W = P->W;
	ZZ[0] = P->Z.w64;
	for (i = 2; i <= 16; i ++) {
		if ((i & 1) == 0) {
			j = i >> 1;
			T.X = win[j - 1].X;
			T.W = win[j - 1].W;
			T.Z.w64 = ZZ[j - 1];
			CN(double)(&T, &T);
			win[i - 1].X = T.X;
			win[i - 1].W = T.W;
			ZZ[i - 1] = T.Z.w64;
		} else {
			j = i - 1;
			T.X = win[j - 1].X;
			T.W = win[j - 1].W;
			T.Z.w64 = ZZ[j - 1];
			CN(add)(&T, &T, P);
			win[i - 1].X = T.X;
			win[i - 1].W = T.W;
			ZZ[i - 1] = T.Z.w64;
		}
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
	 * inverses of the Z coordinates in ZZ. If the initial point was
	 * the neutral, then this sets everything to zero, which is
	 * correct (in the affine point structure, W is allowed to be 0).
	 */
	for (i = 1; i <= 16; i ++) {
		gf zi2;

		gf_sqr(&zi2, &ZZ[i - 1]);
		gf_mul(&win[i - 1].X.w64, &win[i - 1].X.w64, &zi2);
		gf_mul(&win[i - 1].W.w64, &win[i - 1].W.w64, &ZZ[i - 1]);
	}
}

/*
 * Lookup an affine point among 16 values (constant-time).
 * Lookup index is between 0 and 16 (inclusive). The provided array
 * is supposed to hold 1*Q, 2*Q,... 16*Q, in that order, for some
 * point Q. If the index is 0, this returns the neutral; otherwise,
 * this returns index*Q.
 */
static void
window_lookup_16_affine(CN(point_affine) *P,
	const CN(point_affine) *win, size_t index)
{
	uint64_t mf, u;

	/*
	 * Set mf to -1 if index == 0, 0 otherwise.
	 */
	mf = (uint64_t)index - 1;
	mf = (uint64_t)(*(int64_t *)&mf >> 31);

	/*
	 * Set P to the all-zeros. This is a valid affine representation
	 * of the neutral point.
	 */
	P->X.w64.v0 = 0;
	P->X.w64.v1 = 0;
	P->X.w64.v2 = 0;
	P->X.w64.v3 = 0;
	P->W.w64.v0 = 0;
	P->W.w64.v1 = 0;
	P->W.w64.v2 = 0;
	P->W.w64.v3 = 0;

	for (u = 0; u < 16; u ++) {
		uint64_t m;

		/*
		 * m will be -1 for the first point for which index <= u+1,
		 * i.e. such that index - u - 2 < 0.
		 */
		m = (uint64_t)index - u - 2;
		m = (uint64_t)(*(int64_t *)&m >> 31);
		m &= ~mf;
		mf |= m;

		P->X.w64.v0 |= m & win[u].X.w64.v0;
		P->X.w64.v1 |= m & win[u].X.w64.v1;
		P->X.w64.v2 |= m & win[u].X.w64.v2;
		P->X.w64.v3 |= m & win[u].X.w64.v3;
		P->W.w64.v0 |= m & win[u].W.w64.v0;
		P->W.w64.v1 |= m & win[u].W.w64.v1;
		P->W.w64.v2 |= m & win[u].W.w64.v2;
		P->W.w64.v3 |= m & win[u].W.w64.v3;
	}
}

/*
 * Scalar recoding with a 5-bit window: for a 256-bit scalar s, this
 * function computes a 52-digit integer such that all digits are in
 * the -15..+16 range, and s = \sum_i d[i]*2^(5*i). The top digit
 * (index 51) is always in the 0..+2 range.
 *
 * Digits are encoded in sign+mantissa format: sign bit is bit 7 in the
 * byte (1 for negative, 0 for positive).
 */
static void
recode5(uint8_t *sd, const uint8_t *s)
{
	int i;
	unsigned cc;
	unsigned acc;
	int acc_len;

	acc = 0;
	acc_len = 0;
	cc = 0;
	for (i = 0; i < 51; i ++) {
		unsigned b, m;

		/*
		 * Get next 5 bits.
		 */
		if (acc_len < 5) {
			acc |= (unsigned)*s ++ << acc_len;
			acc_len += 8;
		}

		/*
		 * Compute digit value, unsigned, with carry.
		 * Value is in 0..32.
		 */
		b = (acc & 31) + cc;
		acc >>= 5;
		acc_len -= 5;

		/*
		 * If value is in 0..16, then it is correct as is,
		 * and next carry is 0. Otherwise, replace b with b-32
		 * (which will be negative) and set carry to 1. Since
		 * we encode in sign+mantissa, we actually need 32-b.
		 */
		m = (16 - b) >> 8;
		b ^= m & (b ^ ((32 - b) | 0x80));
		cc = m & 1;

		sd[i] = (uint8_t)b;
	}

	/*
	 * Last digit is nonnegative. 'acc' contains the top bit.
	 */
	sd[51] = (uint8_t)(acc + cc);
}
