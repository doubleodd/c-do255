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
 * functions for window creation and lookups. It works with any finite
 * field implementation with 32-bit limbs. A 4-bit window is used because
 * 32-bit archs may be small microcontrollers, which usually have very
 * little available RAM.
 */

/*
 * Fill win[i] with (i+1)*P in affine coordinates, for i = 0..7.
 */
static void
window_fill_8_affine(CN(point_affine) *win, const CN(point) *P)
{
	CN(point) T;
	gf ZZ[8], MZ[8];
	int i, j;

	/*
	 * Compute point multiples; we store the Z coordinates in a
	 * separate array.
	 */
	win[0].X = P->X;
	win[0].W = P->W;
	ZZ[0] = P->Z.w32;
	for (i = 2; i <= 8; i ++) {
		if ((i & 1) == 0) {
			j = i >> 1;
			T.X = win[j - 1].X;
			T.W = win[j - 1].W;
			T.Z.w32 = ZZ[j - 1];
			CN(double)(&T, &T);
			win[i - 1].X = T.X;
			win[i - 1].W = T.W;
			ZZ[i - 1] = T.Z.w32;
		} else {
			j = i - 1;
			T.X = win[j - 1].X;
			T.W = win[j - 1].W;
			T.Z.w32 = ZZ[j - 1];
			CN(add)(&T, &T, P);
			win[i - 1].X = T.X;
			win[i - 1].W = T.W;
			ZZ[i - 1] = T.Z.w32;
		}
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
	 * inverses of the Z coordinates in ZZ. If the initial point was
	 * the neutral, then this sets everything to zero, which is
	 * correct (in the affine point structure, W is allowed to be 0).
	 */
	for (i = 1; i <= 8; i ++) {
		gf zi2;

		gf_sqr(&zi2, &ZZ[i - 1]);
		gf_mul(&win[i - 1].X.w32, &win[i - 1].X.w32, &zi2);
		gf_mul(&win[i - 1].W.w32, &win[i - 1].W.w32, &ZZ[i - 1]);
	}
}

/*
 * Lookup an affine point among 8 values (constant-time).
 * Lookup index is between 0 and 8 (inclusive). The provided array
 * is supposed to hold 1*Q, 2*Q,... 8*Q, in that order, for some
 * point Q. If the index is 0, this returns the neutral; otherwise,
 * this returns index*Q.
 */
static void
window_lookup_8_affine(CN(point_affine) *P,
	const CN(point_affine) *win, size_t index)
{
	uint32_t mf, u;

	/*
	 * Set mf to -1 if index == 0, 0 otherwise.
	 */
	mf = (uint32_t)index - 1;
	mf = (uint32_t)(*(int32_t *)&mf >> 31);

	/*
	 * Set P to the all-zeros. This is a valid affine representation
	 * of the neutral point.
	 */
	memset(P, 0, sizeof *P);

	for (u = 0; u < 8; u ++) {
		uint32_t m;
		int i;

		/*
		 * m will be -1 for the first point for which index <= u+1,
		 * i.e. such that index - u - 2 < 0.
		 */
		m = (uint32_t)index - u - 2;
		m = (uint32_t)(*(int32_t *)&m >> 31);
		m &= ~mf;
		mf |= m;

		for (i = 0; i < 8; i ++) {
			P->X.w32.v[i] |= m & win[u].X.w32.v[i];
			P->W.w32.v[i] |= m & win[u].W.w32.v[i];
		}
	}
}

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
	 * Decode and reduce the scalar modulo r to ensure that it
	 * fits on 255 bits, so that we get no extra carry at the end of
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

/* see do255.h */
void
CN(mulgen)(CN(point) *P3, const void *scalar)
{
	CN(point) P;
	CN(point_affine) Qa;
	int i;
	uint8_t sd[64];
	uint32_t qz;

	/*
	 * Recode the scalar.
	 */
	recode4(sd, scalar);

	/*
	 * We split the digits into four chunks of 16, corresponding to
	 * our four precomputed windows. First batch of lookups is
	 * specialized.
	 */

	/*
	 * Top digit of the full scalar is nonnegative; however, it can
	 * be zero, which requires a specific treatment.
	 */
	window_lookup_8_affine(&Qa, window_G192, sd[63]);
	qz = gf_iszero(&Qa.X.w32);
	P.X = Qa.X;
	P.W = Qa.W;
	memset(&P.Z.w32, 0, sizeof P.Z.w32);
	P.W.w32.v[0] |= qz;
	P.Z.w32.v[0] = 1 - qz;

	/*
	 * Lookups and additions for the top digits of the three other
	 * chunks.
	 */
	window_lookup_8_affine(&Qa, window_G, sd[15] & 15);
	gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[15] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	window_lookup_8_affine(&Qa, window_G64, sd[31] & 15);
	gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[31] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	window_lookup_8_affine(&Qa, window_G128, sd[47] & 15);
	gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[47] >> 7);
	CN(add_mixed)(&P, &P, &Qa);

	for (i = 14; i >= 0; i --) {
		CN(double_x)(&P, &P, 4);

		window_lookup_8_affine(&Qa, window_G, sd[i] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_8_affine(&Qa, window_G64, sd[i + 16] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[i + 16] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_8_affine(&Qa, window_G128, sd[i + 32] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[i + 32] >> 7);
		CN(add_mixed)(&P, &P, &Qa);

		window_lookup_8_affine(&Qa, window_G192, sd[i + 48] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd[i + 48] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}

/*
 * NAF4 recoding, producing 'num' digits out of the provided 32-bit word.
 * Output contains unprocessed bits, with carries added in.
 */
static uint32_t
recode_NAF4_word(int8_t *rc, uint32_t x, int num)
{
	int i;

	for (i = 0; i < num; i ++) {
		/*
		 * We use a branchless algorithm to avoid misprediction
		 * penalties. Use of NAF4 is inherently non-constant-time.
		 *
		 * If x is even, then next digit is a zero.
		 * Otherwise:
		 *  - if the four low bits are in the 1..7 range, then
		 *    this value is the next digit;
		 *  - otherwise, the four low bits are in 9..15, and
		 *    we subtract 16 to make it a negative digit in the
		 *    -7..-1 range; this implies an extra +16 to add to
		 *    the x word (carry).
		 *  Either way, the four low bits of x are then cleared.
		 *
		 * Since x is then even in all cases, we divide it by 2.
		 */
		uint32_t m, t, c;

		m = -(uint32_t)(x & 1);
		t = x & m & (uint32_t)15;
		c = (t & (uint32_t)8) << 1;
		x = (x - t) + c;
		rc[i] = (int8_t)((int)t - (int)c);
		x >>= 1;
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
	 * must call recode_NAF4_word() nine times. We use eight 29-bit chunks
	 * and one final 24-bit chunk.
	 */
	x = c->v[0] & 0x1FFFFFFF;
	x = recode_NAF4_word(rc, x, 29);
	x += ((c->v[0] >> 29) | (c->v[1] << 3)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 29, x, 29);
	x += ((c->v[1] >> 26) | (c->v[2] << 6)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 58, x, 29);
	x += ((c->v[2] >> 23) | (c->v[3] << 9)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 87, x, 29);
	x += ((c->v[3] >> 20) | (c->v[4] << 12)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 116, x, 29);
	x += ((c->v[4] >> 17) | (c->v[5] << 15)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 145, x, 29);
	x += ((c->v[5] >> 14) | (c->v[6] << 18)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 174, x, 29);
	x += ((c->v[6] >> 11) | (c->v[7] << 21)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 203, x, 29);
	x += c->v[7] >> 8;
	x = recode_NAF4_word(rc + 232, x, 24);

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
	 * must call recode_NAF4_word() five times. We use four 29-bit chunks
	 * and one final 12-bit chunk.
	 */
	x = c->v[0] & 0x1FFFFFFF;
	x = recode_NAF4_word(rc, x, 29);
	x += ((c->v[0] >> 29) | (c->v[1] << 3)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 29, x, 29);
	x += ((c->v[1] >> 26) | (c->v[2] << 6)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 58, x, 29);
	x += ((c->v[2] >> 23) | (c->v[3] << 9)) & 0x1FFFFFFF;
	x = recode_NAF4_word(rc + 87, x, 29);
	x += c->v[3] >> 20;
	x = recode_NAF4_word(rc + 116, x, 12);

	return x;
}
