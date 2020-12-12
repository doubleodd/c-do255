/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mul() with a 5-bit window; it works with
 * any finite field implementation with 64-bit limbs.
 */

/*
 * do255e endomorphism:
 * ====================
 *
 * We use one of the cases described by Gallant, Lambert and Vanstone in
 * their 2001 article:
 *   https://www.iacr.org/archive/crypto2001/21390189.pdf
 *
 * Modulus is p = 2^255 - 18651. Curve has order 2*r, with:
 *   r = 2^254 - 131528281291764213006042413802501683931
 *
 * Let eta = sqrt(-1) mod p. There are two such roots, we use this one:
 * 7656063742463026568679823572395325799027601838558345258426535816504372595438
 *
 * Let phi(x, w) = (-x, -eta*w)
 * (In (x, u) coordinates, this is phi(x, u) = (-x, eta*u))
 * phi() is an endomorphism over our group G or order r. For any group
 * element P, phi(P) = mu*P for a constant mu which is a square root of -1
 * modulo r. With our choice of eta, we have the following mu:
 * 23076176648693837106500022901799924463072024427516564762134831823525232195341
 *
 * Let k be a scalar (integer modulo r). We decompose k into two half-size
 * values k0 and k1 such that k = k0 + k1*mu mod r.
 *
 * Since r = 1 mod 4, it can be written as the sum of two squares. Let u
 * and v such that r = u^2 + v^2; the choice is unique up to permutation and
 * change of sign; these values can be found by using Lagrange's algorithm
 * on the lattice ((mu, 1), (r, 0)). We choose the following:
 *
 *   u =  34978546233976132960203755786038370577
 *   v = 166506827525740345966246169588540045182
 *
 * Since (u/v)^2 = -1 mod r, value mu is equal to either u/v or v/u (mod r).
 * With our choices, mu = u/v mod r.
 *
 * It can be verified that:
 *   r = u^2 + v^2
 *   v + mu*u = 0 mod r
 *   -u + mu*v = 0 mod r
 *
 * Given k, we compute integers c and d as:
 *   c = round(k*v / r)
 *   d = round(k*u / r)
 * Note that c and d are nonnegative. Moreover, c and d fit on 127 bits each.
 *
 * We then have:
 *   k0 = k - d*u - c*v
 *   k1 = d*v - c*u
 *
 * It can be shown (see GLV article) that k0^2 + k1^2 <= u^2 + v^2. Since
 * u^2 + v^2 = r < 2^254, this implies that |k0| < 2^127 and |k1| < 2^127.
 * Thus, k0 and k1 (which are signed integers) can fit on 128 bits each,
 * including their sign bit.
 *
 *
 * Rounded division:
 * =================
 *
 * To compute c and d, we need to do a rounded division. We use the
 * fact that r = 2^254 - r0 with r0 < 2^127.
 *
 * Suppose that we have x = k*u or k*v, and we want to compute y = round(x/r).
 * Since r is odd, we have:
 *   y = round(x/r) = floor((x + (r-1)/2) / r)
 * Let z = x + (r-1)/2. We can split z at the 254-bit index:
 *   z = z0 + 2^254*z1
 * with 0 <= z0 < 2^254, and z1 >= 0. Since k < r and given the values of u
 * and v, the maximum value for z1 is about 2^126.97, i.e it fits on 127 bits.
 *
 * We thus have:
 *   z1*r = z1*2^254 - z1*r0 < z1*2^254 < z
 * and
 *   (z1+2)*r = z1^254 - z1*r0 + 2*2^254 - 2*r0
 *            = (z1+1)*2^254 + (2^254 - (z1+2)*r0)
 * Since (z1+2)*r0 < 2^254, we have (z1+2)*r > (z1+1)*2^254 > z.
 *
 * It follows that the rounded division result is necessarily either z1
 * or z1+1. We can thus compute that rounded division with the following
 * algorithm:
 *
 *   Input: integer x such that 0 <= x <= (r-1)*max(u,v)
 *   Output: round(x / r)
 *     1. z <- x + (r-1)/2
 *     2. y <- floor(z / 2^254) + 1
 *     3. if y*r > z, then y <- y-1
 *     4. return y
 *
 * The most expensive operation is the product y*r. However, we only need
 * the sign of z - y*r. We can do that computation as follows:
 *    z - y*r = (y-1)*2^254 + z0 - y*2^254 + y*r0
 *            = z0 + y*r0 - 2^254
 * We thus need to subtract 1 from y if and only if z0 + y*r0 is strictly
 * lower than 2^254. y and r0 both fit on 127 bits each, and z0 is less
 * than 2^254; we can thus do that computation over 255 bits.
 */

static const gf ETA = {
	0xD99E0F1BAA938AEE,
	0xA60D864FB30E6336,
	0xE414983FE53688E3,
	0x10ED2DB33C69B85F
};

static const gf MINUS_ETA = {
	0x2661F0E4556C2C37,
	0x59F279B04CF19CC9,
	0x1BEB67C01AC9771C,
	0x6F12D24CC39647A0
};

/* r */
static const i256 R = {
	0x1F52C8AE74D84525,
	0x9D0C930F54078C53,
	0xFFFFFFFFFFFFFFFF,
	0x3FFFFFFFFFFFFFFF
};

/* (r - 1)/2 */
static const i256 HR = {
	0x8FA964573A6C2292,
	0xCE864987AA03C629,
	0xFFFFFFFFFFFFFFFF,
	0x1FFFFFFFFFFFFFFF
};

/* u */
static const i128 eU = {
	0x2ACCF9DEC93F6111,
	0x1A509F7A53C2C6E6
};

/* v */
static const i128 eV = {
	0x0B7A31305466F77E,
	0x7D440C6AFFBB3A93
};

/*
 * Input:
 *   0 <= k < r
 *   e < 2^127 - 2
 * Output:
 *   d = round(k*e / r)
 */
static void
mul_divr_rounded(i128 *d, const i256 *k, const i128 *e)
{
	i384 z;
	i256 t;
	i128 y;
	unsigned long long lo;
	unsigned char cc;

	/* z <- k*e */
	mul256x128(&z, k, e);

	/* z <- z + (r-1)/2 */
	cc = _addcarry_u64(0, HR.v0, z.v0, (unsigned long long *)&z.v0);
	cc = _addcarry_u64(cc, HR.v1, z.v1, (unsigned long long *)&z.v1);
	cc = _addcarry_u64(cc, HR.v2, z.v2, (unsigned long long *)&z.v2);
	cc = _addcarry_u64(cc, HR.v3, z.v3, (unsigned long long *)&z.v3);
	cc = _addcarry_u64(cc, 0, z.v4, (unsigned long long *)&z.v4);
	(void)_addcarry_u64(cc, 0, z.v5, (unsigned long long *)&z.v5);

	/* y <- floor(z / 2^254) + 1 */
	y.v0 = (z.v3 >> 62) | (z.v4 << 2);
	y.v1 = (z.v4 >> 62) | (z.v5 << 2);
	cc = _addcarry_u64(0, y.v0, 1, (unsigned long long *)&y.v0);
	(void)_addcarry_u64(cc, y.v1, 0, (unsigned long long *)&y.v1);

	/* t <- y*r0 */
	mul128x128(&t, &y, &R0);

	/* t <- t + z0
	   We are only interested in the high limb. */
	z.v3 = z.v3 & 0x3FFFFFFFFFFFFFFF;
	cc = _addcarry_u64(0, z.v0, t.v0, &lo);
	cc = _addcarry_u64(cc, z.v1, t.v1, &lo);
	cc = _addcarry_u64(cc, z.v2, t.v2, &lo);
	(void)_addcarry_u64(cc, z.v3, t.v3, &lo);

	/* The high limb is in 'lo' and it is lower than 2^63. If
	   it is lower than 2^62, then y is too large and we must
	   decrement it; otherwise, we keep it unchanged. */
	cc = _subborrow_u64(0, y.v0, (1 - (lo >> 62)), 
		(unsigned long long *)&d->v0);
	(void)_subborrow_u64(cc, y.v1, 0,
		(unsigned long long *)&d->v1);
}

/*
 * Split scalar k (256 bits) into k0 and k1 (128 bits each, signed),
 * such that k = k0 + k1*mu mod r.
 */
static void
split_scalar(i128 *k0, i128 *k1, const uint8_t *k)
{
	i256 t;
	i128 c, d, f;
	unsigned char cc;

	/*
	 * Decode the scalar and reduce it modulo r.
	 */
	i256_decode(&t, k);
	modr_reduce256_partial(&t, &t, 0);
	modr_reduce256_finish(&t, &t);

	/*
	 * c = round(k*v / r)
	 * d = round(k*u / r)
	 */
	mul_divr_rounded(&c, &t, &eV);
	mul_divr_rounded(&d, &t, &eU);

	/*
	 * k0 = k - d*u - c*v
	 * k1 = d*v - c*u
	 */
	mul128x128trunc(&f, &d, &eU);
	cc = _subborrow_u64(0, t.v0, f.v0, (unsigned long long *)&k0->v0);
	(void)_subborrow_u64(cc, t.v1, f.v1, (unsigned long long *)&k0->v1);
	mul128x128trunc(&f, &c, &eV);
	cc = _subborrow_u64(0, k0->v0, f.v0, (unsigned long long *)&k0->v0);
	(void)_subborrow_u64(cc, k0->v1, f.v1, (unsigned long long *)&k0->v1);

	mul128x128trunc(k1, &d, &eV);
	mul128x128trunc(&f, &c, &eU);
	cc = _subborrow_u64(0, k1->v0, f.v0, (unsigned long long *)&k1->v0);
	(void)_subborrow_u64(cc, k1->v1, f.v1, (unsigned long long *)&k1->v1);
}

/*
 * Recode a signed 128-bit integer with a 5-bit window. The scalar is
 * first replaced with its absolute value x (its original sign is
 * returned: 1 if x was negative, 0 otherwise). The function then
 * computes 26 digits such that all digits are in the -15..+16 range,
 * and x = \sum_i sd[i]*2^(5*i). The top digit (sd[25]) is in the 0..+4
 * range.
 *
 * Digits are encoded in sign+mantissa format: sign bit is bit 7 in the
 * byte (1 for negative, 0 for positive).
 */
static uint64_t
recode5_small(uint8_t *sd, const i128 *s)
{
	i128 x;
	unsigned long long sg, t;
	unsigned char cc;
	unsigned db, b, m;
	int i;

	/* x = abs(s) */
	sg = s->v1 >> 63;
	cc = _addcarry_u64(0, s->v0 ^ -sg, sg, (unsigned long long *)&x.v0);
	(void)_addcarry_u64(cc, s->v1 ^ -sg, 0, (unsigned long long *)&x.v1);

	/* First 12 digits come from the low limb. */
	db = 0;
	t = x.v0;
	for (i = 0; i < 12; i ++) {
		b = ((unsigned)t & 0x1F) + db;
		m = (16 - b) >> 8;
		b ^= m & (b ^ ((32 - b) | 0x80));
		db = m & 1;
		sd[i] = (uint8_t)b;
		t >>= 5;
	}

	/* Get more bits from high limb for next 12 digits. */
	t |= x.v1 << 4;
	for (i = 12; i < 24; i ++) {
		b = ((unsigned)t & 0x1F) + db;
		m = (16 - b) >> 8;
		b ^= m & (b ^ ((32 - b) | 0x80));
		db = m & 1;
		sd[i] = (uint8_t)b;
		t >>= 5;
	}

	/* Last two digits. */
	t = x.v1 >> 56;
	b = ((unsigned)t & 0x1F) + db;
	m = (16 - b) >> 8;
	b ^= m & (b ^ ((32 - b) | 0x80));
	db = m & 1;
	sd[24] = (uint8_t)b;
	t >>= 5;
	sd[25] = (uint8_t)((unsigned)t + db);

	return sg;
}

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point_affine) win0[16], win1[16], Qa;
	CN(point) P;
	i128 k0, k1;
	uint8_t sd0[26], sd1[26];
	uint64_t sg, qz;
	int i;

	/*
	 * Split input scalar into k0 and k1.
	 */
	split_scalar(&k0, &k1, scalar);

	/*
	 * Recode k0. Fill P with a copy of P1 or -P1, depending on
	 * whether k0 is negative or not, then compute a window out
	 * of it.
	 */
	sg = recode5_small(sd0, &k0);
	P = *P1;
	gf_condneg(&P.W.w64, &P.W.w64, sg);
	window_fill_16_affine(win0, &P);

	/*
	 * Recode k1; check whether the sign of k1 differs from that of k0.
	 */
	sg ^= recode5_small(sd1, &k1);

	/*
	 * Apply endomorphism on all points in the first window (win0) to
	 * get the points of the second window (win1). If the signs of
	 * k0 and k1 differ, then points in win1 must be negated.
	 */
	for (i = 0; i < 16; i ++) {
		gf_neg(&win1[i].X.w64, &win0[i].X.w64);
		gf_mul_inline(&win1[i].W.w64, &win0[i].W.w64, &MINUS_ETA);
		gf_condneg(&win1[i].W.w64, &win1[i].W.w64, sg);
	}

	/*
	 * Lookup points corresponding to the top digits, and add them
	 * into P. Take care that if the first looked-up point is the
	 * neutral, then we get affine neutral point (0, 0) and we must
	 * then set W to 1 and Z to 0; otherwise, we use the retrieved W
	 * and set Z to 1.
	 * Top digits in sd0 and sd1 are always nonnnegative.
	 */
	window_lookup_16_affine(&Qa, win0, sd0[25]);
	qz = gf_iszero(&Qa.X.w64);
	P.X = Qa.X;
	P.W = Qa.W;
	P.W.w64.v0 |= qz;
	P.Z.w64.v0 = 1 - qz;
	P.Z.w64.v1 = 0;
	P.Z.w64.v2 = 0;
	P.Z.w64.v3 = 0;
	window_lookup_16_affine(&Qa, win1, sd1[25]);
	CN(add_mixed)(&P, &P, &Qa);

	/*
	 * Process other digits from top to bottom. For each digit:
	 *  - multiply current value by 32 (5 successive doublings);
	 *  - lookup point from first window; negate it if the digit is
	 *    negative;
	 *  - add point to current value;
	 *  - lookup point from second window; negate it if the digit is
	 *    negative;
	 *  - add point to current value.
	 */
	for (i = 24; i >= 0; i --) {
		CN(double_x)(&P, &P, 5);
		window_lookup_16_affine(&Qa, win0, sd0[i] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd0[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
		window_lookup_16_affine(&Qa, win1, sd1[i] & 31);
		gf_condneg(&Qa.W.w64, &Qa.W.w64, sd1[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}

#if 0
/*
 * A version without normalization of the window points to affine.
 * On large x86 (Coffee Lake core), this appears to be substantially
 * slower than the version with affine points (93k cycles instead
 * of 75k!). The exact reason is not known. A generic point addition
 * is 8M+6S, but a mixed addition is 8M+3S; since there are 51
 * point additions to perform per point multiplication, the use of
 * mixed addition should save up to 153S, which is lower than the
 * cost of the inversion used in the window normalization (about 180S,
 * and that's not counting the extra multiplications). Thus, the
 * generic non-affine loop below _should_ be faster. But this is not
 * the case in practive.
 *
 * One theory is that using larger points increases register pressure
 * on lookups, and crosses a threshold somewhere in the compiler
 * optimization routines.
 */

/*
 * Lookup a point among 16 values (constant-time).
 * Lookup index is between 0 and 16 (inclusive). The provided array
 * is supposed to hold 1*Q, 2*Q,... 16*Q, in that order, for some
 * point Q. If the index is 0, this returns the neutral; otherwise,
 * this returns index*Q.
 */
static void
window_lookup_16(CN(point) *P, const CN(point) *win, size_t index)
{
	uint64_t mf, u;

	/*
	 * Set mf to -1 if index == 0, 0 otherwise.
	 */
	mf = (uint64_t)index - 1;
	mf = (uint64_t)(*(int64_t *)&mf >> 31);

	/*
	 * Set P to the all-zeros, except that we set W to a non-zero
	 * value in case the index is 0 (if the index is 0, then this
	 * initial value will be untouched, and we will thus return a
	 * valid representation of the neutral).
	 */
	P->X.w64.v0 = 0;
	P->X.w64.v1 = 0;
	P->X.w64.v2 = 0;
	P->X.w64.v3 = 0;
	P->W.w64.v0 = mf;
	P->W.w64.v1 = 0;
	P->W.w64.v2 = 0;
	P->W.w64.v3 = 0;
	P->Z.w64.v0 = 0;
	P->Z.w64.v1 = 0;
	P->Z.w64.v2 = 0;
	P->Z.w64.v3 = 0;

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
		P->Z.w64.v0 |= m & win[u].Z.w64.v0;
		P->Z.w64.v1 |= m & win[u].Z.w64.v1;
		P->Z.w64.v2 |= m & win[u].Z.w64.v2;
		P->Z.w64.v3 |= m & win[u].Z.w64.v3;
	}
}

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point) win0[16], win1[16];
	CN(point) P, Q;
	i128 k0, k1;
	uint8_t sd0[26], sd1[26];
	uint64_t sg;
	int i;

	/*
	 * Split input scalar into k0 and k1.
	 */
	split_scalar(&k0, &k1, scalar);

	/*
	 * Recode k0. Fill P with a copy of P1 or -P1, depending on
	 * whether k0 is negative or not, then compute a window out
	 * of it.
	 */
	sg = recode5_small(sd0, &k0);
	P = *P1;
	gf_condneg(&P.W.w64, &P.W.w64, sg);
	win0[0] = P;
	CN(double)(&win0[1], &win0[0]);
	for (i = 3; i <= 15; i += 2) {
		CN(add)(&win0[i - 1], &win[i - 2], &P);
		CN(double)(&win0[i], &win0[((i + 1) >> 1) - 1]);
	}

	/*
	 * Recode k1; check whether the sign of k1 differs from that of k0.
	 */
	sg ^= recode5_small(sd1, &k1);

	/*
	 * Apply endomorphism on all points in the first window (win0) to
	 * get the points of the second window (win1). If the signs of
	 * k0 and k1 differ, then points in win1 must be negated.
	 * Points are negated by changing the sign of Z, which is equivalent
	 * to changing the sign of W.
	 */
	for (i = 0; i < 16; i ++) {
		gf_neg(&win1[i].X.w64, &win0[i].X.w64);
		gf_mul_inline(&win1[i].W.w64, &win0[i].W.w64, &MINUS_ETA);
		gf_condneg(&win1[i].Z.w64, &win0[i].Z.w64, sg);
	}

	/*
	 * Lookup points corresponding to the top digits, and add them
	 * into P. Take care that if the first looked-up point is the
	 * neutral, then we get affine neutral point (0, 0) and we must
	 * then set W to 1 and Z to 0; otherwise, we use the retrieved W
	 * and set Z to 1.
	 * Top digits in sd0 and sd1 are always nonnnegative.
	 */
	window_lookup_16(&P, win0, sd0[25]);
	window_lookup_16(&Q, win1, sd1[25]);
	CN(add)(&P, &P, &Q);

	/*
	 * Process other digits from top to bottom. For each digit:
	 *  - multiply current value by 32 (5 successive doublings);
	 *  - lookup point from first window; negate it if the digit is
	 *    negative;
	 *  - add point to current value;
	 *  - lookup point from second window; negate it if the digit is
	 *    negative;
	 *  - add point to current value.
	 */
	for (i = 24; i >= 0; i --) {
		CN(double_x)(&P, &P, 5);
		window_lookup_16(&Q, win0, sd0[i] & 31);
		gf_condneg(&Q.W.w64, &Q.W.w64, sd0[i] >> 7);
		CN(add)(&P, &P, &Q);
		window_lookup_16(&Q, win1, sd1[i] & 31);
		gf_condneg(&Q.W.w64, &Q.W.w64, sd1[i] >> 7);
		CN(add)(&P, &P, &Q);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}
#endif

/*
 * Lookup an affine point among 16 values (constant-time), in (x,u)
 * coordinates.
 *
 * Lookup index is between 0 and 16 (inclusive). The provided array
 * is supposed to hold 1*Q, 2*Q,... 16*Q, in that order, for some
 * point Q. If the index is 0, this returns the neutral; otherwise,
 * this returns index*Q.
 */
static void
window_lookup_16_affine_xu(CN(point_affine_xu) *P,
	const CN(point_affine_xu) *win, size_t index)
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
	P->U.w64.v0 = 0;
	P->U.w64.v1 = 0;
	P->U.w64.v2 = 0;
	P->U.w64.v3 = 0;

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
		P->U.w64.v0 |= m & win[u].U.w64.v0;
		P->U.w64.v1 |= m & win[u].U.w64.v1;
		P->U.w64.v2 |= m & win[u].U.w64.v2;
		P->U.w64.v3 |= m & win[u].U.w64.v3;
	}
}

/* see do255.h */
void
CN(mulgen)(CN(point) *P3, const void *scalar)
{
	CN(point_xu) P;
	CN(point_affine_xu) Qa;
	int i;
	uint8_t sd[52];

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
	 * Top digit of the full scalar is nonnegative.
	 */
	window_lookup_16_affine_xu(&Qa, window_G195_xu, sd[51]);
	P.X = Qa.X;
	P.U = Qa.U;
	P.Z.w64.v0 = 1;
	P.Z.w64.v1 = 0;
	P.Z.w64.v2 = 0;
	P.Z.w64.v3 = 0;
	P.T.w64.v0 = 1;
	P.T.w64.v1 = 0;
	P.T.w64.v2 = 0;
	P.T.w64.v3 = 0;

	/*
	 * Lookups and additions for the top digits of the three other
	 * chunks.
	 */
	window_lookup_16_affine_xu(&Qa, window_G_xu, sd[12] & 31);
	gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[12] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	window_lookup_16_affine_xu(&Qa, window_G65_xu, sd[25] & 31);
	gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[25] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	window_lookup_16_affine_xu(&Qa, window_G130_xu, sd[38] & 31);
	gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[38] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	for (i = 11; i >= 0; i --) {
		CN(double_x_xu)(&P, &P, 5);

		window_lookup_16_affine_xu(&Qa,
			window_G_xu, sd[i] & 31);
		gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[i] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_16_affine_xu(&Qa,
			window_G65_xu, sd[i + 13] & 31);
		gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[i + 13] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_16_affine_xu(&Qa,
			window_G130_xu, sd[i + 26] & 31);
		gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[i + 26] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_16_affine_xu(&Qa,
			window_G195_xu, sd[i + 39] & 31);
		gf_condneg(&Qa.U.w64, &Qa.U.w64, sd[i + 39] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);
	}

	/*
	 * Return the result in Jacobian (x,w) coordinates.
	 *   X3 = X*Z*U^2
	 *   W3 = Z*T   (necessarily non-zero)
	 *   Z3 = Z*U
	 */
	gf_mul(&P3->X.w64, &P.X.w64, &P.U.w64);
	gf_mul(&P3->W.w64, &P.Z.w64, &P.T.w64);
	gf_mul(&P3->Z.w64, &P.Z.w64, &P.U.w64);
	gf_mul(&P3->X.w64, &P3->X.w64, &P3->Z.w64);
}
