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
 * This file implements CURVE_mul() with a 4-bit window; it works with
 * any finite field implementation with 32-bit limbs. A 4-bit window
 * is used because 32-bit archs may be small microcontrollers, which
 * usually have very little available RAM; similarly, a single window
 * is used, the endomorphism being applied dynamically on looked-up
 * points.
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
 * Thus, k0 and k1 (which are signed integer) can fit on 128 bits each,
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

static const gf ETA = { {
	0xAA938AEE, 0xD99E0F1B, 0xB30E6336, 0xA60D864F,
	0xE53688E3, 0xE414983F, 0x3C69B85F, 0x10ED2DB3
} };

static const gf MINUS_ETA = { {
	0x556C2C37, 0x2661F0E4, 0x4CF19CC9, 0x59F279B0,
	0x1AC9771C, 0x1BEB67C0, 0xC39647A0, 0x6F12D24C
} };

/* (r - 1)/2 */
static const i256 HR = { {
	0x3A6C2292, 0x8FA96457, 0xAA03C629, 0xCE864987,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x1FFFFFFF
} };

/* u */
static const i128 eU = { {
	0xC93F6111, 0x2ACCF9DE, 0x53C2C6E6, 0x1A509F7A
} };

/* v */
static const i128 eV = { {
	0x5466F77E, 0x0B7A3130, 0xFFBB3A93, 0x7D440C6A
} };

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
	uint32_t w;
	unsigned char cc;
	int i;

	/* z <- k*e */
	mul256x128(&z, k, e);

	/* z <- z + (r-1)/2 */
	cc = _addcarry_u32(0, HR.v[0], z.v[0], &z.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, HR.v[i], z.v[i], &z.v[i]);
	}
	for (i = 8; i < 12; i ++) {
		cc = _addcarry_u32(cc, 0, z.v[i], &z.v[i]);
	}

	/* y <- floor(z / 2^254) + 1 */
	y.v[0] = (z.v[7] >> 30) | (z.v[8] << 2);
	y.v[1] = (z.v[8] >> 30) | (z.v[9] << 2);
	y.v[2] = (z.v[9] >> 30) | (z.v[10] << 2);
	y.v[3] = (z.v[10] >> 30) | (z.v[11] << 2);
	cc = _addcarry_u32(0, y.v[0], 1, &y.v[0]);
	cc = _addcarry_u32(cc, y.v[1], 0, &y.v[1]);
	cc = _addcarry_u32(cc, y.v[2], 0, &y.v[2]);
	(void)_addcarry_u32(cc, y.v[3], 0, &y.v[3]);

	/* t <- y*r0 */
	mul128x128(&t, &y, &R0);

	/* t <- t + z0
	   We are only interested in the high limb. */
	z.v[7] = z.v[7] & 0x3FFFFFFF;
	cc = _addcarry_u32(0, z.v[0], t.v[0], &w);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, z.v[i], t.v[i], &w);
	}

	/* The high limb is in 'w' and it is lower than 2^31. If
	   it is lower than 2^30, then y is too large and we must
	   decrement it; otherwise, we keep it unchanged. */
	cc = _subborrow_u32(0, y.v[0], 1 - (w >> 30), &d->v[0]);
	cc = _subborrow_u32(cc, y.v[1], 0, &d->v[1]);
	cc = _subborrow_u32(cc, y.v[2], 0, &d->v[2]);
	(void)_subborrow_u32(cc, y.v[3], 0, &d->v[3]);
}

/*
 * Split scalar k (255 bits) into k0 and k1 (128 bits each, signed),
 * such that k = k0 + k1*mu mod r.
 */
static void
split_scalar(i128 *k0, i128 *k1, const void *k)
{
	i256 t;
	i128 c, d, s;
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
	mul128x128trunc(&s, &d, &eU);
	cc = _subborrow_u32(0, t.v[0], s.v[0], &k0->v[0]);
	cc = _subborrow_u32(cc, t.v[1], s.v[1], &k0->v[1]);
	cc = _subborrow_u32(cc, t.v[2], s.v[2], &k0->v[2]);
	(void)_subborrow_u32(cc, t.v[3], s.v[3], &k0->v[3]);
	mul128x128trunc(&s, &c, &eV);
	cc = _subborrow_u32(0, k0->v[0], s.v[0], &k0->v[0]);
	cc = _subborrow_u32(cc, k0->v[1], s.v[1], &k0->v[1]);
	cc = _subborrow_u32(cc, k0->v[2], s.v[2], &k0->v[2]);
	(void)_subborrow_u32(cc, k0->v[3], s.v[3], &k0->v[3]);

	mul128x128trunc(k1, &d, &eV);
	mul128x128trunc(&s, &c, &eU);
	cc = _subborrow_u32(0, k1->v[0], s.v[0], &k1->v[0]);
	cc = _subborrow_u32(cc, k1->v[1], s.v[1], &k1->v[1]);
	cc = _subborrow_u32(cc, k1->v[2], s.v[2], &k1->v[2]);
	(void)_subborrow_u32(cc, k1->v[3], s.v[3], &k1->v[3]);
}

/*
 * Recode a signed 128-bit integer with a 4-bit window. The scalar is
 * first replaced with its absolute value x (its original sign is
 * returned: 1 if x was negative, 0 otherwise). The function then
 * computes 32 digits such that all digits are in the -7..+8 range,
 * and x = \sum_i sd[i]*2^(4*i). The top digit (sd[31]) is in the 0..+8
 * range.
 *
 * Digits are encoded in sign+mantissa format: sign bit is bit 7 in the
 * byte (1 for negative, 0 for positive).
 */
static uint32_t
recode4_small(uint8_t *sd, const i128 *s)
{
	i128 x;
	uint32_t sg, t;
	unsigned char cc;
	unsigned db, b, m;
	int i, j;

	/* x = abs(s) */
	sg = s->v[3] >> 31;
	cc = _addcarry_u32(0, s->v[0] ^ -sg, sg, &x.v[0]);
	cc = _addcarry_u32(cc, s->v[1] ^ -sg, 0, &x.v[1]);
	cc = _addcarry_u32(cc, s->v[2] ^ -sg, 0, &x.v[2]);
	(void)_addcarry_u32(cc, s->v[3] ^ -sg, 0, &x.v[3]);

	/*
	 * Compute all digits. Each limb yields exactly 8 digits.
	 */
	db = 0;
	for (j = 0; j < 4; j ++) {
		t = x.v[j];
		for (i = 0; i < 8; i ++) {
			b = ((unsigned)t & 0x0F) + db;
			m = (8 - b) >> 8;
			b ^= m & (b ^ ((16 - b) | 0x80));
			db = m & 1;
			sd[8 * j + i] = (uint8_t)b;
			t >>= 4;
		}
	}

	return sg;
}

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point_affine) win0[8], Qa;
	CN(point) P;
	i128 k0, k1;
	uint8_t sd0[32], sd1[32];
	uint32_t sg, qz;
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
	sg = recode4_small(sd0, &k0);
	P = *P1;
	gf_condneg(&P.W.w32, &P.W.w32, sg);
	window_fill_8_affine(win0, &P);

	/*
	 * Recode k1; check whether the sign of k1 differs from that of k0.
	 */
	sg ^= recode4_small(sd1, &k1);

	/*
	 * We do NOT apply the endomorphism on all points of win0 in
	 * order to compute a window for the second point: for the 32-bit
	 * implementation, we suppose that RAM is a scarce resource (as
	 * is normally the case on microcontrollers) and we instead
	 * dynamically apply the endomorphism for each relevant lookup.
	 * Extra runtime cost is low (32 field multiplications).
	 */

	/*
	 * Lookup points corresponding to the top digits, and add them
	 * into P. Take care that if the first looked-up point is the
	 * neutral, then we get affine neutral point (0, 0) and we must
	 * then set W to 1 and Z to 0; otherwise, we use the retrieved W
	 * and set Z to 1.
	 */
	window_lookup_8_affine(&Qa, win0, sd0[31]);
	qz = gf_iszero(&Qa.X.w32);
	P.X = Qa.X;
	P.W = Qa.W;
	P.W.w32.v[0] |= qz;
	P.Z.w32.v[0] = 1 - qz;
	for (i = 1; i < 8; i ++) {
		P.Z.w32.v[i] = 0;
	}
	/* Second lookup uses win0, then applies the endomorphism. */
	window_lookup_8_affine(&Qa, win0, sd1[31]);
	gf_neg(&Qa.X.w32, &Qa.X.w32);
	gf_mul_inline(&Qa.W.w32, &Qa.W.w32, &MINUS_ETA);
	gf_condneg(&Qa.W.w32, &Qa.W.w32, sg);
	CN(add_mixed)(&P, &P, &Qa);

	/*
	 * Process other digits from top to bottom. For each digit:
	 *  - multiply current value by 16 (4 successive doublings);
	 *  - lookup point from first window; negate it if the digit is
	 *    negative;
	 *  - add point to current value;
	 *  - lookup point from second window; negate it if the digit is
	 *    negative;
	 *  - add point to current value.
	 */
	for (i = 30; i >= 0; i --) {
		CN(double_x)(&P, &P, 4);
		window_lookup_8_affine(&Qa, win0, sd0[i] & 15);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sd0[i] >> 7);
		CN(add_mixed)(&P, &P, &Qa);
		window_lookup_8_affine(&Qa, win0, sd1[i] & 15);
		gf_neg(&Qa.X.w32, &Qa.X.w32);
		gf_mul_inline(&Qa.W.w32, &Qa.W.w32, &MINUS_ETA);
		gf_condneg(&Qa.W.w32, &Qa.W.w32, sg ^ (sd1[i] >> 7));
		CN(add_mixed)(&P, &P, &Qa);
	}

	/*
	 * Return the result.
	 */
	*P3 = P;
}
