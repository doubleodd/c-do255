/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - included <string.h> (for memset())
 *  - defined gf and operations
 *  - defined curve basic operations
 *  - defined curve multiplication core operations
 *  - defined CURVE to the curve name
 *
 * This file implements CURVE_mul() with a 4-bit window; it works with
 * ARM assembly implementations of finite fields with 32-bit limbs.
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
 * NOTE: in (x,u) coordinates, the endomorphism is:
 *     phi(x, u) = (-x, u/(-eta))
 * but since eta^2 = -1, we have 1/(-eta) = eta, so this also works:
 *     phi(x, u) = (-x, eta*u)
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

/*
 * Split scalar k (255 bits) into k0 anmd k1 (128 bits each, signed),
 * such that k = k0 + k1*mu mod r.
 */
void CN(split_scalar)(i128 *k0, i128 *k1, const void *k);
#define split_scalar   CN(split_scalar)

/*
 * Given scalar k (255 bits, top bit of last byte MUST be zero), split
 * it into k0 and k1 for use with the endomorphism, and recode |k0|
 * and |k1| into sd0[] and sd1[], for use with a 4-bit window.
 *
 * A recoded value consists of 32 digits, each in the -7..+8 range.
 * Negative values are encoded in sign + mantissa (top bit of the byte is
 * the sign). Upper digits (index 31) is necessarily in the 0..+8 range.
 *
 * Returned value contains the signs of k0 and k1:
 *   0   k0 and k1 are nonnegative
 *   1   k0 is negative
 *   2   k1 is negative
 *   3   both k0 and k1 are negative
 */
uint32_t CN(split_recode4_scalar)(uint8_t *sd0, uint8_t *sd1, const void *k);
#define split_recode4_scalar   CN(split_recode4_scalar)

/*
 * Fill win[i] with (i+1)*P in affine (x,u) coordinates, for i = 0..7
 * (the source point P is in Jacobian (x,w) coordinates).
 * (implemented in assembly)
 */
void CN(window_fill_8_to_xu_affine)(
	CN(point_affine_xu) *win, const CN(point) *P);
#define window_fill_8_to_xu_affine   CN(window_fill_8_to_xu_affine)

/*
 * Lookup for affine points in (x,u) coordinates is identical, since
 * the in-memory representations are identical.
 */
static inline void
window_lookup_8_affine_xu(
	CN(point_affine_xu) *P, const CN(point_affine_xu) *win, size_t index)
{
	window_lookup_8_affine((void *)P, (const void *)win, index);
}

/* see do255.h */
void
CN(mulgen)(CN(point) *P3, const void *scalar)
{
	CN(point_xu) P;
	CN(point_affine_xu) Qa;
	int i;
	uint8_t sd[64];

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
	 * Top digit of the full scalar is nonnegative.
	 */
	window_lookup_8_affine_xu(&Qa, window_G192_xu, sd[63]);
	P.X = Qa.X;
	P.U = Qa.U;
	P.Z.w32 = GF_ONE;
	P.T.w32 = GF_ONE;

	/*
	 * Lookups and additions for the top digits of the three other
	 * chunks.
	 */
	window_lookup_8_affine_xu(&Qa, window_G_xu, sd[15] & 15);
	gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[15] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	window_lookup_8_affine_xu(&Qa, window_G64_xu, sd[31] & 15);
	gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[31] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	window_lookup_8_affine_xu(&Qa, window_G128_xu, sd[47] & 15);
	gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[47] >> 7);
	CN(add_mixed_xu)(&P, &P, &Qa);

	for (i = 14; i >= 0; i --) {
		CN(double_x_xu)(&P, &P, 4);

		window_lookup_8_affine_xu(&Qa, window_G_xu, sd[i] & 15);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[i] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_8_affine_xu(&Qa, window_G64_xu, sd[i + 16] & 15);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[i + 16] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_8_affine_xu(&Qa, window_G128_xu, sd[i + 32] & 15);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[i + 32] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);

		window_lookup_8_affine_xu(&Qa, window_G192_xu, sd[i + 48] & 15);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sd[i + 48] >> 7);
		CN(add_mixed_xu)(&P, &P, &Qa);
	}

	/*
	 * Return the result in Jacobian (x,u) coordinates.
	 *   X3 = X*Z*U^2
	 *   W3 = Z*T   (necessarily non-zero)
	 *   Z3 = Z*U
	 */
	gf_mul(&P3->X.w32, &P.X.w32, &P.U.w32);
	gf_mul(&P3->W.w32, &P.Z.w32, &P.T.w32);
	gf_mul(&P3->Z.w32, &P.Z.w32, &P.U.w32);
	gf_mul(&P3->X.w32, &P3->X.w32, &P3->Z.w32);
}

/* see do255.h */
void
CN(mul)(CN(point) *P3, const CN(point) *P1, const void *scalar)
{
	CN(point_affine_xu) win0[8], Qa;
	union {
		CN(point) j;
		CN(point_xu) f;
	} P;
	uint8_t sd0[32], sd1[32];
	uint32_t sg;
	int i;

	/*
	 * Split input scalar into k0 and k1, recoded into sd0 and sd1.
	 */
	sg = split_recode4_scalar(sd0, sd1, scalar);

	/*
	 * Fill P with a copy of P1 or -P1, depending on the sign of k0,
	 * then use it to initialize the window.
	 */
	P.j.X = P1->X;
	gf_condneg(&P.j.W.w32, &P1->W.w32, sg & 1);
	P.j.Z = P1->Z;
	window_fill_8_to_xu_affine(win0, &P.j);

	/*
	 * Set sg to 0 if k0 and k1 have the same sign, 1 otherwise.
	 */
	sg = (sg ^ (sg >> 1)) & 1;

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
	 * Lookup points corresponding to the top digits, and add them
	 * into P.
	 */
	window_lookup_8_affine_xu(&Qa, win0, sd0[31]);
	P.f.X = Qa.X;
	P.f.Z.w32 = GF_ONE;
	P.f.U = Qa.U;
	P.f.T.w32 = GF_ONE;

	/* Second lookup uses win0, then applies the endomorphism. */
	window_lookup_8_affine_xu(&Qa, win0, sd1[31]);
	gf_neg(&Qa.X.w32, &Qa.X.w32);
	gf_mul_inline(&Qa.U.w32, &Qa.U.w32, &ETA);
	gf_condneg(&Qa.U.w32, &Qa.U.w32, sg);
	CN(add_mixed_xu)(&P.f, &P.f, &Qa);

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
		CN(double_x_xu)(&P.f, &P.f, 4);
		window_lookup_8_affine_xu(&Qa, win0, sd0[i] & 15);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sd0[i] >> 7);
		CN(add_mixed_xu)(&P.f, &P.f, &Qa);
		window_lookup_8_affine_xu(&Qa, win0, sd1[i] & 15);
		gf_neg(&Qa.X.w32, &Qa.X.w32);
		gf_mul_inline(&Qa.U.w32, &Qa.U.w32, &ETA);
		gf_condneg(&Qa.U.w32, &Qa.U.w32, sg ^ (sd1[i] >> 7));
		CN(add_mixed_xu)(&P.f, &P.f, &Qa);
	}

	/*
	 * Return the result in Jacobian (x,u) coordinates.
	 *   X3 = X*Z*U^2
	 *   W3 = Z*T   (necessarily non-zero)
	 *   Z3 = Z*U
	 */
	gf_mul(&P3->X.w32, &P.f.X.w32, &P.f.U.w32);
	gf_mul(&P3->W.w32, &P.f.Z.w32, &P.f.T.w32);
	gf_mul(&P3->Z.w32, &P.f.Z.w32, &P.f.U.w32);
	gf_mul(&P3->X.w32, &P3->X.w32, &P3->Z.w32);
}
