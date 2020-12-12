/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations, but NOT gf_sqrt() or gf_issquare()
 *
 * This file is for all implementations of do255e that use 32-bit limbs.
 * It defines:
 *  - gf_sqrt()
 *  - gf_issquare()
 */

/*
 * Square root computation. Returned value is 1 on success (value was a
 * quadratic residue), 0 on failure (value was not a quadratic residue).
 * On success, the returned square root is the one whose least
 * significant bit (as an integer in the 0..p-1 range) is zero. If a
 * failure is reported, then the value written to *d is zero.
 *
 * If d == NULL, the quadratic residue status is still computed and
 * returned.
 */
UNUSED
static uint32_t
gf_sqrt(gf *d, const gf *a)
{
	/*
	 * Since p = 5 mod 8, we use Atkin's algorithm:
	 *   b <- (2*a)^((p-5)/8)
	 *   c <- 2*a*b^2
	 *   return a*b*(c - 1)
	 */
	gf b, c, e, x, x2, x96, y;
	uint32_t qr;
	int i;

	if (d == NULL) {
		return 1 - ((uint32_t)gf_legendre(a) >> 31);
	}

	/* e <- 2*a */
	gf_mul2(&e, a);

	/* Raise e to the power (p-5)/8. This is the expensive step.
	   Sequence below does it in 251 squarings and 13 extra
	   multiplications.
	   (p-5)/8 = (2^240-1)*2^12 + (2^2-1)*2^9 + (2^3-1)*2^5 + 2^2
	   */

	/* x2 <- e^3 */
	gf_sqr(&x2, &e);
	gf_mul(&x2, &x2, &e);

	/* x <- e^(2^4-1) */
	gf_sqr_x(&x, &x2, 2);
	gf_mul(&x, &x, &x2);

	/* x <- e^(2^8-1) */
	gf_sqr_x(&y, &x, 4);
	gf_mul(&x, &y, &x);

	/* x <- e^(2^16-1) */
	gf_sqr_x(&y, &x, 8);
	gf_mul(&x, &y, &x);

	/* x <- e^(2^48-1) */
	gf_sqr_x(&y, &x, 16);
	gf_mul(&y, &y, &x);
	gf_sqr_x(&y, &y, 16);
	gf_mul(&x, &y, &x);

	/* x96 <- e^(2^96-1) */
	gf_sqr_x(&y, &x, 48);
	gf_mul(&x96, &y, &x);

	/* x <- e^(2^240-1) */
	gf_sqr_x(&y, &x96, 96);
	gf_mul(&y, &y, &x96);
	gf_sqr_x(&y, &y, 48);
	gf_mul(&x, &y, &x);

	/* x <- e^((p-5)/8) */
	gf_sqr_x(&x, &x, 3);
	gf_mul(&x, &x, &x2);
	gf_sqr_x(&x, &x, 2);
	gf_mul(&x, &x, &e);
	gf_sqr_x(&x, &x, 2);
	gf_mul(&x, &x, &x2);
	gf_sqr_x(&x, &x, 3);
	gf_mul(&x, &x, &e);
	gf_sqr_x(&b, &x, 2);

	/* We now have b = (2*a)^((p-5)/8). */

	/* c <- 2*a*b^2 */
	gf_sqr(&c, &b);
	gf_mul(&c, &c, &e);

	/* x <- a*b*(c - 1) */
	gf_sub(&x, &c, &GF_ONE);
	gf_mul(&x, &x, a);
	gf_mul(&x, &x, &b);

	/* Normalize and adjust the "sign" if needed. */
	gf_normalize(&x, &x);
	gf_condneg(&x, &x, x.v[0] & 1);

	/*
	 * We now have a potential square root in c. We must check that
	 * it is indeed a square root of a (the source value was not
	 * necessarily a quadratic residue).
	 */

	/* Verify the result; clear it if input was not a square. */
	gf_sqr(&y, &x);
	qr = gf_eq(&y, a);
	for (i = 0; i < 8; i ++) {
		x.v[i] &= -qr;
	}

	/* Return value. */
	*d = x;
	return qr;
}

/*
 * Get quadratic residue status. Returned value is 1 for a quadratic
 * residue, 0 otherwise.
 */
UNUSED
static uint32_t
gf_issquare(const gf *a)
{
	return 1 - ((uint32_t)gf_legendre(a) >> 31);
}
