/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations, but NOT gf_sqrt() or gf_issquare()
 *
 * This file is for all implementations of do255s that use 64-bit limbs.
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
static uint64_t
gf_sqrt(gf *d, const gf *a)
{
	/*
	 * Since p = 3 mod 4, we can compute a potential square root by
	 * raising to power (p+1)/4. Sequence below does it in 252
	 * squarings and 12 extra multiplications. An extra squaring
	 * is used to verify that a square root is indeed obtained.
	 */
	gf x, x2, y;
	uint64_t qr;

	if (d == NULL) {
		return 1 - ((uint64_t)gf_legendre(a) >> 63);
	}

	/* x2 <- a^3 */
	gf_sqr(&x2, a);
	gf_mul(&x2, &x2, a);

	/* x <- a^(2^3-1) */
	gf_sqr(&x, &x2);
	gf_mul(&x, &x, a);

	/* x <- a^(2^9-1) */
	gf_sqr_x(&y, &x, 3);
	gf_mul(&y, &y, &x);
	gf_sqr_x(&y, &y, 3);
	gf_mul(&x, &y, &x);

	/* x <- a^(2^27-1) */
	gf_sqr_x(&y, &x, 9);
	gf_mul(&y, &y, &x);
	gf_sqr_x(&y, &y, 9);
	gf_mul(&x, &y, &x);

	/* x <- a^(2^81-1) */
	gf_sqr_x(&y, &x, 27);
	gf_mul(&y, &y, &x);
	gf_sqr_x(&y, &y, 27);
	gf_mul(&x, &y, &x);

	/* x <- a^(2^243-1) */
	gf_sqr_x(&y, &x, 81);
	gf_mul(&y, &y, &x);
	gf_sqr_x(&y, &y, 81);
	gf_mul(&x, &y, &x);

	/* d <- a^(2^253 - 1024 + 35) */
	gf_sqr_x(&x, &x, 5);
	gf_mul(&x, &x, a);
	gf_sqr_x(&x, &x, 5);
	gf_mul(&x, &x, &x2);

	/* Normalize and adjust the "sign" if needed. */
	gf_normalize(&x, &x);
	gf_condneg(&x, &x, x.v0 & 1);

	/* Verify the result; clear it if input was not a square. */
	gf_sqr(&y, &x);
	qr = gf_eq(&y, a);
	x.v0 &= -qr;
	x.v1 &= -qr;
	x.v2 &= -qr;
	x.v3 &= -qr;

	/* Return value. */
	*d = x;
	return qr;
}

/*
 * Get quadratic residue status. Returned value is 1 for a quadratic
 * residue, 0 otherwise.
 */
UNUSED
static uint64_t
gf_issquare(const gf *a)
{
	return 1 - ((uint64_t)gf_legendre(a) >> 63);
}
