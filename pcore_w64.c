/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined macros CURVE and CN()
 *  - defined gf and operations, including gf_sqrt() and gf_issquare()
 *  - defined CURVE_A and CURVE_4B (to curve parameters a and 4*b)
 *  - defined the conventional curve generator
 *
 * This file is for all implementations that use 64-bit limbs. It defines:
 *  - CURVE_neutral
 *  - CURVE_decode()
 *  - CURVE_encode()
 *  - CURVE_is_neutral()
 *  - CURVE_eq()
 */

/* see do255.h */
const CN(point) CN(neutral) = {
	{ { 0, 0, 0, 0 } },
	{ { 1, 0, 0, 0 } },
	{ { 0, 0, 0, 0 } }
};

/* see do255.h */
int
CN(decode)(CN(point) *P, const void *src)
{
	uint64_t r, qr, zz;
	gf x, w, d;

	/* Decode w. */
	r = gf_decode(&w, src);

	/*
	 * If value is zero, decoding succeeded but the square root below
	 * will fail, but we still want to report a success.
	 */
	zz = r & gf_iszero(&w);

	/* x <- w^2 - a */
	gf_sqr(&x, &w);
	gf_sub(&x, &x, &CURVE_A);

	/* d <- sqrt((w^2 - a)^2 - 4*b) */
	gf_sqr(&d, &x);
	gf_sub(&d, &d, &CURVE_4B);
	if (P == NULL) {
		/* Fast path: if we just want to know whether the point
		   is valid, then we do not have to compute the square
		   root. */
		uint64_t r3;

		r3 = (uint64_t)gf_legendre(&d);
		return (int)((r & ~(r3 >> 1)) | zz);
	}
	r &= gf_sqrt(&d, &d);

	/* x <- ((w^2 - a) + d)/2 */
	gf_add(&x, &x, &d);
	gf_half(&x, &x);

	/* If x is a square, then we must use the other solution,
	   i.e. ((w^2 - a) - d)/2, which we obtain by subtracting d. */
	qr = gf_issquare(&x);
	d.v0 &= -qr;
	d.v1 &= -qr;
	d.v2 &= -qr;
	d.v3 &= -qr;
	gf_sub(&x, &x, &d);

	/* If decoding failed, or ((w^2 - a)^2 - 4*b) was not a square,
	   then we clamp the returned value to X = 0, W = 1, Z = 0 (the
	   neutral point). */
	P->X.w64.v0 = x.v0 & -r;
	P->X.w64.v1 = x.v1 & -r;
	P->X.w64.v2 = x.v2 & -r;
	P->X.w64.v3 = x.v3 & -r;
	P->W.w64.v0 = (w.v0 & -r) | (1 - r);
	P->W.w64.v1 = w.v1 & -r;
	P->W.w64.v2 = w.v2 & -r;
	P->W.w64.v3 = w.v3 & -r;
	P->Z.w64.v0 = r;
	P->Z.w64.v1 = 0;
	P->Z.w64.v2 = 0;
	P->Z.w64.v3 = 0;

	/* If value was zero, the square root failed and the point has
	   been set to neutral, which is what we want, but we also
	   want to report it as a success. */
	r |= zz;

	return (int)r;
}

/* see do255.h */
void
CN(encode)(void *dst, const CN(point) *P)
{
	gf t;

	gf_inv(&t, &P->Z.w64);
	gf_mul(&t, &t, &P->W.w64);
	gf_encode(dst, &t);
}

/* see do255.h */
void
CN(encode_squared_w)(void *dst, const CN(point) *P)
{
	gf t;

	gf_inv(&t, &P->Z.w64);
	gf_mul(&t, &t, &P->W.w64);
	gf_sqr(&t, &t);
	gf_encode(dst, &t);
}

/* see do255.h */
int
CN(is_neutral)(const CN(point) *P)
{
	return (int)gf_iszero(&P->Z.w64);
}

/* see do255.h */
int
CN(eq)(const CN(point) *P1, const CN(point) *P2)
{
	gf t1, t2;

	gf_mul(&t1, &P1->W.w64, &P2->Z.w64);
	gf_mul(&t2, &P1->Z.w64, &P2->W.w64);
	return (int)gf_eq(&t1, &t2);
}

/* see do255.h */
void
CN(neg)(CN(point) *P3, const CN(point) *P1)
{
	P3->X = P1->X;
	gf_neg(&P3->W.w64, &P1->W.w64);
	P3->Z = P1->Z;
}
