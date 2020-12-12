/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements scalar operations modulo r for curve do255s; it
 * works with any implementation with 32-bit limbs.
 */

/*
 * r = 2^254 + r0, with r0 = 56904135270672826811114353017034461895.
 */
static const i128 R0 = { {
	0x396152C7, 0xDCF2AC65, 0x912B7F03, 0x2ACF567A
} };
#define R_top   ((uint32_t)0x40000000)

/*
 * 1/2 mod r = (r+1)/2
 */
static const i256 Rhf = { {
	0x9CB0A964, 0xEE795632, 0x4895BF81, 0x1567AB3D,
	0x00000000, 0x00000000, 0x00000000, 0x20000000
} };

/*
 * 4*r0
 */
static const i128 R0_x4 = { {
	0xE5854B1C, 0x73CAB194, 0x44ADFC0F, 0xAB3D59EA
} };

/*
 * Given input 'a' (up to 2^286-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and thus is less than 2*r. The high
 * bits of 'a' are provided as extra parameter ah.
 */
static void
modr_reduce256_partial(i256 *d, const i256 *a, uint32_t ah)
{
	i256 t;
	uint32_t u[5], x;
	unsigned char cc;
	int i;

	/*
	 * Truncate the source to 254 bits, and apply reduction (with
	 * 2^254 = -r0 mod r) for the top bits.
	 */
	ah = (ah << 2) | (a->v[7] >> 30);
	memcpy(&t.v[0], &a->v[0], 7 * sizeof(uint32_t));
	t.v[7] = a->v[7] & 0x3FFFFFFF;
	x = 0;
	for (i = 0; i < 4; i ++) {
		uint64_t z;

		z = (uint64_t)ah * (uint64_t)R0.v[i] + (uint64_t)x;
		u[i] = (uint32_t)z;
		x = (uint32_t)(z >> 32);
	}
	u[4] = x;

	cc = _subborrow_u32(0, t.v[0], u[0], &t.v[0]);
	for (i = 1; i < 5; i ++) {
		cc = _subborrow_u32(cc, t.v[i], u[i], &t.v[i]);
	}
	for (i = 5; i < 8; i ++) {
		cc = _subborrow_u32(cc, t.v[i], 0, &t.v[i]);
	}

	/*
	 * We may get a borrow here, in which case the value is negative.
	 * But since we subtracted an integer lower than 2^192, adding r
	 * once will be enough to get back to the expected range.
	 */
	x = -(uint32_t)cc;
	cc = _addcarry_u32(0, t.v[0], x & R0.v[0], &d->v[0]);
	for (i = 1; i < 4; i ++) {
		cc = _addcarry_u32(cc, t.v[i], x & R0.v[i], &d->v[i]);
	}
	for (i = 4; i < 7; i ++) {
		cc = _addcarry_u32(cc, t.v[i], 0, &d->v[i]);
	}
	(void)_addcarry_u32(cc, t.v[7], x & R_top, &d->v[7]);
}

/*
 * Given a partially reduced input (less than 2*r), finish reduction
 * (conditional subtraction of r).
 */
static void
modr_reduce256_finish(i256 *d, const i256 *a)
{
	uint32_t t[8], m;
	unsigned char cc;
	int i;

	cc = _subborrow_u32(0, a->v[0], R0.v[0], &t[0]);
	for (i = 1; i < 4; i ++) {
		cc = _subborrow_u32(cc, a->v[i], R0.v[i], &t[i]);
	}
	for (i = 4; i < 7; i ++) {
		cc = _subborrow_u32(cc, a->v[i], 0, &t[i]);
	}
	cc = _subborrow_u32(cc, a->v[7], R_top, &t[7]);

	/*
	 * If there was no borrow, then we keep the subtraction result;
	 * otherwise, we use the source value.
	 */
	m = -(uint32_t)cc;
	for (i = 0; i < 8; i ++) {
		d->v[i] = t[i] ^ (m & (t[i] ^ a->v[i]));
	}
}

/*
 * Given input 'a' (up to 2^384-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits.
 */
static void
modr_reduce384_partial(i256 *d, const i384 *a)
{
	i128 ah;
	i256 t;
	uint32_t t8;
	unsigned char cc;
	int i;

	/*
	 * Multiply the high third (ah) by 4*r0.
	 */
	memcpy(&ah.v[0], &a->v[8], 4 * sizeof(uint32_t));
	mul128x128(&t, &ah, &R0_x4);

	/*
	 * Subtract 4*ah*r0 from the low part of a, then add back 4*r.
	 * Since 4*r0 =~ 2^127.42, the result is necessarily positive;
	 * it may range up to 4*r + 2^256 - 1, which may be slightly
	 * above 2^257.
	 */
	cc = _subborrow_u32(0, a->v[0], t.v[0], &t.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u32(cc, a->v[i], t.v[i], &t.v[i]);
	}
	t8 = -(uint32_t)cc;
	cc = _addcarry_u32(0, t.v[0], R0_x4.v[0], &t.v[0]);
	for (i = 1; i < 4; i ++) {
		cc = _addcarry_u32(cc, t.v[i], R0_x4.v[i], &t.v[i]);
	}
	for (i = 4; i < 8; i ++) {
		cc = _addcarry_u32(cc, t.v[i], 0, &t.v[i]);
	}
	(void)_addcarry_u32(cc, t8, 1, &t8);

	/*
	 * Perform partial reduction.
	 */
	modr_reduce256_partial(d, &t, t8);
}

/* see do255.h */
int
do255s_scalar_is_reduced(const void *a)
{
	i256 t;
	unsigned char cc;
	uint32_t x;

	i256_decode(&t, a);
	cc = _subborrow_u32(0, t.v[0], R0.v[0], &x);
	cc = _subborrow_u32(cc, t.v[1], R0.v[1], &x);
	cc = _subborrow_u32(cc, t.v[2], R0.v[2], &x);
	cc = _subborrow_u32(cc, t.v[3], R0.v[3], &x);
	cc = _subborrow_u32(cc, t.v[4], 0, &x);
	cc = _subborrow_u32(cc, t.v[5], 0, &x);
	cc = _subborrow_u32(cc, t.v[6], 0, &x);
	cc = _subborrow_u32(cc, t.v[7], R_top, &x);
	return cc;
}

/* see do255.h */
void
do255s_scalar_sub(void *d, const void *a, const void *b)
{
	i256 ta, tb, td;
	uint32_t t8, m;
	unsigned char cc;
	int i;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	cc = _subborrow_u32(0, ta.v[0], tb.v[0], &td.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u32(cc, ta.v[i], tb.v[i], &td.v[i]);
	}

	/*
	 * If subtraction yielded a negative value, then we add 4*r. Note
	 * that 2^256 < 4*r < 2^257; thus, we always get a nonnegative
	 * value that fits on 257 bits.
	 */
	m = -(uint32_t)cc;
	cc = _addcarry_u32(0, td.v[0], m & R0_x4.v[0], &td.v[0]);
	for (i = 1; i < 4; i ++) {
		cc = _addcarry_u32(cc, td.v[i], m & R0_x4.v[i], &td.v[i]);
	}
	for (i = 4; i < 8; i ++) {
		cc = _addcarry_u32(cc, td.v[i], 0, &td.v[i]);
	}
	t8 = cc;

	/*
	 * Reduce modulo r.
	 */
	modr_reduce256_partial(&td, &td, t8);
	modr_reduce256_finish(&td, &td);
	i256_encode(d, &td);
}

/* see do255.h */
void
do255s_scalar_half(void *d, const void *a)
{
	i256 x;
	unsigned char cc;
	uint32_t m;
	int i;

	i256_decode(&x, a);
	m = -(x.v[0] & 1);

	/*
	 * Right-shift value. Result is lower than 2^255.
	 */
	for (i = 0; i < 7; i ++) {
		x.v[i] = (x.v[i] >> 1) | (x.v[i + 1] << 31);
	}
	x.v[7] = (x.v[7] >> 1);

	/*
	 * Add (r+1)/2 if the value was odd. This cannot overflow
	 * since r < 2^255.
	 */
	cc = _addcarry_u32(0, x.v[0], m & Rhf.v[0], &x.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, x.v[i], m & Rhf.v[i], &x.v[i]);
	}

	/*
	 * Ensure a reduced value.
	 */
	modr_reduce256_partial(&x, &x, 0);
	modr_reduce256_finish(&x, &x);
	i256_encode(d, &x);
}
