/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements scalar operations modulo r for curve do255e; it
 * works with any implementation with 32-bit limbs.
 */

/*
 * r = 2^254 - r0, with r0 = 131528281291764213006042413802501683931.
 */
static const i128 R0 = { {
	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
} };
static const i256 R = { {
	0x74D84525, 0x1F52C8AE, 0x54078C53, 0x9D0C930F,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x3FFFFFFF
} };

/*
 * 1/2 mod r = (r+1)/2
 */
static const i256 Rhf = { {
	0x3A6C2293, 0x8FA96457, 0xAA03C629, 0xCE864987,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x1FFFFFFF
} };

/*
 * 8*r modulo 2^256.
 */
static const i256 R_x8 = { {
	0xA6C22928, 0xFA964573, 0xA03C6298, 0xE864987A,
	0xFFFFFFFC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF
} };

/*
 * Given input 'a' (up to 2^258-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and is less than 2^254 + 2^131. The
 * high bits of 'a' are provided as extra parameter ah.
 */
static void
modr_reduce256_partial(i256 *d, const i256 *a, uint32_t ah)
{
	uint32_t t[8], u[5], x;
	unsigned char cc;
	int i;

	/*
	 * Truncate the source to 254 bits, and apply reduction (with
	 * 2^254 = r0 mod r) for the four top bits.
	 */
	ah = (ah << 2) | (a->v[7] >> 30);
	memcpy(t, a->v, 7 * sizeof(uint32_t));
	t[7] = a->v[7] & 0x3FFFFFFF;
	x = 0;
	for (i = 0; i < 4; i ++) {
		uint64_t z;

		z = (uint64_t)ah * (uint64_t)R0.v[i] + (uint64_t)x;
		u[i] = (uint32_t)z;
		x = (uint32_t)(z >> 32);
	}
	u[4] = x;

	cc = _addcarry_u32(0, t[0], u[0], &d->v[0]);
	for (i = 1; i < 5; i ++) {
		cc = _addcarry_u32(cc, t[i], u[i], &d->v[i]);
	}
	for (i = 5; i < 8; i ++) {
		cc = _addcarry_u32(cc, t[i], 0, &d->v[i]);
	}
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

	/*
	 * Subtracting r is equivalent to adding r0, then subtracting
	 * 2^254. Since source value was less than 2*r, adding r0 yields
	 * a value that still fits on 255 bits. We can thus use the
	 * high bit of the 256-bit intermediate result to check whether
	 * the subtraction of r yielded a negative result or not.
	 */
	cc = _addcarry_u32(0, a->v[0], R0.v[0], &t[0]);
	for (i = 1; i < 4; i ++) {
		cc = _addcarry_u32(cc, a->v[i], R0.v[i], &t[i]);
	}
	for (i = 4; i < 8; i ++) {
		cc = _addcarry_u32(cc, a->v[i], 0, &t[i]);
	}
	t[7] -= 0x40000000;

	/*
	 * If the subtraction result is not negative, then we keep it;
	 * otherwise, we use the source value.
	 */
	m = -(t[7] >> 31);
	for (i = 0; i < 8; i ++) {
		d->v[i] = t[i] ^ (m & (t[i] ^ a->v[i]));
	}
}

/*
 * Given input 'a' (up to 2^384-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and is less than 2^254 + 2^131.
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
	 * Multiply the high third (ah) by r0.
	 */
	memcpy(&ah.v[0], &a->v[8], 4 * sizeof(uint32_t));
	mul128x128(&t, &ah, &R0);

	/*
	 * Compute 4*ah*r0 and add to the low part of a. Since
	 * 4*r0 =~ 2^128.63, the result fits on 258 bits.
	 */
	t8 = t.v[7] >> 30;
	for (i = 7; i >= 1; i --) {
		t.v[i] = (t.v[i] << 2) | (t.v[i - 1] >> 30);
	}
	t.v[0] = t.v[0] << 2;
	cc = _addcarry_u32(0, t.v[0], a->v[0], &t.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, t.v[i], a->v[i], &t.v[i]);
	}
	(void)_addcarry_u32(cc, t8, 0, &t8);

	/*
	 * Perform partial reduction.
	 */
	modr_reduce256_partial(d, &t, t8);
}

/* see do255.h */
int
do255e_scalar_is_reduced(const void *a)
{
	i256 t;
	unsigned char cc;
	uint32_t x;
	int i;

	i256_decode(&t, a);
	cc = _subborrow_u32(0, t.v[0], R.v[0], &x);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u32(cc, t.v[i], R.v[i], &x);
	}
	return cc;
}

/* see do255.h */
void
do255e_scalar_sub(void *d, const void *a, const void *b)
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
	 * If subtraction yielded a negative value, then we add 8*r. Note
	 * that 2^256 < 8*r < 2^257; thus, we always get a nonnegative
	 * value that fits on 257 bits.
	 */
	m = -(uint32_t)cc;
	cc = _addcarry_u32(0, td.v[0], m & R_x8.v[0], &td.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, td.v[i], m & R_x8.v[i], &td.v[i]);
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
do255e_scalar_half(void *d, const void *a)
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
