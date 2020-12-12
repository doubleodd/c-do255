/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements scalar operations modulo r for curve do255e; it
 * works with any implementation with 64-bit limbs.
 */

/*
 * r = 2^254 - r0, with r0 = 131528281291764213006042413802501683931.
 */
#define R0_lo   ((uint64_t)0xE0AD37518B27BADB)
#define R0_hi   ((uint64_t)0x62F36CF0ABF873AC)
static const i128 R0 = {
	R0_lo, R0_hi
};

#define R_0   ((uint64_t)0x1F52C8AE74D84525)
#define R_1   ((uint64_t)0x9D0C930F54078C53)
#define R_2   ((uint64_t)0xFFFFFFFFFFFFFFFF)
#define R_3   ((uint64_t)0x3FFFFFFFFFFFFFFF)

/*
 * 1/2 modulo r  (equal to (r+1)/2).
 */
#define Rhf_0   ((uint64_t)0x8FA964573A6C2293)
#define Rhf_1   ((uint64_t)0xCE864987AA03C629)
#define Rhf_2   ((uint64_t)0xFFFFFFFFFFFFFFFF)
#define Rhf_3   ((uint64_t)0x1FFFFFFFFFFFFFFF)

/*
 * 8*r modulo 2^256.
 */
#define R_x8_0    ((uint64_t)0xFA964573A6C22928)
#define R_x8_1    ((uint64_t)0xE864987AA03C6298)
#define R_x8_2    ((uint64_t)0xFFFFFFFFFFFFFFFC)
#define R_x8_3    ((uint64_t)0xFFFFFFFFFFFFFFFF)

/*
 * Given input 'a' (up to 2^258-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and is less than 2^254 + 2^131. The
 * high bits of 'a' are provided as extra parameter ah.
 */
static void
modr_reduce256_partial(i256 *d, const i256 *a, uint64_t ah)
{
	unsigned long long t0, t1, t2, t3, u0, u1, u2, x;
	unsigned char cc;

	/*
	 * Truncate the source to 254 bits, and apply reduction (with
	 * 2^254 = r0 mod r) for the four top bits.
	 */
	ah = (ah << 2) | (a->v3 >> 62);
	t0 = a->v0;
	t1 = a->v1;
	t2 = a->v2;
	t3 = a->v3 & 0x3FFFFFFFFFFFFFFF;
	UMUL64(u0, u1, ah, R0_lo);
	UMUL64(x, u2, ah, R0_hi);
	cc = _addcarry_u64(0, u1, x, &u1);
	(void)_addcarry_u64(cc, u2, 0, &u2);

	cc = _addcarry_u64(0, t0, u0, (unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, t1, u1, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, t2, u2, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, t3, 0, (unsigned long long *)&d->v3);
}

/*
 * Given a partially reduced input (less than 2*r), finish reduction
 * (conditional subtraction of r).
 */
static void
modr_reduce256_finish(i256 *d, const i256 *a)
{
	unsigned long long t0, t1, t2, t3, m;
	unsigned char cc;

	/*
	 * Subtracting r is equivalent to adding r0, then subtracting
	 * 2^254. Since source value was less than 2*r, adding r0 yields
	 * a value that still fits on 255 bits. We can thus use the
	 * high bit of the 256-bit intermediate result to check whether
	 * the subtraction of r yielded a negative result or not.
	 */
	cc = _addcarry_u64(0, a->v0, R0_lo, &t0);
	cc = _addcarry_u64(cc, a->v1, R0_hi, &t1);
	cc = _addcarry_u64(cc, a->v2, 0, &t2);
	(void)_addcarry_u64(cc, a->v3, 0, &t3);
	t3 -= 0x4000000000000000;

	/*
	 * If the subtraction result is not negative, then we keep it;
	 * otherwise, we use the source value.
	 */
	m = -(t3 >> 63);
	d->v0 = t0 ^ (m & (t0 ^ a->v0));
	d->v1 = t1 ^ (m & (t1 ^ a->v1));
	d->v2 = t2 ^ (m & (t2 ^ a->v2));
	d->v3 = t3 ^ (m & (t3 ^ a->v3));
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
	unsigned long long t4;
	unsigned char cc;

	/*
	 * Multiply the high third (ah) by r0.
	 */
	ah.v0 = a->v4;
	ah.v1 = a->v5;
	mul128x128(&t, &ah, &R0);

	/*
	 * Compute 4*ah*r0 and add to the low part of a. Since
	 * 4*r0 =~ 2^128.63, the result fits on 258 bits.
	 */
	t4 = t.v3 >> 62;
	t.v3 = (t.v3 << 2) | (t.v2 >> 62);
	t.v2 = (t.v2 << 2) | (t.v1 >> 62);
	t.v1 = (t.v1 << 2) | (t.v0 >> 62);
	t.v0 = t.v0 << 2;
	cc = _addcarry_u64(0, t.v0, a->v0, (unsigned long long *)&t.v0);
	cc = _addcarry_u64(cc, t.v1, a->v1, (unsigned long long *)&t.v1);
	cc = _addcarry_u64(cc, t.v2, a->v2, (unsigned long long *)&t.v2);
	cc = _addcarry_u64(cc, t.v3, a->v3, (unsigned long long *)&t.v3);
	(void)_addcarry_u64(cc, t4, 0, &t4);

	/*
	 * Perform partial reduction.
	 */
	modr_reduce256_partial(d, &t, t4);
}

/* see do255.h */
int
do255e_scalar_is_reduced(const void *a)
{
	i256 t;
	unsigned char cc;
	unsigned long long x;

	i256_decode(&t, a);
	cc = _subborrow_u64(0, t.v0, R_0, &x);
	cc = _subborrow_u64(cc, t.v1, R_1, &x);
	cc = _subborrow_u64(cc, t.v2, R_2, &x);
	cc = _subborrow_u64(cc, t.v3, R_3, &x);
	return cc;
}

/* see do255.h */
void
do255e_scalar_sub(void *d, const void *a, const void *b)
{
	i256 ta, tb, td;
	uint64_t t4, m;
	unsigned char cc;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	cc = _subborrow_u64(0, ta.v0, tb.v0, (unsigned long long *)&td.v0);
	cc = _subborrow_u64(cc, ta.v1, tb.v1, (unsigned long long *)&td.v1);
	cc = _subborrow_u64(cc, ta.v2, tb.v2, (unsigned long long *)&td.v2);
	cc = _subborrow_u64(cc, ta.v3, tb.v3, (unsigned long long *)&td.v3);

	/*
	 * If subtraction yielded a negative value, then we add 8*r. Note
	 * that 2^256 < 8*r < 2^257; thus, we always get a nonnegative
	 * value that fits on 257 bits.
	 */
	m = -(uint64_t)cc;
	cc = _addcarry_u64(0, td.v0, m & R_x8_0, (unsigned long long *)&td.v0);
	cc = _addcarry_u64(cc, td.v1, m & R_x8_1, (unsigned long long *)&td.v1);
	cc = _addcarry_u64(cc, td.v2, m & R_x8_2, (unsigned long long *)&td.v2);
	cc = _addcarry_u64(cc, td.v3, m & R_x8_3, (unsigned long long *)&td.v3);
	t4 = cc;

	/*
	 * Reduce modulo r.
	 */
	modr_reduce256_partial(&td, &td, t4);
	modr_reduce256_finish(&td, &td);
	i256_encode(d, &td);
}

/* see do255.h */
void
do255e_scalar_half(void *d, const void *a)
{
	i256 x;
	unsigned char cc;
	uint64_t m;

	i256_decode(&x, a);
	m = -(x.v0 & 1);

	/*
	 * Right-shift value. Result is lower than 2^255.
	 */
	x.v0 = (x.v0 >> 1) | (x.v1 << 63);
	x.v1 = (x.v1 >> 1) | (x.v2 << 63);
	x.v2 = (x.v2 >> 1) | (x.v3 << 63);
	x.v3 = (x.v3 >> 1);

	/*
	 * Add (r+1)/2 if the value was odd. This cannot overflow
	 * since r < 2^255.
	 */
	cc = _addcarry_u64(0, x.v0, m & Rhf_0, (unsigned long long *)&x.v0);
	cc = _addcarry_u64(cc, x.v1, m & Rhf_1, (unsigned long long *)&x.v1);
	cc = _addcarry_u64(cc, x.v2, m & Rhf_2, (unsigned long long *)&x.v2);
	(void)_addcarry_u64(cc, x.v3, m & Rhf_3, (unsigned long long *)&x.v3);

	/*
	 * Ensure a reduced value.
	 */
	modr_reduce256_partial(&x, &x, 0);
	modr_reduce256_finish(&x, &x);
	i256_encode(d, &x);
}
