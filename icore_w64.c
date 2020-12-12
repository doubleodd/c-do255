/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements some operations on plain integers, used in
 * various places for computations related to scalars. It works with any
 * implementation with 64-bit limbs.
 */

/*
 * Computations on scalars are supposed not to be critical for
 * performance. Therefore, implementations here are perfunctory and
 * favour code compacity over full unrolling.
 */

/*
 * Types for integers that fit on two, four, six or eight 64-bit limbs,
 * respectively. Depending on context, these types may be used with
 * signed or unsigned interpretation.
 */
typedef struct {
	uint64_t v0, v1;
} i128;
typedef struct {
	uint64_t v0, v1, v2, v3;
} i256;
typedef struct {
	uint64_t v0, v1, v2, v3, v4, v5;
} i384;
typedef struct {
	uint64_t v0, v1, v2, v3, v4, v5, v6, v7;
} i512;

/* Forward declaration of functions for reduction modulo r. These
   functions are defined in curve-specific files. */
static void modr_reduce256_partial(i256 *d, const i256 *a, uint64_t ah);
static void modr_reduce256_finish(i256 *d, const i256 *a);
static void modr_reduce384_partial(i256 *d, const i384 *a);

/*
 * Decode a 16-byte integer.
 */
UNUSED
static void
i128_decode(i128 *d, const void *a)
{
	const uint8_t *buf;

	buf = a;
	d->v0 = dec64le(buf);
	d->v1 = dec64le(buf + 8);
}

/*
 * Decode a 32-byte integer.
 */
UNUSED
static void
i256_decode(i256 *d, const void *a)
{
	const uint8_t *buf;

	buf = a;
	d->v0 = dec64le(buf);
	d->v1 = dec64le(buf + 8);
	d->v2 = dec64le(buf + 16);
	d->v3 = dec64le(buf + 24);
}

/*
 * Encode a 32-byte integer.
 */
UNUSED
static void
i256_encode(void *d, const i256 *a)
{
	uint8_t *buf;

	buf = d;
	enc64le(buf, a->v0);
	enc64le(buf + 8, a->v1);
	enc64le(buf + 16, a->v2);
	enc64le(buf + 24, a->v3);
}

/*
 * Multiply two 128-bit integers, result is truncated to 128 bits.
 */
UNUSED
static void
mul128x128trunc(i128 *d, const i128 *a, const i128 *b)
{
	uint64_t t;

	t = a->v0 * b->v1 + a->v1 * b->v0;
	UMUL64(d->v0, d->v1, a->v0, b->v0);
	d->v1 += t;
}

/*
 * Multiply two 128-bit integers, result is 256 bits.
 */
UNUSED
static void
mul128x128(i256 *d, const i128 *a, const i128 *b)
{
	uint64_t lo, hi;
	unsigned char cc;

	UMUL64(d->v0, d->v1, a->v0, b->v0);
	UMUL64(d->v2, d->v3, a->v1, b->v1);
	UMUL64(lo, hi, a->v0, b->v1);
	cc = _addcarry_u64(0, d->v1, lo, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d->v2, hi, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d->v3, 0, (unsigned long long *)&d->v3);
	UMUL64(lo, hi, a->v1, b->v0);
	cc = _addcarry_u64(0, d->v1, lo, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d->v2, hi, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d->v3, 0, (unsigned long long *)&d->v3);
}

/*
 * Multiply a 256-bit integer by a 128-bit integer, result is 384 bits.
 */
UNUSED
static void
mul256x128(i384 *d, const i256 *a, const i128 *b)
{
	i128 al, ah;
	i256 dl, dh;
	unsigned char cc;

	al.v0 = a->v0;
	al.v1 = a->v1;
	mul128x128(&dl, &al, b);
	ah.v0 = a->v2;
	ah.v1 = a->v3;
	mul128x128(&dh, &ah, b);
	d->v0 = dl.v0;
	d->v1 = dl.v1;
	cc = _addcarry_u64(0, dl.v2, dh.v0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, dl.v3, dh.v1, (unsigned long long *)&d->v3);
	cc = _addcarry_u64(cc, 0, dh.v2, (unsigned long long *)&d->v4);
	(void)_addcarry_u64(cc, 0, dh.v3, (unsigned long long *)&d->v5);
}

/*
 * Multiply two 256-bit integers together, result is 512 bits.
 */
UNUSED
static void
mul256x256(i512 *d, const i256 *a, const i256 *b)
{
	i128 al, ah;
	i384 dl, dh;
	unsigned char cc;

	al.v0 = a->v0;
	al.v1 = a->v1;
	mul256x128(&dl, b, &al);
	ah.v0 = a->v2;
	ah.v1 = a->v3;
	mul256x128(&dh, b, &ah);
	d->v0 = dl.v0;
	d->v1 = dl.v1;
	cc = _addcarry_u64(0, dl.v2, dh.v0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, dl.v3, dh.v1, (unsigned long long *)&d->v3);
	cc = _addcarry_u64(cc, dl.v4, dh.v2, (unsigned long long *)&d->v4);
	cc = _addcarry_u64(cc, dl.v5, dh.v3, (unsigned long long *)&d->v5);
	cc = _addcarry_u64(cc, 0, dh.v4, (unsigned long long *)&d->v6);
	(void)_addcarry_u64(cc, 0, dh.v5, (unsigned long long *)&d->v7);
}

/*
 * Multiplication modulo r, with one 256-bit operand and one 128-bit
 * operand. Input operands can use their full range; output is fully
 * reduced.
 */
UNUSED
static void
modr_mul256x128(i256 *d, const i256 *a, const i128 *b)
{
	i256 t;
	i384 e;

	mul256x128(&e, a, b);
	modr_reduce384_partial(&t, &e);
	modr_reduce256_finish(d, &t);
}

/*
 * Multiplication modulo r, with two 256-bit operands.
 * Input operands can use their full range; output is partially reduced.
 */
UNUSED
static void
modr_mul256x256(i256 *d, const i256 *a, const i256 *b)
{
	i256 t;
	i384 e;
	i512 x;

	mul256x256(&x, a, b);
	e.v0 = x.v2;
	e.v1 = x.v3;
	e.v2 = x.v4;
	e.v3 = x.v5;
	e.v4 = x.v6;
	e.v5 = x.v7;
	modr_reduce384_partial(&t, &e);
	e.v0 = x.v0;
	e.v1 = x.v1;
	e.v2 = t.v0;
	e.v3 = t.v1;
	e.v4 = t.v2;
	e.v5 = t.v3;
	modr_reduce384_partial(&t, &e);
	modr_reduce256_finish(d, &t);
}

/* see do255.h */
void
CN(scalar_reduce)(void *d, const void *a, size_t a_len)
{
	const uint8_t *buf;
	size_t k;
	i256 t;

	/*
	 * If the source value is shorter than 32 bytes, then we can
	 * just use it, it's already reduced.
	 */
	if (a_len < 32) {
		memmove(d, a, a_len);
		memset((unsigned char *)d + a_len, 0, 32 - a_len);
		return;
	}

	/*
	 * Decode high bytes; we use as many bytes as possible, but no
	 * more than 32, and such that the number of undecoded bytes is
	 * a multiple of 16. Also ensure that these initial contents are
	 * partially reduced.
	 */
	buf = a;
	k = a_len & 31;
	if (k == 0) {
		i256_decode(&t, buf + a_len - 32);
		buf += a_len - 32;
		modr_reduce256_partial(&t, &t, 0);
	} else if (k == 16) {
		i256_decode(&t, buf + a_len - 32);
		buf += a_len - 32;
	} else {
		uint8_t tmp[32];

		if (k < 16) {
			k += 16;
		}
		memcpy(tmp, buf + a_len - k, k);
		memset(tmp + k, 0, 32 - k);
		i256_decode(&t, tmp);
		buf += a_len - k;
	}

	/*
	 * For each subsequent chunk of 16 bytes, shift and inject the
	 * new chunk, then reduce (partially).
	 */
	while (buf != a) {
		i384 t2;

		buf -= 16;
		t2.v0 = dec64le(buf);
		t2.v1 = dec64le(buf + 8);
		t2.v2 = t.v0;
		t2.v3 = t.v1;
		t2.v4 = t.v2;
		t2.v5 = t.v3;
		modr_reduce384_partial(&t, &t2);
	}

	/*
	 * Value is almost reduced; it fits on 255 bits.
	 * A single conditional subtraction yields the correct result.
	 */
	modr_reduce256_finish(&t, &t);
	i256_encode(d, &t);
}

/* see do255.h */
void
CN(scalar_add)(void *d, const void *a, const void *b)
{
	i256 ta, tb, td;
	uint64_t t4;
	unsigned char cc;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	cc = _addcarry_u64(0, ta.v0, tb.v0, (unsigned long long *)&td.v0);
	cc = _addcarry_u64(cc, ta.v1, tb.v1, (unsigned long long *)&td.v1);
	cc = _addcarry_u64(cc, ta.v2, tb.v2, (unsigned long long *)&td.v2);
	cc = _addcarry_u64(cc, ta.v3, tb.v3, (unsigned long long *)&td.v3);
	t4 = cc;
	modr_reduce256_partial(&td, &td, t4);
	modr_reduce256_finish(&td, &td);
	i256_encode(d, &td);
}

/* see do255.h */
void
CN(scalar_neg)(void *d, const void *a)
{
	static const uint8_t zero[32] = { 0 };

	CN(scalar_sub)(d, zero, a);
}

/* see do255.h */
void
CN(scalar_mul)(void *d, const void *a, const void *b)
{
	i256 ta, tb, td;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	modr_mul256x256(&td, &ta, &tb);
	i256_encode(d, &td);
}
