/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements some operations on plain integers, used in
 * various places for computations related to scalars. It works with any
 * implementation with 32-bit limbs.
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
	uint32_t v[4];
} i128;
typedef struct {
	uint32_t v[8];
} i256;
typedef struct {
	uint32_t v[12];
} i384;
typedef struct {
	uint32_t v[16];
} i512;

/* Forward declaration of functions for reduction modulo r. These
   functions are defined in curve-specific files. */
static void modr_reduce256_partial(i256 *d, const i256 *a, uint32_t ah);
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
	int i;

	buf = a;
	for (i = 0; i < 4; i ++) {
		d->v[i] = dec32le(buf + 4 * i);
	}
}

/*
 * Decode a 32-byte integer.
 */
UNUSED
static void
i256_decode(i256 *d, const void *a)
{
	const uint8_t *buf;
	int i;

	buf = a;
	for (i = 0; i < 8; i ++) {
		d->v[i] = dec32le(buf + 4 * i);
	}
}

/*
 * Encode a 32-byte integer.
 */
UNUSED
static void
i256_encode(void *d, const i256 *a)
{
	uint8_t *buf;
	int i;

	buf = d;
	for (i = 0; i < 8; i ++) {
		enc32le(buf + 4 * i, a->v[i]);
	}
}

/*
 * Multiply two 128-bit integers, result is truncated to 128 bits.
 */
UNUSED
static void
mul128x128trunc(i128 *d, const i128 *a, const i128 *b)
{
	int i, j;
	uint32_t f, g;
	uint64_t z;
	i128 t;

	f = b->v[0];
	g = 0;
	for (i = 0; i < 4; i ++) {
		z = (uint64_t)f * (uint64_t)a->v[i] + (uint64_t)g;
		t.v[i] = (uint32_t)z;
		g = (uint32_t)(z >> 32);
	}
	for (j = 1; j < 4; j ++) {
		f = b->v[j];
		g = 0;
		for (i = 0; i < (4 - j); i ++) {
			z = (uint64_t)f * (uint64_t)a->v[i] + (uint64_t)g;
			z += (uint64_t)t.v[i + j];
			t.v[i + j] = (uint32_t)z;
			g = (uint32_t)(z >> 32);
		}
	}
	*d = t;
}

/*
 * Multiply two 128-bit integers, result is 256 bits.
 */
UNUSED
static void
mul128x128(i256 *d, const i128 *a, const i128 *b)
{
	int i, j;
	uint32_t f, g;
	uint64_t z;

	f = b->v[0];
	g = 0;
	for (i = 0; i < 4; i ++) {
		z = (uint64_t)f * (uint64_t)a->v[i] + (uint64_t)g;
		d->v[i] = (uint32_t)z;
		g = (uint32_t)(z >> 32);
	}
	d->v[4] = g;
	for (j = 1; j < 4; j ++) {
		f = b->v[j];
		g = 0;
		for (i = 0; i < 4; i ++) {
			z = (uint64_t)f * (uint64_t)a->v[i] + (uint64_t)g;
			z += (uint64_t)d->v[i + j];
			d->v[i + j] = (uint32_t)z;
			g = (uint32_t)(z >> 32);
		}
		d->v[j + 4] = g;
	}
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
	int i;

	memcpy(&al.v[0], &a->v[0], 4 * sizeof(uint32_t));
	mul128x128(&dl, &al, b);
	memcpy(&ah.v[0], &a->v[4], 4 * sizeof(uint32_t));
	mul128x128(&dh, &ah, b);
	memcpy(&d->v[0], &dl.v[0], 4 * sizeof(uint32_t));
	cc = _addcarry_u32(0, dl.v[4], dh.v[0], &d->v[4]);
	for (i = 1; i < 4; i ++) {
		cc = _addcarry_u32(cc, dl.v[4 + i], dh.v[i], &d->v[4 + i]);
	}
	for (i = 4; i < 8; i ++) {
		cc = _addcarry_u32(cc, 0, dh.v[i], &d->v[4 + i]);
	}
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
	int i;

	memcpy(&al.v[0], &a->v[0], 4 * sizeof(uint32_t));
	mul256x128(&dl, b, &al);
	memcpy(&ah.v[0], &a->v[4], 4 * sizeof(uint32_t));
	mul256x128(&dh, b, &ah);
	memcpy(&d->v[0], &dl.v[0], 4 * sizeof(uint32_t));
	cc = _addcarry_u32(0, dl.v[4], dh.v[0], &d->v[4]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dl.v[4 + i], dh.v[i], &d->v[4 + i]);
	}
	for (i = 8; i < 12; i ++) {
		cc = _addcarry_u32(cc, 0, dh.v[i], &d->v[4 + i]);
	}
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
	memcpy(&e.v[0], &x.v[4], 12 * sizeof(uint32_t));
	modr_reduce384_partial(&t, &e);
	memcpy(&e.v[0], &x.v[0], 4 * sizeof(uint32_t));
	memcpy(&e.v[4], &t.v[0], 8 * sizeof(uint32_t));
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
		t2.v[0] = dec32le(buf);
		t2.v[1] = dec32le(buf + 4);
		t2.v[2] = dec32le(buf + 8);
		t2.v[3] = dec32le(buf + 12);
		memcpy(&t2.v[4], &t.v[0], 8 * sizeof(uint32_t));
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
	uint32_t t4;
	unsigned char cc;
	int i;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	cc = _addcarry_u32(0, ta.v[0], tb.v[0], &td.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, ta.v[i], tb.v[i], &td.v[i]);
	}
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
