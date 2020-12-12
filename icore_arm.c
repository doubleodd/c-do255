/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements some operations on plain integers, used in
 * various places for computations related to scalars. It is meant for
 * 32-bit ARM platforms.
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

/*
 * Perform partial reduction of an integer slightly larger than 256 bits.
 * The extra bits are provided in 'ah'. Result is partially reduced.
 * (implemented in assembly)
 */
void CN(modr_reduce256_partial)(i256 *d, const i256 *a, uint32_t ah);
#define modr_reduce256_partial   CN(modr_reduce256_partial)

/*
 * Finish reduction modulo r of a partially reduced value.
 * (implemented in assembly)
 */
void CN(modr_reduce256_finish)(i256 *d, const i256 *a);
#define modr_reduce256_finish   CN(modr_reduce256_finish)

/*
 * Reduce a 384-bit integer modulo r. Result is partially reduced.
 * (implemented in assembly)
 */
void CN(modr_reduce384_partial)(i256 *d, const i384 *a);
#define modr_reduce384_partial   CN(modr_reduce384_partial)

/*
 * Decode a 16-byte integer.
 */
UNUSED
static inline void
i128_decode(i128 *d, const void *a)
{
	memcpy(d, a, 16);
}

/*
 * Decode a 32-byte integer.
 */
UNUSED
static inline void
i256_decode(i256 *d, const void *a)
{
	memcpy(d, a, 32);
}

/*
 * Encode a 32-byte integer.
 */
UNUSED
static inline void
i256_encode(void *d, const i256 *a)
{
	memcpy(d, a, 32);
}

/*
 * Multiply two 128-bit integers, result is truncated to 128 bits.
 * (implemented in assembly)
 */
void CN(mul128x128trunc)(i128 *d, const i128 *a, const i128 *b);
#define mul128x128trunc   CN(mul128x128trunc)

/*
 * Multiply two 128-bit integers, result is 256 bits.
 * (implemented in assembly)
 */
void CN(mul128x128)(i256 *d, const i128 *a, const i128 *b);
#define mul128x128   CN(mul128x128)

/*
 * Multiply a 256-bit integer by a 128-bit integer, result is 384 bits.
 * (implemented in assembly)
 */
void CN(mul256x128)(i384 *d, const i256 *a, const i128 *b);
#define mul256x128   CN(mul256x128)

/*
 * Multiply two 256-bit integers together, result is 512 bits.
 * (implemented in assembly)
 */
void CN(mul256x256)(i512 *d, const i256 *a, const i256 *b);
#define mul256x256   CN(mul256x256)

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
		memset(d + a_len, 0, 32 - a_len);
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
		memcpy(&t2.v[0], buf, 4 * sizeof(uint32_t));
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

/*
 * 256-bit addition; the carry is returned.
 * (implemented in assembly)
 */
uint32_t CN(i256_add)(i256 *d, const i256 *a, const i256 *b);
#define i256_add   CN(i256_add)

/*
 * 256-bit subtraction; the borrow is returned (1 if borrow, 0 otherwise).
 * (implemented in assembly)
 */
uint32_t CN(i256_sub)(i256 *d, const i256 *a, const i256 *b);
#define i256_sub   CN(i256_sub)

/* see do255.h */
void
CN(scalar_add)(void *d, const void *a, const void *b)
{
	i256 ta, tb, td;
	uint32_t t4;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	t4 = i256_add(&td, &ta, &tb);
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
