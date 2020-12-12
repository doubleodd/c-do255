/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined CURVE to the curve name
 *
 * This file implements scalar operations modulo r for curve do255e; it
 * works with 32-bit ARM implementations.
 */

/*
 * Given input 'a' (up to 2^258-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and is less than 2^254 + 2^131. The
 * high bits of 'a' are provided as extra parameter ah.
 * (implemented in assembly)
 */
void CN(modr_reduce256_partial)(i256 *d, const i256 *a, uint32_t ah);
#define modr_reduce256_partial   CN(modr_reduce256_partial)

/*
 * Given a partially reduced input (less than 2*r), finish reduction
 * (conditional subtraction of r).
 * (implemented in assembly)
 */
void CN(modr_reduce256_finish)(i256 *d, const i256 *a);
#define modr_reduce256_finish   CN(modr_reduce256_finish)

/*
 * Given input 'a' (up to 2^384-1), perform a partial reduction modulo r;
 * output (into 'd') fits on 255 bits and is less than 2^254 + 2^131.
 * (implemented in assembly)
 */
void CN(modr_reduce384_partial)(i256 *d, const i384 *a);
#define modr_reduce384_partial   CN(modr_reduce384_partial)

/*
 * 8*r modulo 2^256.
 */
static const i256 R_x8 = { {
	0xA6C22928, 0xFA964573, 0xA03C6298, 0xE864987A,
	0xFFFFFFFC, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF
} };

/* do255e_scalar_is_reduced() is implemented in assembly */

/* see do255.h */
void
do255e_scalar_sub(void *d, const void *a, const void *b)
{
	i256 ta, tb;
	uint32_t t8, m, c;
	int i;

	i256_decode(&ta, a);
	i256_decode(&tb, b);
	c = i256_sub(&tb, &ta, &tb);

	/*
	 * If subtraction yielded a negative value, then we add 8*r. Note
	 * that 2^256 < 8*r < 2^257; thus, we always get a nonnegative
	 * value that fits on 257 bits.
	 */
	m = -c;
	for (i = 0; i < 8; i ++) {
		ta.v[i] = m & R_x8.v[i];
	}
	t8 = i256_add(&tb, &tb, &ta);

	/*
	 * Reduce modulo r.
	 */
	modr_reduce256_partial(&tb, &tb, t8);
	modr_reduce256_finish(&tb, &tb);
	i256_encode(d, &tb);
}

/*
 * 1/2 mod r = (r+1)/2
 */
static const i256 Rhf = { {
	0x3A6C2293, 0x8FA96457, 0xAA03C629, 0xCE864987,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x1FFFFFFF
} };

/* see do255.h */
void
do255e_scalar_half(void *d, const void *a)
{
	i256 x, t;
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
	for (i = 0; i < 8; i ++) {
		t.v[i] = m & Rhf.v[i];
	}
	i256_add(&x, &x, &t);

	/*
	 * Ensure a reduced value.
	 */
	modr_reduce256_partial(&x, &x, 0);
	modr_reduce256_finish(&x, &x);
	i256_encode(d, &x);
}
