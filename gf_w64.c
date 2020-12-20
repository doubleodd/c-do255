/*
 * This file implements elementary operations in finite field GF(2^255-MQ),
 * with no inline assembly. The _addcarry_u64() and _subborrow_u64()
 * intrinsics are used to get 64-bit additions and subtractions with
 * carries; for 64x64->128 multiplications, the 'unsigned __int128' type
 * or the _umul128() intrinsic are used, depending on the local compiler
 * (MSVC does not support 'unsigned __int128').
 *
 * It is not meant to be compiled by itself, but included in an outer
 * file. The outer file must provide the following:
 *
 *    MQ            macro that evaluates to a literal constant
 *    GF_INVT508    value 1/2^508 mod 2^255-MQ
 *
 * Moreover, the common header "do255.h" must have already been included.
 * Value GF_INVT508 is a 'struct do255_int256_w64' and the value is
 * represented in four 64-bit limbs (v0 is lowest limb, v3 is highest).
 */

#include <immintrin.h>

typedef struct do255_int256_w64 gf;

/*
 * Constants below need not be modified if MQ is changed.
 */
static const gf GF_ZERO = { 0, 0, 0, 0 };
static const gf GF_ONE = { 1, 0, 0, 0 };
static const gf GF_TWO = { 2, 0, 0, 0 };
static const gf GF_SEVEN = { 7, 0, 0, 0 };
static const gf GF_P = {
	(uint64_t)-MQ, (uint64_t)-1, (uint64_t)-1, (uint64_t)-1 >> 1
};

/*
 * Helpful macro to convert expressions to strings, for inline assembly.
 */
#define ST(x)    ST_(x)
#define ST_(x)   #x

#if defined __GNUC__ || defined __clang__
#define UMUL64(lo, hi, x, y)   do { \
		unsigned __int128 umul64_tmp; \
		umul64_tmp = (unsigned __int128)(x) * (unsigned __int128)(y); \
		(lo) = (uint64_t)umul64_tmp; \
		(hi) = (uint64_t)(umul64_tmp >> 64); \
	} while (0)
#define UMUL64x2(lo, hi, x1, y1, x2, y2)   do { \
		unsigned __int128 umul64_tmp; \
		umul64_tmp = (unsigned __int128)(x1) * (unsigned __int128)(y1) \
			+ (unsigned __int128)(x2) * (unsigned __int128)(y2); \
		(lo) = (uint64_t)umul64_tmp; \
		(hi) = (uint64_t)(umul64_tmp >> 64); \
	} while (0)
#define UMUL64x2_ADD(lo, hi, x1, y1, x2, y2, z3)   do { \
		unsigned __int128 umul64_tmp; \
		umul64_tmp = (unsigned __int128)(x1) * (unsigned __int128)(y1) \
			+ (unsigned __int128)(x2) * (unsigned __int128)(y2) \
			+ (unsigned __int128)(z3); \
		(lo) = (uint64_t)umul64_tmp; \
		(hi) = (uint64_t)(umul64_tmp >> 64); \
	} while (0)
#else
#define UMUL64(lo, hi, x, y)   do { \
		uint64_t umul64_hi; \
		(lo) = _umul128((x), (y), &umul64_hi); \
		(hi) = umul64_hi; \
	} while (0)
#define UMUL64x2(lo, hi, x1, y1, x2, y2)   do { \
		uint64_t umul64_lo1, umul64_lo2, umul64_hi1, umul64_hi2; \
		unsigned char umul64_cc; \
		umul64_lo1 = _umul128((x1), (y1), &umul64_hi1); \
		umul64_lo2 = _umul128((x2), (y2), &umul64_hi2); \
		umul64_cc = _addcarry_u64(0, \
			umul64_lo1, umul64_lo2, &umul64_lo1); \
		(void)_addcarry_u64(umul64_cc, \
			umul64_hi1, umul64_hi2, &umul64_hi1); \
		(lo) = umul64_lo1; \
		(hi) = umul64_hi1; \
	} while (0)
#define UMUL64x2_ADD(lo, hi, x1, y1, x2, y2, z3)   do { \
		uint64_t umul64_lo1, umul64_lo2, umul64_hi1, umul64_hi2; \
		unsigned char umul64_cc; \
		umul64_lo1 = _umul128((x1), (y1), &umul64_hi1); \
		umul64_lo2 = _umul128((x2), (y2), &umul64_hi2); \
		umul64_cc = _addcarry_u64(0, \
			umul64_lo1, umul64_lo2, &umul64_lo1); \
		(void)_addcarry_u64(umul64_cc, \
			umul64_hi1, umul64_hi2, &umul64_hi1); \
		umul64_cc = _addcarry_u64(0, umul64_lo1, (z3), &umul64_lo1); \
		(void)_addcarry_u64(umul64_cc, umul64_hi1, 0, &umul64_hi1); \
		(lo) = umul64_lo1; \
		(hi) = umul64_hi1; \
	} while (0)
#endif

/*
 * A field element is represented as four limbs, in base 2^64. Operands
 * and result may be up to 2^256-1.
 * A _normalized_ value is in 0..p-1. gf_normalize() ensures normalization;
 * it is called before encoding, and also when starting inversion.
 */

/* d <- a + b */
UNUSED
static void
gf_add(gf *d, const gf *a, const gf *b)
{
	unsigned long long d0, d1, d2, d3;
	unsigned char cc;

	/*
	 * Compute sum over 256 bits + carry.
	 */
	cc = _addcarry_u64(0, a->v0, b->v0, &d0);
	cc = _addcarry_u64(cc, a->v1, b->v1, &d1);
	cc = _addcarry_u64(cc, a->v2, b->v2, &d2);
	cc = _addcarry_u64(cc, a->v3, b->v3, &d3);

	/*
	 * If there is a carry, subtract 2*p, i.e. add 2*MQ and
	 * (implicitly) subtract 2^256.
	 */
	cc = _addcarry_u64(0, d0, -(unsigned long long)cc & (2 * MQ), &d0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);

	/*
	 * Usually, there won't be an extra carry here, and the
	 * implicit subtraction of 2^256 is enough. However, if the
	 * the original sum was at least 2^257 - 2*MQ, then this
	 * will show up as an extra carry here. In that case, the
	 * low limb (d0) is necessarily less than 2*MQ, and thus adding
	 * back 2*MQ to it will not trigger any carry.
	 */
	(void)_addcarry_u64(0, d0, -(unsigned long long)cc & (2 * MQ),
		(unsigned long long *)&d->v0);
}

/* d <- a - b */
UNUSED
static void
gf_sub(gf *d, const gf *a, const gf *b)
{
	unsigned long long d0, d1, d2, d3, e;
	unsigned char cc;

	/*
	 * Compute subtraction over 256 bits.
	 */
	cc = _subborrow_u64(0, a->v0, b->v0, &d0);
	cc = _subborrow_u64(cc, a->v1, b->v1, &d1);
	cc = _subborrow_u64(cc, a->v2, b->v2, &d2);
	cc = _subborrow_u64(cc, a->v3, b->v3, &d3);

	/*
	 * If the result is negative, add back 2*p = 2^256 - 2*MQ (which
	 * we do by subtracting 2*MQ).
	 */
	e = -(unsigned long long)cc;
	cc = _subborrow_u64(0, d0, e & (2 * MQ), &d0);
	cc = _subborrow_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _subborrow_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _subborrow_u64(cc, d3, 0, (unsigned long long *)&d->v3);

	/*
	 * If there is still a carry, then we must add 2*p again. In that
	 * case, subtracting 2*MQ from the low limb won't trigger a carry.
	 */
	d->v0 = d0 - (-(unsigned long long)cc & (2 * MQ));
}

/* d <- -a */
UNUSED
static void
gf_neg(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, e;
	unsigned char cc;

	/*
	 * Compute 2*p - a over 256 bits.
	 */
	cc = _subborrow_u64(0, (unsigned long long)-(2 * MQ), a->v0, &d0);
	cc = _subborrow_u64(cc, (unsigned long long)-1, a->v1, &d1);
	cc = _subborrow_u64(cc, (unsigned long long)-1, a->v2, &d2);
	cc = _subborrow_u64(cc, (unsigned long long)-1, a->v3, &d3);

	/*
	 * If the result is negative, add back p = 2^255 - MQ.
	 */
	e = -(unsigned long long)cc;
	cc = _addcarry_u64(0, d0, e & -MQ, (unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, d1, e, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, e, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d3, e >> 1, (unsigned long long *)&d->v3);
}

/* If ctl = 1, then d <- -a
   If ctl = 0, then d is unchanged
   ctl MUST be 0 or 1 */
UNUSED
static void
gf_condneg(gf *d, const gf *a, uint64_t ctl)
{
	gf t;

	gf_neg(&t, a);
	d->v0 = a->v0 ^ (-ctl & (a->v0 ^ t.v0));
	d->v1 = a->v1 ^ (-ctl & (a->v1 ^ t.v1));
	d->v2 = a->v2 ^ (-ctl & (a->v2 ^ t.v2));
	d->v3 = a->v3 ^ (-ctl & (a->v3 ^ t.v3));
}

/* d <- a - b - c */
UNUSED
static void
gf_sub2(gf *d, const gf *a, const gf *b, const gf *c)
{
	/* TODO: see if this can be optimized a bit */
	gf t;

	gf_sub(&t, a, b);
	gf_sub(d, &t, c);
}

/* d <- a/2 */
UNUSED
static void
gf_half(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, tt;
	unsigned char cc;

	d0 = (a->v0 >> 1) | (a->v1 << 63);
	d1 = (a->v1 >> 1) | (a->v2 << 63);
	d2 = (a->v2 >> 1) | (a->v3 << 63);
	d3 = (a->v3 >> 1);
	tt = -(a->v0 & 1);

	cc = _addcarry_u64(0, d0, tt & (uint64_t)-((MQ - 1) >> 1),
		(unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, d1, tt, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, tt, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d3, tt >> 2, (unsigned long long *)&d->v3);
}

/* If ctl = 1, then d <- a
   If ctl = 0, then d <- b
   ctl MUST be 0 or 1 */
UNUSED
static inline void
gf_sel2(gf *d, const gf *a, const gf *b, uint64_t ctl)
{
	d->v0 = b->v0 ^ (-ctl & (a->v0 ^ b->v0));
	d->v1 = b->v1 ^ (-ctl & (a->v1 ^ b->v1));
	d->v2 = b->v2 ^ (-ctl & (a->v2 ^ b->v2));
	d->v3 = b->v3 ^ (-ctl & (a->v3 ^ b->v3));
}

/*
 * Set d to a, b or c, depending on flags:
 *  - if fa = 1, then set d to a
 *  - otherwise, if fb = 1, then set d to b
 *  - otherwise, set d to c
 * Flags fa and fb MUST have value 0 or 1. If both are 1, fa takes
 * precedence (c is set to a).
 */
UNUSED
static void
gf_sel3(gf *d, const gf *a, const gf *b, const gf *c, uint64_t fa, uint64_t fb)
{
	uint64_t ma, mb, mc;

	ma = -fa;
	mb = -fb & ~ma;
	mc = ~(ma | mb);
	d->v0 = (a->v0 & ma) | (b->v0 & mb) | (c->v0 & mc);
	d->v1 = (a->v1 & ma) | (b->v1 & mb) | (c->v1 & mc);
	d->v2 = (a->v2 & ma) | (b->v2 & mb) | (c->v2 & mc);
	d->v3 = (a->v3 & ma) | (b->v3 & mb) | (c->v3 & mc);
}

/* d <- 2*a */
UNUSED
static void
gf_mul2(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, tt;
	unsigned char cc;

	d0 = (a->v0 << 1);
	d1 = (a->v1 << 1) | (a->v0 >> 63);
	d2 = (a->v2 << 1) | (a->v1 >> 63);
	d3 = (a->v3 << 1) | (a->v2 >> 63);

	tt = (*(int64_t *)&a->v3 >> 63) & (uint64_t)(2 * MQ);
	cc = _addcarry_u64(0, d0, tt, &d0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);
	d->v0 = d0 + (-(unsigned long long)cc & (2 * MQ));
}

/* d <- 4*a */
UNUSED
static void
gf_mul4(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, tt;
	unsigned char cc;

	d0 = (a->v0 << 2);
	d1 = (a->v1 << 2) | (a->v0 >> 62);
	d2 = (a->v2 << 2) | (a->v1 >> 62);
	d3 = (a->v3 << 2) | (a->v2 >> 62);

	tt = (a->v3 >> 62) * (uint64_t)(2 * MQ);
	cc = _addcarry_u64(0, d0, tt, &d0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);
	d->v0 = d0 + (-(unsigned long long)cc & (2 * MQ));
}

/* d <- 8*a */
UNUSED
static void
gf_mul8(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, tt;
	unsigned char cc;

	d0 = (a->v0 << 3);
	d1 = (a->v1 << 3) | (a->v0 >> 61);
	d2 = (a->v2 << 3) | (a->v1 >> 61);
	d3 = (a->v3 << 3) | (a->v2 >> 61);

	tt = (a->v3 >> 61) * (uint64_t)(2 * MQ);
	cc = _addcarry_u64(0, d0, tt, &d0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);
	d->v0 = d0 + (-(unsigned long long)cc & (2 * MQ));
}

/* d <- 32*a */
UNUSED
static void
gf_mul32(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, tt;
	unsigned char cc;

	d0 = (a->v0 << 5);
	d1 = (a->v1 << 5) | (a->v0 >> 59);
	d2 = (a->v2 << 5) | (a->v1 >> 59);
	d3 = (a->v3 << 5) | (a->v2 >> 59);

	tt = (a->v3 >> 59) * (uint64_t)(2 * MQ);
	cc = _addcarry_u64(0, d0, tt, &d0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	cc = _addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);
	d->v0 = d0 + (-(unsigned long long)cc & (2 * MQ));
}

/* d <- a*b, with b < 2^32 */
UNUSED
static void
gf_mul_small(gf *d, const gf *a, uint64_t b)
{
	unsigned long long d0, d1, d2, d3, d4, lo, hi;
	unsigned char cc;

	/* Multiplication into d0..d4 */
	UMUL64(d0, d1, a->v0, b);
	UMUL64(d2, d3, a->v2, b);
	UMUL64(lo, hi, a->v1, b);
	cc = _addcarry_u64(0, d1, lo, &d1);
	cc = _addcarry_u64(cc, d2, hi, &d2);
	UMUL64(lo, d4, a->v3, b);
	cc = _addcarry_u64(cc, d3, lo, &d3);
	(void)_addcarry_u64(cc, d4, 0, &d4);

	/* Fold upper bits. */
	d4 = (d4 << 1) | (d3 >> 63);
	d3 &= 0x7FFFFFFFFFFFFFFF;
	d4 *= MQ;
	cc = _addcarry_u64(0, d0, d4, (unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, d1, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, 0, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d3, 0, (unsigned long long *)&d->v3);
}

/* d <- a*b  (always inlined) */
FORCE_INLINE
static inline void
gf_mul_inline(gf *d, const gf *a, const gf *b)
{
	unsigned long long e0, e1, e2, e3, e4, e5, e6, e7;
	unsigned long long h0, h1, h2, h3;
	unsigned long long lo, hi, lo2, hi2;
	unsigned char cc;

	UMUL64(e0, e1, a->v0, b->v0);
	UMUL64(e2, e3, a->v1, b->v1);
	UMUL64(e4, e5, a->v2, b->v2);
	UMUL64(e6, e7, a->v3, b->v3);

	UMUL64(lo, hi, a->v0, b->v1);
	cc = _addcarry_u64(0, e1, lo, &e1);
	cc = _addcarry_u64(cc, e2, hi, &e2);
	UMUL64(lo, hi, a->v0, b->v3);
	cc = _addcarry_u64(cc, e3, lo, &e3);
	cc = _addcarry_u64(cc, e4, hi, &e4);
	UMUL64(lo, hi, a->v2, b->v3);
	cc = _addcarry_u64(cc, e5, lo, &e5);
	cc = _addcarry_u64(cc, e6, hi, &e6);
	(void)_addcarry_u64(cc, e7, 0, &e7);

	UMUL64(lo, hi, a->v1, b->v0);
	cc = _addcarry_u64(0, e1, lo, &e1);
	cc = _addcarry_u64(cc, e2, hi, &e2);
	UMUL64(lo, hi, a->v3, b->v0);
	cc = _addcarry_u64(cc, e3, lo, &e3);
	cc = _addcarry_u64(cc, e4, hi, &e4);
	UMUL64(lo, hi, a->v3, b->v2);
	cc = _addcarry_u64(cc, e5, lo, &e5);
	cc = _addcarry_u64(cc, e6, hi, &e6);
	(void)_addcarry_u64(cc, e7, 0, &e7);

	UMUL64(lo, hi, a->v0, b->v2);
	cc = _addcarry_u64(0, e2, lo, &e2);
	cc = _addcarry_u64(cc, e3, hi, &e3);
	UMUL64(lo, hi, a->v1, b->v3);
	cc = _addcarry_u64(cc, e4, lo, &e4);
	cc = _addcarry_u64(cc, e5, hi, &e5);
	cc = _addcarry_u64(cc, e6, 0, &e6);
	(void)_addcarry_u64(cc, e7, 0, &e7);

	UMUL64(lo, hi, a->v2, b->v0);
	cc = _addcarry_u64(0, e2, lo, &e2);
	cc = _addcarry_u64(cc, e3, hi, &e3);
	UMUL64(lo, hi, a->v3, b->v1);
	cc = _addcarry_u64(cc, e4, lo, &e4);
	cc = _addcarry_u64(cc, e5, hi, &e5);
	cc = _addcarry_u64(cc, e6, 0, &e6);
	(void)_addcarry_u64(cc, e7, 0, &e7);

	UMUL64(lo, hi, a->v1, b->v2);
	UMUL64(lo2, hi2, a->v2, b->v1);
	cc = _addcarry_u64(0, lo, lo2, &lo);
	cc = _addcarry_u64(cc, hi, hi2, &hi);
	(void)_addcarry_u64(cc, 0, 0, &hi2);
	cc = _addcarry_u64(0, e3, lo, &e3);
	cc = _addcarry_u64(cc, e4, hi, &e4);
	cc = _addcarry_u64(cc, e5, hi2, &e5);
	cc = _addcarry_u64(cc, e6, 0, &e6);
	(void)_addcarry_u64(cc, e7, 0, &e7);

	UMUL64(lo, h0, e4, 2 * MQ);
	cc = _addcarry_u64(0, e0, lo, &e0);
	UMUL64(lo, h1, e5, 2 * MQ);
	cc = _addcarry_u64(cc, e1, lo, &e1);
	UMUL64(lo, h2, e6, 2 * MQ);
	cc = _addcarry_u64(cc, e2, lo, &e2);
	UMUL64(lo, h3, e7, 2 * MQ);
	cc = _addcarry_u64(cc, e3, lo, &e3);
	(void)_addcarry_u64(cc, 0, h3, &h3);

	h3 = (h3 << 1) | (e3 >> 63);
	e3 &= 0x7FFFFFFFFFFFFFFF;
	cc = _addcarry_u64(0, e0, h3 * MQ, &e0);
	cc = _addcarry_u64(cc, e1, h0, &e1);
	cc = _addcarry_u64(cc, e2, h1, &e2);
	(void)_addcarry_u64(cc, e3, h2, &e3);

	d->v0 = e0;
	d->v1 = e1;
	d->v2 = e2;
	d->v3 = e3;
}

/* d <- a*b  (never inlined) */
NO_INLINE UNUSED
static void
gf_mul(gf *d, const gf *a, const gf *b)
{
	gf_mul_inline(d, a, b);
}

/* d <- a^2  (always inlined) */
FORCE_INLINE
static inline void
gf_sqr_inline(gf *d, const gf *a)
{
	unsigned long long e0, e1, e2, e3, e4, e5, e6, e7;
	unsigned long long h0, h1, h2, h3;
	unsigned long long lo, hi;
	unsigned char cc;

	UMUL64(e1, e2, a->v0, a->v1);
	UMUL64(e3, e4, a->v0, a->v3);
	UMUL64(e5, e6, a->v2, a->v3);
	UMUL64(lo, hi, a->v0, a->v2);
	cc = _addcarry_u64(0, e2, lo, &e2);
	cc = _addcarry_u64(cc, e3, hi, &e3);
	UMUL64(lo, hi, a->v1, a->v3);
	cc = _addcarry_u64(cc, e4, lo, &e4);
	cc = _addcarry_u64(cc, e5, hi, &e5);
	(void)_addcarry_u64(cc, e6, 0, &e6);
	UMUL64(lo, hi, a->v1, a->v2);
	cc = _addcarry_u64(0, e3, lo, &e3);
	cc = _addcarry_u64(cc, e4, hi, &e4);
	cc = _addcarry_u64(cc, e5, 0, &e5);
	(void)_addcarry_u64(cc, e6, 0, &e6);

	/* There cannot be extra carry here because the partial sum is
	   necessarily lower than 2^448 at this point. */

	e7 = e6 >> 63;
	e6 = (e6 << 1) | (e5 >> 63);
	e5 = (e5 << 1) | (e4 >> 63);
	e4 = (e4 << 1) | (e3 >> 63);
	e3 = (e3 << 1) | (e2 >> 63);
	e2 = (e2 << 1) | (e1 >> 63);
	e1 = e1 << 1;

	UMUL64(e0, hi, a->v0, a->v0);
	cc = _addcarry_u64(0, e1, hi, &e1);
	UMUL64(lo, hi, a->v1, a->v1);
	cc = _addcarry_u64(cc, e2, lo, &e2);
	cc = _addcarry_u64(cc, e3, hi, &e3);
	UMUL64(lo, hi, a->v2, a->v2);
	cc = _addcarry_u64(cc, e4, lo, &e4);
	cc = _addcarry_u64(cc, e5, hi, &e5);
	UMUL64(lo, hi, a->v3, a->v3);
	cc = _addcarry_u64(cc, e6, lo, &e6);
	(void)_addcarry_u64(cc, e7, hi, &e7);

	/* Reduction */

	UMUL64(lo, h0, e4, 2 * MQ);
	cc = _addcarry_u64(0, e0, lo, &e0);
	UMUL64(lo, h1, e5, 2 * MQ);
	cc = _addcarry_u64(cc, e1, lo, &e1);
	UMUL64(lo, h2, e6, 2 * MQ);
	cc = _addcarry_u64(cc, e2, lo, &e2);
	UMUL64(lo, h3, e7, 2 * MQ);
	cc = _addcarry_u64(cc, e3, lo, &e3);
	(void)_addcarry_u64(cc, 0, h3, &h3);

	h3 = (h3 << 1) | (e3 >> 63);
	e3 &= 0x7FFFFFFFFFFFFFFF;
	cc = _addcarry_u64(0, e0, h3 * MQ, &e0);
	cc = _addcarry_u64(cc, e1, h0, &e1);
	cc = _addcarry_u64(cc, e2, h1, &e2);
	(void)_addcarry_u64(cc, e3, h2, &e3);

	d->v0 = e0;
	d->v1 = e1;
	d->v2 = e2;
	d->v3 = e3;
}

/* d <- a^2  (never inlined) */
NO_INLINE UNUSED
static void
gf_sqr(gf *d, const gf *a)
{
	gf_sqr_inline(d, a);
}

/* d <- a^(2^num)  (repeated squarings, always inlined) */
FORCE_INLINE
static void
gf_sqr_x_inline(gf *d, const gf *a, long num)
{
	gf t;

	t = *a;
	do {
		gf_sqr_inline(&t, &t);
	} while (-- num > 0);
	*d = t;
}

/* d <- a^(2^num)  (repeated squarings, never inlined) */
NO_INLINE UNUSED
static void
gf_sqr_x(gf *d, const gf *a, long num)
{
	gf_sqr_x_inline(d, a, num);
}

/*
 * Normalize a value to 0..p-1.
 */
UNUSED
static void
gf_normalize(gf *d, const gf *a)
{
	unsigned long long d0, d1, d2, d3, e;
	unsigned char cc;

	/*
	 * If top bit is set, propagate it.
	 */
	e = -(unsigned long long)(a->v3 >> 63);
	cc = _addcarry_u64(0, a->v0, e & MQ, &d0);
	cc = _addcarry_u64(cc, a->v1, 0, &d1);
	cc = _addcarry_u64(cc, a->v2, 0, &d2);
	(void)_addcarry_u64(cc, a->v3, e << 63, &d3);

	/*
	 * Value is now at most 2^255+MQ-1. Subtract p, and add it back
	 * if the result is negative.
	 */
	cc = _subborrow_u64(0, d0, (unsigned long long)-MQ, &d0);
	cc = _subborrow_u64(cc, d1, (unsigned long long)-1, &d1);
	cc = _subborrow_u64(cc, d2, (unsigned long long)-1, &d2);
	cc = _subborrow_u64(cc, d3, (unsigned long long)-1 >> 1, &d3);

	e = -(unsigned long long)cc;
	cc = _addcarry_u64(0, d0, e & -MQ, (unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, d1, e, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2, e, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d3, e >> 1, (unsigned long long *)&d->v3);
}

/*
 * Compare 'a' with zero; returned value is 1 if the value is zero (as
 * a field element), 0 otherwise.
 */
UNUSED
static uint64_t
gf_iszero(const gf *a)
{
	unsigned long long t0, t1, t2;

	/*
	 * Since values are over 256 bits, there are three possible
	 * representations for 0: 0, p and 2*p.
	 */
	t0 = a->v0 | a->v1 | a->v2 | a->v3;
	t1 = (a->v0 + MQ) | ~a->v1 | ~a->v2 | (a->v3 ^ 0x7FFFFFFFFFFFFFFF);
	t2 = (a->v0 + (2 * MQ)) | ~a->v1 | ~a->v2 | ~a->v3;

	/*
	 * Value is zero if and only if one of t0, t1 or t2 is zero.
	 */
	return 1 - (((t0 | -t0) & (t1 | -t1) & (t2 | -t2)) >> 63);
}

/*
 * Compare two values as field elements; returned value is 1 on equality,
 * 0 otherwise.
 */
UNUSED
static uint64_t
gf_eq(const gf *a, const gf *b)
{
	gf t;

	gf_sub(&t, a, b);
	return gf_iszero(&t);
}

/*
 * Compute (a*f+v*g)/2^31. Parameters f and g are provided with an
 * unsigned type, but they are signed integers in the -2^31..+2^31 range.
 * Values a, b and d are not field elements, but 255-bit integers (value
 * is in 0..2^255-1, top bit is zero). The division by 2^31 is assumed to
 * be exact (low 31 bits of a*f+v*g are dropped). The result is assumed to
 * fit in 256 bits in signed two's complement notation (truncation is
 * applied on higher bits).
 *
 * If the result turns out to be negative, then it is negated. Returned
 * value is 1 if the result was negated, 0 otherwise.
 */
static uint64_t
s256_lin_div31_abs(gf *d, const gf *a, const gf *b,
	unsigned long long f, unsigned long long g)
{
	gf ta, tb;
	unsigned long long sf, sg, d0, d1, d2, d3, t;
	unsigned char cc;

	/*
	 * If f < 0, replace f with -f but keep the sign in sf.
	 * Similarly for g.
	 */
	sf = f >> 63;
	f = (f ^ -sf) + sf;
	sg = g >> 63;
	g = (g ^ -sg) + sg;

	/*
	 * Apply signs sf and sg to a and b, respectively.
	 */
	cc = _addcarry_u64(0, a->v0 ^ -sf, sf, (unsigned long long *)&ta.v0);
	cc = _addcarry_u64(cc, a->v1 ^ -sf, 0, (unsigned long long *)&ta.v1);
	cc = _addcarry_u64(cc, a->v2 ^ -sf, 0, (unsigned long long *)&ta.v2);
	cc = _addcarry_u64(cc, a->v3 ^ -sf, 0, (unsigned long long *)&ta.v3);

	cc = _addcarry_u64(0, b->v0 ^ -sg, sg, (unsigned long long *)&tb.v0);
	cc = _addcarry_u64(cc, b->v1 ^ -sg, 0, (unsigned long long *)&tb.v1);
	cc = _addcarry_u64(cc, b->v2 ^ -sg, 0, (unsigned long long *)&tb.v2);
	cc = _addcarry_u64(cc, b->v3 ^ -sg, 0, (unsigned long long *)&tb.v3);

	/*
	 * Now that f and g are nonnegative, compute a*f+b*g into
	 * d0:d1:d2:d3:t. Since f and g are at most 2^31, we can
	 * add two 128-bit products with no overflow (they are actually
	 * 95 bits each at most).
	 */
	UMUL64x2(d0, t, ta.v0, f, tb.v0, g);
	UMUL64x2_ADD(d1, t, ta.v1, f, tb.v1, g, t);
	UMUL64x2_ADD(d2, t, ta.v2, f, tb.v2, g, t);
	UMUL64x2_ADD(d3, t, ta.v3, f, tb.v3, g, t);

	/*
	 * Don't forget the signs: if a < 0, then the result is
	 * overestimated by 2^256*f; similarly, if b < 0, then the
	 * result is overestimated by 2^256*g. We thus must subtract
	 * 2^256*(sa*f+sb*g), where sa and sb are the signs of a and b,
	 * respectively (1 if negative, 0 otherwise).
	 */
	t -= -(unsigned long long)(ta.v3 >> 63) & f;
	t -= -(unsigned long long)(tb.v3 >> 63) & g;

	/*
	 * Apply the shift.
	 */
	d0 = (d0 >> 31) | (d1 << 33);
	d1 = (d1 >> 31) | (d2 << 33);
	d2 = (d2 >> 31) | (d3 << 33);
	d3 = (d3 >> 31) | (t << 33);

	/*
	 * Perform conditional negation, if the result is negative.
	 */
	t >>= 63;
	cc = _addcarry_u64(0, d0 ^ -t, t, (unsigned long long *)&d->v0);
	cc = _addcarry_u64(cc, d1 ^ -t, 0, (unsigned long long *)&d->v1);
	cc = _addcarry_u64(cc, d2 ^ -t, 0, (unsigned long long *)&d->v2);
	(void)_addcarry_u64(cc, d3 ^ -t, 0, (unsigned long long *)&d->v3);

	return t;
}

/*
 * Compute u*f+v*g (modulo p). Parameters f and g are provided with
 * an unsigned type, but they are signed integers in the -2^62..+2^62 range.
 */
static void
gf_lin(gf *d, const gf *u, const gf *v,
	unsigned long long f, unsigned long long g)
{
	gf tu, tv;
	unsigned long long sf, sg, d0, d1, d2, d3, t, lo, hi;
	unsigned char cc;

	/*
	 * If f < 0, replace f with -f but keep the sign in sf.
	 * Similarly for g.
	 */
	sf = f >> 63;
	f = (f ^ -sf) + sf;
	sg = g >> 63;
	g = (g ^ -sg) + sg;

	/*
	 * Apply signs sf and sg to u and v.
	 */
	gf_condneg(&tu, u, sf);
	gf_condneg(&tv, v, sg);

	/*
	 * Since f <= 2^62 and g <= 2^62, we can add 128-bit products: they
	 * fit on 126 bits each.
	 */
	UMUL64x2(d0, t, tu.v0, f, tv.v0, g);
	UMUL64x2_ADD(d1, t, tu.v1, f, tv.v1, g, t);
	UMUL64x2_ADD(d2, t, tu.v2, f, tv.v2, g, t);
	UMUL64x2_ADD(d3, t, tu.v3, f, tv.v3, g, t);

	/*
	 * Upper word t can be up to 63 bits.
	 */
	UMUL64(lo, hi, t, 2 * MQ);
	cc = _addcarry_u64(0, d0, lo, &d0);
	cc = _addcarry_u64(cc, d1, hi, &d1);
	cc = _addcarry_u64(cc, d2, 0, &d2);
	cc = _addcarry_u64(cc, d3, 0, &d3);

	/*
	 * There could be a carry, but in that case, the folding won't
	 * propagate beyond the second limb.
	 */
	cc = _addcarry_u64(0, d0, -(unsigned long long)cc & (2 * MQ), &d0);
	(void)_addcarry_u64(cc, d1, 0, &d1);

	d->v0 = d0;
	d->v1 = d1;
	d->v2 = d2;
	d->v3 = d3;
}

/*
 * Inversion in the field: d <- 1/y
 * If y = 0, then d is set to zero.
 * Returned value is 1 if the value was invertible, 0 otherwise.
 */
UNUSED
static uint64_t
gf_inv(gf *d, const gf *y)
{
	gf a, b, u, v;
	unsigned long long f0, f1, g0, g1, xa, xb, fg0, fg1;
	unsigned long long nega, negb;
	int i, j;

	/*
	 * Extended binary GCD:
	 *
	 *   a <- y
	 *   b <- q
	 *   u <- 1
	 *   v <- 0
	 *
	 * a and b are nonnnegative integers (in the 0..q range). u
	 * and v are integers modulo q.
	 *
	 * Invariants:
	 *    a = y*u mod q
	 *    b = y*v mod q
	 *    b is always odd
	 *
	 * At each step:
	 *    if a is even, then:
	 *        a <- a/2, u <- u/2 mod q
	 *    else:
	 *        if a < b:
	 *            (a, u, b, v) <- (b, v, a, u)
	 *        a <- (a-b)/2, u <- (u-v)/2 mod q
	 *
	 * At one point, value a reaches 0; it will then stay there
	 * (all subsequent steps will only keep a at zero). The value
	 * of b is then GCD(y, p), i.e. 1 if y is invertible; value v
	 * is then the inverse of y modulo p.
	 *
	 * If y = 0, then a is initially at zero and stays there, and
	 * v is unchanged. The function then returns 0 as "inverse" of
	 * zero, which is the desired behaviour; therefore, no corrective
	 * action is necessary.
	 *
	 * It can be seen that each step will decrease the size of one of
	 * a and b by at least 1 bit; thus, starting with two values of
	 * at most 255 bits each, 509 iterations are always sufficient.
	 *
	 *
	 * In practice, we optimize this code in the following way:
	 *  - We do iterations by group of 31.
	 *  - In each group, we use _approximations_ of a and b that
	 *    fit on 64 bits each:
	 *      Let n = max(len(a), len(b)).
	 *      If n <= 64, then xa = na and xb = nb.
	 *      Otherwise:
	 *         xa = (a mod 2^31) + 2^31*floor(a / 2^(n - 33))
	 *         xb = (b mod 2^31) + 2^31*floor(b / 2^(n - 33))
	 *    I.e. we keep the same low 31 bits, but we remove all the
	 *    "middle" bits to just keep the higher bits. At least one
	 *    of the two values xa and xb will have maximum length (64 bits)
	 *    if either of a or b exceeds 64 bits.
	 *  - We record which subtract and swap we did into update
	 *    factors, which we apply en masse to a, b, u and v after
	 *    the 31 iterations.
	 *
	 * Since we kept the correct low 31 bits, all the "subtract or
	 * not" decisions are correct, but the "swap or not swap" might
	 * be wrong, since these use the difference between the two
	 * values, and we only have approximations thereof. A consequence
	 * is that after the update, a and/or b may be negative. If a < 0,
	 * we negate it, and also negate u (which we do by negating the
	 * update factors f0 and g0 before applying them to compute the
	 * new value of u); similary for b and v.
	 *
	 * It can be shown that the 31 iterations, along with the
	 * conditional negation, ensure that len(a)+len(b) is still
	 * reduced by at least 31 bits (compared with the classic binary
	 * GCD, the use of the approximations may make the code "miss one
	 * bit", but the conditional subtraction regains it). Thus,
	 * 509 iterations are sufficient in total. As explained later on,
	 * we can skip the final iteration as well.
	 */

	gf_normalize(&a, y);
	b = GF_P;
	u = GF_ONE;
	v = GF_ZERO;

	/*
	 * Generic loop first does 15*31 = 465 iterations.
	 */
	for (i = 0; i < 15; i ++) {
		unsigned long long m1, m2, m3, tnz1, tnz2, tnz3;
		unsigned long long tnzm, tnza, tnzb, snza, snzb;
		unsigned long long s, sm;
		gf na, nb, nu, nv;

		/*
		 * Get approximations of a and b over 64 bits:
		 *  - If len(a) <= 64 and len(b) <= 64, then we just
		 *    use the value (low limb).
		 *  - Otherwise, with n = max(len(a), len(b)), we use:
		 *       (a mod 2^31) + 2^31*(floor(a / 2^(n-33)))
		 *       (b mod 2^31) + 2^31*(floor(b / 2^(n-33)))
		 * I.e. we remove the "middle bits".
		 */
		m3 = a.v3 | b.v3;
		m2 = a.v2 | b.v2;
		m1 = a.v1 | b.v1;
		tnz3 = -((m3 | -m3) >> 63);
		tnz2 = -((m2 | -m2) >> 63) & ~tnz3;
		tnz1 = -((m1 | -m1) >> 63) & ~tnz3 & ~tnz2;
		tnzm = (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1);
		tnza = (a.v3 & tnz3) | (a.v2 & tnz2) | (a.v1 & tnz1);
		tnzb = (b.v3 & tnz3) | (b.v2 & tnz2) | (b.v1 & tnz1);
		snza = (a.v2 & tnz3) | (a.v1 & tnz2) | (a.v0 & tnz1);
		snzb = (b.v2 & tnz3) | (b.v1 & tnz2) | (b.v0 & tnz1);

		/*
		 * If both len(a) <= 64 and len(b) <= 64, then:
		 *    tnzm = 0
		 *    tnza = 0, snza = 0, tnzb = 0, snzb = 0
		 * Otherwise:
		 *    tnzm != 0, length yields value of n
		 *    tnza contains the top limb of a, snza the second limb
		 *    tnzb contains the top limb of b, snzb the second limb
		 *
		 * We count the number of leading zero bits in tnzm:
		 *  - If s <= 31, then the top 31 bits can be extracted
		 *    from tnza and tnzb alone.
		 *  - If 32 <= s <= 63, then we need some bits from snza
		 *    as well.
		 *
		 * We rely on the fact shifts don't reveal the shift count
		 * through side channels. This would not have been true on
		 * the Pentium IV, but it is true on all known x86 CPU that
		 * have 64-bit support and implement the LZCNT opcode.
		 */
		s = _lzcnt_u64(tnzm);
		sm = -(unsigned long long)(s >> 5);
		tnza ^= sm & (tnza ^ ((tnza << 32) | (snza >> 32)));
		tnzb ^= sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32)));
		s &= 31;
		tnza <<= s;
		tnzb <<= s;

		/*
		 * At this point:
		 *  - If len(a) <= 64 and len(b) <= 64, then:
		 *       tnza = 0
		 *       tnzb = 0
		 *       tnz1 = tnz2 = tnz3 = 0
		 *  - Otherwise, we need to use the top 33 bits of tnza
		 *    and tnzb in combination with the low 31 bits of
		 *    a.v0 and b.v0, respectively.
		 */
		tnza |= a.v0 & ~(tnz1 | tnz2 | tnz3);
		tnzb |= b.v0 & ~(tnz1 | tnz2 | tnz3);
		xa = (a.v0 & 0x7FFFFFFF) | (tnza & 0xFFFFFFFF80000000);
		xb = (b.v0 & 0x7FFFFFFF) | (tnzb & 0xFFFFFFFF80000000);

		/*
		 * We can now run the binary GCD on xa and xb for 31
		 * rounds.
		 */
		fg0 = (uint64_t)1;
		fg1 = (uint64_t)1 << 32;
		for (j = 0; j < 31; j ++) {
			uint64_t a_odd, swap;
			unsigned long long t;

			a_odd = -(uint64_t)(xa & 1);
			swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
			t = swap & (xa ^ xb);
			xa ^= t;
			xb ^= t;
			t = swap & (fg0 ^ fg1);
			fg0 ^= t;
			fg1 ^= t;
			xa -= a_odd & xb;
			fg0 -= a_odd & fg1;
			xa >>= 1;
			fg1 <<= 1;
		}
		fg0 += 0x7FFFFFFF7FFFFFFF;
		fg1 += 0x7FFFFFFF7FFFFFFF;
		f0 = (fg0 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
		g0 = (fg0 >> 32) - (uint64_t)0x7FFFFFFF;
		f1 = (fg1 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
		g1 = (fg1 >> 32) - (uint64_t)0x7FFFFFFF;

		/*
		 * We now need to propagate updates to a, b, u and v.
		 */
		nega = s256_lin_div31_abs(&na, &a, &b, f0, g0);
		negb = s256_lin_div31_abs(&nb, &a, &b, f1, g1);
		f0 = (f0 ^ -nega) + nega;
		g0 = (g0 ^ -nega) + nega;
		f1 = (f1 ^ -negb) + negb;
		g1 = (g1 ^ -negb) + negb;
		gf_lin(&nu, &u, &v, f0, g0);
		gf_lin(&nv, &u, &v, f1, g1);
		a = na;
		b = nb;
		u = nu;
		v = nv;
	}

	/*
	 * At that point, if y is invertible, then the final GCD is 1,
	 * and len(a) + len(b) <= 45, so it is known that the values
	 * fully fit in a single register each. We can do the remaining
	 * 44 iterations in one go (they are exact, no approximation
	 * here). In fact, we can content ourselves with 43 iterations,
	 * because when arriving at the last iteration, we know that a = 0
	 * or 1 and b = 1 (the final state of the algorithm is that a = 0
	 * and b is the GCD of y and q), so there would be no swap. Since
	 * we only care about the update factors f1 and g1, we can simply
	 * avoid the final iteration.
	 *
	 * The update values f1 and g1, for v, will be up to 2^43 (in
	 * absolute value) but this is supported by gf_lin().
	 *
	 * If y is zero, then none of the iterations changes anything
	 * to b and v, and therefore v = 0 at this point, which is the
	 * value we want to obtain in that case.
	 */
	xa = a.v0;
	xb = b.v0;

	/*
	 * We first do the first 31 iterations with two update factors
	 * paired in the same variable.
	 */
	fg0 = (uint64_t)1;
	fg1 = (uint64_t)1 << 32;
	for (j = 0; j < 31; j ++) {
		uint64_t a_odd, swap;
		unsigned long long t;

		a_odd = -(uint64_t)(xa & 1);
		swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
		t = swap & (xa ^ xb);
		xa ^= t;
		xb ^= t;
		t = swap & (fg0 ^ fg1);
		fg0 ^= t;
		fg1 ^= t;
		xa -= a_odd & xb;
		fg0 -= a_odd & fg1;
		xa >>= 1;
		fg1 <<= 1;
	}

	/*
	 * Pairing the update factors in the same variable works only for
	 * up to 31 iterations; we split them into separate variables
	 * to finish the work with 12 extra iterations.
	 */
	fg0 += 0x7FFFFFFF7FFFFFFF;
	f0 = (fg0 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
	g0 = (fg0 >> 32) - (uint64_t)0x7FFFFFFF;
	fg1 += 0x7FFFFFFF7FFFFFFF;
	f1 = (fg1 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
	g1 = (fg1 >> 32) - (uint64_t)0x7FFFFFFF;
	for (j = 0; j < 12; j ++) {
		uint64_t a_odd, swap;
		unsigned long long t;

		a_odd = -(uint64_t)(xa & 1);
		swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
		t = swap & (xa ^ xb);
		xa ^= t;
		xb ^= t;
		t = swap & (f0 ^ f1);
		f0 ^= t;
		f1 ^= t;
		t = swap & (g0 ^ g1);
		g0 ^= t;
		g1 ^= t;
		xa -= a_odd & xb;
		f0 -= a_odd & f1;
		g0 -= a_odd & g1;
		xa >>= 1;
		f1 <<= 1;
		g1 <<= 1;
	}

	gf_lin(&v, &u, &v, f1, g1);

	/*
	 * We did 31*15+43 = 508 iterations, and each injected a factor 2,
	 * thus we must divide by 2^508 (mod q).
	 *
	 * Result is correct if source operand was invertible, i.e.
	 * distinct from zero (since all non-zero values are invertible
	 * modulo a prime integer); the inverse is then also non-zero.
	 * If the source was zero, then the result is zero as well. We
	 * can thus test d instead of a.
	 */
	gf_mul_inline(d, &v, &GF_INVT508);
	return gf_iszero(d) ^ 1;
}

/*
 * Legendre symbol computation. Return value:
 *   1   if y != 0 and is a quadratic residue
 *  -1   if y != 0 and is not a quadratic residue
 *   0   if y == 0
 */
UNUSED
static int64_t
gf_legendre(const gf *y)
{
	/*
	 * Algorithm is basically the same as the binary GCD, with the
	 * following differences:
	 *
	 *  - We do not keep track or update the 'u' and 'v' values.
	 *  - In the inner loop, we need to access the three bottom bits
	 *    of 'b'; the last two iterations do not necessarily have
	 *    these values available, so we must compute the updated a
	 *    and b (low words only) at some point to get the correct bits.
	 */
	gf a, b;
	unsigned long long f0, f1, g0, g1, xa, xb, fg0, fg1, ls, a0, b0;
	unsigned long long nega;
	int i, j;

	gf_normalize(&a, y);
	b = GF_P;
	ls = 0;

	/*
	 * Generic loop first does 15*31 = 465 iterations.
	 */
	for (i = 0; i < 15; i ++) {
		unsigned long long m1, m2, m3, tnz1, tnz2, tnz3;
		unsigned long long tnzm, tnza, tnzb, snza, snzb;
		unsigned long long s, sm;
		gf na, nb;

		/*
		 * Get approximations of a and b over 64 bits:
		 *  - If len(a) <= 64 and len(b) <= 64, then we just
		 *    use the value (low limb).
		 *  - Otherwise, with n = max(len(a), len(b)), we use:
		 *       (a mod 2^31) + 2^31*(floor(a / 2^(n-33)))
		 *       (b mod 2^31) + 2^31*(floor(b / 2^(n-33)))
		 * I.e. we remove the "middle bits".
		 */
		m3 = a.v3 | b.v3;
		m2 = a.v2 | b.v2;
		m1 = a.v1 | b.v1;
		tnz3 = -((m3 | -m3) >> 63);
		tnz2 = -((m2 | -m2) >> 63) & ~tnz3;
		tnz1 = -((m1 | -m1) >> 63) & ~tnz3 & ~tnz2;
		tnzm = (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1);
		tnza = (a.v3 & tnz3) | (a.v2 & tnz2) | (a.v1 & tnz1);
		tnzb = (b.v3 & tnz3) | (b.v2 & tnz2) | (b.v1 & tnz1);
		snza = (a.v2 & tnz3) | (a.v1 & tnz2) | (a.v0 & tnz1);
		snzb = (b.v2 & tnz3) | (b.v1 & tnz2) | (b.v0 & tnz1);

		/*
		 * If both len(a) <= 64 and len(b) <= 64, then:
		 *    tnzm = 0
		 *    tnza = 0, snza = 0, tnzb = 0, snzb = 0
		 *    tnzm = 0
		 * Otherwise:
		 *    tnzm != 0, length yields value of n
		 *    tnza contains the top limb of a, snza the second limb
		 *    tnzb contains the top limb of b, snzb the second limb
		 *
		 * We count the number of leading zero bits in tnzm:
		 *  - If s <= 31, then the top 31 bits can be extracted
		 *    from tnza and tnzb alone.
		 *  - If 32 <= s <= 63, then we need some bits from snza
		 *    as well.
		 *
		 * We rely on the fact shifts don't reveal the shift count
		 * through side channels. This would not have been true on
		 * the Pentium IV, but it is true on all known x86 CPU that
		 * have 64-bit support and implement the LZCNT opcode.
		 */
		s = _lzcnt_u64(tnzm);
		sm = -((unsigned long long)(31 - s) >> 63);
		tnza ^= sm & (tnza ^ ((tnza << 32) | (snza >> 32)));
		tnzb ^= sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32)));
		s -= 32 & sm;
		tnza <<= s;
		tnzb <<= s;

		/*
		 * At this point:
		 *  - If len(a) <= 64 and len(b) <= 64, then:
		 *       tnza = 0
		 *       tnzb = 0
		 *       tnz1 = tnz2 = tnz3 = 0
		 *  - Otherwise, we need to use the top 33 bits of tnza
		 *    and tnzb in combination with the low 31 bits of
		 *    a.v0 and b.v0, respectively.
		 */
		tnza |= a.v0 & ~(tnz1 | tnz2 | tnz3);
		tnzb |= b.v0 & ~(tnz1 | tnz2 | tnz3);
		xa = (a.v0 & 0x7FFFFFFF) | (tnza & 0xFFFFFFFF80000000);
		xb = (b.v0 & 0x7FFFFFFF) | (tnzb & 0xFFFFFFFF80000000);

		/*
		 * We can now run the binary GCD on xa and xb for 29
		 * rounds.
		 */
		fg0 = (uint64_t)1;
		fg1 = (uint64_t)1 << 32;
		for (j = 0; j < 29; j ++) {
			uint64_t a_odd, swap;
			unsigned long long t;

			a_odd = -(uint64_t)(xa & 1);
			swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
			ls += swap & ((xa & xb) >> 1);
			t = swap & (xa ^ xb);
			xa ^= t;
			xb ^= t;
			t = swap & (fg0 ^ fg1);
			fg0 ^= t;
			fg1 ^= t;
			xa -= a_odd & xb;
			fg0 -= a_odd & fg1;
			xa >>= 1;
			fg1 <<= 1;
			ls += (xb + 2) >> 2;
		}

		/*
		 * For the last two iterations, we need to recompute the
		 * updated a and b (low words only) to get their low bits.
		 */
		f0 = ((fg0 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF)
			- (uint64_t)0x7FFFFFFF;
		g0 = ((fg0 + 0x7FFFFFFF7FFFFFFF) >> 32)
			- (uint64_t)0x7FFFFFFF;
		f1 = ((fg1 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF)
			- (uint64_t)0x7FFFFFFF;
		g1 = ((fg1 + 0x7FFFFFFF7FFFFFFF) >> 32)
			- (uint64_t)0x7FFFFFFF;
		a0 = (a.v0 * f0 + b.v0 * g0) >> 29;
		b0 = (a.v0 * f1 + b.v0 * g1) >> 29;
		for (j = 0; j < 2; j ++) {
			uint64_t a_odd, swap;
			unsigned long long t;

			a_odd = -(uint64_t)(xa & 1);
			swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
			ls += swap & ((a0 & b0) >> 1);
			t = swap & (xa ^ xb);
			xa ^= t;
			xb ^= t;
			t = swap & (fg0 ^ fg1);
			fg0 ^= t;
			fg1 ^= t;
			t = swap & (a0 ^ b0);
			a0 ^= t;
			b0 ^= t;
			xa -= a_odd & xb;
			fg0 -= a_odd & fg1;
			a0 -= a_odd & b0;
			xa >>= 1;
			fg1 <<= 1;
			a0 >>= 1;
			ls += (b0 + 2) >> 2;
		}

		fg0 += 0x7FFFFFFF7FFFFFFF;
		fg1 += 0x7FFFFFFF7FFFFFFF;
		f0 = (fg0 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
		g0 = (fg0 >> 32) - (uint64_t)0x7FFFFFFF;
		f1 = (fg1 & 0xFFFFFFFF) - (uint64_t)0x7FFFFFFF;
		g1 = (fg1 >> 32) - (uint64_t)0x7FFFFFFF;

		/*
		 * We now need to propagate updates to a and b.
		 */
		nega = s256_lin_div31_abs(&na, &a, &b, f0, g0);
		s256_lin_div31_abs(&nb, &a, &b, f1, g1);
		ls += nega & (nb.v0 >> 1);
		a = na;
		b = nb;
	}

	/*
	 * At that point, if y is invertible, then the final GCD is 1,
	 * and len(a) + len(b) <= 45, so it is known that the values
	 * fully fit in a single register each. We can do the remaining
	 * 44 iterations in one go. We do not need to compute update
	 * factors. Since values are exact, we always have the proper
	 * bit values. We can even stop at 43 iterations, because at the
	 * last iteration, b = 1 and a = 0 or 1; in either case, no
	 * modification of the Legendre symbol will happen.
	 */
	xa = a.v0;
	xb = b.v0;
	for (j = 0; j < 43; j ++) {
		uint64_t a_odd, swap;
		unsigned long long t;

		a_odd = -(uint64_t)(xa & 1);
		swap = a_odd & -(uint64_t)_subborrow_u64(0, xa, xb, &t);
		ls += swap & ((xa & xb) >> 1);
		t = swap & (xa ^ xb);
		xa ^= t;
		xb ^= t;
		xa -= a_odd & xb;
		xa >>= 1;
		ls += (xb + 2) >> 2;
	}

	/*
	 * At this point, if y != 0, then the low bit of ls contains the
	 * QR status (0 = square, 1 = non-square). If y == 0, then we
	 * replace the output value with zero, as per the API.
	 */
	return (int64_t)((uint64_t)(1 - ((int64_t)(ls & 1) << 1))
		& (gf_iszero(y) - 1));
}

static unsigned long long
dec64le(const uint8_t *buf)
{
	return (unsigned long long)buf[0]
		| ((unsigned long long)buf[1] << 8)
		| ((unsigned long long)buf[2] << 16)
		| ((unsigned long long)buf[3] << 24)
		| ((unsigned long long)buf[4] << 32)
		| ((unsigned long long)buf[5] << 40)
		| ((unsigned long long)buf[6] << 48)
		| ((unsigned long long)buf[7] << 56);
}

static void
enc64le(uint8_t *buf, unsigned long long x)
{
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
	buf[4] = (uint8_t)(x >> 32);
	buf[5] = (uint8_t)(x >> 40);
	buf[6] = (uint8_t)(x >> 48);
	buf[7] = (uint8_t)(x >> 56);
}

/*
 * Decode a value from 32 bytes. Decoding always succeeds. Returned value
 * is 1 if the value was normalized (in the 0..p-1 range), 0 otherwise.
 */
UNUSED
static uint64_t
gf_decode(gf *d, const void *src)
{
	const uint8_t *buf;
	unsigned long long t;
	unsigned char cc;

	buf = src;
	d->v0 = dec64le(buf);
	d->v1 = dec64le(buf + 8);
	d->v2 = dec64le(buf + 16);
	d->v3 = dec64le(buf + 24);
	cc = _subborrow_u64(0, d->v0, GF_P.v0, &t);
	cc = _subborrow_u64(cc, d->v1, GF_P.v1, &t);
	cc = _subborrow_u64(cc, d->v2, GF_P.v2, &t);
	cc = _subborrow_u64(cc, d->v3, GF_P.v3, &t);
	return (uint64_t)cc;
}

/*
 * Decode a value from bytes, with modular reduction.
 */
UNUSED
static void
gf_decode_reduce(gf *d, const void *src, size_t len)
{
	const uint8_t *buf;
	unsigned long long d0, d1, d2, d3;

	buf = src;
	if (len == 0) {
		d->v0 = 0;
		d->v1 = 0;
		d->v2 = 0;
		d->v3 = 0;
		return;
	} else if ((len & 31) != 0) {
		uint8_t tmp[32];
		size_t n;

		n = len & 31;
		len -= n;
		memcpy(tmp, buf + len, n);
		memset(tmp + n, 0, (sizeof tmp) - n);
		d0 = dec64le(tmp);
		d1 = dec64le(tmp + 8);
		d2 = dec64le(tmp + 16);
		d3 = dec64le(tmp + 24);
	} else {
		len -= 32;
		d0 = dec64le(buf + len);
		d1 = dec64le(buf + len + 8);
		d2 = dec64le(buf + len + 16);
		d3 = dec64le(buf + len + 24);
	}

	while (len > 0) {
		unsigned long long e0, e1, e2, e3;
		unsigned long long h0, h1, h2, h3;
		unsigned char cc;

		len -= 32;
		e0 = dec64le(buf + len);
		e1 = dec64le(buf + len + 8);
		e2 = dec64le(buf + len + 16);
		e3 = dec64le(buf + len + 24);
		UMUL64(d0, h0, d0, 2 * MQ);
		cc = _addcarry_u64(0, d0, e0, &d0);
		UMUL64(d1, h1, d1, 2 * MQ);
		cc = _addcarry_u64(cc, d1, e1, &d1);
		UMUL64(d2, h2, d2, 2 * MQ);
		cc = _addcarry_u64(cc, d2, e2, &d2);
		UMUL64(d3, h3, d3, 2 * MQ);
		cc = _addcarry_u64(cc, d3, e3, &d3);
		(void)_addcarry_u64(cc, 0, h3, &h3);

		h3 = (h3 << 1) | (d3 >> 63);
		d3 &= 0x7FFFFFFFFFFFFFFF;
		cc = _addcarry_u64(0, d0, h3 * MQ, &d0);
		cc = _addcarry_u64(cc, d1, h0, &d1);
		cc = _addcarry_u64(cc, d2, h1, &d2);
		(void)_addcarry_u64(cc, d3, h2, &d3);
	}

	d->v0 = d0;
	d->v1 = d1;
	d->v2 = d2;
	d->v3 = d3;
}

/*
 * Encode a value into 32 bytes. Encoding is always normalized; a
 * consequence is that the highest bit of the last byte is always 0.
 */
UNUSED
static void
gf_encode(void *dst, const gf *a)
{
	gf t;
	uint8_t *buf;

	gf_normalize(&t, a);
	buf = dst;
	enc64le(buf, t.v0);
	enc64le(buf + 8, t.v1);
	enc64le(buf + 16, t.v2);
	enc64le(buf + 24, t.v3);
}
