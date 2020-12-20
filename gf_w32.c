/*
 * This file implements elementary operations in finite field GF(2^255-MQ),
 * with only portable C code (though intrinsics are used for addition
 * and subtraction with carries).
 *
 * It is not meant to be compiled by itself, but included in an outer
 * file. The outer file must provide the following:
 *
 *    MQ            macro that evaluates to a literal constant
 *    GF_INVT508    value 1/2^508 mod 2^255-MQ
 *
 * Moreover, the common header "do255.h" must have already been included.
 * Value GF_INVT508 is a 'struct do255_int256_w32' and the value is
 * represented in four 32-bit limbs (v[0] is lowest limb, v[7] is highest).
 */

#include <immintrin.h>

typedef struct do255_int256_w32 gf;

/*
 * Constants below need not be modified if MQ is changed.
 */
static const gf GF_ZERO = { { 0, 0, 0, 0, 0, 0, 0, 0 } };
static const gf GF_ONE = { { 1, 0, 0, 0, 0, 0, 0, 0 } };
static const gf GF_TWO = { { 2, 0, 0, 0, 0, 0, 0, 0 } };
static const gf GF_SEVEN = { { 7, 0, 0, 0, 0, 0, 0, 0 } };
static const gf GF_P = { {
	(uint32_t)-MQ, (uint32_t)-1, (uint32_t)-1, (uint32_t)-1,
	(uint32_t)-1, (uint32_t)-1, (uint32_t)-1, (uint32_t)-1 >> 1
} };

/*
 * Helpful macro to convert expressions to strings, for inline assembly.
 */
#define ST(x)    ST_(x)
#define ST_(x)   #x

/*
 * A field element is represented as eight limbs, in base 2^32. Operands
 * and result may be up to 2^256-1.
 * A _normalized_ value is in 0..p-1. gf_normalize() ensures normalization;
 * it is called before encoding, and also when starting inversion.
 */

/* d <- a + b */
UNUSED
static void
gf_add(gf *d, const gf *a, const gf *b)
{
	int i;
	uint32_t dw[8];
	unsigned char cc;

	/*
	 * Compute sum over 256 bits + carry.
	 */
	cc = 0;
	for (i = 0; i < 8; i ++) {
		cc = _addcarry_u32(cc, a->v[i], b->v[i], &dw[i]);
	}

	/*
	 * If there is a carry, subtract 2*p, i.e. add 2*MQ and
	 * (implicitly) subtract 2^256.
	 */
	cc = _addcarry_u32(0, dw[0], -(uint32_t)cc & (2 * MQ), &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}

	/*
	 * Usually, there won't be an extra carry here, and the
	 * implicit subtraction of 2^256 is enough. However, if the
	 * the original sum was at least 2^257 - 2*MQ, then this
	 * will show up as an extra carry here. In that case, the
	 * low limb (dw[0]) is necessarily less than 2*MQ, and thus adding
	 * back 2*MQ to it will not trigger any carry.
	 */
	(void)_addcarry_u32(0, dw[0], -(uint32_t)cc & (2 * MQ), &d->v[0]);
}

/* d <- a - b */
UNUSED
static void
gf_sub(gf *d, const gf *a, const gf *b)
{
	int i;
	uint32_t dw[8];
	unsigned char cc;

	/*
	 * Compute subtraction over 256 bits.
	 */
	cc = 0;
	for (i = 0; i < 8; i ++) {
		cc = _subborrow_u32(cc, a->v[i], b->v[i], &dw[i]);
	}

	/*
	 * If the result is negative, add back 2*p = 2^256 - 2*MQ (which
	 * we do by subtracting 2*MQ).
	 */
	cc = _subborrow_u32(0, dw[0], -(uint32_t)cc & (2 * MQ), &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u32(cc, dw[i], 0, &d->v[i]);
	}

	/*
	 * If there is still a carry, then we must add 2*p again. In that
	 * case, subtracting 2*MQ from the low limb won't trigger a carry.
	 */
	d->v[0] = dw[0] - (-(uint32_t)cc & (2 * MQ));
}

/* d <- -a */
UNUSED
static void
gf_neg(gf *d, const gf *a)
{
	int i;
	uint32_t dw[8], e;
	unsigned char cc;

	/*
	 * Compute 2*p - a over 256 bits.
	 */
	cc = _subborrow_u32(0, (uint32_t)-(2 * MQ), a->v[0], &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u32(cc, (uint32_t)-1, a->v[i], &dw[i]);
	}

	/*
	 * If the result is negative, add back p = 2^255 - MQ.
	 */
	e = -(uint32_t)cc;
	cc = _addcarry_u32(0, dw[0], e & -MQ, &d->v[0]);
	for (i = 1; i < 7; i ++) {
		cc = _addcarry_u32(cc, dw[i], e, &d->v[i]);
	}
	(void)_addcarry_u32(cc, dw[7], e >> 1, &d->v[7]);
}

/* If ctl = 1, then d <- -a
   If ctl = 0, then d <- a
   ctl MUST be 0 or 1 */
UNUSED
static void
gf_condneg(gf *d, const gf *a, uint32_t ctl)
{
	gf t;
	int i;

	gf_neg(&t, a);
	for (i = 0; i < 8; i ++) {
		d->v[i] = a->v[i] ^ (-ctl & (a->v[i] ^ t.v[i]));
	}
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
	int i;
	uint32_t dw[8], tt;
	unsigned char cc;

	for (i = 0; i < 7; i ++) {
		dw[i] = (a->v[i] >> 1) | (a->v[i + 1] << 31);
	}
	dw[7] = a->v[7] >> 1;
	tt = -(a->v[0] & 1);

	cc = _addcarry_u32(0, dw[0], tt & (uint32_t)-((MQ - 1) >> 1), &d->v[0]);
	for (i = 1; i < 7; i ++) {
		cc = _addcarry_u32(cc, dw[i], tt, &d->v[i]);
	}
	(void)_addcarry_u32(cc, dw[7], tt >> 2, &d->v[7]);
}

/* If ctl = 1, then d <- a
   If ctl = 0, then d <- b
   ctl MUST be 0 or 1 */
UNUSED
static inline void
gf_sel2(gf *d, const gf *a, const gf *b, uint32_t ctl)
{
	int i;

	for (i = 0; i < 8; i ++) {
		d->v[i] = b->v[i] ^ (-ctl & (a->v[i] ^ b->v[i]));
	}
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
gf_sel3(gf *d, const gf *a, const gf *b, const gf *c, uint32_t fa, uint32_t fb)
{
	int i;
	uint32_t ma, mb, mc;

	ma = -fa;
	mb = -fb & ~ma;
	mc = ~(ma | mb);
	for (i = 0; i < 8; i ++) {
		d->v[i] = (a->v[i] & ma) | (b->v[i] & mb) | (c->v[i] & mc);
	}
}

/* d <- 2*a */
UNUSED
static void
gf_mul2(gf *d, const gf *a)
{
	int i;
	uint32_t dw[8], tt;
	unsigned char cc;

	dw[0] = a->v[0] << 1;
	for (i = 1; i < 8; i ++) {
		dw[i] = (a->v[i] << 1) | (a->v[i - 1] >> 31);
	}

	tt = (*(int32_t *)&a->v[7] >> 31) & (uint32_t)(2 * MQ);
	cc = _addcarry_u32(0, dw[0], tt, &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}
	d->v[0] = dw[0] + (-(uint32_t)cc & (2 * MQ));
}

/* d <- 4*a */
UNUSED
static void
gf_mul4(gf *d, const gf *a)
{
	int i;
	uint32_t dw[8], tt;
	unsigned char cc;

	dw[0] = a->v[0] << 2;
	for (i = 1; i < 8; i ++) {
		dw[i] = (a->v[i] << 2) | (a->v[i - 1] >> 30);
	}

	tt = (a->v[7] >> 30) * (uint32_t)(2 * MQ);
	cc = _addcarry_u32(0, dw[0], tt, &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}
	d->v[0] = dw[0] + (-(uint32_t)cc & (2 * MQ));
}

/* d <- 8*a */
UNUSED
static void
gf_mul8(gf *d, const gf *a)
{
	int i;
	uint32_t dw[8], tt;
	unsigned char cc;

	dw[0] = a->v[0] << 3;
	for (i = 1; i < 8; i ++) {
		dw[i] = (a->v[i] << 3) | (a->v[i - 1] >> 29);
	}

	tt = (a->v[7] >> 29) * (uint32_t)(2 * MQ);
	cc = _addcarry_u32(0, dw[0], tt, &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}
	d->v[0] = dw[0] + (-(uint32_t)cc & (2 * MQ));
}

/* d <- 32*a */
UNUSED
static void
gf_mul32(gf *d, const gf *a)
{
	int i;
	uint32_t dw[8], tt;
	unsigned char cc;

	dw[0] = a->v[0] << 5;
	for (i = 1; i < 8; i ++) {
		dw[i] = (a->v[i] << 5) | (a->v[i - 1] >> 27);
	}

	tt = (a->v[7] >> 27) * (uint32_t)(2 * MQ);
	cc = _addcarry_u32(0, dw[0], tt, &dw[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}
	d->v[0] = dw[0] + (-(uint32_t)cc & (2 * MQ));
}

/* d <- a*b, with b < 2^32 */
UNUSED
static void
gf_mul_small(gf *d, const gf *a, uint32_t b)
{
	int i;
	uint64_t z;
	uint32_t g;
	unsigned char cc;

	g = 0;
	for (i = 0; i < 8; i ++) {
		z = (uint64_t)a->v[i] * (uint64_t)b + (uint64_t)g;
		d->v[i] = (uint32_t)z;
		g = (uint32_t)(z >> 32);
	}
	g = (g << 1) | (d->v[7] >> 31);
	d->v[7] &= 0x7FFFFFFF;
	z = (uint64_t)g * (uint64_t)MQ;
	cc = _addcarry_u32(0, d->v[0], (uint32_t)z, &d->v[0]);
	cc = _addcarry_u32(cc, d->v[1], (uint32_t)(z >> 32), &d->v[1]);
	for (i = 2; i < 8; i ++) {
		cc = _addcarry_u32(cc, d->v[i], 0, &d->v[i]);
	}
}

/* d <- a*b  (always inlined) */
FORCE_INLINE
static inline void
gf_mul_inline(gf *d, const gf *a, const gf *b)
{
	int i, j;
	uint32_t e[16], f, g;
	uint64_t z;
	unsigned char cc;

	f = a->v[0];
	g = 0;
	for (j = 0; j < 8; j ++) {
		z = (uint64_t)f * (uint64_t)b->v[j] + (uint64_t)g;
		e[j] = (uint32_t)z;
		g = (uint32_t)(z >> 32);
	}
	e[8] = g;

	for (i = 1; i < 8; i ++) {
		f = a->v[i];
		g = 0;
		for (j = 0; j < 8; j ++) {
			z = (uint64_t)f * (uint64_t)b->v[j] + (uint64_t)g;
			z += (uint64_t)e[i + j];
			e[i + j] = (uint32_t)z;
			g = (uint32_t)(z >> 32);
		}
		e[i + 8] = g;
	}

	g = 0;
	cc = 0;
	for (i = 0; i < 8; i ++) {
		z = (uint64_t)e[i + 8] * (uint64_t)(2 * MQ) + (uint64_t)g;
		g = (uint32_t)(z >> 32);
		cc = _addcarry_u32(cc, e[i], (uint32_t)z, &e[i]);
	}
	g += cc;

	z = (uint64_t)g * (uint64_t)(2 * MQ);
	cc = _addcarry_u32(0, e[0], (uint32_t)z, &e[0]);
	cc = _addcarry_u32(cc, e[1], (uint32_t)(z >> 32), &e[1]);
	for (i = 2; i < 8; i ++) {
		cc = _addcarry_u32(cc, e[i], 0, &d->v[i]);
	}

	cc = _addcarry_u32(0, e[0], -(uint32_t)cc & (2 * MQ), &d->v[0]);
	(void)_addcarry_u32(cc, e[1], 0, &d->v[1]);
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
	// TODO: optimize this
	gf_mul_inline(d, a, a);
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
	int i;
	uint32_t dw[8], e;
	unsigned char cc;

	/*
	 * If top bit is set, propagate it.
	 */
	e = (uint32_t)(*(int32_t *)&a->v[7] >> 31);
	cc = _addcarry_u32(0, a->v[0], e & MQ, &dw[0]);
	for (i = 1; i < 7; i ++) {
		cc = _addcarry_u32(cc, a->v[i], 0, &dw[i]);
	}
	(void)_addcarry_u32(cc, a->v[7], e << 31, &dw[7]);

	/*
	 * Value is now at most 2^255+MQ-1. Subtract p, and add it back
	 * if the result is negative.
	 */
	cc = _subborrow_u32(0, dw[0], (uint32_t)-MQ, &dw[0]);
	for (i = 1; i < 7; i ++) {
		cc = _subborrow_u32(cc, dw[i], (uint32_t)-1, &dw[i]);
	}
	cc = _subborrow_u32(cc, dw[7], (uint32_t)-1 >> 1, &dw[7]);

	e = -(uint32_t)cc;
	cc = _addcarry_u32(0, dw[0], e & -MQ, &d->v[0]);
	for (i = 1; i < 7; i ++) {
		cc = _addcarry_u32(cc, dw[i], e, &d->v[i]);
	}
	(void)_addcarry_u32(cc, dw[7], e >> 1, &d->v[7]);
}

/*
 * Compare 'a' with zero; returned value is 1 if the value is zero (as
 * a field element), 0 otherwise.
 */
UNUSED
static uint32_t
gf_iszero(const gf *a)
{
	int i;
	uint32_t t0, t1, t2;

	/*
	 * Since values are over 256 bits, there are three possible
	 * representations for 0: 0, p and 2*p.
	 */
	t0 = a->v[0];
	t1 = a->v[0] + MQ;
	t2 = a->v[0] + 2 * MQ;
	for (i = 1; i < 7; i ++) {
		t0 |= a->v[i];
		t1 |= ~a->v[i];
		t2 |= ~a->v[i];
	}
	t0 |= a->v[7];
	t1 |= a->v[7] ^ 0x7FFFFFFF;
	t2 |= ~a->v[7];

	/*
	 * Value is zero if and only if one of t0, t1 or t2 is zero.
	 */
	return 1 - (((t0 | -t0) & (t1 | -t1) & (t2 | -t2)) >> 31);
}

/*
 * Compare two values as field elements; returned value is 1 on equality,
 * 0 otherwise.
 */
UNUSED
static uint32_t
gf_eq(const gf *a, const gf *b)
{
	gf t;

	gf_sub(&t, a, b);
	return gf_iszero(&t);
}

/*
 * Compute (a*f+v*g)/2^15. Parameters f and g are provided with an
 * unsigned type, but they are signed integers in the -2^15..+2^15 range.
 * Values a, b and d are not field elements, but 255-bit integers (value
 * is in 0..2^255-1, top bit is zero). The division by 2^15 is assumed to
 * be exact (low 15 bits of a*f+v*g are dropped). The result is assumed to
 * fit in 256 bits in signed two's complement notation (truncation is
 * applied on higher bits).
 *
 * If the result turns out to be negative, then it is negated. Returned
 * value is 1 if the result was negated, 0 otherwise.
 */
static uint32_t
s256_lin_div15_abs(gf *d, const gf *a, const gf *b, uint32_t f, uint32_t g)
{
	gf ta, tb;
	int i;
	uint32_t sf, sg, dw[8], t;
	uint64_t z;
	unsigned char cc;

	/*
	 * If f < 0, replace f with -f but keep the sign in sf.
	 * Similarly for g.
	 */
	sf = f >> 31;
	f = (f ^ -sf) + sf;
	sg = g >> 31;
	g = (g ^ -sg) + sg;

	/*
	 * Apply signs sf and sg to a and b, respectively.
	 */
	cc = _addcarry_u32(0, a->v[0] ^ -sf, sf, &ta.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, a->v[i] ^ -sf, 0, &ta.v[i]);
	}
	cc = _addcarry_u32(0, b->v[0] ^ -sg, sg, &tb.v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, b->v[i] ^ -sg, 0, &tb.v[i]);
	}

	/*
	 * Now that f and g are nonnegative, compute a*f+b*g into
	 * d0:d1:d2:d3:t. Since f and g are at most 2^15, we can
	 * add two 64-bit products with no overflow (they are actually
	 * 47 bits each at most).
	 */
	z = (uint64_t)ta.v[0] * (uint64_t)f + (uint64_t)tb.v[0] * (uint64_t)g;
	dw[0] = (uint32_t)z;
	t = (uint32_t)(z >> 32);
	for (i = 1; i < 8; i ++) {
		z = (uint64_t)ta.v[i] * (uint64_t)f
			+ (uint64_t)tb.v[i] * (uint64_t)g
			+ (uint64_t)t;
		dw[i] = (uint32_t)z;
		t = (uint32_t)(z >> 32);
	}

	/*
	 * Don't forget the signs: if a < 0, then the result is
	 * overestimated by 2^256*f; similarly, if b < 0, then the
	 * result is overestimated by 2^256*g. We thus must subtract
	 * 2^256*(sa*f+sb*g), where sa and sb are the signs of a and b,
	 * respectively (1 if negative, 0 otherwise).
	 */
	t -= -(uint32_t)(ta.v[7] >> 31) & f;
	t -= -(uint32_t)(tb.v[7] >> 31) & g;

	/*
	 * Apply the shift.
	 */
	for (i = 0; i < 7; i ++) {
		dw[i] = (dw[i] >> 15) | (dw[i + 1] << 17);
	}
	dw[7] = (dw[7] >> 15) | (t << 17);

	/*
	 * Perform conditional negation, if the result is negative.
	 */
	t >>= 31;
	cc = _addcarry_u32(0, dw[0] ^ -t, t, &d->v[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i] ^ -t, 0, &d->v[i]);
	}

	return t;
}

/*
 * Compute u*f+v*g (modulo p). Parameters f and g are provided with
 * an unsigned type, but they are signed integers in the -2^30..+2^30 range.
 */
static void
gf_lin(gf *d, const gf *u, const gf *v,
	uint32_t f, uint32_t g)
{
	gf tu, tv;
	int i;
	uint32_t sf, sg, dw[8], t;
	uint64_t z;
	unsigned char cc;

	/*
	 * If f < 0, replace f with -f but keep the sign in sf.
	 * Similarly for g.
	 */
	sf = f >> 31;
	f = (f ^ -sf) + sf;
	sg = g >> 31;
	g = (g ^ -sg) + sg;

	/*
	 * Apply signs sf and sg to u and v.
	 */
	gf_condneg(&tu, u, sf);
	gf_condneg(&tv, v, sg);

	/*
	 * Since f <= 2^30 and g <= 2^30, we can add 64-bit products: they
	 * fit on 61 bits each.
	 */
	z = (uint64_t)tu.v[0] * (uint64_t)f + (uint64_t)tv.v[0] * (uint64_t)g;
	dw[0] = (uint32_t)z;
	t = (uint32_t)(z >> 32);
	for (i = 1; i < 8; i ++) {
		z = (uint64_t)tu.v[i] * (uint64_t)f
			+ (uint64_t)tv.v[i] * (uint64_t)g
			+ (uint64_t)t;
		dw[i] = (uint32_t)z;
		t = (uint32_t)(z >> 32);
	}

	/*
	 * Upper word t can be up to 31 bits.
	 */
	z = (uint64_t)t * (uint64_t)(2 * MQ);
	cc = _addcarry_u32(0, dw[0], (uint32_t)z, &dw[0]);
	cc = _addcarry_u32(cc, dw[1], (uint32_t)(z >> 32), &dw[1]);
	for (i = 2; i < 8; i ++) {
		cc = _addcarry_u32(cc, dw[i], 0, &d->v[i]);
	}
	cc = _addcarry_u32(0, dw[0], -(uint32_t)cc & (2 * MQ), &d->v[0]);
	(void)_addcarry_u32(cc, dw[1], 0, &d->v[1]);
}

/*
 * Count leading zeros in provided value.
 */
static inline uint32_t
count_leading_zeros(uint32_t x)
{
	uint32_t r, c;

	r = 0;
	c = -(((x >> 16) - 1) >> 31);
	r += c & 16;
	x ^= c & (x ^ (x << 16));
	c = -(((x >> 24) - 1) >> 31);
	r += c & 8;
	x ^= c & (x ^ (x << 8));
	c = -(((x >> 28) - 1) >> 31);
	r += c & 4;
	x ^= c & (x ^ (x << 4));
	c = -(((x >> 30) - 1) >> 31);
	r += c & 2;
	x ^= c & (x ^ (x << 2));
	c = -(((x >> 31) - 1) >> 31);
	r += c & 1;
	x ^= c & (x ^ (x << 1));
	r += 1 - ((x | -x) >> 31);
	return r;
}

/*
 * Inversion in the field: d <- 1/y
 * If y = 0, then d is set to zero.
 * Returned value is 1 if the value was invertible, 0 otherwise.
 */
UNUSED
static uint32_t
gf_inv(gf *d, const gf *y)
{
	gf a, b, u, v;
	uint32_t f0, f1, g0, g1, xa, xb, fg0, fg1;
	uint32_t nega, negb;
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
	 *  - We do iterations by group of 15.
	 *  - In each group, we use _approximations_ of a and b that
	 *    fit on 32 bits each:
	 *      Let n = max(len(a), len(b)).
	 *      If n <= 32, then xa = na and xb = nb.
	 *      Otherwise:
	 *         xa = (a mod 2^15) + 2^15*floor(a / 2^(n - 17))
	 *         xb = (b mod 2^15) + 2^15*floor(b / 2^(n - 17))
	 *    I.e. we keep the same low 15 bits, but we remove all the
	 *    "middle" bits to just keep the higher bits. At least one
	 *    of the two values xa and xb will have maximum length (32 bits)
	 *    if either of a or b exceeds 32 bits.
	 *  - We record which subtract and swap we did into update
	 *    factors, which we apply en masse to a, b, u and v after
	 *    the 15 iterations.
	 *
	 * Since we kept the correct low 15 bits, all the "subtract or
	 * not" decisions are correct, but the "swap or not swap" might
	 * be wrong, since these use the difference between the two
	 * values, and we only have approximations thereof. A consequence
	 * is that after the update, a and/or b may be negative. If a < 0,
	 * we negate it, and also negate u (which we do by negating the
	 * update factors f0 and g0 before applying them to compute the
	 * new value of u); similary for b and v.
	 *
	 * It can be shown that the 15 iterations, along with the
	 * conditional negation, ensure that len(a)+len(b) is still
	 * reduced by at least 15 bits (compared with the classic binary
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
	 * Generic loop first does 32*15 = 480 iterations.
	 */
	for (i = 0; i < 32; i ++) {
		uint32_t c_hi, c_lo, a_hi, a_lo, b_hi, b_lo, s;
		gf na, nb, nu, nv;

		/*
		 * Get approximations of a and b over 32 bits:
		 *  - If len(a) <= 32 and len(b) <= 32, then we just
		 *    use the value (low limb).
		 *  - Otherwise, with n = max(len(a), len(b)), we use:
		 *       (a mod 2^15) + 2^15*(floor(a / 2^(n-17)))
		 *       (b mod 2^15) + 2^15*(floor(b / 2^(n-17)))
		 * I.e. we remove the "middle bits".
		 */
		c_hi = (uint32_t)-1;
		c_lo = (uint32_t)-1;
		a_hi = 0;
		a_lo = 0;
		b_hi = 0;
		b_lo = 0;
		for (j = 7; j >= 0; j --) {
			uint32_t aw, bw, mw;

			aw = a.v[j];
			bw = b.v[j];
			a_hi ^= (a_hi ^ aw) & c_hi;
			a_lo ^= (a_lo ^ aw) & c_lo;
			b_hi ^= (b_hi ^ bw) & c_hi;
			b_lo ^= (b_lo ^ bw) & c_lo;
			c_lo = c_hi;
			mw = aw | bw;
			c_hi &= ((mw | -mw) >> 31) - (uint32_t)1;
		}

		/*
		 * If c_lo = 0, then we grabbed two words for a and b.
		 * If c_lo != 0 but c_hi = 0, then we grabbed one word
		 * (in a_hi / b_hi), which means that both values are
		 * at most 32 bits.
		 * It is not possible that c_hi != 0 because b != 0,
		 * and thus we must have encountered a non-zero word.
		 *
		 * Note that a_hi and b_hi cannot be both zero; therefore,
		 * count_leading_zeros() will return a value in the 0..31
		 * range.
		 *
		 * We rely on the assumption that left/right shifts will
		 * not reveal the shift count through timing (this is
		 * a valid assumption on ARM CPUs, including the small
		 * ones from the Cortex M line; on x86, this is true as
		 * well, except for the infamous Pentium IV).
		 */
		s = count_leading_zeros(a_hi | b_hi);
		xa = (a_hi << s) | ((a_lo >> 1) >> (31 - s));
		xb = (b_hi << s) | ((b_lo >> 1) >> (31 - s));
		xa = (xa & 0xFFFF8000) | (a.v[0] & 0x00007FFF);
		xb = (xb & 0xFFFF8000) | (b.v[0] & 0x00007FFF);

		/*
		 * If c_lo != 0, then the computed values for xa and xb
		 * should be ignored, since both values would fit on one
		 * word only.
		 */
		xa ^= c_lo & (xa ^ a.v[0]);
		xb ^= c_lo & (xb ^ b.v[0]);

		/*
		 * We can now run the binary GCD on xa and xb for 15
		 * rounds.
		 */
		fg0 = (uint32_t)1;
		fg1 = (uint32_t)1 << 16;
		for (j = 0; j < 15; j ++) {
			uint32_t a_odd, swap, t;

			a_odd = -(uint32_t)(xa & 1);
			swap = a_odd & -(uint32_t)_subborrow_u32(0, xa, xb, &t);
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
		fg0 += 0x7FFF7FFF;
		fg1 += 0x7FFF7FFF;
		f0 = (fg0 & 0xFFFF) - (uint32_t)0x7FFF;
		g0 = (fg0 >> 16) - (uint32_t)0x7FFF;
		f1 = (fg1 & 0xFFFF) - (uint32_t)0x7FFF;
		g1 = (fg1 >> 16) - (uint32_t)0x7FFF;

		/*
		 * We now need to propagate updates to a, b, u and v.
		 */
		nega = s256_lin_div15_abs(&na, &a, &b, f0, g0);
		negb = s256_lin_div15_abs(&nb, &a, &b, f1, g1);
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
	 * and len(a) + len(b) <= 30, so it is known that the values
	 * fully fit in a single register each. We can do the remaining
	 * 29 iterations in one go (they are exact, no approximation
	 * here). In fact, we can content ourselves with 28 iterations,
	 * because when arriving at the last iteration, we know that a = 0
	 * or 1 and b = 1 (the final state of the algorithm is that a = 0
	 * and b is the GCD of y and q), so there would be no swap. Since
	 * we only care about the update factors f1 and g1, we can simply
	 * avoid the final iteration.
	 *
	 * The update values f1 and g1, for v, will be up to 2^28 (in
	 * absolute value) but this is supported by gf_lin().
	 *
	 * If y is zero, then none of the iterations changes anything
	 * to b and v, and therefore v = 0 at this point, which is the
	 * value we want to obtain in that case.
	 */
	xa = a.v[0];
	xb = b.v[0];

	f0 = 1;
	g0 = 0;
	f1 = 0;
	g1 = 1;
	for (j = 0; j < 28; j ++) {
		uint32_t a_odd, swap, t;

		a_odd = -(uint32_t)(xa & 1);
		swap = a_odd & -(uint32_t)_subborrow_u32(0, xa, xb, &t);
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
	 * We did 32*15+28 = 508 iterations, and each injected a factor 2,
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
 *  -1   if y != 0 and is a not a quadratic residue
 *   0   if y == 0
 */
UNUSED
static int32_t
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
	uint32_t f0, f1, g0, g1, xa, xb, fg0, fg1, ls, a0, b0;
	uint32_t nega;
	int i, j;

	gf_normalize(&a, y);
	b = GF_P;
	ls = 0;

	/*
	 * Generic loop first does 32*15 = 480 iterations.
	 */
	for (i = 0; i < 32; i ++) {
		uint32_t c_hi, c_lo, a_hi, a_lo, b_hi, b_lo, s;
		gf na, nb;

		/*
		 * Get approximations of a and b over 32 bits:
		 *  - If len(a) <= 32 and len(b) <= 32, then we just
		 *    use the value (low limb).
		 *  - Otherwise, with n = max(len(a), len(b)), we use:
		 *       (a mod 2^15) + 2^15*(floor(a / 2^(n-17)))
		 *       (b mod 2^15) + 2^15*(floor(b / 2^(n-17)))
		 * I.e. we remove the "middle bits".
		 */
		c_hi = (uint32_t)-1;
		c_lo = (uint32_t)-1;
		a_hi = 0;
		a_lo = 0;
		b_hi = 0;
		b_lo = 0;
		for (j = 7; j >= 0; j --) {
			uint32_t aw, bw, mw;

			aw = a.v[j];
			bw = b.v[j];
			a_hi ^= (a_hi ^ aw) & c_hi;
			a_lo ^= (a_lo ^ aw) & c_lo;
			b_hi ^= (b_hi ^ bw) & c_hi;
			b_lo ^= (b_lo ^ bw) & c_lo;
			c_lo = c_hi;
			mw = aw | bw;
			c_hi &= ((mw | -mw) >> 31) - (uint32_t)1;
		}

		/*
		 * If c_lo = 0, then we grabbed two words for a and b.
		 * If c_lo != 0 but c_hi = 0, then we grabbed one word
		 * (in a_hi / b_hi), which means that both values are
		 * at most 32 bits.
		 * It is not possible that c_hi != 0 because b != 0,
		 * and thus we must have encountered a non-zero word.
		 *
		 * Note that a_hi and b_hi cannot be both zero; therefore,
		 * count_leading_zeros() will return a value in the 0..31
		 * range.
		 *
		 * We rely on the assumption that left/right shifts will
		 * not reveal the shift count through timing (this is
		 * a valid assumption on ARM CPUs, including the small
		 * ones from the Cortex M line; on x86, this is true as
		 * well, except for the infamous Pentium IV).
		 */
		s = count_leading_zeros(a_hi | b_hi);
		xa = (a_hi << s) | ((a_lo >> 1) >> (31 - s));
		xb = (b_hi << s) | ((b_lo >> 1) >> (31 - s));
		xa = (xa & 0xFFFF8000) | (a.v[0] & 0x00007FFF);
		xb = (xb & 0xFFFF8000) | (b.v[0] & 0x00007FFF);

		/*
		 * If c_lo != 0, then the computed values for xa and xb
		 * should be ignored, since both values would fit on one
		 * word only.
		 */
		xa ^= c_lo & (xa ^ a.v[0]);
		xb ^= c_lo & (xb ^ b.v[0]);

		/*
		 * We can now run the binary GCD on xa and xb for 13
		 * rounds.
		 */
		fg0 = (uint32_t)1;
		fg1 = (uint32_t)1 << 16;
		for (j = 0; j < 13; j ++) {
			uint32_t a_odd, swap, t;

			a_odd = -(uint32_t)(xa & 1);
			swap = a_odd & -(uint32_t)_subborrow_u32(0, xa, xb, &t);
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
		f0 = ((fg0 + 0x7FFF7FFF) & 0xFFFF) - (uint32_t)0x7FFF;
		g0 = ((fg0 + 0x7FFF7FFF) >> 16) - (uint32_t)0x7FFF;
		f1 = ((fg1 + 0x7FFF7FFF) & 0xFFFF) - (uint32_t)0x7FFF;
		g1 = ((fg1 + 0x7FFF7FFF) >> 16) - (uint32_t)0x7FFF;
		a0 = (a.v[0] * f0 + b.v[0] * g0) >> 13;
		b0 = (a.v[0] * f1 + b.v[0] * g1) >> 13;
		for (j = 0; j < 2; j ++) {
			uint32_t a_odd, swap, t;

			a_odd = -(uint32_t)(xa & 1);
			swap = a_odd & -(uint32_t)_subborrow_u32(0, xa, xb, &t);
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

		fg0 += 0x7FFF7FFF;
		fg1 += 0x7FFF7FFF;
		f0 = (fg0 & 0xFFFF) - (uint32_t)0x7FFF;
		g0 = (fg0 >> 16) - (uint32_t)0x7FFF;
		f1 = (fg1 & 0xFFFF) - (uint32_t)0x7FFF;
		g1 = (fg1 >> 16) - (uint32_t)0x7FFF;

		/*
		 * We now need to propagate updates to a and b.
		 */
		nega = s256_lin_div15_abs(&na, &a, &b, f0, g0);
		s256_lin_div15_abs(&nb, &a, &b, f1, g1);
		ls += nega & (nb.v[0] >> 1);
		a = na;
		b = nb;
	}

	/*
	 * At that point, if y is invertible, then the final GCD is 1,
	 * and len(a) + len(b) <= 30, so it is known that the values
	 * fully fit in a single register each. We can do the remaining
	 * 29 iterations in one go. We do not need to compute update
	 * factors. Since values are exact, we always have the proper
	 * bit values. We can even stop at 28 iterations, because at the
	 * last iteration, b = 1 and a = 0 or 1; in either case, no
	 * modification of the Legendre symbol will happen.
	 */
	xa = a.v[0];
	xb = b.v[0];
	for (j = 0; j < 28; j ++) {
		uint32_t a_odd, swap, t;

		a_odd = -(uint32_t)(xa & 1);
		swap = a_odd & -(uint32_t)_subborrow_u32(0, xa, xb, &t);
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
	return (int32_t)((uint32_t)(1 - ((int32_t)(ls & 1) << 1))
		& (gf_iszero(y) - 1));
}

static uint32_t
dec32le(const uint8_t *buf)
{
	return (uint32_t)buf[0]
		| ((uint32_t)buf[1] << 8)
		| ((uint32_t)buf[2] << 16)
		| ((uint32_t)buf[3] << 24);
}

static void
enc32le(uint8_t *buf, uint32_t x)
{
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
}

/*
 * Decode a value from 32 bytes. Decoding always succeeds. Returned value
 * is 1 if the value was normalized (in the 0..p-1 range), 0 otherwise.
 */
UNUSED
static uint32_t
gf_decode(gf *d, const void *src)
{
	const uint8_t *buf;
	int i;
	uint32_t t;
	unsigned char cc;

	buf = src;
	for (i = 0; i < 8; i ++) {
		d->v[i] = dec32le(buf + 4 * i);
	}
	cc = 0;
	for (i = 0; i < 8; i ++) {
		cc = _subborrow_u32(cc, d->v[i], GF_P.v[i], &t);
	}
	return (uint32_t)cc;
}

/*
 * Decode a value from bytes, with modular reduction.
 */
UNUSED
static void
gf_decode_reduce(gf *d, const void *src, size_t len)
{
	const uint8_t *buf;
	uint32_t dd[8];
	int i;

	buf = src;
	if (len == 0) {
		*d = GF_ZERO;
		return;
	} else if ((len & 31) != 0) {
		uint8_t tmp[32];
		size_t n;

		n = len & 31;
		len -= n;
		memcpy(tmp, buf + len, n);
		memset(tmp + n, 0, (sizeof tmp) - n);
		for (i = 0; i < 8; i ++) {
			dd[i] = dec32le(tmp + 4 * i);
		}
	} else {
		len -= 32;
		for (i = 0; i < 8; i ++) {
			dd[i] = dec32le(buf + len + 4 * i);
		}
	}

	while (len > 0) {
		uint32_t e[8], g;
		uint64_t z;
		unsigned char cc;

		len -= 32;
		for (i = 0; i < 8; i ++) {
			e[i] = dec32le(buf + len + 4 * i);
		}

		g = 0;
		cc = 0;
		for (i = 0; i < 8; i ++) {
			z = (uint64_t)dd[i] * (uint64_t)(2 * MQ) + (uint64_t)g;
			g = (uint32_t)(z >> 32);
			cc = _addcarry_u32(cc, e[i], (uint32_t)z, &dd[i]);
		}
		g += cc;

		g = (g << 1) | (dd[7] >> 31);
		dd[7] &= 0x7FFFFFFF;
		z = (uint64_t)g * (uint64_t)MQ;
		cc = _addcarry_u32(0, dd[0], (uint32_t)z, &dd[0]);
		cc = _addcarry_u32(cc, dd[1], (uint32_t)(z >> 32), &dd[1]);
		for (i = 2; i < 8; i ++) {
			cc = _addcarry_u32(cc, dd[i], 0, &dd[i]);
		}
	}

	for (i = 0; i < 8; i ++) {
		d->v[i] = dd[i];
	}
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
	int i;

	gf_normalize(&t, a);
	buf = dst;
	for (i = 0; i < 8; i ++) {
		enc32le(buf + 4 * i, t.v[i]);
	}
}
