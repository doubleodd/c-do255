/*
 * This file implements elementary operations in finite field GF(2^255-MQ),
 * for ARM Cortex M0 and M0+.
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
void CN(gf_add)(gf *d, const gf *a, const gf *b);
#define gf_add   CN(gf_add)

/* d <- a - b */
void CN(gf_sub)(gf *d, const gf *a, const gf *b);
#define gf_sub   CN(gf_sub)

/* d <- -a */
void CN(gf_neg)(gf *d, const gf *a);
#define gf_neg   CN(gf_neg)

/* If ctl = 1, then d <- -a
   If ctl = 0, then d <- a
   ctl MUST be 0 or 1 */
void CN(gf_condneg)(gf *d, const gf *a, uint32_t ctl);
#define gf_condneg   CN(gf_condneg)

/* d <- a - b - c */
void CN(gf_sub2)(gf *d, const gf *a, const gf *b, const gf *c);
#define gf_sub2   CN(gf_sub2)

/* d <- a/2 */
void CN(gf_half)(gf *d, const gf *a);
#define gf_half   CN(gf_half)

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
 * precedence (d is set to a).
 */
void CN(gf_sel3)(gf *d, const gf *a, const gf *b, const gf *c,
	uint32_t fa, uint32_t fb);
#define gf_sel3   CN(gf_sel3)

/* d <- 2*a */
void CN(gf_mul2)(gf *d, const gf *a);
#define gf_mul2   CN(gf_mul2)

/* d <- 4*a */
void CN(gf_mul4)(gf *d, const gf *a);
#define gf_mul4   CN(gf_mul4)

/* d <- 8*a */
void CN(gf_mul8)(gf *d, const gf *a);
#define gf_mul8   CN(gf_mul8)

/* d <- 32*a */
void CN(gf_mul32)(gf *d, const gf *a);
#define gf_mul32   CN(gf_mul32)

/* d <- a*b, with b < 2^16 */
void CN(gf_mul_small)(gf *d, const gf *a, uint32_t b);
#define gf_mul_small   CN(gf_mul_small)

/* d <- a*b  (never actually inlined) */
void CN(gf_mul)(gf *d, const gf *a, const gf *b);
#define gf_mul          CN(gf_mul)
#define gf_mul_inline   gf_mul

/* d <- a^2  (never actually inlined) */
void CN(gf_sqr)(gf *d, const gf *a);
#define gf_sqr          CN(gf_sqr)
#define gf_sqr_inline   gf_sqr

/* d <- a^(2^num)  (repeated squarings, never actually inlined) */
void CN(gf_sqr_x)(gf *d, const gf *a, unsigned num);
#define gf_sqr_x          CN(gf_sqr_x)
#define gf_sqr_x_inline   gf_sqr_x

/*
 * Normalize a value to 0..p-1.
 */
void CN(gf_normalize)(gf *d, const gf *a);
#define gf_normalize   CN(gf_normalize)

/*
 * Compare 'a' with zero; returned value is 1 if the value is zero (as
 * a field element), 0 otherwise.
 */
uint32_t CN(gf_iszero)(const gf *a);
#define gf_iszero   CN(gf_iszero)

/*
 * Compare two values as field elements; returned value is 1 on equality,
 * 0 otherwise.
 */
uint32_t CN(gf_eq)(const gf *a, const gf *b);
#define gf_eq   CN(gf_eq)

/*
 * Inversion in the field: d <- 1/y
 * If y = 0, then d is set to zero.
 * Returned value is 1 if the value was invertible, 0 otherwise.
 */
uint32_t CN(gf_inv)(gf *d, const gf *y);
#define gf_inv   CN(gf_inv)

/*
 * Legendre symbol computation. Return value:
 *   1   if y != 0 and is a quadratic residue
 *  -1   if y != 0 and is not a quadratic residue
 *   0   if y == 0
 */
int32_t CN(gf_legendre)(const gf *y);
#define gf_legendre   CN(gf_legendre)

/*
 * Decode a value from 32 bytes. Decoding always succeeds. Returned value
 * is 1 if the value was normalized (in the 0..p-1 range), 0 otherwise.
 */
uint32_t CN(gf_decode)(gf *d, const void *src);
#define gf_decode   CN(gf_decode)

/*
 * Decode a value from bytes, with modular reduction.
 */
uint32_t CN(gf_decode_reduce)(gf *d, const void *src, size_t len);
#define gf_decode_reduce   CN(gf_decode_reduce)

/*
 * Encode a value into 32 bytes. Encoding is always normalized; a
 * consequence is that the highest bit of the last byte is always 0.
 */
void CN(gf_encode)(void *dst, const gf *d);
#define gf_encode   CN(gf_encode)
