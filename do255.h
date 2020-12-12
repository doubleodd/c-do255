#ifndef DO255_H__
#define DO255_H__

#include <stddef.h>
#include <stdint.h>

/* ==================================================================== */

/*
 * Generic type for a 256-bit integer. This can be used for signed
 * integers (-2^255 to +2^255-1) or unsigned integers (0 to 2^256-1).
 * Each implementation will use one of the three variants, depending on
 * internal word size. Contents are not defined by this header; there
 * is no attempt at compatibility between variants with different word
 * sizes.
 */
typedef union {
	struct do255_int256_w64 { uint64_t v0, v1, v2, v3; } w64;
	struct do255_int256_w32 { uint32_t v[8]; } w32;
	struct do255_int256_w16 { uint16_t v[16]; } w16;
} do255_int256;

/*
 * For each curve, the curve order is 2*r for a prime r, but we always
 * stay on points which are NOT of r-torsion. If N is the point (0,0)
 * (unique point of order 2), then the group law (which will hereafter
 * call addition) between elements P1 and P2 is equal to the sum P1+P2+N
 * over the full curve. The neutral is N; there is no "point at
 * infinity".
 *
 * Curve do255s:
 *    Field: integers modulo p = 2^255 - 3957
 *    Equation: y^2 = x*(x^2 - x + 1/2)
 *    Order: 2*r, with r = 2^254 + 56904135270672826811114353017034461895
 *
 * Curve do255e:
 *    Field: integers modulo p = 2^255 - 18651
 *    Equation: y^2 = x*(x^2 - 2)
 *    Order: 2*r, with r = 2^254 - 131528281291764213006042413802501683931
 *
 * Jacobian (x,w) coordinates: (X:W:Z)
 *    x = X/Z^2
 *    w = y/x = W/Z
 *    W != 0
 *    neutral: N = (0:W:0) for any W != 0
 *
 * Conventional generator is point G = (x,w) with w = 1/u, where u is
 * the smallest possible non-zero integer (in the 1..p-1 range).
 *
 * All points can be encoded into a sequence of 32 bytes, and decoded from
 * a sequence of 32 bytes. Encoding is canonical and verified: there is a
 * unique way to encode a point (group element), and the decoding function
 * checks that the obtained point would indeed encode to exactly the source
 * bytes. In the encoding format, the highest bit of the last byte is
 * always zero; thus, this bit can be used to store other data, as long as
 * it is cleared before calling the point decoding function. The neutral
 * element encodes as an all-zero sequence.
 *
 * Unless explicitly noted otherwise:
 *  - All functions are complete: input points can be the neutral, and
 *    can be equal to each other.
 *  - All functions are constant-time: execution time and memory access
 *    pattern are independent of the input or output values (exception:
 *    the *_verify_helper_vartime() functions are not constant-time, as
 *    they are meant to be used on public values only).
 *  - Overlaps are allowed between inputs and outputs; e.g. the point
 *    addition functions can use the same structure as destination and
 *    as one or both sources (however, _partial_ overlaps, in which
 *    an input and an output structure are not at the same address
 *    but still share some storage bytes, are forbidden).
 *
 * In general, the right-to-left convention is used: output parameters
 * come before input parameters (similar to the C function memcpy()).
 */

/*
 * Point encoding structures.
 */
typedef struct {
	do255_int256 X, W, Z;
} do255s_point;
typedef struct {
	do255_int256 X, W, Z;
} do255e_point;

/*
 * A statically allocated neutral element.
 */
extern const do255s_point do255s_neutral;
extern const do255e_point do255e_neutral;

/*
 * A statically allocated conventional group generator.
 */
extern const do255s_point do255s_generator;
extern const do255e_point do255e_generator;

/*
 * Decode a point from 32 bytes. Returned value is 1 on success, 0 on
 * failure. A failure is reported for both invalid and non-canonical
 * encodings. If a failure is reported, then the point 'P' is set to
 * the neutral.
 *
 * All 256 bits of the input are used. For a canonical encoding, the
 * highest bit of the last byte MUST be zero. If that bit is reused
 * for storing some other data, then it MUST be cleared before calling
 * this function.
 *
 * Note that a canonical encoding of the neutral (all 32 bytes are zero)
 * will decode successfully (this function will return 1).
 *
 * If P == NULL, then the decoded point is not returned; the function
 * will still verify that the input is a canonical encoding of a proper
 * point.
 */
int do255s_decode(do255s_point *P, const void *src);
int do255e_decode(do255e_point *P, const void *src);

/*
 * Encode a point into 32 bytes. This function can only produce a
 * canonical encoding.
 */
void do255s_encode(void *dst, const do255s_point *P);
void do255e_encode(void *dst, const do255e_point *P);

/*
 * Get the squared w coordinate of a point. That value is normalized
 * (into 0..p-1 range) and encoded over 32 bytes, with unsigned
 * little-endian convention. This function is meant to support ECDH.
 */
void do255s_encode_squared_w(void *dst, const do255s_point *P);
void do255e_encode_squared_w(void *dst, const do255e_point *P);

/*
 * Compare a point with the neutral element. Returned value is 1 if the
 * point is the neutral, 0 otherwise.
 */
int do255s_is_neutral(const do255s_point *P);
int do255e_is_neutral(const do255e_point *P);

/*
 * Compare two points. Returned value is 1 if the points are equal to
 * each other, 0 otherwise.
 */
int do255s_eq(const do255s_point *P1, const do255s_point *P2);
int do255e_eq(const do255e_point *P1, const do255e_point *P2);

/*
 * Add two points together (P3 <- P1 + P2).
 */
void do255s_add(do255s_point *P3,
	const do255s_point *P1, const do255s_point *P2);
void do255e_add(do255e_point *P3,
	const do255e_point *P1, const do255e_point *P2);

/*
 * Double a point (P3 <- 2*P1).
 */
void do255s_double(do255s_point *P3, const do255s_point *P1);
void do255e_double(do255e_point *P3, const do255e_point *P1);

/*
 * Double a point n times (P3 <- (2^n)*P1). This is equivalent to
 * calling the do255*_double() function n times, but it may be slightly
 * faster in some implementations.
 *
 * The function is not constant-time relatively to n (the number of
 * point doublings): processing time will be proportional to the value
 * of n.
 */
void do255s_double_x(do255s_point *P3, const do255s_point *P1, unsigned n);
void do255e_double_x(do255e_point *P3, const do255e_point *P1, unsigned n);

/*
 * Negate a point (P3 <- -P1).
 */
void do255s_neg(do255s_point *P3, const do255s_point *P1);
void do255e_neg(do255e_point *P3, const do255e_point *P1);

/*
 * Multiply a point by a scalar (P3 <- scalar*P1).
 *
 * The scalar is provided as a sequence of bytes in unsigned
 * little-endian convention, over exactly 32 bytes. The scalar value
 * may range up to 2^256-1 (inclusive).
 *
 * Note that the final point P3 may be the neutral, if P1 is itself
 * the neutral, and/or the scalar is a multiple of r.
 */
void do255s_mul(do255s_point *P3, const do255s_point *P1, const void *scalar);
void do255e_mul(do255e_point *P3, const do255e_point *P1, const void *scalar);

/*
 * Multiply the conventional generator by a scalar (P3 <- scalar*G).
 *
 * The scalar is provided as a sequence of bytes in unsigned
 * little-endian convention, over exactly 32 bytes. The scalar value
 * may range up to 2^256-1 (inclusive).
 *
 * Note that the final point P3 may be the neutral if the scalar is a
 * multiple of r.
 *
 * These functions are usually faster than calling do255*_mul() over a
 * point structure set to a copy of the generator.
 */
void do255s_mulgen(do255s_point *P3, const void *scalar);
void do255e_mulgen(do255e_point *P3, const void *scalar);

/*
 * Signature verification helper: given scalars k0 and k1 and points P
 * and R, verify that k0*G + k1*P = R (with G being the conventional
 * generator element). Scalars k0 and k1 are encoded over 256 bits (32
 * bytes, unsigned little-endian convention) and may be up to 2^256-1.
 * Point R is provided encoded as R_enc (32 bytes).
 *
 * Return value is 1 if the equation is fulfilled, 0 otherwise.
 *
 * This function is not necessarily constant-time; side channels may
 * leak information about the points P and R, the scalars k0 and k1, and
 * the overall outcome. This is appropriate for signature verification,
 * where signature and public key and considered to be public data.
 */
int do255e_verify_helper_vartime(const void *k0,
	const do255e_point *P, const void *k1, const void *R_enc);
int do255s_verify_helper_vartime(const void *k0,
	const do255s_point *P, const void *k1, const void *R_enc);

/*
 * Map a source value (arbitrary sequence of 'len' bytes) onto a point.
 * The mapping is not one-way and not uniform, but can be used to
 * implement a hash-to-curve process by doing the following:
 *  - Hash the input data into 512 bits (64 bytes) with a suitable
 *    hash function (e.g. SHA-512 or SHAKE256).
 *  - Split the 64 bytes into two halves of 32 bytes.
 *  - Map each half onto a point.
 *  - Add the two points together.
 */
void do255e_map_to_curve(do255e_point *P, const void *src, size_t len);
void do255s_map_to_curve(do255s_point *P, const void *src, size_t len);

/*
 * Scalar check: the input value is a 256-bit integer (unsigned
 * little-endian convention). Returned value is 1 if that integer
 * is less than r, 0 otherwise.
 */
int do255e_scalar_is_reduced(const void *a);
int do255s_scalar_is_reduced(const void *a);

/*
 * Scalar modular reduction: the input unsigned integer 'a' (in
 * little-endian order, over 'a_len' bytes) is reduced from r (the
 * curve subgroup order) and encoded over exactly 32 bytes into 'd'
 * (unsigned little-endian encoding). This function is constant-time
 * with regard to the integer value (but the number of source bytes
 * 'a_len' may leak).
 * Input and output buffers may overlap freely.
 */
void do255e_scalar_reduce(void *d, const void *a, size_t a_len);
void do255s_scalar_reduce(void *d, const void *a, size_t a_len);

/*
 * Scalar modular addition: given two input scalars 'a' and 'b'
 * (unsigned little-endian order, 32 bytes each), their sum modulo
 * r is written into 'd' (exactly 32 bytes). The two input values
 * may be arbitrary in their range (0 to 2^256-1); the output is
 * reduced modulo r (value is in 0..r-1).
 * Input and output buffers may overlap freely.
 */
void do255e_scalar_add(void *d, const void *a, const void *b);
void do255s_scalar_add(void *d, const void *a, const void *b);

/*
 * Scalar modular subtraction: given two input scalars 'a' and 'b'
 * (unsigned little-endian order, 32 bytes each), their difference
 * modulo r (i.e. a-b mod r) is written into 'd' (exactly 32 bytes).
 * The two input values may be arbitrary in their range (0 to 2^256-1);
 * the output is reduced modulo r (value is in 0..r-1).
 * Input and output buffers may overlap freely.
 */
void do255e_scalar_sub(void *d, const void *a, const void *b);
void do255s_scalar_sub(void *d, const void *a, const void *b);

/*
 * Scalar modular negation: given input 'a' (unsigned little-endian order,
 * 32 bytes), -a mod r is written into 'd' (exactly 32 bytes). The
 * input value may be arbitrary in its range (0 to 2^256-1); the output is
 * reduced modulo r (value is in 0..r-1).
 */
void do255e_scalar_neg(void *d, const void *a);
void do255s_scalar_neg(void *d, const void *a);

/*
 * Scalar modular halving: given input 'a' (unsigned little-endian order,
 * 32 bytes), a/2 mod r is written into 'd' (exactly 32 bytes). The
 * input value may be arbitrary in its range (0 to 2^256-1); the output is
 * reduced modulo r (value is in 0..r-1).
 */
void do255e_scalar_half(void *d, const void *a);
void do255s_scalar_half(void *d, const void *a);

/*
 * Scalar modular multiplication: given two input scalars 'a' and 'b'
 * (unsigned little-endian order, 32 bytes each), their product
 * modulo r (i.e. a*b mod r) is written into 'd' (exactly 32 bytes).
 * The two input values may be arbitrary in their range (0 to 2^256-1);
 * the output is reduced modulo r (value is in 0..r-1).
 * Input and output buffers may overlap freely.
 */
void do255e_scalar_mul(void *d, const void *a, const void *b);
void do255s_scalar_mul(void *d, const void *a, const void *b);

/* ==================================================================== */

#endif
