/*
 * This file implements Lagrange's algorithm for reduction of a lattice
 * basis in dimension 2 (optimized). It is used for signature verification
 * in do255s (curve do255e instead uses the curve endomorphism to achieve
 * a similar optimization with lower overhead). This implementation is
 * for ARM architectures.
 *
 * This file is not meant to be compiled by itself, but included in an outer
 * file.
 */

/*
 * Apply Lagrange's algorithm on the provided scalar k. Result is
 * encoded over two 17-byte arrays c0 and c1 (signed integers,
 * absolute value is lower than 2^128, little-endian encoding).
 * (implemented in assembly)
 */
void CN(reduce_basis_vartime)(uint8_t *c0, uint8_t *c1, const uint8_t *k);
#define reduce_basis_vartime   CN(reduce_basis_vartime)
