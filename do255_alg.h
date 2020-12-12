#ifndef DO255_ALG_H__
#define DO255_ALG_H__

#include <stddef.h>
#include <stdint.h>
#include "sha3.h"

/* ==================================================================== */
/*
 * High-level cryptographic operations using curves do255e and do255s.
 *
 * Private keys are non-zero integers modulo the curve subgroup order r.
 * All functions below accept arbitrary 256-bit integers as input; the
 * key pair generation functions produce only fully reduced values, i.e.
 * in the 1 to r-1 range. Encoding uses unsigned little-endian convention.
 * Since r < 2^255, this implies that when a fully reduced private key is
 * encoded, the top bit of the last byte is 0 (for curve do255e, where r
 * is slightly below 2^254, the two top bits of the last byte are 0).
 *
 * Public keys are encoded curve points. A curve point is encoded over
 * exactly 32 bytes. The top bit of the last byte is always 0 when
 * encoding; that bit is verified when decoding (i.e. if it is not 0,
 * then decoding fails). The neutral element in the curve subgroup we
 @ are using is not accepted as a public key; key exchange and signature
 * verification functions explicitly reject such a public key.
 *
 * In this API, private keys and public keys are wrapped into structures
 * with a single field called 'b', that contains the byte representation
 * of the private or public key, as explained above. The structures are
 * meant to provide some type checking, to avoid inadvertently using a
 * public key where a private key is expected, or vice versa. Input/output
 * (writing and reading keys) is meant to be handled by the application,
 * by reading from and writing to the b[] field of the relevant structures.
 *
 * There is no type-based distinction between keys used for signatures,
 * and keys used for key exchange. However, it is NOT RECOMMENDED to use
 * the same private key for both key exchange and signatures.
 *
 *
 * A signature consists of the concatenation of the encoding of a curve
 * point, and an integer modulo r. The same encoding rules as above are
 * used, for these two elements, except that:
 *  - The point is allowed to be the neutral element (encoded as 32 bytes
 *    of value 0x00).
 *  - The integer modulo r is allowed to be zero (again encoded as 32 bytes
 *    of value 0x00).
 * Total signature size is 64 bytes. As explained above, 2 or 3 specific
 * bits in the signature (depending on the used curve) are always 0, and
 * could be reused for storing other values, provided that they are set
 * back to 0 before trying to verify the signature. A dedicated wrapper
 * structure for signatures is provided.
 *
 * Signatures are not malleable (a given signature can be encoded in
 * only one way that the verifying function will accept; making a
 * distinct but equivalent signature on the same data, relatively to the
 * same public key, requires knowledge of the private key). Signature
 * generation is by default deterministic (no reliance on a random
 * source) but can be made non-deterministic by providing a non-empty
 * extra seed. Signature security is maintained regardless of whether
 * the seed is present or not, and regardless of the quality of the
 * randomness in the seed if present. Injecting a varying seed (even if
 * not random) can provide additional security benefits against some
 * classes of physical attacks, in particular fault attacks.
 *
 *
 * All functions are constant-time, except signature verification, which
 * is explicitly tagged "vartime" (and normally works only over public
 * data).
 */

/*
 * Private key structure. The b[] field contains the encoded private key.
 */
typedef struct {
	unsigned char b[32];
} do255e_private_key;
typedef struct {
	unsigned char b[32];
} do255s_private_key;

/*
 * Public key structure. The b[] field contains the encoded public key.
 */
typedef struct {
	unsigned char b[32];
} do255e_public_key;
typedef struct {
	unsigned char b[32];
} do255s_public_key;

/*
 * Signature structure. The b[] field contains the encoded signature.
 */
typedef struct {
	unsigned char b[64];
} do255e_signature;
typedef struct {
	unsigned char b[64];
} do255s_signature;

/*
 * Key pair generation.
 *
 * A private/public key pair is produced, using the provided random
 * generator as source of randomness. That generator is a SHAKE context
 * which must have already been initialized, filled with enough entropy,
 * and "flipped" into output generation mode. The private key is written
 * in *sk, and the public key in *pk. Both have size exactly 32
 * bytes. If pk is NULL, then only the private key is produced.
 */
void do255e_keygen(shake_context *rng,
	do255e_private_key *sk, do255e_public_key *pk);
void do255s_keygen(shake_context *rng,
	do255s_private_key *sk, do255s_public_key *pk);

/*
 * Rebuild the public key (pk) from the private key (sk). Both keys have
 * size exactly 32 bytes.
 *
 * Normally, when storing a private key, the public key should be stored
 * along with it. If storage capacity is scarce in a given application,
 * it is possible to store only the private key, and use these functions
 * to recompute the public key dynamically when needed. Recomputing the
 * public key from the private key has a cost similar to full key pair
 * generation, or to signature generation.
 */
void do255e_make_public(do255e_public_key *pk, const do255e_private_key *sk);
void do255s_make_public(do255s_public_key *pk, const do255s_private_key *sk);

/*
 * Verify a public key. If the provided public key (pk) is the canonical
 * encoding of a valid non-neutral group point, then this function
 * returns 1; otherwise, it returns 0.
 * Note that other functions in this API that consume public keys
 * inherently verify public key validity; it is thus not necessary to
 * call this function prior to using an incoming public key in other
 * functions such as do255e_key_exchange() or do255e_verify_vartime().
 */
int do255e_check_public(const do255e_public_key *pk);
int do255s_check_public(const do255s_public_key *pk);

/*
 * Compute a key exchange. The private key sk is combined with the public
 * key from the peer (pk_peer) to produce the shared secret. An internal
 * key derivation function (SHAKE256) is used to produce the shared
 * secret, which can have arbitrary length (secret_len, expressed in bytes),
 * and is suitable for use in other cryptographic algorithms.
 *
 * The public key (pk_self) corresponding to the private sk must also be
 * provided: it is not part of curve computation itself, but is used in
 * the key derivation. This "binds" the output to the exact points that
 * were used, which is convenient in some security proofs. This function
 * DOES NOT verify that pk_self corresponds to sk or even is a valid
 * public key.
 *
 * Returned value is 1 if the exchange succeeded, 0 otherwise. An error
 * status (0) is returned if the public key from the peer (pk_peer) is
 * not a valid point, or is the subgroup neutral point, or if our own
 * private key (sk) is invalid (i.e. is equal to 0 modulo r). When such
 * an error case is reached, the generated secret is deterministically
 * derived from sk, pk_self and pk_peer in such a way that the resulting
 * secret value is unpredictable by outsiders who do not know sk.
 *
 * This function is constant-time with regard to all parameters, and also
 * to the returned status: no timing-based side channel may reveal
 * whether the call succeeded or failed.
 */
int do255e_key_exchange(void *secret, size_t secret_len,
	const do255e_private_key *sk,
	const do255e_public_key *pk_self, const do255e_public_key *pk_peer);
int do255s_key_exchange(void *secret, size_t secret_len,
	const do255s_private_key *sk,
	const do255s_public_key *pk_self, const do255s_public_key *pk_peer);

/*
 * Sign some data with the provided private key (sk). The signature is
 * written in sig[] and has length exactly 64 bytes; it is the
 * concatenation (in that order) of an encoded curve point, and an
 * encoded fully reduced scalar value.
 *
 * The public key (pk) corresponding to the private key is also provided;
 * it is used in the internal hashing process. This function DOES NOT
 * verify that the public key matches the private key or even is a valid
 * public key.
 *
 * The data to sign is normally a hash value. An identifier for the
 * used hash value is provided as a character string (hash_oid); it
 * should be the decimal-dotted representation of the object identifier
 * (OID) allocated to the hash function. Some standard OID strings are
 * defined below for the SHA-2 and SHA-3 functions. The hash value itself
 * is given as hv[], of size hv_len (in bytes). The used hash function
 * must provide collision resistance with sufficient security.
 *
 * It is possible to set hash_oid to NULL (or an empty string), in case
 * the provided hash value is the raw data itself, not a hash thereof.
 * There is no limit on the size of that data. Take care that using raw
 * data directly, while safe, can make verification harder for
 * constrained systems, because the verifier will need to know the
 * public and the signature value before starting to process the data
 * itself; thus, not hashing the data tends to hinder streamed
 * processing of data, which is a problem for verifier with little
 * available RAM.
 *
 * An extra seed[] (of size seed_len bytes) can be provided to randomize
 * signatures. There is no hard requirements on that seed value; it can
 * have arbitrary length and needs not be unpredictable. The seed may be
 * a timestamp or a monotonous counter. The seed may even be absent (by
 * setting seed_len to 0, in which case the 'seed' pointer is ignored).
 * If there is no seed, then signatures are deterministic (you would
 * always get the same signature value for the same data and private
 * key). Deterministic signatures are safe. Adding a non-repeating seed
 * MAY provide additional resistance against some classes of physical
 * attacks (especially fault attacks).
 */
void do255e_sign(do255e_signature *sig,
	const do255e_private_key *sk, const do255e_public_key *pk,
	const char *hash_oid, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len);
void do255s_sign(do255s_signature *sig,
	const do255s_private_key *sk, const do255s_public_key *pk,
	const char *hash_oid, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len);

/*
 * Verify a signature with the provided public key. The signed data
 * is provided as a hash value (hv[], of size hv_len) and identifier
 * for the hash function (hash_oid). hash_oid may be NULL if the
 * signed data is raw; see above for details.
 *
 * Returned value is 1 on success (signature is valid), 0 otherwise
 * (signature and/or public key is invalid, or does not match the
 * data).
 *
 * This function is not constant-time: it assumes that the public key,
 * signature and signed data are public.
 */
int do255e_verify_vartime(
	const do255e_signature *sig, const do255e_public_key *pk,
	const char *hash_oid, const void *hv, size_t hv_len);
int do255s_verify_vartime(
	const do255s_signature *sig, const do255s_public_key *pk,
	const char *hash_oid, const void *hv, size_t hv_len);

/* Hash function identifier: SHA-224 */
#define DO255_OID_SHA224        "2.16.840.1.101.3.4.2.4"

/* Hash function identifier: SHA-256 */
#define DO255_OID_SHA256        "2.16.840.1.101.3.4.2.1"

/* Hash function identifier: SHA-384 */
#define DO255_OID_SHA384        "2.16.840.1.101.3.4.2.2"

/* Hash function identifier: SHA-512 */
#define DO255_OID_SHA512        "2.16.840.1.101.3.4.2.3"

/* Hash function identifier: SHA-512-224 */
#define DO255_OID_SHA512_224    "2.16.840.1.101.3.4.2.5"

/* Hash function identifier: SHA-512-256 */
#define DO255_OID_SHA512_256    "2.16.840.1.101.3.4.2.6"

/* Hash function identifier: SHA3-224 */
#define DO255_OID_SHA3_224      "2.16.840.1.101.3.4.2.7"

/* Hash function identifier: SHA3-256 */
#define DO255_OID_SHA3_256      "2.16.840.1.101.3.4.2.8"

/* Hash function identifier: SHA3-384 */
#define DO255_OID_SHA3_384      "2.16.840.1.101.3.4.2.9"

/* Hash function identifier: SHA3-512 */
#define DO255_OID_SHA3_512      "2.16.840.1.101.3.4.2.10"

#endif
