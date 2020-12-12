/*
 * This file is not meant to be compiled by itself, but is included
 * from either alg_do255e.c or alg_do255s.c.
 */

/* see do255_alg.h */
int
CN(check_public)(const CN(public_key) *pk)
{
	unsigned r, t;
	int i;

	/*
	 * Public key is valid if and only if it can decode as a
	 * non-zero field element w such that (w^2 - a)^2 - 4*b is
	 * a quadratic residue. CURVE_decode() will do all that,
	 * except that it returns 1 for an all-zero input, but we
	 * defined that public keys cannot be the group neutral.
	 */
	r = (unsigned)CN(decode)(NULL, pk->b);
	t = 0;
	for (i = 0; i < 32; i ++) {
		t |= pk->b[i];
	}
	r &= ~((t - 1) >> 8);
	return (int)r;
}

/* see do255_alg.h */
void
CN(keygen)(shake_context *rng, CN(private_key) *sk, CN(public_key) *pk)
{
	/*
	 * If the RNG produces an output which happens to be equal to
	 * zero modulo r, we try again. This is overwhelmingly
	 * improbable, and nobody knows an input to SHAKE that would
	 * lead to such an outcome; for all practical purposes, only a
	 * single iteration of the loop below is executed.
	 *
	 * We use only 32 bytes of RNG output because the group order r
	 * is very close to 2^254 (the difference is less than 2^127 for
	 * both do255e and do255s curves); thus, the bias is negligible.
	 */
	for (;;) {
		unsigned t;
		int i;

		shake_extract(rng, sk->b, 32);
		CN(scalar_reduce)(sk->b, sk->b, 32);
		t = 0;
		for (i = 0; i < 32; i ++) {
			t |= sk->b[i];
		}
		if (t == 0) {
			continue;
		}

		if (pk != NULL) {
			CN(point) Q;

			CN(mulgen)(&Q, sk->b);
			CN(encode)(pk->b, &Q);
		}
		break;
	}
}

/* see do255_alg.h */
void
CN(make_public)(CN(public_key) *pk, const CN(private_key) *sk)
{
	CN(point) Q;

	CN(mulgen)(&Q, sk->b);
	CN(encode)(pk->b, &Q);
}

/*
 * Return 1 if a < b, 0 otherwise. This uses unsigned little-endian
 * convention, over exactly 32 bytes.
 */
static unsigned
i256_is_lower(const unsigned char *a, const unsigned char *b)
{
	int i;
	unsigned cc;

	cc = 0;
	for (i = 0; i < 32; i ++) {
		unsigned x;

		x = (unsigned)a[i] - (unsigned)b[i] - cc;
		cc = (x >> 8) & 1;
	}
	return cc;
}

/*
 * If ctl == 1, then write a[] into d[]; otherwise, if ctl == 0, then
 * write b[] into d[].
 * All arrays have size 32 bytes. ctl MUST be 0 or 1. This is constant-time.
 */
static void
csel(unsigned char *d,
	const unsigned char *a, const unsigned char *b, unsigned ctl)
{
	int i;

	ctl = -ctl;
	for (i = 0; i < 32; i ++) {
		d[i] = (a[i] & ctl) | (b[i] & ~ctl);
	}
}

/* see do255_alg.h */
int
CN(key_exchange)(void *secret, size_t secret_len,
	const CN(private_key) *sk,
	const CN(public_key) *pk_self, const CN(public_key) *pk_peer)
{
	CN(point) P;
	shake_context sc;
	uint32_t f;
	unsigned char pm[32], fx;
	int i;
	unsigned ctl;

	CN(decode)(&P, pk_peer->b);
	CN(mul)(&P, &P, sk->b);

	/*
	 * The two failure conditions are an invalid input point (pk_peer),
	 * a neutral input point, and an invalid private key (sk). In all
	 * cases, this leads to point P being the neutral. If there was
	 * no failure, then P cannot be the neutral. Therefore, we only
	 * need to check here whether the neutral point was obtained.
	 */
	f = CN(is_neutral)(&P);

	/*
	 * The pre-master secret is the encoding of the square of the w
	 * coordinate of the obtained point; a square is used to make it
	 * simpler to support alternate ladder implementations of the
	 * point multiplication.
	 * On failure, we replace the pre-master secret with the private
	 * key itself.
	 */
	CN(encode_squared_w)(pm, &P);
	for (i = 0; i < 32; i ++) {
		pm[i] ^= -f & (pm[i] ^ sk->b[i]);
	}

	/*
	 * Find out whether pk_self comes before pk_peer in reverse
	 * lexicographic order (equivalently, whether pk_self is lower
	 * than pk_peer when both are interpreted as 32-byte integers
	 * in unsigned little-endian convention).
	 */
	ctl = i256_is_lower(pk_self->b, pk_peer->b);

	/*
	 * Key derivation: SHAKE256(DOM, f, pm, pk1, pk2)
	 * with:
	 *    DOM = doman-separation string ("do255e-ecdh:" or "do255s-ecdh:")
	 *    f = success/failure (one byte, 0x00 on success, 0xFF on failure)
	 *    pk1 = first public key
	 *    pk2 = second public key
	 * Public keys are ordered numerically with unsigned little-endian
	 * convention.
	 */
	shake_init(&sc, 256);
	shake_inject(&sc, DOM_ECDH, strlen(DOM_ECDH));
	fx = (unsigned char)-f;
	shake_inject(&sc, &fx, 1);
	shake_inject(&sc, pm, 32);
	csel(pm, pk_self->b, pk_peer->b, ctl);
	shake_inject(&sc, pm, 32);
	csel(pm, pk_peer->b, pk_self->b, ctl);
	shake_inject(&sc, pm, 32);
	shake_flip(&sc);
	shake_extract(&sc, secret, secret_len);
	return (int)(1 - f);
}

/*
 * Compute the "challenge" part of the signature process. It is a hash
 * value (SHAKE256 output, 32 bytes) computed over the concatenation of:
 *  - a domain separation string
 *  - the "commitment" (encoded R point, first half of the signature value)
 *  - the public key of the signer
 *  - the signed hash value, with an explicit hash function label
 */
static void
make_e(shake_context *sc, unsigned char *e, const unsigned char *R_enc,
	const CN(public_key) *pk,
	const char *hash_oid, const void *hv, size_t hv_len)
{
	shake_init(sc, 256);
	shake_inject(sc, DOM_SIGN_E, strlen(DOM_SIGN_E));
	shake_inject(sc, R_enc, 32);
	shake_inject(sc, pk->b, 32);
	shake_inject(sc, hash_oid, strlen(hash_oid));
	shake_inject(sc, ":", 1);
	shake_inject(sc, hv, hv_len);
	shake_flip(sc);
	shake_extract(sc, e, 32);
}

/* see do255_alg.h */
void
CN(sign)(CN(signature) *sig,
	const CN(private_key) *sk, const CN(public_key) *pk,
	const char *hash_oid, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len)
{
	shake_context sc;
	unsigned char tmp[8], k[32], e[32];
	CN(point) R;
	int i;

	/*
	 * Harmonize the label for non-hashed data.
	 */
	if (hash_oid == NULL) {
		hash_oid = "";
	}

	/*
	 * Make the secret k value. How we do that is not important for
	 * interoperability, but crucial for security. We use SHAKE256
	 * over the concatenation of:
	 *  - a domain separation string
	 *  - the private key
	 *  - the extra seed, if provided (with an explicit length header)
	 *  - the signed data (with label)
	 * We then extract 32 bytes and reduce that modulo r. Since r is
	 * close to 2^254 (difference is less than 2^127 for both curves),
	 * the bias is negligible.
	 */
	shake_init(&sc, 256);
	shake_inject(&sc, DOM_SIGN_K, strlen(DOM_SIGN_K));
	shake_inject(&sc, sk->b, 32);
	for (i = 0; i < 8; i ++) {
		tmp[i] = (unsigned char)((uint64_t)seed_len >> (8 * i));
	}
	shake_inject(&sc, tmp, 8);
	shake_inject(&sc, seed, seed_len);
	shake_inject(&sc, hash_oid, strlen(hash_oid));
	shake_inject(&sc, ":", 1);
	shake_inject(&sc, hv, hv_len);
	shake_flip(&sc);
	shake_extract(&sc, k, 32);
	CN(scalar_reduce)(k, k, 32);

	/*
	 * Note that k may theoretically be zero (with negligible
	 * probability). This is not a problem for our signatures, since
	 * all relevant internal functions work well with the neutral
	 * (we explicitly reject the neutral for public keys only because
	 * we want to ensure contributory behaviour for key exchange,
	 * and we don't want distinct rules for signature public keys).
	 */

	/*
	 * R = k*G.
	 */
	CN(mulgen)(&R, k);
	CN(encode)(sig->b, &R);

	/*
	 * Compute challenge e.
	 */
	make_e(&sc, e, sig->b, pk, hash_oid, hv, hv_len);

	/*
	 * Answer to challenge: s = k + e*sk.
	 * (This implicitly reduces e modulo r.)
	 */
	CN(scalar_mul)(e, e, sk->b);
	CN(scalar_add)(sig->b + 32, k, e);
}

/* see do255_alg.h */
int
CN(verify_vartime)(const CN(signature) *sig, const CN(public_key) *pk,
	const char *hash_oid, const void *hv, size_t hv_len)
{
	shake_context sc;
	unsigned char e[32];
	const unsigned char *s;
	CN(point) Q;

	static const unsigned char zero[32] = { 0 };

	if (hash_oid == NULL) {
		hash_oid = "";
	}

	/*
	 * Check that s is in the correct range; otherwise, the signature
	 * is invalid.
	 */
	s = sig->b + 32;
	if (!CN(scalar_is_reduced)(s)) {
		return 0;
	}

	/*
	 * Decode the signer public key. If it is invalid, then verification
	 * should fail.
	 */
	if (!CN(decode)(&Q, pk->b) || CN(is_neutral)(&Q)) {
		return 0;
	}

	/*
	 * Recompute challenge e then negate it.
	 */
	make_e(&sc, e, sig->b, pk, hash_oid, hv, hv_len);
	CN(scalar_sub)(e, zero, e);

	/*
	 * Signature is valid if and only if s*G = k*G + e*sk*G,
	 * i.e. s*G - e*pk = R.
	 */
	return CN(verify_helper_vartime)(s, &Q, e, sig->b);
}
