/*
 * This file implements Lagrange's algorithm for reduction of a lattice
 * basis in dimension 2 (optimized). It is used for signature verification
 * in do255s (curve do255e instead uses the curve endomorphism to achieve
 * a similar optimization with lower overhead). This implementation is
 * for architectures with 64-bit limbs.
 *
 * This file is not meant to be compiled by itself, but included in an outer
 * file.
 */

#include <immintrin.h>

static void
array_mul256x256(uint64_t *d, const uint64_t *a, const uint64_t *b)
{
	i256 ta, tb;
	i512 td;

	ta.v0 = a[0];
	ta.v1 = a[1];
	ta.v2 = a[2];
	ta.v3 = a[3];
	tb.v0 = b[0];
	tb.v1 = b[1];
	tb.v2 = b[2];
	tb.v3 = b[3];
	mul256x256(&td, &ta, &tb);
	d[0] = td.v0;
	d[1] = td.v1;
	d[2] = td.v2;
	d[3] = td.v3;
	d[4] = td.v4;
	d[5] = td.v5;
	d[6] = td.v6;
	d[7] = td.v7;
}

/*
 * Return 1 if a < b, 0 otherwise. Input values are nonnegative.
 */
static uint64_t
cmp512lt(const uint64_t *a, const uint64_t *b)
{
	unsigned char cc;
	unsigned long long t;
	int i;

	cc = _subborrow_u64(0, a[0], b[0], &t);
	for (i = 1; i < 8; i ++) {
		cc = _subborrow_u64(cc, a[i], b[i], &t);
	}
	return cc;
}

/*
 * Get the bitlength of a 512-bit signed integer (possibly negative).
 */
static int
bitlength512(const uint64_t *a)
{
	int i;
	uint64_t m;

	m = (uint64_t)(*(int64_t *)&a[7] >> 63);
	for (i = 7; i >= 0; i --) {
		uint64_t aw;

		aw = a[i] ^ m;
		if (aw != 0) {
			return 64 * i + 64 - (int)_lzcnt_u64(aw);
		}
	}
	return 0;
}

/*
 * Add b*2^s to a (512-bit integers).
 */
static void
add_lshift_512(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 8; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 8; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 8; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 8; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (512-bit integers).
 */
static void
sub_lshift_512(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 8; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 8; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 8; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 8; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Return 1 if a < b, 0 otherwise. Input values are nonnegative.
 */
static uint64_t
cmp384lt(const uint64_t *a, const uint64_t *b)
{
	unsigned char cc;
	unsigned long long t;
	int i;

	cc = _subborrow_u64(0, a[0], b[0], &t);
	for (i = 1; i < 6; i ++) {
		cc = _subborrow_u64(cc, a[i], b[i], &t);
	}
	return cc;
}

/*
 * Get the bitlength of a 384-bit signed integer (possibly negative).
 */
static int
bitlength384(const uint64_t *a)
{
	int i;
	uint64_t m;

	m = (uint64_t)(*(int64_t *)&a[5] >> 63);
	for (i = 5; i >= 0; i --) {
		uint64_t aw;

		aw = a[i] ^ m;
		if (aw != 0) {
			return 64 * i + 64 - (int)_lzcnt_u64(aw);
		}
	}
	return 0;
}

/*
 * Add b*2^s to a (384-bit integers).
 */
static void
add_lshift_384(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 6; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 6; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 6; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 6; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (384-bit integers).
 */
static void
sub_lshift_384(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 6; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 6; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 6; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 6; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Add b*2^s to a (192-bit integers).
 */
static void
add_lshift_192(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		if (s >= 192) {
			return;
		}
		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 3; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 3; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 3; i ++) {
				cc = _addcarry_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _addcarry_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 3; i ++) {
				cc = _addcarry_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (192-bit integers).
 */
static void
sub_lshift_192(uint64_t *d, const uint64_t *a, const uint64_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 64) {
		int i, j;
		unsigned char cc;

		if (s >= 192) {
			return;
		}
		j = s >> 6;
		s &= 63;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u64(0, a[j], b[0],
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 3; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i - j],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[j], b[0] << s,
				(unsigned long long *)&d[j]);
			for (i = j + 1; i < 3; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u64(0, a[0], b[0],
				(unsigned long long *)&d[0]);
			for (i = 1; i < 3; i ++) {
				cc = _subborrow_u64(cc, a[i], b[i],
					(unsigned long long *)&d[i]);
			}
		} else {
			cc = _subborrow_u64(0, a[0], b[0] << s,
				(unsigned long long *)&d[0]);
			for (i = 1; i < 3; i ++) {
				cc = _subborrow_u64(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (64 - s)),
					(unsigned long long *)&d[i]);
			}
		}
	}
}

/*
 * Apply Lagrange's algorithm on the provided scalar k. Result is
 * encoded over two 17-byte arrays c0 and c1 (signed integers,
 * absolute value is lower than 2^128, little-endian encoding).
 */
static void
reduce_basis_vartime(uint8_t *c0, uint8_t *c1, const uint8_t *k)
{
	/*
	 * Algorithm:
	 *
	 *   Init:
	 *     We first reduce k down to 255 bits at most.
	 *     u = [r, 0]
	 *     v = [k, 1]
	 *     nu = r^2
	 *     nv = k^2 + 1
	 *     sp = r*k
	 *
	 *   Loop:
	 *     - if nu < nv then:
	 *          (u, v) <- (v, u)
	 *          (nu, nv) <- (nv, nu)
	 *     - if bitlength(nv) <= 255 then:
	 *          return (v0, v1)
	 *     - s <- max(0, bitlength(sp) - bitlength(nv))
	 *     - if sp > 0 then:
	 *          u <- u - lshift(v, s)
	 *          nu <- nu + lshift(nv, 2*s) - lshift(sp, s+1)
	 *          sp <- sp - lshift(nv, s)
	 *       else:
	 *          u <- u + lshift(v, s)
	 *          nu <- nu + lshift(nv, 2*s) + lshift(sp, s+1)
	 *          sp <- sp + lshift(nv, s)
	 *
	 * TODO: group updates to u and v by recording shifts and
	 * subtractions in update factors, and applying them when the
	 * total shift count may make factors no longer fit in registers.
	 */

	/* r */
	static const uint64_t R1[] = {
		0xDCF2AC65396152C7, 0x2ACF567A912B7F03,
		0x0000000000000000, 0x4000000000000000
	};

	/* r^2 */
	static const uint64_t RR[] = {
		0xA31F34E2739216B1, 0x86A297C9835B5211,
		0x95DCE66BF04303AD, 0x8728B04D2F0F9E3C,
		0xEE7956329CB0A963, 0x1567AB3D4895BF81,
		0x0000000000000000, 0x1000000000000000
	};

	i256 ivk;
	uint64_t vk[4];
	uint64_t nu[8], nv[8], sp[8];
	uint64_t u0[3], u1[3];
	uint64_t v0[3], v1[3];
	int i;
	unsigned char cc;

	/* Decode scalar reduce it modulo r. */
	i256_decode(&ivk, k);
	modr_reduce256_partial(&ivk, &ivk, 0);
	modr_reduce256_finish(&ivk, &ivk);
	vk[0] = ivk.v0;
	vk[1] = ivk.v1;
	vk[2] = ivk.v2;
	vk[3] = ivk.v3;

	/* u <- (r, 0) */
	memcpy(u0, R1, sizeof u0);
	memset(u1, 0, sizeof u1);

	/* v <- (k, 1) */
	memcpy(v0, vk, sizeof v0);
	memset(v1, 0, sizeof v1);
	v1[0] = 1;

	/* nu <- r^2 */
	memcpy(nu, RR, sizeof RR);

	/* nv <- k^2 + 1 */
	array_mul256x256(nv, vk, vk);
	cc = _addcarry_u64(0, nv[0], 1, (unsigned long long *)&nv[0]);
	for (i = 1; i < 8; i ++) {
		cc = _addcarry_u64(cc, nv[i], 0, (unsigned long long *)&nv[i]);
	}

	/* sp <- u*v */
	array_mul256x256(sp, vk, R1);

	/* Main algorithm loop. */
	for (;;) {
		uint32_t bl_nv, bl_sp, s;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		if (cmp512lt(nu, nv)) {
			uint64_t tt;

			for (i = 0; i < 3; i ++) {
				tt = u0[i];
				u0[i] = v0[i];
				v0[i] = tt;
				tt = u1[i];
				u1[i] = v1[i];
				v1[i] = tt;
			}
			for (i = 0; i < 8; i ++) {
				tt = nu[i];
				nu[i] = nv[i];
				nv[i] = tt;
			}
		}

		/*
		 * nu is the larger of the two; if it fits on 6 words
		 * (including room for a sign bit) then we switch
		 * to the second loop, which considers only 384-bit
		 * values.
		 */
		if ((nu[6] | nu[7]) == 0 && (nu[5] >> 63) == 0) {
			break;
		}

		/*
		 * Get bit lengths; if v is small enough, return.
		 * We know that we can get ||v||^2 down to about 1.075*r.
		 * Since r is approximately 2^254, this means that we
		 * can always have nv <= 2^255. This guarantees that
		 * both v0 and v1 will at most 2^127.5 in absolute value.
		 */
		bl_nv = bitlength512(nv);
		bl_sp = bitlength512(sp);
		if (bl_nv <= 255) {
			enc64le(c0, v0[0]);
			enc64le(c0 + 8, v0[1]);
			c0[16] = (uint8_t)v0[2];
			enc64le(c1, v1[0]);
			enc64le(c1 + 8, v1[1]);
			c1[16] = (uint8_t)v1[2];
			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * Subtract or add, depending on sign of sp.
		 */
		if ((sp[7] >> 63) == 0) {
			sub_lshift_192(u0, u0, v0, s);
			sub_lshift_192(u1, u1, v1, s);
			add_lshift_512(nu, nu, nv, 2 * s);
			sub_lshift_512(nu, nu, sp, s + 1);
			sub_lshift_512(sp, sp, nv, s);
		} else {
			add_lshift_192(u0, u0, v0, s);
			add_lshift_192(u1, u1, v1, s);
			add_lshift_512(nu, nu, nv, 2 * s);
			add_lshift_512(nu, nu, sp, s + 1);
			add_lshift_512(sp, sp, nv, s);
		}
	}

	/* Second-step algorithm loop, when values have shrunk below
	   384 bits. */
	for (;;) {
		uint32_t bl_nv, bl_sp, s;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		if (cmp384lt(nu, nv)) {
			uint64_t tt;

			for (i = 0; i < 3; i ++) {
				tt = u0[i];
				u0[i] = v0[i];
				v0[i] = tt;
				tt = u1[i];
				u1[i] = v1[i];
				v1[i] = tt;
			}
			for (i = 0; i < 6; i ++) {
				tt = nu[i];
				nu[i] = nv[i];
				nv[i] = tt;
			}
		}

		/*
		 * Get bit lengths; if v is small enough, return.
		 * We know that we can get ||v||^2 down to about 1.075*r.
		 * Since r is approximately 2^254, this means that we
		 * can always have nv <= 2^255. This guarantees that
		 * both v0 and v1 will at most 2^127.5 in absolute value.
		 */
		bl_nv = bitlength384(nv);
		bl_sp = bitlength384(sp);
		if (bl_nv <= 255) {
			enc64le(c0, v0[0]);
			enc64le(c0 + 8, v0[1]);
			c0[16] = (uint8_t)v0[2];
			enc64le(c1, v1[0]);
			enc64le(c1 + 8, v1[1]);
			c1[16] = (uint8_t)v1[2];
			return;
		}

		/*
		 * Compute shift amount.
		 */
		s = bl_sp - bl_nv;
		s &= ~-(s >> 31);

		/*
		 * Subtract or add, depending on sign of sp.
		 */
		if ((sp[5] >> 63) == 0) {
			sub_lshift_192(u0, u0, v0, s);
			sub_lshift_192(u1, u1, v1, s);
			add_lshift_384(nu, nu, nv, 2 * s);
			sub_lshift_384(nu, nu, sp, s + 1);
			sub_lshift_384(sp, sp, nv, s);
		} else {
			add_lshift_192(u0, u0, v0, s);
			add_lshift_192(u1, u1, v1, s);
			add_lshift_384(nu, nu, nv, 2 * s);
			add_lshift_384(nu, nu, sp, s + 1);
			add_lshift_384(sp, sp, nv, s);
		}
	}
}
