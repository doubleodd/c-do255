/*
 * This file implements Lagrange's algorithm for reduction of a lattice
 * basis in dimension 2 (optimized). It is used for signature verification
 * in do255s (curve do255e instead uses the curve endomorphism to achieve
 * a similar optimization with lower overhead). This implementation is
 * for architectures with 32-bit limbs.
 *
 * This file is not meant to be compiled by itself, but included in an outer
 * file.
 */

#include <immintrin.h>

static void
array_mul256x256(uint32_t *d, const uint32_t *a, const uint32_t *b)
{
	i256 ta, tb;
	i512 td;

	memcpy(&ta.v[0], a, 8 * sizeof(uint32_t));
	memcpy(&tb.v[0], b, 8 * sizeof(uint32_t));
	mul256x256(&td, &ta, &tb);
	memcpy(d, &td.v[0], 16 * sizeof(uint32_t));
}

/*
 * Return 1 if a < b, 0 otherwise. Input values are nonnegative.
 */
static uint32_t
cmp512lt(const uint32_t *a, const uint32_t *b)
{
	int i;

	for (i = 15; i >= 0; i --) {
		if (a[i] != b[i]) {
			return a[i] < b[i];
		}
	}
	return 0;
}

/*
 * Get the bitlength of a 512-bit signed integer (possibly negative).
 */
static int
bitlength512(const uint32_t *a)
{
	int i;
	uint32_t m;

	m = (uint32_t)(*(int32_t *)&a[15] >> 31);
	for (i = 15; i >= 0; i --) {
		uint32_t aw;

		aw = a[i] ^ m;
		if (aw != 0) {
			return 32 * i + 32 - _lzcnt_u32(aw);
		}
	}
	return 0;
}

/*
 * Add b*2^s to a (512-bit integers).
 */
static void
add_lshift_512(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 16; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 16; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 16; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 16; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (512-bit integers).
 */
static void
sub_lshift_512(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 16; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 16; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 16; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 16; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
			}
		}
	}
}

/*
 * Return 1 if a < b, 0 otherwise. Input values are nonnegative.
 */
static uint32_t
cmp384lt(const uint32_t *a, const uint32_t *b)
{
	int i;

	for (i = 11; i >= 0; i --) {
		if (a[i] != b[i]) {
			return a[i] < b[i];
		}
	}
	return 0;
}

/*
 * Get the bitlength of a 384-bit signed integer (possibly negative).
 */
static int
bitlength384(const uint32_t *a)
{
	int i;
	uint32_t m;

	m = (uint32_t)(*(int32_t *)&a[11] >> 31);
	for (i = 11; i >= 0; i --) {
		uint32_t aw;

		aw = a[i] ^ m;
		if (aw != 0) {
			return 32 * i + 32 - _lzcnt_u32(aw);
		}
	}
	return 0;
}

/*
 * Add b*2^s to a (384-bit integers).
 */
static void
add_lshift_384(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 12; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 12; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 12; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 12; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (384-bit integers).
 */
static void
sub_lshift_384(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 12; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 12; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 12; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 12; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
			}
		}
	}
}

/*
 * Add b*2^s to a (160-bit integers).
 */
static void
add_lshift_160(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		if (s >= 160) {
			return;
		}
		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _addcarry_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 5; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 5; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _addcarry_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 5; i ++) {
				cc = _addcarry_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _addcarry_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 5; i ++) {
				cc = _addcarry_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
			}
		}
	}
}

/*
 * Subtract b*2^s from a (160-bit integers).
 */
static void
sub_lshift_160(uint32_t *d, const uint32_t *a, const uint32_t *b, unsigned s)
{
	/*
	 * We assume large shifts are rare.
	 */
	if (s >= 32) {
		int i, j;
		unsigned char cc;

		if (s >= 160) {
			return;
		}
		j = s >> 5;
		s &= 31;
		memmove(d, a, j * sizeof *a);
		if (s == 0) {
			cc = _subborrow_u32(0, a[j], b[0], &d[j]);
			for (i = j + 1; i < 5; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i - j], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[j], b[0] << s, &d[j]);
			for (i = j + 1; i < 5; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i - j] << s)
						| (b[i - j - 1] >> (32 - s)),
					&d[i]);
			}
		}
	} else {
		int i;
		unsigned char cc;

		if (s == 0) {
			cc = _subborrow_u32(0, a[0], b[0], &d[0]);
			for (i = 1; i < 5; i ++) {
				cc = _subborrow_u32(cc, a[i], b[i], &d[i]);
			}
		} else {
			cc = _subborrow_u32(0, a[0], b[0] << s, &d[0]);
			for (i = 1; i < 5; i ++) {
				cc = _subborrow_u32(cc, a[i],
					(b[i] << s) | (b[i - 1] >> (32 - s)),
					&d[i]);
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
	 * We group updates to u and v by recording shifts and subtractions
	 * in update factors, and applying them when the total shift
	 * count may make factors no longer fit in registers.
	 */

	/* r */
	static const uint32_t R_1[] = {
		0x396152C7, 0xDCF2AC65, 0x912B7F03, 0x2ACF567A,
		0x00000000, 0x00000000, 0x00000000, 0x40000000
	};

	/* r^2 */
	static const uint32_t RR[] = {
		0x739216B1, 0xA31F34E2, 0x835B5211, 0x86A297C9,
		0xF04303AD, 0x95DCE66B, 0x2F0F9E3C, 0x8728B04D,
		0x9CB0A963, 0xEE795632, 0x4895BF81, 0x1567AB3D,
		0x00000000, 0x00000000, 0x00000000, 0x10000000
	};

	i256 vk;
	uint32_t nu[16], nv[16], sp[16];
	uint32_t u0[5], u1[5];
	uint32_t v0[5], v1[5];
	int i;
	unsigned char cc;

	/* Decode scalar and reduce it modulo r. */
	i256_decode(&vk, k);
	modr_reduce256_partial(&vk, &vk, 0);
	modr_reduce256_finish(&vk, &vk);

	/* u <- (r, 0) */
	memcpy(u0, R_1, sizeof u0);
	memset(u1, 0, sizeof u1);

	/* v <- (k, 1) */
	memcpy(v0, vk.v, sizeof v0);
	memset(v1, 0, sizeof v1);
	v1[0] = 1;

	/* nu <- r^2 */
	memcpy(nu, RR, sizeof RR);

	/* nv <- k^2 + 1 */
	array_mul256x256(nv, vk.v, vk.v);
	cc = _addcarry_u32(0, nv[0], 1, &nv[0]);
	for (i = 1; i < 16; i ++) {
		cc = _addcarry_u32(cc, nv[i], 0, &nv[i]);
	}

	/* sp <- u*v */
	array_mul256x256(sp, vk.v, R_1);

	/* Main algorithm loop. */
	for (;;) {
		uint32_t bl_nv, bl_sp, s;

		/*
		 * If nu < nv, then swap(u,v) and swap(nu,nv).
		 */
		if (cmp512lt(nu, nv)) {
			uint32_t tt;

			for (i = 0; i < 5; i ++) {
				tt = u0[i];
				u0[i] = v0[i];
				v0[i] = tt;
				tt = u1[i];
				u1[i] = v1[i];
				v1[i] = tt;
			}
			for (i = 0; i < 16; i ++) {
				tt = nu[i];
				nu[i] = nv[i];
				nv[i] = tt;
			}
		}

		/*
		 * nu is the larger of the two; if it fits on 12 words
		 * (including room for a sign bit) then we switch
		 * to the second loop, which considers only 384-bit
		 * values.
		 */
		if ((nu[12] | nu[13] | nu[14] | nu[15]) == 0
			&& (nu[11] >> 31) == 0)
		{
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
			for (i = 0; i < 4; i ++) {
				enc32le(c0 + 4 * i, v0[i]);
			}
			c0[16] = (uint8_t)v0[4];
			for (i = 0; i < 4; i ++) {
				enc32le(c1 + 4 * i, v1[i]);
			}
			c1[16] = (uint8_t)v1[4];
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
		if ((sp[15] >> 31) == 0) {
			sub_lshift_160(u0, u0, v0, s);
			sub_lshift_160(u1, u1, v1, s);
			add_lshift_512(nu, nu, nv, 2 * s);
			sub_lshift_512(nu, nu, sp, s + 1);
			sub_lshift_512(sp, sp, nv, s);
		} else {
			add_lshift_160(u0, u0, v0, s);
			add_lshift_160(u1, u1, v1, s);
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
			uint32_t tt;

			for (i = 0; i < 5; i ++) {
				tt = u0[i];
				u0[i] = v0[i];
				v0[i] = tt;
				tt = u1[i];
				u1[i] = v1[i];
				v1[i] = tt;
			}
			for (i = 0; i < 12; i ++) {
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
			for (i = 0; i < 4; i ++) {
				enc32le(c0 + 4 * i, v0[i]);
			}
			c0[16] = (uint8_t)v0[4];
			for (i = 0; i < 4; i ++) {
				enc32le(c1 + 4 * i, v1[i]);
			}
			c1[16] = (uint8_t)v1[4];
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
		if ((sp[11] >> 31) == 0) {
			sub_lshift_160(u0, u0, v0, s);
			sub_lshift_160(u1, u1, v1, s);
			add_lshift_384(nu, nu, nv, 2 * s);
			sub_lshift_384(nu, nu, sp, s + 1);
			sub_lshift_384(sp, sp, nv, s);
		} else {
			add_lshift_160(u0, u0, v0, s);
			add_lshift_160(u1, u1, v1, s);
			add_lshift_384(nu, nu, nv, 2 * s);
			add_lshift_384(nu, nu, sp, s + 1);
			add_lshift_384(sp, sp, nv, s);
		}
	}
}
