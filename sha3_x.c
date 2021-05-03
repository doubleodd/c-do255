/*
 * SHA3 and SHAKE implementation, using an external process_block()
 * function (for assembly implementations).
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "sha3.h"

#define process_block   do255_sha3_process_block
void process_block(uint64_t *A);

/* see sha3.h */
void
shake_init(shake_context *sc, unsigned size)
{
	sc->rate = 200 - (size_t)(size >> 2);
	sc->dptr = 0;
	memset(sc->A, 0, sizeof sc->A);
}

/* see sha3.h */
void
shake_inject(shake_context *sc, const void *in, size_t len)
{
	size_t dptr, rate;
	const uint8_t *buf;

	dptr = sc->dptr;
	rate = sc->rate;
	buf = in;
	while (len > 0) {
		size_t clen, u;

		clen = rate - dptr;
		if (clen > len) {
			clen = len;
		}
		for (u = 0; u < clen; u ++) {
			size_t v;

			v = u + dptr;
			sc->A[v >> 3] ^= (uint64_t)buf[u] << ((v & 7) << 3);
		}
		dptr += clen;
		buf += clen;
		len -= clen;
		if (dptr == rate) {
			process_block(sc->A);
			dptr = 0;
		}
	}
	sc->dptr = dptr;
}

/* see sha3.h */
void
shake_flip(shake_context *sc)
{
	/*
	 * We apply padding and pre-XOR the value into the state. We
	 * set dptr to the end of the buffer, so that first call to
	 * shake_extract() will process the block.
	 */
	unsigned v;

	v = (unsigned)sc->dptr;
	sc->A[v >> 3] ^= (uint64_t)0x1F << ((v & 7) << 3);
	v = (unsigned)(sc->rate - 1);
	sc->A[v >> 3] ^= (uint64_t)0x80 << ((v & 7) << 3);
	sc->dptr = sc->rate;
}

/* see sha3.h */
void
shake_extract(shake_context *sc, void *out, size_t len)
{
	size_t dptr, rate;
	uint8_t *buf;

	dptr = sc->dptr;
	rate = sc->rate;
	buf = out;
	while (len > 0) {
		size_t clen;

		if (dptr == rate) {
			process_block(sc->A);
			dptr = 0;
		}
		clen = rate - dptr;
		if (clen > len) {
			clen = len;
		}
		len -= clen;
		while (clen -- > 0) {
			*buf ++ = (uint8_t)
				(sc->A[dptr >> 3] >> ((dptr & 7) << 3));
			dptr ++;
		}
	}
	sc->dptr = dptr;
}

/* see sha3.h */
void
sha3_init(sha3_context *sc, unsigned size)
{
	shake_init(sc, size);
}

/* see sha3.h */
void
sha3_update(sha3_context *sc, const void *in, size_t len)
{
	shake_inject(sc, in, len);
}

/* see sha3.h */
void
sha3_close(sha3_context *sc, void *out)
{
	unsigned v;
	uint8_t *buf;
	size_t u, len;

	/*
	 * Apply padding. It differs from the SHAKE padding in that
	 * we append '01', not '1111'.
	 */
	v = (unsigned)sc->dptr;
	sc->A[v >> 3] ^= (uint64_t)0x06 << ((v & 7) << 3);
	v = (unsigned)(sc->rate - 1);
	sc->A[v >> 3] ^= (uint64_t)0x80 << ((v & 7) << 3);

	/*
	 * Process the padded block.
	 */
	process_block(sc->A);

	/*
	 * Write output. Output length (in bytes) is obtained from the rate.
	 */
	buf = out;
	len = (200 - sc->rate) >> 1;
	for (u = 0; u < len; u ++) {
		buf[u] = (uint8_t)(sc->A[u >> 3] >> ((u & 7) << 3));
	}
}
