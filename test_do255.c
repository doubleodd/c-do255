#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "do255.h"
#include "do255_alg.h"
#include "sha3.h"

#ifndef DO_BENCH86
#if defined __i386__ || defined _M_IX86 || defined __x86_64__ || defined _M_X64
#define DO_BENCH86   1
#else
#define DO_BENCH86   0
#endif
#endif

#if DO_BENCH86
/* We use __rdtsc() for benchmark purposes. */
#include <immintrin.h>
#endif

static size_t
hextobin(void *dst, size_t max_len, const char *src)
{
	uint8_t *buf;
	size_t u;
	unsigned acc;
	int z;

	buf = dst;
	u = 0;
	acc = 0;
	z = 0;
	while (*src) {
		int c;

		c = *src ++;
		if (c >= '0' && c <= '9') {
			c -= '0';
		} else if (c >= 'A' && c <= 'F') {
			c -= 'A' - 10;
		} else if (c >= 'a' && c <= 'f') {
			c -= 'a' - 10;
		} else if (c == ' ' || c == ':') {
			continue;
		} else {
			fprintf(stderr, "Invalid hex character: %c\n", c);
			exit(EXIT_FAILURE);
		}
		if (z) {
			if (u >= max_len) {
				fprintf(stderr, "hextobin: overflow\n");
				exit(EXIT_FAILURE);
			}
			buf[u ++] = (uint8_t)((acc << 4) + (unsigned)c);
		} else {
			acc = (unsigned)c;
		}
		z = !z;
	}
	if (z) {
		fprintf(stderr, "hextobin: half final byte\n");
		exit(EXIT_FAILURE);
	}
	return u;
}

#define HEXTOBIN_LEN(dst, len, src)   do { \
		size_t p_hextobin_t = (len); \
		if (hextobin((dst), p_hextobin_t, (src)) != p_hextobin_t) { \
			fprintf(stderr, "Unexpected string length\n"); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)

#define HEXTOBIN(dst, src)   HEXTOBIN_LEN(dst, sizeof (dst), src)

static void
check_equals(const void *a1, const void *a2, size_t len, const char *msg)
{
	const uint8_t *b1, *b2;
	size_t u;

	if (memcmp(a1, a2, len) == 0) {
		return;
	}
	fprintf(stderr, "ERR: %s\n", msg);
	b1 = a1;
	b2 = a2;
	fprintf(stderr, "a1 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b1[u]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "a2 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b2[u]);
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void
scal_reduce_finish(uint8_t *d, const uint8_t *r)
{
	for (;;) {
		uint8_t tmp[32];
		unsigned cc;
		int j;

		cc = 0;
		for (j = 0; j < 32; j ++) {
			unsigned m;

			m = d[j] - r[j] - cc;
			tmp[j] = (uint8_t)m;
			cc = (m >> 8) & 1;
		}
		if (cc != 0) {
			return;
		}
		memcpy(d, tmp, 32);
	}
}

static void
scal_reduce(uint8_t *d, const uint8_t *a, size_t a_len, const uint8_t *r)
{
	uint8_t tmp[32];
	int i, j;

	memset(tmp, 0, sizeof tmp);
	for (i = (int)a_len * 8 - 1; i >= 0; i --) {
		for (j = 31; j > 0; j --) {
			tmp[j] = (tmp[j] << 1) | (tmp[j - 1] >> 7);
		}
		tmp[0] = (tmp[0] << 1) | ((a[i >> 3] >> (i & 7)) & 1);
		while (tmp[31] >= 0x80) {
			unsigned cc;

			cc = 0;
			for (j = 0; j < 32; j ++) {
				unsigned m;

				m = tmp[j] - r[j] - cc;
				tmp[j] = (uint8_t)m;
				cc = (m >> 8) & 1;
			}
		}
	}
	scal_reduce_finish(tmp, r);
	memcpy(d, tmp, 32);
}

static void
scal_add(uint8_t *d, const uint8_t *a, const uint8_t *b, const uint8_t *r)
{
	int i;
	unsigned cc;

	cc = 0;
	for (i = 0; i < 32; i ++) {
		unsigned m;

		m = a[i] + b[i] + cc;
		d[i] = (uint8_t)m;
		cc = m >> 8;
	}
	while (cc != 0) {
		cc = 0;
		for (i = 0; i < 32; i ++) {
			unsigned m;

			m = d[i] - r[i] - cc;
			d[i] = (uint8_t)m;
			cc = (m >> 8) & 1;
		}
		cc = 1 - cc;
	}
	scal_reduce_finish(d, r);
}

static void
scal_sub(uint8_t *d, const uint8_t *a, const uint8_t *b, const uint8_t *r)
{
	int i;
	unsigned cc;

	cc = 0;
	for (i = 0; i < 32; i ++) {
		unsigned m;

		m = a[i] - b[i] - cc;
		d[i] = (uint8_t)m;
		cc = (m >> 8) & 1;
	}
	while (cc != 0) {
		cc = 0;
		for (i = 0; i < 32; i ++) {
			unsigned m;

			m = d[i] + r[i] + cc;
			d[i] = (uint8_t)m;
			cc = m >> 8;
		}
		cc = 1 - cc;
	}
	scal_reduce_finish(d, r);
}

static void
scal_half(uint8_t *d, const uint8_t *a, const uint8_t *r)
{
	int i;

	memmove(d, a, 32);
	scal_reduce_finish(d, r);
	if ((d[0] & 1) != 0) {
		unsigned cc;

		cc = 0;
		for (i = 0; i < 32; i ++) {
			unsigned m;

			m = d[i] + r[i] + cc;
			d[i] = (uint8_t)m;
			cc = m >> 8;
		}
	}
	for (i = 0; i < 31; i ++) {
		d[i] = (uint8_t)((d[i] >> 1) | (d[i + 1] << 7));
	}
	d[31] = d[31] >> 1;
}

static void
scal_mul(uint8_t *d, const uint8_t *a, const uint8_t *b, const uint8_t *r)
{
	int i, j;
	uint8_t t[64];

	memset(t, 0, sizeof t);
	for (i = 0; i < 32; i ++) {
		unsigned cc;

		cc = 0;
		for (j = 0; j < 32; j ++) {
			unsigned m;

			m = a[i] * b[j] + t[i + j] + cc;
			t[i + j] = (uint8_t)m;
			cc = m >> 8;
		}
		t[i + 32] = (uint8_t)cc;
	}
	scal_reduce(d, t, 64, r);
}

static const uint8_t DO255E_R[] = {
	0x25, 0x45, 0xD8, 0x74, 0xAE, 0xC8, 0x52, 0x1F,
	0x53, 0x8C, 0x07, 0x54, 0x0F, 0x93, 0x0C, 0x9D,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0x3F
};

static const uint8_t DO255S_R[] = {
	0xC7, 0x52, 0x61, 0x39, 0x65, 0xAC, 0xF2, 0xDC,
	0x03, 0x7F, 0x2B, 0x91, 0x7A, 0x56, 0xCF, 0x2A,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40
};

static void
test_do255e_scalar(void)
{
	shake_context rng;
	int i;

	printf("Test do255e scalar: ");
	fflush(stdout);

	shake_init(&rng, 256);
	shake_inject(&rng, "test do255e_scalar", 18);
	shake_flip(&rng);

	for (i = 0; i < 1000; i ++) {
		uint8_t bb[70], a[32], b[32], c[32], d[32];
		int j;
		unsigned w;
		static const uint8_t zero[32] = { 0 };

		shake_extract(&rng, bb, sizeof bb);
		for (j = 0; j <= (int)sizeof bb; j ++) {
			do255e_scalar_reduce(a, bb, j);
			scal_reduce(b, bb, j, DO255E_R);
			check_equals(a, b, 32, "scalar_reduce");
		}
		memcpy(a, bb, sizeof a);
		memcpy(b, bb + 32, sizeof b);

		do255e_scalar_add(c, a, b);
		scal_add(d, a, b, DO255E_R);
		check_equals(c, d, 32, "scalar_add");

		do255e_scalar_sub(c, a, b);
		scal_sub(d, a, b, DO255E_R);
		check_equals(c, d, 32, "scalar_sub");

		do255e_scalar_neg(c, a);
		scal_sub(d, zero, a, DO255E_R);
		check_equals(c, d, 32, "scalar_neg");

		do255e_scalar_half(c, a);
		scal_half(d, a, DO255E_R);
		check_equals(c, d, 32, "scalar_half");

		do255e_scalar_mul(c, a, b);
		scal_mul(d, a, b, DO255E_R);
		check_equals(c, d, 32, "scalar_mul");

		memcpy(a, DO255E_R, 32);
		w = (unsigned)a[0] + ((unsigned)a[1] << 8);
		w -= 500;
		w += i;
		a[0] = (uint8_t)w;
		a[1] = (uint8_t)(w >> 8);
		if (i < 500) {
			if (!do255e_scalar_is_reduced(a)) {
				fprintf(stderr, "test_reduced (%d)\n", i);
				exit(EXIT_FAILURE);
			}
		} else {
			if (do255e_scalar_is_reduced(a)) {
				fprintf(stderr, "test_reduced (%d)\n", i);
				exit(EXIT_FAILURE);
			}
		}

		if (i % 20 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255s_scalar(void)
{
	shake_context rng;
	int i;

	printf("Test do255s scalar: ");
	fflush(stdout);

	shake_init(&rng, 256);
	shake_inject(&rng, "test do255s_scalar", 18);
	shake_flip(&rng);

	for (i = 0; i < 1000; i ++) {
		uint8_t bb[70], a[32], b[32], c[32], d[32];
		int j;
		unsigned w;
		static const uint8_t zero[32] = { 0 };

		shake_extract(&rng, bb, sizeof bb);
		for (j = 0; j <= (int)sizeof bb; j ++) {
			do255s_scalar_reduce(a, bb, j);
			scal_reduce(b, bb, j, DO255S_R);
			check_equals(a, b, 32, "scalar_reduce");
		}
		memcpy(a, bb, sizeof a);
		memcpy(b, bb + 32, sizeof b);

		do255s_scalar_add(c, a, b);
		scal_add(d, a, b, DO255S_R);
		check_equals(c, d, 32, "scalar_add");

		do255s_scalar_sub(c, a, b);
		scal_sub(d, a, b, DO255S_R);
		check_equals(c, d, 32, "scalar_sub");

		do255s_scalar_neg(c, a);
		scal_sub(d, zero, a, DO255S_R);
		check_equals(c, d, 32, "scalar_neg");

		do255s_scalar_half(c, a);
		scal_half(d, a, DO255S_R);
		check_equals(c, d, 32, "scalar_half");

		do255s_scalar_mul(c, a, b);
		scal_mul(d, a, b, DO255S_R);
		check_equals(c, d, 32, "scalar_mul");

		memcpy(a, DO255S_R, 32);
		w = (unsigned)a[0] + ((unsigned)a[1] << 8);
		w -= 500;
		w += i;
		a[0] = (uint8_t)w;
		a[1] = (uint8_t)(w >> 8);
		if (i < 500) {
			if (!do255s_scalar_is_reduced(a)) {
				fprintf(stderr, "test_reduced (%d)\n", i);
				exit(EXIT_FAILURE);
			}
		} else {
			if (do255s_scalar_is_reduced(a)) {
				fprintf(stderr, "test_reduced (%d)\n", i);
				exit(EXIT_FAILURE);
			}
		}

		if (i % 20 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

/* These values can all be decoded. */
static const char *KAT_DO255S_DECODE_OK[] = {
	"0000000000000000000000000000000000000000000000000000000000000000",
	"c591e486e9f0f87589ba486b51629c358cddf4bd8c50db5c53d5631a3f785074",
	"3280b9a23784102fb56338df91bc417efa9d04a4718a36a55552d0e55a8d6a25",
	"5357ab7306780b1b2fac2f107ff1effe3f3c39ad3c9fa4e77702b88059d52f25",
	"ed12d7281a7638d58234df32ec1c4985d1ae099c5f7ffe066df6693c5b242f04",
	"252f89df45487c2d67b066bd2bb5ecad1e8f6970a63ce456fe2ba430ea665a0d",
	"5cd741e897f34f6769f49e6d192e091b021394b7c29fcd0f92211dc473b22a23",
	"736f644265f9d4ab1d3520c004ef6f46c6ac48eb36ca216c0a3f35a4eb55d347",
	"b58cc302de47ce281bd63edfd6464ff56b0742d1b99c540b63e7427d823c277c",
	"54f9b832bc6ba95f300235f0593d38564b9041df76b27c76c4bbb170c5892404",
	"ded481afe5d5ee743be61c06f6306dc6885bb1241022d7feed4f7654eafb5c50",
	"468e63bd5efe7dd6b88855804c465c94a929c6051bb411dc5259e02ff6d3e072",
	"b3e8e5972e218178596832d80b946f89da310666fb3751160ea7fc8a6b2de163",
	"03b720f4d6b9e76485eeb8fc82c183679ee25d0d69d4511c91d742ea60391a74",
	"83ca018c52abee2354aed77f55d344dce48c9a30408e9531904999a7b28bec0c",
	"f57f1503088132b922c7888b742e82565e4826e01968c035fde90d11b456353c",
	"b32fdcf6a19cc515f6b819c431ba5d2f3b80a111ea3dc4adff749d802060d640",
	"e38a8462fb16cd79cebe97c582294e9019a801b21c28578c036d35dae7b6e11a",
	"096ef5f368322cce5a902cc9e316007352519a52c37c2ddf7e59b7f1697f0f00",
	"b89e8f053f295e983f03bddd2761faf59eba3c1f69071697a388da00b1457e3f",
	"6982b80b863d40f177f0c48402cb246ef85644354eb84003f84d9ffe4754592a",
	NULL
};

/* These values can all be decoded. */
static const char *KAT_DO255E_DECODE_OK[] = {
	"0000000000000000000000000000000000000000000000000000000000000000",
	"5863d164f563f2fe6b720174c440d41666397b4aaf6d2243c57f5ac77b6f5361",
	"81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73",
	"eb01802d010a2e2dd6af94ee08ce72f89b9a7996f197174e3aac3086d45c5c6b",
	"8b38c5596c5f334db20143d106813e8b47573452098af9dca9fb35d3a0f9c77d",
	"17f96c91779fa4e8de9efe799dc513bb3830a9295b28f7ec78d1f651eeeeb861",
	"901818436d2422c5c4b9120506922ca9209c422b3465acf4b499e416679c076a",
	"946b683865a5d8daf4c643bf378ba87175e333ea53f5c28803e7da1e8e37c87f",
	"061958de1ff3560bfc352b6b1fa4e4169ff443f540f4a3356274efe9a504ee37",
	"add42bb95400b3d2a07f6e096bc1a762578083c38c67d6a053648973e950cc62",
	"c081d24b00b5e41ae485411795241b46799d4c8e6c9c92c372e220021628423d",
	"e9e9b4b160a1c391b0abbc321bdfb81168960249c785e812910bc7e110017065",
	"b1463c47e5d9015a6df289dd3c72b2f8f92b619683561aa4b518ea2c9346a06d",
	"e66400cdde69bd0f4740c18dc734b5a9c8a73a731a9fc2a2d2bc2d82607e8a64",
	"ed31ae129c1d51524ac74b0e880f9c4e899c9d01246fd1e14e6647c6b07e4f4f",
	"402e84ac250247b758aea5bcf350342d351f4c1d136be05bbaf3b67536e87f00",
	"04f27a42e40833d894f41830a2bb69bb704e0b964a2ba6241b916bc5b2f9e442",
	"d611039fb328c085ba83e560a1f4b6fb302a2a19e27afcf79b9eded0cb91761a",
	"c370faf1b9aa1ea5768ae36bd1712d7673cad20416d5325b5135b3acde3ad25a",
	"a0e880d8fdd28fa5d90fed2c78454c7c1a7ee7254e4166a626e9ddbe8aac9d54",
	"3a9418d0fdc91071c347682c78cc113b6c3c5a8e6790817e54ad798e958e1f31",
	NULL
};

static const char *KAT_DO255S_DECODE_BAD[] = {
	/* These values cannot be decoded (w is out of range). */
	"3a92541c88f5ca5620985ec23400c1149969685ff403f9ff8fa63349a9f518a9",
	"7dcebc2ba5e16650eef4f69d4fb444786858d42576e2fcbb867938df7451e4c6",
	"f60268231297120bced5af37a3513caa0e2b5733a30cd04c3034084e5e0f5acc",
	"10f225adbaa42e20723ee481a584269c0b573603b65565ffa8f39d8c96237683",
	"02f8239691757d594850446ffb6248b889af05fb9580789fdb0d827faef2bbd9",
	"8cfcce45e6979c13b7a20129f49bf857dac14fb6dd6311237111c278747860af",
	"07321e9fa5affc3dae33da1af30ee0eb05db0fbbccee8bc990fb2cf0c44ec2d2",
	"65f3c2f7ee8c16f053c047d10b185507cea824a35f9b38ed9db54f3bd2616c83",
	"ff57953b5b8497d3266be9e27b05e974a8dec70452ddc421df115dbc0670b7a0",
	"a0741124d8c6cd1aa30528026bed28b9ee00096d57aa12048f0f7b67dd6e0c93",
	"d0745a09c13838989c83ec2ad350fdb94ea9517d20527354067ad1f43b37fef1",
	"36c2472d2a704844e64024dfcb1077fa715bf394994915143157f90337a3bdb8",
	"830271462ef4b2430012cb6b089a2d087d2b982e91ecb194d4997e643cad0487",
	"9632aa4a38a6dddd15d8f7877ba46ccfbe6b80f7254c008a973cd0f3ade8eaf5",
	"f3db83e51c406741d20710fe4aa5f04e06c0ed4794b21d6931946007c68fcbb9",
	"dc8c9ca8b34e132d1ac61694d8e457f296cd21f3d75f640ede46267f3b7c80bf",
	"922c234934c93ae72825d2167430d5939305fea8019d63c014527e7f2727a2c2",
	"ea6fbe2ec7177b20ad1b7bab4c61bdff817e9afb9dc9b7f7ce8ff4df2ae721cf",
	"4b0c025553ed615139140e7852682cd87fda491eec9ffafa5e515da9e71a6ad6",
	"dff93a41bef66c4725e63b20634a758c61a0fff9d92e5874fbed4d44cb317fc5",

	/* These values cannot be decoded (w matches no point). */
	"01e244ae27c9bf0406559c04ab7b25587f56ff0b79968ed05cee35ac7e527e31",
	"0904f4a4622446cd7f453e13968f22aaf1841fe17461e0998753e8a6c3f1db1b",
	"abf5bb06f321a68875a285c8ae8a1f10352c566a9a347c8758b385b05151c147",
	"e0102e5cb6cf5e1f4b9d027d147d7cfc161ae6441f33f55a0acc2737d3117229",
	"2d81ecbfec945c26df0d4cbe401f55fd550bea5ca5facf26f592e18f3e398957",
	"7d781c87ecd3ae16b13d72ced0bf9fad2ab00014e1c5035f95cb32bb8d295e60",
	"3cbba04f2ec5f32c2fc68183de992759d0b6600983ea5e836c190bbaddef1b5b",
	"f53d4f896d6cb6c0ac1f36939df609419c90c27b0533cee21d0893173be4a111",
	"b2ec464f0ce61e3015e68e2f8c8ca55da784b605694a740391f0b6d405756b63",
	"e8e297612bb58431e08168f6b27739461cc5d8e5b8e02edc7dcf0209440bdb62",
	"28c94ceb6c5b83b802cfbd4ae3e04455061269a78b09299c2767c90ed2aea714",
	"5dc5a2ee12aaaf093fa2c1fa845fa45431ad630760424ce24c4e15951fb69645",
	"e1434027d7502a74417c663183df9a570746c873e02d0ec2cbf7a33d5024fa33",
	"4120309222b712ffe3390b263c018eb3957154337578affef09b8dc918c57e14",
	"f045fcd44e2c4582bc929cac99b8fe62013a48b7a15a8c7102b283eb18ec8068",
	"548982ec9b6a9d396bf36c3d7df1fb24975cdb750b7c316f218bbc388aebc732",
	"3309ba8cad6523b7984f19342e5f63b83b83964a798a1388ecbf2673aa15fb40",
	"60db10a949d2e9ed81a4fe4b6b0d6e4454ac420f22e7dfcaf4503aefd7113542",
	"b93e2061042a1e65640bc8632e880a565e459ee4ac4755705cf9f0769cf58c57",
	"80493859a49bf1ebc8a3e7195157b2b98a0ed135d77959ef8a94a64b709f7628",

	NULL
};

static const char *KAT_DO255E_DECODE_BAD[] = {
	/* These values cannot be decoded (w is out of range). */
	"08761e618ba81e1cb145e406e8c6d63b115d2ca460a585a85275ffed14ec1ea9",
	"6743c991104ac09dd95f81607645aa69892ac9c79e186fb666c201b465f40e94",
	"435f64012646f31d93edc5f20a4bd03c51b30cd314f18ac23ad74779e2190eb1",
	"59a13eaefc3e1043f748adef252abf956ca58287be56a2a8990fbe77c6790689",
	"d9001e28d9468785a32fa452a99a5fa617e3cc25dac80c7b4dc668883ee19eef",
	"4b80e321fc3b29f3dee7ee0450eb309bbc3f628e0b0e799e1161a3c0f0ff8aa6",
	"632d62209619b3a388d00001fcb25d627f830517814267f3cff381c6d83b1ef2",
	"7ad205df3300a6350b03c8337c89550564fd8c91779e9d0342dbed07346e3a80",
	"a9747cda67e5ad2bfe85253f094c105856c4ca94de7f58ec1ca64353eb4eafa1",
	"0049a4befd5678f03d3b23a986dc3350fee4cf41140cbaf78db10919ecd288a9",
	"45c00fec5f4ee698f8fb70fa6fbb102bd2707878cdc22dbec2a073d7670752ad",
	"f32961c81d13acd47d56e10abd4cf21086ee9ba34cea03384bcdd19b8336cc8d",
	"120106c537b05cd7ec4a24b83e99c399f212776c0e3f43ed923a7cf8370316b1",
	"7160d3059ab1eef1cb3f33d4296e0ca05594f39f59b11e32e1e61ba478ce05ed",
	"d700eecb2c254be73f7d19aae82c883fea2362225769b4cdeba60f63cccc7aed",
	"22c91e4ae9566785b1fa88bcf6d57011300c63af5f4c4f16a6366d57db0973e1",
	"b5ec84834affc9d9026d2051da7b4584fda756e39f1b69df8a36a8b35902e3ee",
	"bdf310230866f6e8ed2497956f148e01f7d712ba46255de4e6aa2eaec6d25adb",
	"0b03203dd0d310b4f203ff7216bb0b666e91c2680ed06aef96c78f7a4e72a5f2",
	"b7dadcce34e857c65acf1e2c3f36376aeb906cf48f37559582d1b6067ea0d5a8",

	/* These values cannot be decoded (w matches no point). */
	"7d901cb1f92cd676a64241f368a415e48c6a8a662566b25a309d8f1823eab363",
	"e4f8ba84d29399e935a8f8dfee14f1ffb5b5bc20c4889ab391c0de6464a5ca46",
	"e3603e6b73b33cb33c9cea3119d0e223520ca31dee0c9238a10f7036b0a62f5a",
	"d0c7fb74ee0c7315b689cf84efff6c4407d86a1620138840512090c26ad58972",
	"630a1f53aa31c5809408d113a3537d2d63815c221378469b309ac4544cd33877",
	"bda70de7c625827b184dc76747e1f52f218259720f8790082c409e66fc771a2e",
	"25128aecc25469119568d04c274e1cb523bf5007fce1487957db96e14954f11d",
	"20c393a3b0f36c5461a5edc04b574105d2d9c8ec615626810b838cc2b77f9a43",
	"808fab937bda301503242d8d6f8da0809b1e9c9d958b46f1c9b2cde0ba365d29",
	"bcc21cf263eb2840bb6e7a974c91ecbfdb6be13f9cf37c749652086abf4c7526",
	"be0e14ae5cec6ec93c83274c83f84dfd9296243da9209151f5e020af6b29103b",
	"ad82e5d3f63b9507e8cb5faa79052259da755a7f1437fc7cb07a5398330b6147",
	"e7fa1c035f07d020e6b91d6b65f08f6385e3b39cd449c7f67729b1e058103038",
	"f855695bf1c0fc5618437e1c10d30f40ed0b9b823d9424d291dc3bf3f0203953",
	"1d41b6b8376bbd4316623b800fa253d674feafe3b2a51bc709b7b3619fd5f025",
	"af249c269234142d7d87ac6a23e2c65ae26d18e5493ae15e3469ee46c73b9125",
	"3044fa93f5d5a1acbb4bbdd050aabc1b15d6c95fe4c0ab74cf7cb0cb8e503d3b",
	"b88a6ffe9e90bd7546a22ea54bef4cde2e5e6cb36bd004fef175bc4567ce3120",
	"e28e882d49e4bd2dc3e3301ffc2e8e234810d067ed7783e3f84515985a65f812",
	"f326a01c98b91363347b6e6816366be025f1437a680b99993f28c0e41f01b323",

	NULL
};

static void
test_do255s_decode(void)
{
	const char *const *s;

	static const uint8_t zero[32] = { 0 };

	printf("Test do255s decode/encode: ");
	fflush(stdout);

	s = KAT_DO255S_DECODE_OK;
	while (*s != NULL) {
		uint8_t e1[32], e2[32];
		do255s_point P;
		do255s_public_key pk;

		HEXTOBIN(e1, *s ++);
		if (!do255s_decode(&P, e1) || !do255s_decode(NULL, e1)) {
			fprintf(stderr, "Decoding failed\n");
			exit(EXIT_FAILURE);
		}
		memcpy(pk.b, e1, 32);
		if (memcmp(pk.b, zero, 32) == 0) {
			if (do255s_check_public(&pk)) {
				fprintf(stderr,
					"All-zeros pubkey not rejected\n");
				exit(EXIT_FAILURE);
			}
		} else {
			if (!do255s_check_public(&pk)) {
				fprintf(stderr,
					"Valid pubkey was rejected\n");
				exit(EXIT_FAILURE);
			}
		}
		do255s_encode(e2, &P);
		check_equals(e1, e2, 32, "decode/encode");
		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	s = KAT_DO255S_DECODE_BAD;
	while (*s != NULL) {
		uint8_t e1[32];
		do255s_point P;

		HEXTOBIN(e1, *s ++);
		if (do255s_decode(&P, e1) || do255s_decode(NULL, e1)) {
			fprintf(stderr, "Decoding should have failed\n");
			exit(EXIT_FAILURE);
		}
		if (!do255s_is_neutral(&P)) {
			fprintf(stderr,
				"Failed decoding does not yield neutral\n");
			exit(EXIT_FAILURE);
		}
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_decode(void)
{
	const char *const *s;

	static const uint8_t zero[32] = { 0 };

	printf("Test do255e decode/encode: ");
	fflush(stdout);

	s = KAT_DO255E_DECODE_OK;
	while (*s != NULL) {
		uint8_t e1[32], e2[32];
		do255e_point P;
		do255e_public_key pk;

		HEXTOBIN(e1, *s ++);
		if (!do255e_decode(&P, e1) || !do255e_decode(NULL, e1)) {
			fprintf(stderr, "Decoding failed\n");
			exit(EXIT_FAILURE);
		}
		memcpy(pk.b, e1, 32);
		if (memcmp(pk.b, zero, 32) == 0) {
			if (do255e_check_public(&pk)) {
				fprintf(stderr,
					"All-zeros pubkey not rejected\n");
				exit(EXIT_FAILURE);
			}
		} else {
			if (!do255e_check_public(&pk)) {
				fprintf(stderr,
					"Valid pubkey was rejected\n");
				exit(EXIT_FAILURE);
			}
		}
		do255e_encode(e2, &P);
		check_equals(e1, e2, 32, "decode/encode");
		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	s = KAT_DO255E_DECODE_BAD;
	while (*s != NULL) {
		uint8_t e1[32];
		do255e_point P;

		HEXTOBIN(e1, *s ++);
		if (do255e_decode(&P, e1) || do255e_decode(NULL, e1)) {
			fprintf(stderr, "Decoding should have failed\n");
			exit(EXIT_FAILURE);
		}
		if (!do255e_is_neutral(&P)) {
			fprintf(stderr,
				"Failed decoding does not yield neutral\n");
			exit(EXIT_FAILURE);
		}
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * Map-to-curve test vectors for Do255s
 * Each group of two values is: input bytes, mapped point
 */
static const char *const KAT_DO255S_POINT_MAP[] = {
	"0100000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"ef8ddec9c3906319fb74f84036829ba34a1f51f7700b78c81bc91389b8adcb29",
	"8c09ce78bc83500190e7059310f5e8637170769d844233bc557e5e2497dcf820",

	"09aa856c6506ed9848570423dc8fbf4e916fe7f36ee895109893895a73415f76",
	"336f8f5872a29e0defd157585246037e0e349e6ef3bd1c990efdeda52d9d181a",

	"2615809de9b07f3c5cf0170ceccdab6843798f628743daab41b235ebd93f0916",
	"f801553749eaab1c17ab414754bdf0214ae534a6254cc3f74767b5c86732b472",

	"0cc3f41eef7051a5bd1b32ff447618375a6efef1adebb9331abf4bd7b172d672",
	"17a43022c6b25554d98355bb937ec9d3179c999c2082c96da20d9cdd44666b00",

	"d706eec3dc6f112f1fa8ee531dc6af144c03a7ef261c117f7d3ead67f4f6a835",
	"702cf2796ae47b5b0892c8eac02567c28fab77f50ae42a01167899aa08ac0646",

	"927fadf4f0357ee03d8f52503692c404c9b3ac05443990619852098ed1ffe741",
	"a6815504de2232cc30e4f68b7329e5824339d47429fcd9ac374c9f312a239108",

	"3ed2a27bb19a5d0e4015e310ea2c7d5be9363d2e3a699c50ff1073d89d984b3e",
	"d30a3eadf2a120cb196da94463c216e0a5b5c239892c4c21caa1995a2c7cb86e",

	"184cffce57049542aab9a816dc806b4419a6e98a538685c627bc1cb86acbebad",
	"74ea8109db1b8c0829bd1fcae8ccdbae2df1683909992c2bbeff1b9338c2c47a",

	"299039d72adb0f9afb61e59c3d7f0fdb102cf152e3cac224922ca5be9087b12e",
	"350b451b72b66f4d0995d468cc51716a36aaa45f4280c6e378de3d277e2e5422",

	"27411f5d3c65de7d6c3a34bfe03985c3b9ab86d67a33af1851f025c4dc2e45c8",
	"8272b971dbdad55b02c5ba6bbffd02c41a9979a7a6373514c10d4a861a81404e",

	"27caf691b8003ef948ffbbf0520bf129d47024e0576ac3f28e73911b75afed33",
	"c88444c61ca57d298a2c40e59a26d2fe7261e12ea8e04a2b19e5f10b5c61a26f",

	"aa09eba29f536707455c34328d74b417f5100f3cb499cb360fc7cc368153d4a1",
	"7a70f2f0076852e9eb0331333021ab2d41eab2dd72ec2e143f25a77796a4f647",

	"8c8f7c49e620854b260b5f5ad5e3d221580b24aaac46377b13024f7d3198d6f2",
	"f406e61169965432bdc314b8dfc35e72ad886c76c72ccdcca8abd227709c2c35",

	"6fddae9d444eaf137aa2831b1713c93dbade201ea2a811a608ae0cfbcacffb6a",
	"a829e7440a42b80665eea59eb1cf5480e348be004a8b4b4b07e493122b0e5513",

	"0193677a09da4f1d819ef8418d15e0cae8eb9a89aafb96f13a262a3e6d950f11",
	"dd027f5ce90e7afa06c890fad453f76005e9cdf827342e91ca7042b3e722120b",

	"1d947f967ddc8fe741eb3533f7f25f78c541f53d30fcb7d6a5b5d8f269020030",
	"8c12d16626cde555cb922dccaa8ff8feed931df1bcce5f4f85a240ae2dd26821",

	"32049bcf13de6ac721a02a4bdc0737abc4c3c0f28b9a6123ccb33a9e7647d136",
	"70b18ef2e9c0a33a8c6a79f6823e405cf8451c293b5999eb7a7c1a67cdfba460",

	"69cf1efafc674d1a41ba28a16b58759bb42d30bc77e8ec5abe19eae43f150085",
	"1234cd73cffb6350e8b816d349d79064ae1d5c01869a08d11384070b53027f75",

	"3453d8c0de0389d588ae9791b6105f17c7407fabd3e2a665194e591600392a64",
	"4fe4ff767a01552863bc22363b2431f19c3321e2bb0323a2c5bca6dc94be560d",

	"530195eb8bf148eecdd97ecebf7580ebc685a2a7a1a7aa54bfa27cf87a203ec4",
	"c26b2c32991126301254bed67ca20674dd8f42151a7b80898eab4b2a03aa712b",

	"5f915e83780fdf6cc484b2d890c2a5f053f2eb845722882c074a271ce3bde224",
	"f81dac9fa4565539ab631a282c86d2cd49b9f18fd18357cf6ec77c93330c6a7b",

	"c412fee4f10c9b8950077e31438f9bfdc644a28cc9f3cfd83ef480a82569a41f",
	"178e85745e49482e6ca1fa345ec910f1e2f8724ea03c5093fa47b6491d54c174",

	"3121bd886339c8004216c982f690b50f4e633d60c9988b669d6d2c1fccab4ddf",
	"30e95dca620183bb105b09b13b8f9e2ed3ad5856824975e9728019e18d57fc1b",

	"14b087ed68b14798d2d1bc4ea65ff28305862b94ca6601d9bb8f807043974e67",
	"38f3860c889afc82974d2aa0bf33afdba9bc4b40cbfe59e331ee94e4dba9167e",

	"b26fcbc480eb4d531fcdc50bed5d292bfca679ffc358655e0135b169a351ddd6",
	"92f06ee33ee6c62b24b8546b9e235bf21dc924e9e2709ce1e912271c3156a237",

	"39f63d9002b06a64fbd6e22e541f6374e6731fc831c1f1c2e46c5391c959b8c6",
	"32c7847e59e8e0c39d5bd9efe2cd656b2d7654c733316dc7e5fa4dd1fccbd314",

	"7b7c740f9fffb909e45fb3d5ef8b705f009590f7cd55a96047a30407f2553e5e",
	"bb9eef67a061a575c3b71a3e29993eb8c523fe9ae11a707a30c048f0fe2d8007",

	"0932c9c6acde82775d0bb254731bb096833362990f68528c9b14f9d076e93c02",
	"906894b1b036386b4ddf2a8a58d45f308c9321109fc12a34341e50f0e80df216",

	"ba4ad4152fd792db8ea1575c87494d23036ee200ceef8d8c1bc99f6df38ab514",
	"dd9ab7a2c8778a7cd89769a71246551ed70f7606fdecd717565184740c1d4f03",

	"5896d16b48209ca651ba8f686688d091072805b245f136871b43bb8167fca960",
	"ed209723b4d1d06ffec847ae89708474dc4fc133fe793614347fab61e67b2c74",

	"ce5fdf706ea9c16e47d4cb9fab4ba6173858acfefec7096669aab4c40b928346",
	"39fbb562fefcbeaa069d6ba0c975e99bc36ea3f017c24b3bcf36fa6ab84e0c17",

	"a5f64a58f56dfb212f63cbde27a246e96f05dcc60687a8eb51bb84200e3be5e9",
	"ff8094c66be12b8008b3360fd7ca828798f497797f879c4b72fb1575e4b85b69",

	"70b2451d5465130a6cf7998d3b752dc8ccd7547598149435e77573a86bff1ea5",
	"9969658848ac05f24f3387ef5f6ffb7e4aa3339171537143a23984205341cb3a",

	"b87e45cd4369dba7cd7b877a20325b027ba87c90d23939844c0c9a194ee39d3a",
	"2e6b5b348687087b30c13bb9ed59fd9e324490ea09b66d1a55a7601a8fb25b77",

	"163eec28be72a61876479a7e5838e438c5d053ac09e4a9b6c4f5e91ece2a9244",
	"06aa16b19bf14ef4607f31f7a782ec38da885f912b24fa19ffc834a16861d805",

	"9e3b7632ae68362ed9c716bd1e52c040315531a4f75c67ce3e3c8253d36c98e5",
	"f42ec2b036a422009273ca5d4edaa6e1a910cf288f6e1435856769feddc99564",

	"f1ec32a37b6a82cbe322e049687af65ab32383e5b4ff1a40b40876c4c9876ad8",
	"ce2661546a9ffe2c7bcfb6782290870b68544b2624a52031c179be82d530497b",

	"8f6b4c01d3090814d79008367f0aec97e27a6aecf751eb2cd678b744018649ab",
	"46bf5c2ad04470fd09f209af2f65f1c8063bf58ae55fbb06758d8c150a78c105",

	"395c4cd8ca231c8a20169c456908ce8b88905e740f389837b33e693ce6d9e1f1",
	"227fcab45ba482a85ac0b4f30b0e1c4fbca733aae08022a017a814cb20601c2d",

	NULL
};

/*
 * Hash-to-curve tests for Do255s
 * For i = 0..99, hash-to-curve using as data the first i bytes
 * of the sequence 00 01 02 03 .. 62  (raw data, no hash function)
 */
static const char *const KAT_DO255S_HASH_TO_CURVE[] = {
	"8c35143563bb29d47d0d6ec499f7a87c7d10d17884d238b36c853bb8c185fe60",
	"59576551095a30f97af32498b5a0e6c186d8ed6fac7b3e6def5065ad3d26646c",
	"6724b037c9c1c71c486ed00623188c3bec884d7f97c0c7d0c54332b173610d20",
	"1a7a33303a6a05050132197e22f4fedfa42c06ed337f83c06fc9218dc0039759",
	"e2826a696c1c9fdebbf40069cd33f1c5dc9537cc01d89286dba57b3b70533f06",
	"23269765cb7facd0de7c614984285154b6a01415ba6d98cc66a66ccfeb730a61",
	"af038217e410f5b9966fead7d3dbaf1884adc303d97563b1ee9716aa465f8901",
	"07a3041a40f33fd15cf5f88a56c54ff4ba8eb121ae8b1571e39a0c64c3fdc243",
	"65940de3733e4c9aecac4e108e0e45cbe94fded31acbcf167539079b5423b407",
	"60767c7eaaba817450aa0c414de0cab0d44bb93279743f571f801ef22d2a191c",
	"9ea34c6ab1af1ccf4bc9fe4a04f7866fadbd31c639faa6495bd451449d022d23",
	"ae727ca185765c1c4a360a4f1b83f5ab320bea3b535b32df79d59f080954ae72",
	"964f944c3464bdca56ea87ec98e6342df8b664a0d30ef5633b10bdf8cfbe363f",
	"d20863ea8cae7d03ac9fe747a8c041a7fe92c04311fead520677278190d79104",
	"e5e36632efb9872cf81bf6fae52068c813b1137f87979bae7b81092fbc952e6f",
	"c124f7e1cc16d942c9187a202016d1a27d957e9e0df9cd93cfc1df49ab811f13",
	"d911125a17a764a5c7011612f3b932215b52c9cd4ce34d3bf92ffd91c6296f41",
	"8e0bdaeb8f8604a1650fcf659eb3a387b61fb6e44138632fe8ad875c7d2e0d4b",
	"a0fc0dc322d579ba389f3a75d3ad7fd59e04304cd4126ad55c8e0c3dceba7f74",
	"7a012a1827936b5a710e3557b50ebebfaaf1e51198a4a144e611ed62cb0e3476",
	"345fc77f3816626db1796a39a2c82757d01fb930cf34f4a336a4cbd50ee3e706",
	"250dd149cc80248f890ca3a73cb98038faf71d1a46cfe5dfdc67866256305a40",
	"86a4760aaaf6f65f3eb4faf7ad8e8a5a29a50779b41c70f19afc1f2c56726272",
	"0d3e23a30aa6b37cc76badcb68e15c603c43a934be686779735d12014af2274f",
	"ee9f9d815ffa0ba67f2fc5a9e13b64b9a58fd5634c86f57936a26b8925915e29",
	"af948113d4226fafa1eb4e3cefbcc47c919904325f68ba2d20e1d1b7fd1e7015",
	"215bb4bcbfb8827aa7416c790ac5965c7b22ca5428a276b414a72b1d561fb80d",
	"b38f1f8f5de6b71732853f786e966486dbb0c5b3f6923f55b4b2ed390c56c518",
	"23794c432a0f61e7b6022ee213ecc72f2fa2f1851e0209e7b85250be111a5d3d",
	"460213c323a4c8f68f826c45e7d80bd26f467d3c35a3446fde031f46cb206216",
	"38e1488c19b78c4c8756bcf4eed42f8a887a7a166d432d1a3264c89f41644e01",
	"d989985ad0c6ba22ceaf73efee43a3585034fa82fe89c66531646b58c36cfc55",
	"5dc6676c53893b75172525e9a87d91baabd790c83ea9a9e1a1c7da67a4bffe29",
	"b388676f85d778ea7faa86d6b73d31d21d4f18bb8be31066b9fa9aa4b7818d72",
	"d14718fc239ed0eea42deb5e68be1ff253d375ab743e91e6076e0a7198355753",
	"c0b9115fff3c70266ea4f999c3cd22ce8789f777f2bf228d9797e2c384501250",
	"4beec4d5a7504f4dfbc96094e9b339276cc40343b7623d8799fef4f213aa3b27",
	"54946aff70e5b27dace1487098df4805bd00804322437fded07698f161c65b6d",
	"63a8d805f7006abcb1c92b7a21578c08ea2855e66324c3cba638d652c44d2e55",
	"32151b640ffeaa0411c2b9ec13174277abc0aaa38115b668bf4a75b932fc5e30",
	"5c57700824d40b0777645528b9a30beb3e483359917a2a659a5c77751d70325a",
	"22841f145dd9275bc0c0b05beb4d3d865412983db93b98e5895619c734002465",
	"d79cd79635edafa2a735e3c005efe6f767ec7244637c264db99e4103d364de64",
	"9231bf2482ced363717ff81b3d3c2e388bd43c7a6cc3adc305a132ae405e2c0a",
	"87bb320a6a669d77455e17c103883fee41a070cd4781e9e2ce22d11a2aed885e",
	"b8f329f2c7bf551948dfdfb8f3990da0be285694e42ccef73bcd6772f55a2803",
	"58a457c3610d720bfa80ff401bbf63bbf1e619196bcbf880280b85763401d479",
	"547ed55a6355a7f8563f11d4217d325d9b37d4380789484563f219d43a6af17b",
	"88dfc63b8fb76966fa0bbf0515cc90b4c9b82a60bc6fd5c4fd55623d2f81b17a",
	"412e9dfd55bda895af4f39da96407f06be635c17637e32b1e0d4fcebda0d255a",
	"3758ac7e0387ed2ca95b942d8ecf3ef6f36cadb2cd48c8cd6f5789cb91e7e72b",
	"d7aa8e7f91639e07614b1b7363e1ded71891f5ea741cd2f43aab0afa2260b65a",
	"37eaa8c86f93d20f95c45187b9606befe0567bc2ad4b3f568b2036d5165a8b6c",
	"d63c099ccb5ed3dcd7747da52f4aae3f75fa868327768483cd6d730a8ae2283e",
	"0a4150f14f2d042491c1df67b1e16b65a681d5fcce69049cd756cffff3698f18",
	"df25f9f2498b180f6e826fc6628112430cd6895d533edac1e58f0300b5ce9e77",
	"f04125d1e3f3f9d3a042cee2a8b58086f76f2ba2f262cfe3210d7f39ddf4483f",
	"0232b4fe7aaeffafe8e806f346d3d930bb3cfa5d43d5d85c1ec01eadc4f9e33d",
	"17c8657d17021ce73427b981a7ac83a149edae0a0c360f31ee4359ae9af7dc5d",
	"15bdf2e72c95d9d99eb6136c918c957c4373cf110633e06a4f8b1e270c0fad00",
	"493a5f2fb1f9e55b6fe76cfc7d0b6b2e76baa9ddcb3ae73e3a03b12b049eb83e",
	"79dd571998c4801ae39941a431e1243f9c02bb530291bb7ff0cabcc3988a6733",
	"02860af5337512aac3f99d8174dd226618d0635c77f1d0c227db4c3a8f08386c",
	"dab3b53ea1486faa3b9881afde91751d50708f8fce2751c3aca51c6e9ed13b44",
	"098e257d4575448f3ab9c334f567e5f8ab9bd9423ff548ad55d53cc78b9e6e46",
	"33868cd10db59c109b30d1f93a7a92b0abeb472b9cf1e8f859b04d849697e728",
	"03da6fe617e7071b43604887c685b2fdcca08dd456a6ca12acb31173702bb643",
	"b4f55d9e129835fab1d9dfd98e356756a5bcae6e978e828024c31ad38d92aa76",
	"8d417184564ed3f60070dfedbacd3ba62664546d53a7b575130b8951a0be5371",
	"dfb363b9c8f034831d0de31e1e2f5ce0a2867d7fde8076a76e1dba34b763053f",
	"fdb60f76c6601507395092af8f9f81ed0ec7e8d0383ba8778dec726b4fea4770",
	"95cf00eb6599914728598671b5c4a3ec247ce0eb0faa212311751e7f27e1aa42",
	"998efe9e82d9cb862324e6c97fd9c7176c4cf09510f78a50062491ffed7a6742",
	"13769fc61477769d5418014a8065173b3f89e804c4f0ca87b71ddb98c605ab1a",
	"0c28957cde7c0219f51cda1872e9e75243157723ea9a261ad6a25988ca80177d",
	"2608b260a1fda0af848408f49ada63415972d771038e532ff7811eb1913f3e4b",
	"dedad45da7a1421f5b1c3dd01bb81906d0a735233cc50d5622121f7f7d962e4c",
	"4b93bb568dd4e2fc1074054c67b82e6b2b67e5986c9b8fbc7b16536dddb60c72",
	"abc4afe954f3f17d9544465f86d02e17b8744ea2e8d84a6a1e02c1559d536a1d",
	"dce51bfc4a032d8849f91ab01fd25722f3cb11e1560fc4dbae6b5fc20ff24670",
	"dbb5b802adb973f25feaceb6ea5a7471ba14389acc6890d821a6e008eba83a20",
	"b17c1d1dc897189e03d2c54aa10acaa76eb8298d9b86b81224fdf969a571cf11",
	"b6b628500dce4a23fabe9a65e4a11d749c1f42f86cc5c9675ba036e9a831ba65",
	"71e5fa4d352e389a3ce2953c4ecd9b718af81150baa5c47c15c0c3164ff7b472",
	"6aae4cf5e9a91f0c78bd00b4a7afc5b50eb61ace9989e3a526a9b4370fc7442f",
	"6af3ae261d3d779c8ea9c92c73f0e51c7d80b8403cf27f506311f4a2cded876f",
	"d5c4a4e295622e36751c553a6e9ca78b56a9a672a525bc03875412af7f9c1658",
	"3c918a93c63911f754c4c553ba92086a8a4031968c80c41d2b23784b8c7a9673",
	"cb4372ffab23a8a130dc7ec042a47888c2642db74e31c4484f332005fb56fc4d",
	"67632c77e3fb561bb377c87123d30a3e67cb64382bef2dea898c3af419237250",
	"e692c9bf8381a9f14b7706c7b54fdf7199b2f92a2f7b79123cfc219643b75b17",
	"4577bb695fdef6ef95670a0b73981a9b34191a34d0817bf36acb9ab2cc675400",
	"41bcdc50b33fa715e357b49b1f8202008cf47ba9634c3089bae3602027bb9f77",
	"c52a328eb1512390c29342530ff8e26b00a2dc57b41b5b6a76837eb7effbd61b",
	"761d6ad42448abf97429f7f06ebf3752960a45bad8e63893eaca8d635ea7f16a",
	"014c4564543df3137cc3d6b4915f6694427b3ddb736d450950d82bb3dd224c1c",
	"4eae98d914d1407a9463f76bc44acb17f3debe57875d692a17ce87c547311260",
	"c9bf5cb469ba38e1cdedf51dd7b78ae95a387560a06b6cbd8649c54c96dd5868",
	"91f20a822a430faab66529e3ff7fa9b05b916950bfbc79f209fbfafdde1ec22e",
	"d8deb1cb4e5597e6f28c279e0db145c2d9241d54ea80795dfe6e5ee50d0dc429",
	NULL
};

static void
test_do255s_map_to_curve(void)
{
	const char *const *s;
	uint8_t data[100];
	int i;

	printf("Test do255s map-to-curve: ");
	fflush(stdout);

	s = KAT_DO255S_POINT_MAP;
	while (*s != NULL) {
		uint8_t bb[32], rr[32], ee[32];
		do255s_point P;

		HEXTOBIN(bb, *s ++);
		HEXTOBIN(rr, *s ++);
		do255s_map_to_curve(&P, bb, sizeof bb);
		do255s_encode(ee, &P);
		check_equals(rr, ee, 32, "map-to-curve");

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	for (i = 0; i < 100; i ++) {
		data[i] = (uint8_t)i;
	}
	for (i = 0; i < 100; i ++) {
		uint8_t rr[32], ee[32], tmp[64];
		do255s_point P1, P2;
		shake_context sc;

		HEXTOBIN(rr, KAT_DO255S_HASH_TO_CURVE[i]);
		shake_init(&sc, 256);
		shake_inject(&sc, "do255s-hash-to-curve::", 22);
		shake_inject(&sc, data, i);
		shake_flip(&sc);
		shake_extract(&sc, tmp, 64);
		do255s_map_to_curve(&P1, tmp, 32);
		do255s_map_to_curve(&P2, tmp + 32, 32);
		do255s_add(&P1, &P1, &P2);
		do255s_encode(ee, &P1);
		check_equals(rr, ee, 32, "hash-to-curve");

		if (i % 10 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * Map-to-curve test vectors for Do255e
 * Each group of two values is: input bytes, mapped point
 */
static const char *const KAT_DO255E_POINT_MAP[] = {
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"7c77393af403852b621bbae660adeb83939dd0a5f50a986438f1cee5aa076f2b",
	"f22b8c41ef41d627ccb0315b63a344b4bbd0c468126d48c5fb464c2a10234e25",

	"2a24450b202dbd95efdf0926740f914b548057b52840dd6954d3361c135ababe",
	"28a3522313142c44e6495452ece93e213b4a794f9238f9bcc80389dead84d955",

	"82a075600b10b133622709204280900fc0386d666676fee3708b129838c6bed8",
	"bdcbf7989c2d5e9f266c766bdeb8290a0e0d0f23faa5d21906d290f541ca8847",

	"e21baddbe4292c21084e50113560b5b6192a34c5a349d55eafc368fd62ff904c",
	"dd3cf979869404d984555aa660a60cc9fbea3c8c7efb6bd67bc366f63acedd38",

	"e6dd4f7067fb4c0cd3a4dcd4e6c48ee0bbfdb4745fbdd548cd16405bac258181",
	"07a37c1c788e206e713df4318daaee2997427f3016849fd09319ccb2584b1117",

	"3a396c8cab2a781765b7ee28485320ee8e662ffad75b6e45ae16eca1cfc9d62b",
	"d64c989d65d58f1d5f3a0d214305246fe91a2573877125a01ff92c2bc465d76f",

	"296e39ee957395e29a2053c7d7222555030500f88b9135cd7d1b48a95e654a93",
	"4927392fa7d750f00559acec28f480ac4937b7d881c59f63170292435cc6ed69",

	"c5562800cc863d83f735c7ac9fbe1c36f2e69a3b5af9ffb8dc4a7e094a771455",
	"5597a65c7391679627df68674303b900e0625a499a822251f271f1a51ee80e5d",

	"b3f09c45abf1c3b8e1198a7572df6dfb4f62df6cdba0f76ec55b225e0fe77ce5",
	"5bd3b7683a0cbbf0b9f5dba9893d74bff57db1ec1ce56071cb14d00d60e0cd58",

	"0a19a23eb556a74203b62a06ed663afacf691ee13c62b384148e56e6c2b2b37b",
	"0e630f3953e0132cfd8db0decdaaa33060948d8e3b36ada436b5078b451ed361",

	"c28a168b21c3a1107533c998cecd45597b03c2be70554013521918256f5baf15",
	"9f780a60fe7435024e455049ca6c4e0154d4848e6172c3fbf1396249b2cc206e",

	"552ab32fc9e4d4226334a31b08c05e5e5c0ddaff8b5f291bc1762a9f8d920ffc",
	"be5de20456aa8904c98ae2dfcab82fd47ab492c59aa3a9af26a9bbf52683205d",

	"cab9ab03b89c099894e759a01fc0faececae6eaa7a0f67cc166bb5a4c422634f",
	"1896d0514ed6f9795175e22b34d19248b117dc684cc3ffeed96dfb28eaeb633d",

	"471e3c38e2550fddc47246ce8583d0a6a2316c7aa207861a357f6e8995c5d3bf",
	"8f6085b268cef2d676ee879f608abd613eb593716d9b82f6e659f95c443cf013",

	"b4d07e65ae3a7b808e619350b7f39a118cf779cd2fdd691668d6cc0d6684cd01",
	"f3ac72e9b70cfb74f1f53d90abe98123fbc8cf19af9d1c04eb4906455f5b3c26",

	"c081b2e772e34361a99001ce607c75271a615588da226a0095ef3f6fca15ea06",
	"7ae24441ae9477d82e4d8bbfd815c120d6def83ab992f6a13e30b47b59587377",

	"e5b8e5d2aa21225f8bf2a219952233b328fb975360c8c7adb1e56335abc6ad41",
	"ae3840621877c76f4572ba238ccc4411a68fe4596d8634b5925bade40357274b",

	"79e264d13a72e2c262000f89429828ed66040cca0a1e7c949270221e12dbb80c",
	"71628b93d86ba6b8a7e72daa6c44ff6df82a80f91c3cab3e6124594cf284396c",

	"16788253ca890ad995049642ce00bf04c1d8083d25d974f7780ba963034d7141",
	"2d77a65d15122d977149e7012bd70899514545f42e9d3585251dd39fb90f222b",

	"0d80f3763d03cf4746439badfea8f28bae79220d636b374b31c401088e83ef8b",
	"1a3e56d52ca4d670dad6de6a34db373cbcc3f0b763638006dedac85882570127",

	"7f0004fc41e569b73041691b7974c51a9b1b04ea8f76fa1caf628e34bfb1f4f6",
	"3c6545119f989679e9e363c322e50828861d9a9e25ff8a21a26aa5bbfb7ec172",

	"88b2470104159787e3db1e024c3db4eba7bb795f347270b29d5ec7f02aad49c8",
	"f6e875b6f81dda89c73b79621ec942756b383db0a510fa931354c65c1d50871a",

	"0262670011cb98108e1c1e776eb7cefdda204c66f5e6bd1664229e61ed08817b",
	"2f307bbc3a3181dc432e869ac40468c587e1f469748892add24bcdbbc9971046",

	"828a4cca2cc9cccf7164aaf20770536c3174921d0b55c650b10faa588df4f2ac",
	"dbfc23a3b67028ff389609d6af5177c29da61b34452c7ccf8d6e6562fdaf2b2a",

	"8d7bd596075453d76c4562bbab98589c8e21f5623d60189335362329d15f9265",
	"f0097dc3b33994fcdaaeb6bf0629a3c6da6adf03c497ffd30768d8065cf07303",

	"a1770a277066fac683ac2c2df85e21e657f12c2e818e54b80e3bb1ce2c8d5191",
	"8e5207791c09593db77edef73767aa2ad9108ddcff5ea0fa0a2e7d711d9b5f67",

	"8bfaf9e763343eecde43e3a22ea33626a4a016d23148c25fe3ff48386f3d8836",
	"0e7847eae28288c147064370beaed29895bfdc22e5b1e9e89e0c7d6069889d21",

	"e7a2946b955bb7be27adcf4c59867c24abd018d46e57eaef53a0d531d1def7c7",
	"66d6c79a8ef0651c10228c80634254d5e948ea301fba87b76c72e9c1a51b6e41",

	"2583f37d0878c0cbf4bb6478a804c9f2f7b779c9c22104a4ebb320a3b501c707",
	"5c90c5abc2bc789c82bca2aa7fa247084a9514d4da8c7db077ed27aa7b105e4d",

	"8381b36d9e53d498d93446895736efd43a9d932201925d7bb851269347736e10",
	"1966062ec5e0de2601e9615a08d2eca562c4b55c44ecaa7abb487cc6b5b2955d",

	"a9f08891e93c5750784f818e885917e05e60877c8f6a80f2f756f68a7c19ea22",
	"c9b1b5412d6741e59aecf34a9b4f3c967a75ea47cdb31953d7a88a793a808525",

	"5672edfc5480a3fe3ecea6d3421efa43bb23f74a549df865e3f31b3ba083be72",
	"9755fe7a49517013d6a0a9c300ee4e9d363df8426a96acd7b1b4245282f95a3f",

	"8ef25822dea237140212c8fa43a6ac5b4c369f032c020b06a9d4dfb337612ada",
	"5d2f36e163362a8ef2050f7a56e5ebcc5cf47cb1841ee89097f0864eb942996e",

	"114606d85fdca97ef2717c6819212ce8823e1ba09157e96f1e10239c90ea731a",
	"085e7474203f9f8d10069edb5cd0eb53d9bcda19402a239ddeec703c9e92f534",

	"2fa6413288fd64a5390a01cc2bde69181113f2b3d4ac87b024577e0e8036f9c5",
	"0e9de1be3279625b8013c4a1646aedf1b0be75bedd4e73f049816dce6651126e",

	"dcdf5d997984a90b61be027b9d0aa5a9e6d4f93c22ce72162fe47c3ff6d92294",
	"a8ba440f9460ab0667b50119ea613300fb941ddc82b291ebefaa5e73a7bdff00",

	"041a24caf2b427154320009bcdb12f8a859dc02f4100cd9f731d5ab995d25c16",
	"1e560340a15e33c02edf628e2b4bf3e33f3ec50c469d88600b31e2aeae67a12e",

	"f1617c38e005c0d8452edcdab243e7da129484e87970c462a1582773e5166874",
	"43886629b6814a7f449b12077afd9e1990d9b50f4488fc6c8e755d669575d220",

	"f6c657965b8a999310799f4c55f9e352f9f78602a4a37c2e937aa117a153f2fe",
	"025834988b1c7719797b8020e7bba4a3b449c7ee3a74119c396c066db8d2992c",

	NULL
};

/*
 * Hash-to-curve tests for Do255e
 * For i = 0..99, hash-to-curve using as data the first i bytes
 * of the sequence 00 01 02 03 .. 62  (raw data, no hash function)
 */
static const char *const KAT_DO255E_HASH_TO_CURVE[] = {
	"861e72edc639eedb213402247d6407ae53560ef27e95c382f49f919bd2178a1d",
	"b72d63aaba58dc2452f076548f568a930650fda789ba831fe6240ad49e760335",
	"740f0b266692be55685a7ccd5128e831d6e12d23047e120dc80cc30393998226",
	"322a80507ea5ed3cd4c75679444de8467a2e2f80abb62793481dfcef5224a56c",
	"5a59995930f15ccb4c1a29755557a32ed272384cc0c99cdf455ad00a72afed1c",
	"57a2f7c94bdba928c49cad61eecf373454ea534fd8fcbb692a634b588f3ed547",
	"945c4fb41de6bd382a190b723bbea250836be9e916516aeed0198d640731ee44",
	"9bf508fe68ded333a854fc2c74a4d69e93df7ae49273e85ce4c9018e14e05c5a",
	"cd1da924eaa026f7a902ae14521174e307b66cd500a965445a07727656136c32",
	"7f3e88c73609d0aa799a197483de06d78f32d77f125dbb7626fef51b419fe072",
	"989d3c16763d1142928ff28234350dc31cdf34fbb411ef5272316f6d609ed007",
	"baed042923e5e0b1221e1c8a2e85819c00a821d6d31aa497f90cb688d4bf600b",
	"46dd3eee793887ec953fbc03299f29ada825bf1d3ded0c15866616d342a32545",
	"0e6e1e7722f1e56eefab4e29819045b0db5f88573fcb3fe2328e6bf8524e6214",
	"d34da5de99e0dcb27081bd490624bc24b2bc708458e9277648e3c10558b98862",
	"8b694aa9a5a1e9072f39d33d59f8fc4bcb23cae237b9c10dbf2c3b850efd364d",
	"b3823f5c9f0dbb47ae3e135bd9e3f1210600404c948b97388930ecb428191925",
	"2a7905eb92f67332ede28a9f9deb85155f3bf760d998c29860b07c69fc51f70c",
	"fe0d9083b7e90f4d17e978c2fb135fec262253b9d7d9be8a39ce70a17fe6812f",
	"0fdb2aae6916c8db6cf279a1293aafa49e44c8c029b7a83938c3cc0de2192c4d",
	"4474184a631e1875a82c8b56d954d684f5205b84aae1a6ea34387fd3fc618965",
	"34c3f83e195520a176e5b466c68526b732e86d4741253b83bf6a3fe5b1207000",
	"c5a7bf1c5c78df74d08f7755544f539db08a4712a749a296fc46f8ed68370636",
	"72577dc246ad6177c0ee965b9fce81d75d92858546204cc14f3122f30bb5a414",
	"63cfff6bfd336867f98f817089aad08a7e5b34ccdaae180939c1e200f736e077",
	"60c17f8e508606a337df2722aa23380ffb8e76d916b7b74bf29061777524b344",
	"1a83af7a615af8b93747a30bd37dec320d343d38bd7270bf499bc0c92c271325",
	"402529801d5d2702c44fcfd99cca07242574ab06465f795bca6ddca0c109784f",
	"0c2c2229b494070e7db740176437bb1213861eb05822b397c2ea7928972d7768",
	"c5b1e745e72bc8fdc521da65d00f1166b249c1921f2d1a295f5505d891ce7109",
	"f734e15a84cdef03602903ee596eb1e95d5a321866c77df49c7d4b1ae25a4e1e",
	"9a9f35fad601974fdb8a7a552aab3288f95eb207d20125f4289c3ffaea146739",
	"cfe2afee945e072176ae09903695920f9cb004df25158e663756d82ccc102a7b",
	"9fcb71820e2341ad3a1882d4857f496a93405920172514c89e0157ffbfb7513b",
	"348645699c431bd7199077c783932674cc2ddbc0a01c7a392558d6640f673216",
	"40d19d9df629e299e755be3ab4b44c54d4dd7f08e65288a2d9e1d26fba10da4b",
	"284bbeeb0cc4bd2d7208e758d902b8d71a01e9f72f3c783a7e8a458d52d78e7f",
	"cf46126f17b648ef1c0f1cb0b4a8fa22aabfbfb5a66008fbcc6719892f32e37c",
	"6bacd78ff0d6dd24df385b9c50ac8e9ad02a5a104d7c9ef74365de204d5e2963",
	"526c3e6a7923481480ce9254d34b3b9b174869ede1f7761cd794df4383d82f66",
	"5de4ef7c6c994d59c7d5af1d212836f7c7ebf4256585f3e47bca7399144cde29",
	"76f6811c3d1d3fc296dfd0f452117af0237cc2bacefee0f13eac2d04e0280c66",
	"1a51b5db42c78b7feb285212f7c2af8510cd400be3982c5f8aaf5a7b4911674f",
	"4c5f58d2b79df33e563fd8a42d62d86dfeb6bbceaf7c3db257be813cfe54f062",
	"6a2220c4a7507c1b89c95c7831382280aa9957b334af424d16f9bfb1ecf9932d",
	"83285ed4e62e64ad6cfda1484c042a352ebc4de94d5ae38f598a5dc049b97b59",
	"35acbee663ddf2202d8b1afed3f283f93c2201a908ebbdf3946d1c21e919ae4c",
	"4efe182b92402fa143ce11ae801c933f1cca888a431b811b6fe97c0e2a9cfb45",
	"4b8e7e48ce7a761befc3ba4df4b253e46c5134a6eddf344f3d44e571d9157a18",
	"fa639dfb30a6063c51925909d9b192283b90b6a78b0707824b5d7ccc3bebdf54",
	"e49bbecd48904995fc67602953da2f61c91b9ad490549ec94768b1249e214110",
	"946fb1f93254f175397e807c553c46f896d765cf6b12e0c80d7e2a4916a7b107",
	"b0936cd7e35593093db4954bffd6c477a2ffaec633430ff6f4e9ecc8b043683a",
	"9ac0ba7cc2542fd589cff9c339f137097f83f42f0f37f270b6cff6a9aaa41b58",
	"a88bc5dc7ea06c141933575c0226ceec5d1681513e5118dfceb090a6338f8e55",
	"2fa434804d5bd27cf6e9c2ccadb61b1162d3006e980fcae5a6af62e3d93b3e25",
	"19d3b3a009c0cbaff7431f8cf3ab4435822e953bc6b39b5082a0c6e9185a862d",
	"c5120a5f6916e859db7edb8c7bf4e0ee6d7718285e5bbb2b8bcb357310a8a246",
	"0a28e36e20cad5379cfe17898fb0191f581c42f7d3caf6c9b9e61bc1aebebb48",
	"3c94d00a09d2d6d171766785749e419486abfdf0a86bcc983c126419cf3a570a",
	"efa54d52539265c3b5eae3320487cae499b990c46ee458610b0e74a34ad5b106",
	"664efc21abfca2baa595a347c0556c8f9ce23de8529b778de3c8d13cdfe44a2b",
	"b8ee466d830b0c629f34ca7f4ceed98c92d8c76fa36a3961d81492b0c3892a3d",
	"bf78e67c6e22723ed6746d0baf083e40d759949c05a96ca2e3a065b32c882b16",
	"4eebf090d62e5eea40e3e452042c7f89e77a3e8eb77ff432f30957f62bb76009",
	"9b8425fc814c17e4102d01389c62f508963da6d98a8ca11d8ea6ddd008ec3e63",
	"bb1fce2d586d25c08496ac15e3fb737440b39583a25583b1683926d7db0cdb23",
	"18826a53178fe649cc877d8e5dcfaa844459b3dcc5760a6b9fa115a6cd3f3458",
	"8b05f1e1fd43b0be9344d5dfbbd41ed73094a070061b70d6b4de999ca871f510",
	"814003828a1f6e5daf54fa34710a935ed078b35ff4b16c16052791ff00b90d7d",
	"1b01be1f45997cbc2ad6ce5d5a85ddf79b899c338c65d21d8c1e4d00a4a1fe70",
	"653044b1a923fac6d5eedd8d86efac47923d7ae11a9615b380ce3e29f896b957",
	"44e0f3e7c29aa910c484aceed07433722de4c3eada4943b40616bb9485c6a107",
	"37550c2b6bf40c052a17c43e981fb0115197f3ce7402ef35f03abf68d5c55261",
	"ba8debe4587d2025a35741c5a89b2d2cee6df51c44b2e8b0b897ee8a981e5d6b",
	"e9f72566ca13f3e19504e205c199c05dc3da7228a08ec26d170abc25f9e47628",
	"891c557e298ff2fb581b2cec8a9eaca30e9b8502550b416ef6ff3c3194ccab02",
	"d356347493b35348983730c43211bfe5c68862e6779097e0ee6e72a2be3e8d72",
	"6942efa7a9c4682050257a9730fce592dc963aec071f5809165858709eaaa739",
	"4b05b689074073976ea72d8be80d7a4f6fc901285b168627cbe3b23c80db4729",
	"5df5049b00dcc0428fded1bd4f1ae66a7e95294d0d1dd19c02ad6b8b4f2f0f64",
	"2d897859114204fb9f8d975c6c084eeb24a7982ed49cc0fb66443b2113d81142",
	"c269f91836a3bc9cfbd199d3ea37cb0eca512cc3c05177e75885285485ba8504",
	"30a21295a0f3ba4bde9967f5493f73d2636736ff9d0857f31e60aebce5e18e6d",
	"12d0e154797f0ae938fd4538d9d2fff33a23e16f69d5b095217006d0e8d5c176",
	"1ab485cb1e96134090d5ee1885aaf5f7d90a5a96bf87225862b9249b4a1a2213",
	"aeaaab64314a6b48bdd173c801cb302025a25396116c633e79e165d1e9ee1e17",
	"cadc4f7262ed89b282a20ca38092504ccf1f8033f65f04fd8bacdd76f2eade52",
	"c54228e6f7fae77745b4ecfbfe323c546b95baa58da671d613c0bcc490292b40",
	"c68b7342e06739eb9d8632b7429342b9beb6a7660923aa6a0e3d992bce5f9a07",
	"7c36b5eb9f039ced61321e96d0ed446be221fefc5f701233366077abffb1743d",
	"3bed3db90e9aaf6adb728063e5518b723f1994a40262cb535bf16df102689f20",
	"9d57edbd744fcbef378eaedf6e51eae974a2e8bbc6b2b8f892938b5f354ab334",
	"519d2ba9129644043c610c6e45928fdb00e529ca1287d5471722f0dfd5a9dc60",
	"dbdcaabb3af4b7acbaf2f0893ee675c296e4ef7253fd80858c56da3890149e59",
	"3007e953e6f76886beafa8ee131d0c74c662461941c8cbcf1ce261c211f9df67",
	"da4c362b1def08eda514ef82e56bf4dbe7bd6af14cbd16414e8942080e900229",
	"bf21a32b133c49b0eb048ebb47c346e2cc4a1d52b4e92825414b2f6cf7185450",
	"8574291eb03ee42ac3212d697552bbe524ee954d6bed5b8844b4651f7929d515",
	"d84c40df970a0e8f7fb74571021823e38d4bcda3a381dff71572fe507ef32d09",
	NULL
};

static void
test_do255e_map_to_curve(void)
{
	const char *const *s;
	uint8_t data[100];
	int i;

	printf("Test do255e map-to-curve: ");
	fflush(stdout);

	s = KAT_DO255E_POINT_MAP;
	while (*s != NULL) {
		uint8_t bb[32], rr[32], ee[32];
		do255e_point P;

		HEXTOBIN(bb, *s ++);
		HEXTOBIN(rr, *s ++);
		do255e_map_to_curve(&P, bb, sizeof bb);
		do255e_encode(ee, &P);
		check_equals(rr, ee, 32, "map-to-curve");

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	for (i = 0; i < 100; i ++) {
		data[i] = (uint8_t)i;
	}
	for (i = 0; i < 100; i ++) {
		uint8_t rr[32], ee[32], tmp[64];
		do255e_point P1, P2;
		shake_context sc;

		HEXTOBIN(rr, KAT_DO255E_HASH_TO_CURVE[i]);
		shake_init(&sc, 256);
		shake_inject(&sc, "do255e-hash-to-curve::", 22);
		shake_inject(&sc, data, i);
		shake_flip(&sc);
		shake_extract(&sc, tmp, 64);
		do255e_map_to_curve(&P1, tmp, 32);
		do255e_map_to_curve(&P2, tmp + 32, 32);
		do255e_add(&P1, &P1, &P2);
		do255e_encode(ee, &P1);
		check_equals(rr, ee, 32, "hash-to-curve");

		if (i % 10 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * Each entry is: P1, P2, P3, P4, P5, P6
 * with:
 *    P3 = P1 + P2
 *    P4 = 2*P1
 *    P5 = P4 + P2 = P3 + P1
 *    P6 = P5 + P2 = P4 + 2*P2
 */
static const char *KAT_DO255S_POINT_ADD[] = {
	"a640cf808233d9ba7d04df6eb7f79bfca282e16788a06268458a23a67c786915",
	"6227d95a323ef8dcec8f963ee75b0c3610aaa85eb141eaf0c9351c3bb7614357",
	"3e662bf0ac88715ece7ae60455d8e3cd62ad6015eb15f3e9bfb1affa13ff0e2a",
	"d96d674f7963b382007d7969307150c4cc95b1063e4b073a0183a7519ab67752",
	"61e387b94adf41923549c7178bb6b8c924423eb96d6d0fee123203400dd0a44c",
	"2b8e33dd62d6f3572ea7c53770089717c25bac5279bf2f43c0fe2f1b0932f27f",

	"1a3a1178a2bd338bdf7ca49fa59ecb9a652bc96c6a48577aa881b51adf4d282e",
	"6040ce129c4c483fb92cf164a63a6a2d6e4848999391d4b94147fb8d9fe59551",
	"4dad0c9677b1fd43be0e79d2e7b6d90e1a5c55b925bbdfb50bd52352bbcfe11a",
	"dd11bd55b0b4079954655b36ffb20a96f1f4f163904daffe14b8e0efd750727a",
	"a250018185e4123367020d813e5b728b55c73c5e1ea7ac2f0116d54ba701ee7d",
	"32b7903197d16ab2fde54a282910aae877b858918c5627e9c4126e4e5dc8c479",

	"4ddc18a7d00a0bdda5a875971e6c731f2a2959130cdec37ee4913e54fd674b1a",
	"13e378c4c73a554ce2b28b9b26cd0f5dcb8fdd73cddbc71194d98da886e2fd69",
	"4b1d837735e6f1785f16ae1112989d538fb11ed2a1d26de666016d65be5a0934",
	"2d324d8c027cb2618369443b1dac74b2fcc99bb5d1cf24082b33fec22d3c906f",
	"97dfecbad133037a8b2848e3c8e29beb398fe331f8ac9c1e5646ef847c379a36",
	"1a5a97a862136f2f36777005d9e82947faa03b8cbc76520ddcd833c0f0300d74",

	"980578f78bf511d15c4bef18686a0ead48cd74fd15631e3300949e64d507865a",
	"e305b704759244b0fddbd82cfb5496b72a7eb8a76943fd73467df6073662d948",
	"7bc92979187737d63345c79420ba1251a011872357e023aac8539b9c0ec81331",
	"3e35c89f8742a9bc007e1fcbf5e736267b8f4fb3b1145b6c397662233e05f145",
	"8931d362632b2c3c19b14bd425d2c91d356573bf544c628e586b909eb849306b",
	"9a5bc0169afdad2c9e22478ce8b0b10dd881b4278e0b19561d42529c9c721432",

	"ef4cf925122ab9574e5c8e96da58adca49faadca63fc5049f8e61c89988b3d27",
	"14d1512620692e704b03839103075cc452aca3ebef3f9ab4dcff7bd0d62ef723",
	"dd6e8f2a38be87cbefece4adbc8fd3dfb0eb901da1adda39318f04333eb49a6f",
	"6af14986a4b3abe3caf085e3e7c373737cafbaffb8bb7f828d6c830806c2125f",
	"4d1c54d59ab47e71bc0e02354f1826bf9090d5d727badbbcb5f4a68f75a9dc13",
	"bf4ca4c88648f0a8bc6d608219da5c32e3c23f0bb0e0c65df40eba1c0f4b4116",

	"5a866e0e5b955a299d93417d16532948b06579c37b22b829f6751fcb1f222f5e",
	"c81e212420d832cae7575262188d2303e7a245c22725ad79d10ceda84ecf6c4a",
	"a5fb129ef16915a67975e42c6eaf9c766a4943d8983b104dfb2095e28daebe43",
	"a4aa6b8bbb5e54d4b59dd4695ba9369ab207981a5bfeb9e10d67630bcae67a3a",
	"44e8e363febe36514664a9f57713a10d5aa31b2acc13a93b7a411d386712146c",
	"e78f9b616bd6438c894cf37740e574163487b266fec78a7b962fc69372ddcb0d",

	"f4e77d81028ad7ea08fbff56163fc8b3b58956f09e0fef5c1dbe58e8b630c277",
	"83478fef1acff3b8e642530c6002e0c4093b39ec127fa5f568c88bd8a2e0fd46",
	"c358d1adb7f48a0c3baf35dc1c12734ae63132071380f26b0237bd0b6527b17a",
	"059a8aba8d5c99995ec568c9b799d4025e973088fa08665b6e844d32d754e027",
	"308e175d7a39cc88a1881c4476501c097cef861b65786c0e890b6fb6cee7890a",
	"9e7ea170a7c12bed18c6ce91455f14e90a310d18a60de655aec9106555683301",

	"d94b17ae0d156cd7701980163adfe43f53a52e4b493c24231aaacd1994ee1e2f",
	"217a0af174b9c48e16c18183323db11e7e152ba85ea64681fd5f3caa489c1470",
	"b538611b042a80842504533c0fc968eed008509773f9b97fa8b93416e8a96925",
	"480817661be3ea8a9f6cbf5213c7ef7f1007b207ec016ca141139f503212c005",
	"1928cebc87cfab8f727ee04c4f296604001f93e9110fd963f4e162f925c97517",
	"4950d17b374f18b9db841928058f58eff7f2d41208cd7066f8a14dab091bfd3d",

	"3451b6465903cd506a60b5d19e828f8fba9a30ae1e4a9a75a58dcced87db0a58",
	"87449490955acf02069a0677d351f2afcb5cc5e8aa2bbb2cd8b2fd53bac7a138",
	"4909536e0e2ce88ab529803022ba132a12ec1e140b59b1bdefc01273b6c00405",
	"7268449992fc72695f3d36908cfb70d9789940f03eaf1491d7f5fbdd93367158",
	"3b930f298efa8eb77007d1817c042bf1d439e5a1a6f9969160f5a17e9115fa29",
	"3a7e6b21aff84f0e95a627f9915a9b6c4393eef8e06c2184d6fd346006e6ce42",

	"01f3d5e31341701444bf1fff853e10ef39cd5e9da454e136190aa0c38370e46b",
	"c20930495669cf2ce86a3870a7341e881dd763b919cd17d9c6dcea5f9b1cef3f",
	"70c85608f178aaf809689cbe0e6be6344bf4c8ebfd04cedd822b809c34606172",
	"05e6d3585f6bc0e5bb6dbc12b6a3544cc16657738407674743a5c7f4d4356212",
	"3e1a83ab6c56128c9c789d65263b757072ab98731245430b979fa5a50479a917",
	"5962a17a505e76348a07a06fea32ddd380cd4efa4fca24bec30b18110281fc43",

	"9ba12f4e0eafb933c308fd5da5f4d3d55cdad1c197d848165d1bdf4a8d01b04b",
	"c48ac3ce9e271636bde728cfbf88bce1c7d4df91068aaabdf124c377f52b8a7c",
	"183a4140f451b54a43df9a6147a4846cb63cf829bb0e87547e082ad49f92f16f",
	"23e790504acaa7e2d5927dfe2cb98eee7ae47c1817c30334d5a5bd3ef4db1d2f",
	"1dd985dc467a5187d2ba6687597153085460c1dee9cff40705db34434c1d4f30",
	"0de054052572e7696e766ab3851e3c4fabf04d9c3fd1e4872f7f5736c1b34426",

	"4c3ffae19afb8620c6e867983a06e3cd3cc238323d93fba4142b5c801faae25d",
	"86a0d63d0cb41621a1f1c56030799ae34a57932e92d5081c3733f9496116101c",
	"a54b96b61371fcf4c3a2fa68c9bcfa804ff28b70d64dee42dd64cefa28849c47",
	"f0c949a2897ff7a9009f7d0f9b6464a82a64f8aaf16a452eefdff25e69503a71",
	"909236068a12858179fe0dbbdb67466b861469c835b3ac1b57eb9229dd083b60",
	"b18c21a34114f914cd8d155e7e6089b3c8b72e262aa83fd388c420a78a056c5b",

	"7c7e0d9f029cebc1efba013595a28b06f8ebcb5286a0d4611d336dfe4fbef73d",
	"ae068c18b240bc106b43e7aec262e441c700296284f83677f2141a6d5d7f6a67",
	"7bbc38a70a3064d8f52d8b6a0bbfa49578a595b0af49c8ddeac09fbdfd3f544d",
	"ea98385f84356007e7ffba03702e12a89086926936506fde777cb197d166bb10",
	"c6099de75d6711fed77a5344e8c7ee5355b608043ce9ee4c0c35e5f1e52b8543",
	"1fdf7f1fa1faf2041d58e757f65f251b03eccc093ab8fb206cbb70e0732b9527",

	"c062787a188868f598801bb2b45244d3a64987bc05aed3b2ba0238cbf1902a0e",
	"6bfaf1582a3518bf52c297cd0419aa4a197972c8f40dd9862008bf596f49d841",
	"8aafb17e7cc0ba1054dd06ffaa81803e20429b2dc9ed1d6594df8f5081add408",
	"9540226944ce86fbca7ec3d0041288c93f8de3ce714742dc9e16a43a631b8650",
	"c16ff9ff1e6bec20ac97addf4fb0027a7dd04841384c4a37a12e03c3083e3179",
	"a4bb87c6c9b78269929d0f24f428f118ec605400a5c26a22ac22b0b15f0efc0c",

	"8b50fbc3ef43eaacfd453ca05c20b74e22ab95db9507f9c9109d8c37af407215",
	"c877e5939c705a57fbaf5f75506a6172bb4a4436b23cf83e80c4a9edf242022f",
	"c8982f4dff39bc51419625d17bf0cce215f520e1a33d8d3b73cc21b1b2bbc658",
	"1e9077eee79ef046056ce2e3057d720f1870c59696f116cd54971152e639876e",
	"13bedf646bef836bf86f4aa03f08b3cd2f2532b45111053e92aba56986accc2f",
	"47cc22449fdc2af3e88904cd7714670505773b3d5483c4d332d25496555bb978",

	"8c691786576a60074ba1156ca5a4d0d7a71133b0c4cb7bbb8c87efbf0811b45e",
	"189d3e7a5ae4e805a4a1b9f89e410a70aad1927e1bce117fc7e365a89206822e",
	"03a508232dc80744e628fc8eca32b7a0abb85aaa7f58e49da8e3843a4ae58a42",
	"a1d95266f2390bcec21ac2ed4eba6ea5a506c9bbdc89c543dc894bfeb4457a77",
	"30a601da257fa562472a2f785f6a9518cf0fdf4de683ca2651ef12df2058bf32",
	"bb652bd7f4e631c461ad7c251a5e638c49b76af05daee1626e55cb7e12f0d146",

	"24156a730debf2b2a0f5efe84a6cfdd1d07ab61b6a16ea56048c99294d254721",
	"0f2722057dbcd9cd3124a52166ec5eb73e01bedadb16dce73860ddfbda1ddc6b",
	"2ff1a64d2ec074fe9d060d11015f8fd835dcc885ed266968c7baf437dbbc3977",
	"88c9a3bf21647f769fabe5440170f43f34e720785ef17a52bb487873894a3130",
	"7c5d1b2632051bc4304d2ad4353f484ad2500c52c2bf5d98e222337b699f6e6e",
	"df6fb6079fae6e04be70a8b0f1ee5ecb89392f14c564b698b96d99f4bbfabb0c",

	"48dff31893d9a18b9d140f545b6eaf96a03e752cdf7c099f4cd4e85dc3190845",
	"4fcb08eca6160c48057efcab290f3232d1f5f50b1b7f9be529f71ec5e40dfa09",
	"0a314c221238d38e37aa1600d319a88c5dda964809ab43d07ba4033f670ad46f",
	"d54657b2b5ba22cdd1c9ebe846ae4a53d6484feb641578dd26bc01c6e0d94409",
	"01816d3e0361f5f0bf007d7bded050b1f923881ed6be118542ee45ba9d95ae07",
	"089170ef5470057a02db477cb9e993eb00563459dd68e5ee4b1eda7c9be09038",

	"e4276557bb7938928a671a1d740a6eb174f46aae97da5b27a71ea3321bf57267",
	"664554c073fb6ba6c95721802949ee9dfb22adee718b5dd1fc1205a14aaacc7b",
	"a28a97509a17d7b05b73c57eb334257b05ef220f009fa9e7ea5d1b1df617486c",
	"2dda9bd368aa4eb61737bf91e9913f0f50d6aa10e441695c81f96306a48e506e",
	"7b128440f290897d78b553f29df73ab4831d4de1f617ef7d9759177e3bb35e17",
	"5af72fe3775742e97db400783f033da4a0620cb8305ffd8570788decc3035f65",

	"c9d8d4d274e82877e09ed0c0b7c32c46779bce29f22b7ca3852bacf297ec2c46",
	"a56593047ac3ccca17975e80655a44f32cebe8272c8d0415e70d65a8b3d69762",
	"ea59dc471f4de1f0937b20d2748c324b8d0e1b5f2592b44c769f44362545b852",
	"68ddfa06d3116251c20b464046e072811beb0351a2107b475a4515605bff9f19",
	"df74d55747c447d7ffa55545207ab3a26654590ad631174171548734572e8215",
	"89c435bd0da90d2007629c58db2f5220226363dc2c324e8242241a14c7f7de10",

	NULL
};

/*
 * Each entry is: P1, P2, P3, P4, P5, P6
 * with:
 *    P3 = P1 + P2
 *    P4 = 2*P1
 *    P5 = P4 + P2 = P3 + P1
 *    P6 = P5 + P2 = P4 + 2*P2
 */
static const char *KAT_DO255E_POINT_ADD[] = {
	"8bf6108baa2a18e584a9d9ddb84d60d1fb7152c0b070dc32211454aa63182028",
	"71cb77fb7b314b5917f984007e4b9e8da84e30df23b6361f3d5748859b3c3d46",
	"05efef11a57d4bddc59507920f67c702744485e17bca3960b9c75a6bd19ac535",
	"32a7f632e99ef0e67af738cbbb35cf6e3191b614840ff5071a67f0cd92a66169",
	"d3e338bd7a3d585b1d77795426a8ddd08c08b92a981087887d29e4643599316c",
	"9d9a190314a3a2f0640faafe0d63821380be7ce4b41bf2cab0bc0f3c1978196f",

	"7827dbeea7d4909ce79a9cca5c1683196353d9e15ce378e3b503c6f8fb937156",
	"a73ad61994300390f82e1a1bc730cfab1b47bb15e9feacf6831bfbc075af2554",
	"abc43d03dc36abdec998f1076a57cf673e9ab6ea064e3787d652fcbc84267f6d",
	"cbb97881af8dd4c407f4bfed6e61f02c74e87a16bcd8046336e324bc8c1cbe3f",
	"6553437175dd8e34e6a977f4abac327624245d7237d57133af319b46c76baf71",
	"b890a1ce8535b4b7191b666aaa8740f29ec858e6b417b287ba786b0eef881f1f",

	"16da527e26afc209d365df1cd4ea3613f2d40272d51d03e16414f662b4645748",
	"ab7791b15b2a06b5aed5ae8ff815db1e6a62495c30b8e4021a0d3296b93c4c41",
	"c439ec6d77de43c6ec146c843bf098b668200abf23831f4bab5aa1c4738e7c19",
	"2f82160afb2b3c6b394ae6861bab9e3077c810c7f136213d6ba94cdfdddc1d42",
	"fc69d563fa3c7b5b00710d7e5b940520e8990474c7bb60bc07d823ef5ae3837f",
	"a900860305b2398ac56d3a01c7823eaa80551b55f47e70eee3b9186be1da5662",

	"949f530930d463015efbd3aac6cee9e11e3522785f889ee5c8163a35fa06a73d",
	"b35f9e085825f3d62ff28700db54be2eb4efe30d4ca65ce16ea1b76b6c44fd2f",
	"06cc9347b401d20480fe0e2e85dcfe5b72bcefe30f5eb1c7cd900f94c084141b",
	"339e4e5c897420df32f25ebe5838b581e02f811685ce80f57c4406def28cd911",
	"efbf2eea194160d60259e400f546028f1de1ce57e41c24070e37855748d9b92a",
	"ef9fe3a03923d1a53a3f72f8f92bb2f538096e5414bacd9ef1b82baf4650c717",

	"522645d388ba218cbf28ab87b1ae6358dd1c609c46e1e174aff0d44c38855b56",
	"32cede21959628e69241141709677523a162fd0f9e7872cb8e4981e618293f69",
	"a813de7a5dea2df59bd5db75322fe62aada06cd0ca6e4fab52c364c3f9a02f35",
	"a1bb202cf9cc7237c8a8529e7512b3fc931e4b83feb4dde7ffaebaef4f5b4d46",
	"013a508f242ad1e205eb9c9437d49a80a355fbc5ecde40d14ed2e37918a2c651",
	"0984e413fd315d3ceedfc33b25d1119c4e5ac3153f481fd6304cfa05b69ca112",

	"0375359b6d6614120e445bfc974b888858fcc20469a1362c386f530bc76f8f07",
	"fb26cc0184231ee374cc0fb319bbebe07727e954fcd5316c9d752b2d5360733b",
	"4928c88e81caf3722037798d20bfbfac13ca3f89faf1fb0275ae334328a0d532",
	"d3caa6985f1d538c973a4c5bb65961f4c6b210fbad79f3f4a7a2875c591bee18",
	"022f0fe0313a6f07246906d42508aaeac38d7b6fc2b15faee197a757ef89a90f",
	"5b1e42fcdae23b9c3c3b6f76b75a151ac634eccc003afaa8fe28951c282de648",

	"4440458c78397f91bb5b8d55e62fa2ca76de17c15d9451605b9f05a482fdbd38",
	"e6ddb801f2dc0b11ad0c0f22216bee8bd32af5a706168e1bb6f57d8bbf492312",
	"5764a6d58f7b80a2f7ef3a5cdb4c15eff9fd2636e25a99db632658f5342fb925",
	"186e3ebd608eeff2980ae50b4fc7704890d83830d7f503d2d3c5f0b0521a4e28",
	"7f21526d9a0e2f70cef68007fc7f424429792c7a10025cfff6761242c0aa6954",
	"18ebdb169f70b5a97846b5534f2d92c5200a9362dc9680a53c80ff29bdef485d",

	"019c01e44f1b3533557ff98e13b51e9f821f672e410db1aa94ea0009e0474a74",
	"804bbe8b9847e5f2b79bbe85704f6db18e67c56d6d9f8e47169d7b923a27d07d",
	"e31925dcb78adc7ff9931dfb10e472096819de19e9f9ecbd9b564fdff1b2bd5e",
	"780068e4858c04e2bf6f5eb4ffa0a6b244440314f926ca3b64502590dfab442e",
	"781bbbb46597f4f8df52549bbcad7121f058ea9b9af2606fbb39f177c5734d69",
	"330839d3cdbc637542af3cb3d883a8c0df814f99f4f65a627ab55bd5c0ca9418",

	"97127912865b7a61972c2a26529afc606792395dad31ccc63e886b3baffdcd07",
	"cd613277de9313f158807e1b5a7aca4d4f88ad8acac8d562db7dfb7406d74970",
	"126b737a30166bad6b2b7067914969b83d299f94e2a0eac859dbe4fa622e6b16",
	"31dd812b5d40f895f411d7a212995dff94c4906a8904149f23293e7d3ca47e2c",
	"1ffab5efdfeb44526bf8a88c22efd5e79288c80f9271f835faa214d64060380b",
	"0a80970ae1707caadd76acd4fc3c1dc9f1f7f9db97542c3836842e27b72c6b1d",

	"d875a8b4ed944713a54541f8be96f6b5d240cb08d02ef13059ca883be716d53c",
	"795671ed1a8167aca295b1ce4f876839488a0515142b18b806671ad8820dd25a",
	"444fb1cd52b2e78f456848213770c024cc63c26dee0c95d780f393f6dd0da65a",
	"3105e2ecf53d159145db28f117415229aaf4f78249cffda30f6ae112b25de554",
	"0a612d272a741a0975064b575252da0f46d8e05f31e0d87610e7e51f376a5509",
	"06b816c9d80bf4dda119218f10d979c88147cf66dbcffc8494fde5e0d4df6314",

	"6d1136b1d187c6375e0330127b89308972cd648ca8d4c0b3b3084faef6a0ae7b",
	"4361d9f6b265d5759ba21f16b7ec7a69408273380e275540b6af742a1708c259",
	"ae8a77a19adf2dbc82c5c04d92c902d04f4a85ee4c80f64053a8843285df800a",
	"6327e44a9b9796328256748870f3f073f73746085f16a97fe053ec34e98aff10",
	"5d7fb07b9d2f53ec351bbd163b0b4aeb0685f1901b249d4d7cef4fc0c956e25d",
	"1884a71c2496a3c9f5898c31cec62c941b7a1bafb584cf8a325ca7823fd85847",

	"91cc0ccdd5d8d2b910458051630886f0673132a6b76a2c23eb645a89f99e3717",
	"9f0a91c033ac72c60d572c1f344ef9fc1e9c0172a5a7e281e6639fd7e4cedd38",
	"0e5fbc3bf39d05a7e147541e1ca171e7a15675b2373b0db070dcbaf80d5dc253",
	"1e3a2d1bb41c39867070d81067434ead00a750ba8627d0d480410e3184731c02",
	"fb58beb3611f9355ce14d9087fc2259928f58b2fd462f1d3836455810985a118",
	"24be6c7b725b3af82bdb75d71e97f105cff5cc5da34be848eefe48c9342ffd31",

	"cf817a9bd2d22183b03272ff3276e42c4634051f1315987dd9e0ae5a445ac617",
	"8ab17277cfed468c8c652b27856d24793590690b40be91a8468dc03abd7ad61c",
	"cd824cc69fb91c47a9ab94058e208eb1733629d3dc3f625435753b1fdbc47563",
	"f18abf4199b14cdb8ad2b326da607592604fd59c9704560ab188ca0dcf8f4370",
	"b244624df14cc53ff5a01bcf25e873c5c0c4dbeee33fbf55d84283f142cbf27b",
	"9c547dad9a911895a33189e6885849baba64ca7c63d9c377e6728a86d901f339",

	"5330654d131b9def2c388564f5f78b059fbb99ba7c6f06741972dcea83e75d67",
	"f0783bcf87b5924f38ab75810a8de6f588f2b222680f07b6f75988f38b1a2d2b",
	"68bdb7f049133af828eab84184a632c9aa0c0ed441584bd18547ff03ba7d6a47",
	"dc32f7696728eb88eba440d2f543b25d4ac9e8433dcab897671a9663d0e1ea08",
	"9051f8a29c4072a28075b1e4e4b44c50d2a502f9c25f05af85609311da7e5d21",
	"9ce241d02fd4ad88677bd5bb887e1d05631e672cdab4f49d390fe38d1b419202",

	"23aa1c0155bb5e9759ed2a0532db1f8e51f89dab6c1e75cc258445d53c227223",
	"e8a502c19b44324b448bf78d402e915000b607adb9b438af5fa869a687774c15",
	"e8e5546ebba258b43bc7f28b8bcb52f9e07ea7971271ee7fcb0ca5d43b824861",
	"2833cad71970a23547cd2cb00914420cef26dc28a1b950fc2d773bec352ba05c",
	"f82abc383a3aa6f4465fb4cafdbe8614a824337a86aa5ba53cefce8820772336",
	"68136ffba26af05621be27b015755821fad002b0f55e586786ebe358592cd702",

	"03148f09d7a2648c00abd3da21dcae3df8473d1b612bfb9a401ffc049b0c8354",
	"80c087a47761b790653e1d86a518115163df2857950391e8af9ca97c5a280f60",
	"3de574edf3a9f2fcbf23390226235a10a1459a298b62adabde41c2cd11396e61",
	"73c30df91adc16faf853ed42b81113caf584d6ae4eea7fa6e9ce6c951e5d130b",
	"f96be16653d7c87ca872603fd1b2fb61419c72beb843a6073d6dd0553b82cc0f",
	"f700f1d273e0d817fd8280ac57d1e502e99889713e0ebbbcbe59bf0b150fab56",

	"4e58f499625b53cc971e69a28145f03295b903fd9f3f84178fc6b8b634334414",
	"f6b39e655d6778231927cb055855d748253bb1b19b4d428096d36f66df9c4617",
	"581e8de1b87d2db07b44deca726b42fe92f124c35b401f3a1eabcd60a95de93a",
	"b5d527e23c48a98f62013d300ab15b9109d7017dbd372fc46aeca2d98009ba6e",
	"9e84aa553f5e341d425ae035a78c4f596675d1a1b6f7a778600bd6198e960e13",
	"635e4679ba93e1ff5aabe350cbfd9f5632c6175d603d3c010d90b4d62cb0bb4b",

	"8de5c537d165cb2d1e19f6935e6b0b4cd0f683434a6dbf7d7d23d3a868b26b27",
	"7db3664ddea75e39c231e003d3a2986e9030c02e796c14886a9152daa135572f",
	"7bf14a8eb8d4f292863b40d7a6fb24517cfd8412f10cd7006d0a8b5270609431",
	"df19b33b430b5b3cce3c3db41bf4726b3256c9e1ff532ad752cc9e24aa08370b",
	"4cab0b529c42814c8d479ab3f658482114ca5a91e0684741985ba945fe8c687d",
	"a281c63147dd8fb29faddfc2bc7c6f37ee37e6feec7c8acf46b7d7608267140f",

	"c024f64a2411676d42122c4b79854ef4e0dc00a22c547227e1557f236971dd5d",
	"50c8657ccabec0d17de623a2c41d34cf88cf9796fc7e76df820612b70c0cdf03",
	"fa7cdfc900600c4f1549a40deaa7bb692621e127ad23ba096175eaa1bcf44730",
	"093b873c1c294550da42e04c77521198d5e8b94e7e6c4bd6b980083f2324c067",
	"677c68c92caf8986d442a2cc3160f3906262c5e7a3467369fbaa1feba4e34033",
	"76bc5021fecd8c27bf08cf9a32ca5a7495b40d97d07d9eda73e7e221c7aaa50f",

	"9c354f93143666937ebed8402584e0a60e9ba754bc4e8ef7a7ca1dc5d38d0c7e",
	"1046ef3436f1cae1294503bd0c107a170206b9eadc25f2c6c7fdb27df5b31b74",
	"bc76ec59ae200a641878917894207f2d59f4045fec7bfdbb5a120f07b289b879",
	"30e61ad45b8c89cbec642f83cb6d7ac0b44f8687396cb39f747a084b57111f02",
	"57e17f471ef8b2d3568d59a1330f0ccbbcc15de23f43d2e6222ab7bfc58dc934",
	"7897bb18009af0f6c0bc83662dc8df9a00e26e030904fdfff13047e9648b3206",

	NULL
};

static void
test_do255s_add(void)
{
	const char *const *s;

	printf("Test do255s add: ");
	fflush(stdout);

	s = KAT_DO255S_POINT_ADD;
	while (*s != NULL) {
		uint8_t e1[32], e2[32], e3[32], e4[32], e5[32], e6[32];
		uint8_t tmp[32];
		do255s_point P1, P2, P3, P4, P5, P6, T, U;

		HEXTOBIN(e1, *s ++);
		HEXTOBIN(e2, *s ++);
		HEXTOBIN(e3, *s ++);
		HEXTOBIN(e4, *s ++);
		HEXTOBIN(e5, *s ++);
		HEXTOBIN(e6, *s ++);
		if (!do255s_decode(&P1, e1)
			|| !do255s_decode(&P2, e2)
			|| !do255s_decode(&P3, e3)
			|| !do255s_decode(&P4, e4)
			|| !do255s_decode(&P5, e5)
			|| !do255s_decode(&P6, e6))
		{
			fprintf(stderr, "Decoding failed\n");
			exit(EXIT_FAILURE);
		}

		/* P1 + N = P1 */
		do255s_add(&T, &P1, &do255s_neutral);
		do255s_encode(tmp, &T);
		check_equals(tmp, e1, 32, "KAT add 1");
		do255s_add(&T, &do255s_neutral, &P1);
		do255s_encode(tmp, &T);
		check_equals(tmp, e1, 32, "KAT add 2");

		/* P3 = P1 + P2 */
		do255s_add(&T, &P1, &P2);
		do255s_encode(tmp, &T);
		check_equals(tmp, e3, 32, "KAT add 3");
		do255s_add(&T, &P2, &P1);
		do255s_encode(tmp, &T);
		check_equals(tmp, e3, 32, "KAT add 4");

		/* P4 = 2*P1 */
		do255s_add(&T, &P1, &P1);
		do255s_encode(tmp, &T);
		check_equals(tmp, e4, 32, "KAT add 5");
		do255s_double(&T, &P1);
		do255s_encode(tmp, &T);
		check_equals(tmp, e4, 32, "KAT add 6");

		/* P5 = 2*P1 + P2 */
		do255s_double(&T, &P1);
		do255s_add(&T, &T, &P2);
		do255s_encode(tmp, &T);
		check_equals(tmp, e5, 32, "KAT add 7");
		do255s_add(&T, &P1, &P2);
		do255s_add(&T, &T, &P1);
		do255s_encode(tmp, &T);
		check_equals(tmp, e5, 32, "KAT add 8");

		/* P6 = 2*(P1 + P2) = 2*P1 + 2*P2 */
		do255s_add(&T, &P1, &P2);
		do255s_double(&T, &T);
		do255s_encode(tmp, &T);
		check_equals(tmp, e6, 32, "KAT add 9");
		do255s_double(&T, &P1);
		do255s_double(&U, &P2);
		do255s_add(&U, &T, &U);
		do255s_encode(tmp, &U);
		check_equals(tmp, e6, 32, "KAT add 10");

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_add(void)
{
	const char *const *s;

	printf("Test do255e add: ");
	fflush(stdout);

	s = KAT_DO255E_POINT_ADD;
	while (*s != NULL) {
		uint8_t e1[32], e2[32], e3[32], e4[32], e5[32], e6[32];
		uint8_t tmp[32];
		do255e_point P1, P2, P3, P4, P5, P6, T, U;

		HEXTOBIN(e1, *s ++);
		HEXTOBIN(e2, *s ++);
		HEXTOBIN(e3, *s ++);
		HEXTOBIN(e4, *s ++);
		HEXTOBIN(e5, *s ++);
		HEXTOBIN(e6, *s ++);
		if (!do255e_decode(&P1, e1)
			|| !do255e_decode(&P2, e2)
			|| !do255e_decode(&P3, e3)
			|| !do255e_decode(&P4, e4)
			|| !do255e_decode(&P5, e5)
			|| !do255e_decode(&P6, e6))
		{
			fprintf(stderr, "Decoding failed\n");
			exit(EXIT_FAILURE);
		}

		/* P1 + N = P1 */
		do255e_add(&T, &P1, &do255e_neutral);
		do255e_encode(tmp, &T);
		check_equals(tmp, e1, 32, "KAT add 1");
		do255e_add(&T, &do255e_neutral, &P1);
		do255e_encode(tmp, &T);
		check_equals(tmp, e1, 32, "KAT add 2");

		/* P3 = P1 + P2 */
		do255e_add(&T, &P1, &P2);
		do255e_encode(tmp, &T);
		check_equals(tmp, e3, 32, "KAT add 3");
		do255e_add(&T, &P2, &P1);
		do255e_encode(tmp, &T);
		check_equals(tmp, e3, 32, "KAT add 4");

		/* P4 = 2*P1 */
		do255e_add(&T, &P1, &P1);
		do255e_encode(tmp, &T);
		check_equals(tmp, e4, 32, "KAT add 5");
		do255e_double(&T, &P1);
		do255e_encode(tmp, &T);
		check_equals(tmp, e4, 32, "KAT add 6");

		/* P5 = 2*P1 + P2 */
		do255e_double(&T, &P1);
		do255e_add(&T, &T, &P2);
		do255e_encode(tmp, &T);
		check_equals(tmp, e5, 32, "KAT add 7");
		do255e_add(&T, &P1, &P2);
		do255e_add(&T, &T, &P1);
		do255e_encode(tmp, &T);
		check_equals(tmp, e5, 32, "KAT add 8");

		/* P6 = 2*(P1 + P2) = 2*P1 + 2*P2 */
		do255e_add(&T, &P1, &P2);
		do255e_double(&T, &T);
		do255e_encode(tmp, &T);
		check_equals(tmp, e6, 32, "KAT add 9");
		do255e_double(&T, &P1);
		do255e_double(&U, &P2);
		do255e_add(&U, &T, &U);
		do255e_encode(tmp, &U);
		check_equals(tmp, e6, 32, "KAT add 10");

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * Each entry is: point, scalar, result.
 */
static const char *KAT_DO255S_POINT_MUL[] = {
	"5dc536d2997635a8a3bbff585cebccb8b05f4689535d464a687daf26b0a85e50",
	"ee5f3f239a06a05bb113cd4f568e4f85ec3a6fd7f394b0d35c30af7ef25cb088",
	"294d9711c65e9d92d6a461390c74af1e6736d9eb1beab7956c95c5f93f7b0016",

	"22cdffbbb8f87841aab2ec1afa2c5fb5cb8f67e286d0fdcaf6124951a21a597a",
	"24c84823084e64146d09b64a58ed5f39256b568834c4215687ac642436a5d17e",
	"1ea8dba51fd35925488a671414cbb276c1e9bae3731967ca63cf53cae6428b0d",

	"2ca39783b2d55263d0b1c8f0cdce65f35536959e0127d8522b9624535d750e7c",
	"0fe8dd3dacef7df60e4ba71de3d35e5e6aeb6fbf7c960f72626fd3aae43c8b88",
	"ba32589f0ec3475f0ec62059e65a370df09ddcca269c0548d14e66967b95f664",

	"b514f908b623f7858ab324cdef9eea187b094192cf3c5c55e0a657cc13faeb14",
	"84d3e49037e45ecde991d0742d2f82c3ba9e2a0f0902fde2596369d9b0e5dd94",
	"8715245a20000eb3af5970f87237c4bba588af8e703143089b6b9f86fa1df51b",

	"73156faf9fbd65fd4af9ad2fec458760cb13ae044adf35bc4212ea0698e46d4e",
	"d350128f36b85d365060770fa450cb4a37a95842b5b84de39528c2a66b5425e6",
	"89daddb7da6f0379478aa027e963d140c6a3696b791140fee1e6b50ded5c1c56",

	"e74f7f5f0deaebbecac8ec772db37e5d495ba594e4f97a42e36ae544d1b6b748",
	"8ee094aed7a63b40981dddc40ec2e8ac9e1b8d3d3e40bdb4fdf6e35ad2afb262",
	"0a5fd0843da111a93b585a877690136cc610a0042d90937cbd5c88a7b9bfd76d",

	"ea9fab33e79dc4c47f83d95eb0f8cc0d486c78d1cb512603cd83aa7f7ecda23a",
	"a8e4a6496f60ea4cc879bfce649ffde1651c6258e6fdb6109dc431faf1f1676a",
	"7fd8e7464f704e0e4476ff6f188ee8b89ba8687d2dcf2097e30926289859c45c",

	"1bcb940392743b88e7bd5162aba895bd215e6726ae560d7c50162ccf60d61229",
	"078d80335d06b2f52557ba6382a32d6531bef339e7b6fd96881a332800f4d85a",
	"270537a047f33d4b82c4d9f12c3c3e1a3ee981dbaa78c7715f5874add68b9f7f",

	"fed66edc33f4bab3130c952ad4593bdf572b2b05ee0c422f4f7681d14549a55b",
	"cce339e0bfd3f9e30f14a44217e1c0bafdcb03120031a9eb234db4b25a6b359e",
	"433eed46a714dee7c7ec6b28fadec8d4ec67d9e0064fd12cc82e11a0b23d815a",

	"95603f9756b59436ba21c374e2ec38d179d719de24f9500bad1b4d06a0f7bd55",
	"0f0653922e02f8b71c7deb00d61dd6bb88ac94d158cac6f6018fef6cc932eb34",
	"87e70397eb175f392cd150dbf05ac8660ef8edc5ce9063fdc9a8d3002c90013d",

	"e01ca2fc36b41a801ac05c3d9f07befd1c08e25ff5a7b5f43956fb49a69b8e76",
	"dda5d59e8d25b5f6819d70bc3121120f2ab43cf61a676adeaa2237d16cbd35f7",
	"b55d2dffb68ee646b2b65bcd7164ee73412d8b2d81c7c53bf0b0f3fdc5367a14",

	"f989d9044394cb23198341262a09b701546acd9a096eebb4ca6c6fce9e790207",
	"e41ea9fdf13d5740295d30b7546a15bf74cc6fa08329839e11d1e902a1e8e338",
	"98a8193f4ac5ef3025992819a146e8d506e9c0896180cc89c9dc9f92f81b3f7e",

	"634e2ec626c1b57aa596122b196cd6cc4bb22ad2ade86703d9a45d24c8f19c1f",
	"55519d76b70e2a4fcf5b7f2161dfbdcb2ec0366f33f614343999bf5a43f07d4c",
	"92b4f42339ae4c0d290c4ffc543cbc13df85de7de8307126a595e6a32ea83349",

	"6046ca0c0c42306a8b600a0ead94364d22249b4a422d363c9858d512f6fc8b6a",
	"7e88fd3c2ff9fe0675e6b0e660c20a6ce0c111f096fe492a411edd9ffeb6cd48",
	"d317b8695e7797f2148695e2235efd7b764cf700d9010e28ccba1383d615864e",

	"dbd2df5aee4e9c100cb43fa8882df37d26bfb1cdfd732f9023cd9b626ac01a5b",
	"6e1aea679e9989af043f8b8ec0501fd62fe5f912d281f162f5702fceb659ba71",
	"a0a0c612850f2f761d47732b12f1b870218f77f01cbd9ceb9e96239e939b1a51",

	"a4b5c4f695e5a7da5081f2d00037220184cbe49cc69f3f0b6a6f9196c55ab924",
	"3f4d9db8fb2b46ea67329dc8110333dc2cc7c4d71a00a2fdda4ccd32fa204e23",
	"354446a1a2e23ece38fab88a85c926bdd4e338107386aba1f76192fc9565c178",

	"12c6501fba2cbbda12625bb61c197132a81fa722ef3a1401b5fca5104fa65e40",
	"1a688a922425c286ef6ff140ae11ca95cc175412233ce0a3abb22500c52c5bf2",
	"01598c40b7500e3c6f406d70932e5e2e116cf1eb93d0dc4f6fda6fa1ac5ad77c",

	"d00c41d984b324d821942f265bcfddef1b92d63d702470f50cd60015c52c012f",
	"2c1c60d6ad54a3a0e7b2f1b2dcf220ff51780217d6a6f1a88e677ef676b6580d",
	"99f8e913e519e0c94dea397dd62b2f80394cb0e3fd8e4c0eb66e040199fee156",

	"9c4401b972d38bc839bcc503950ce8d5cbdd1bd7b389124a21798bfb94dcf659",
	"198c7dabc8da7215fa19efb2d633f75790b02ef2e9b5d15117da2dac07aa8e7e",
	"a969d581f759f9551941c0f6ea801258ee7590ba0dce346d200ebc745a79677a",

	"137c9e0678b34e923b0e96e03d42181812105228f541b2aa2ca047eabf0fb96c",
	"a5f3e423c54bfb63b6db2908dadccb141ec9c8da8f20f635bd7542fe4bc311ba",
	"b9263ce912675c3602b6941f735610bd1ea46749b9f8f86af48adb203d631b6d",

	NULL
};

/*
 * Each entry is: point, scalar, result.
 */
static const char *KAT_DO255E_POINT_MUL[] = {
	"ef599ab0984eba6f1aee6c55415e9c30d1cb856b58492b645cec774b8f5d6e1c",
	"fef20582903e6fac0f8c4c4f09598d4feca23851250200204b3988dc73d3713e",
	"3fd5579d954e26abcf2cc8f5b452c3bc060eb012ecbfb81c8a061bd1ab54cc4d",

	"735b7c942443c90187b57707ef3e70cd6ebee4119e3e27a22eb27e19ba2dea23",
	"4606e325b06211bc19d2ea57e5413f7b7c3c9a2b85787af4164289f80239197d",
	"a7c405102137937c2589cdac9c73ec8afb939ef5e4508d159855d4c4cbfdc10e",

	"1d0128c75bc84f26393d0e8ed60ac65c49823aff69a03a57058f2a9f99d9fc7d",
	"4ff84fe04d1b5e9348138dbc179125c6b8dc8bfe67ffdb759cd317b794e00697",
	"40637e11cd493c67ce9664b91c1587f9199a1bffb62182c1de099c270f78da63",

	"659490c6898a71e9f7ddd11fc3d58a26f3673e4b8ffc73c62384d032393d522f",
	"beb7868f12ab5d17d83db56959f126cdde4d38d27a09dff17d64e6e917ea01b7",
	"6cdf044d38e8ee606e77be414255c4885c584d231f921bb9c04f88e5c4f88226",

	"3bc044686531f9a9475083a43fccb2186cd4b55c7a404aae1528d3912f65c028",
	"376c73d83abb5c746ae5cc34f86d2e782c109f3a1fed8fd0ad65c988e0fbca9a",
	"4f1aac7047637ac387c71bc2024c471e51639a9d9cc9e4aa046dfdd341322360",

	"62ab48031c59fa983709d9f8026c1caabc5d9576574e00e6a6a98d951b8f353b",
	"26e99a6e445c253168a96528b46aac08d94132a1780cb03c6646bc40a5aba447",
	"da26085b3f34ac6927ed9762b2ce5875e974282d8e8288b704c7149643d9a735",

	"6d3e16376b59f7ffc4a348e990e9ad2b66a0b04679b9cdd210bd59073ff55411",
	"3122b512408cb45009b26c266e2e868cc68c6bf5c91b1e292b56643eec5cb7b5",
	"800d80f6386bf5f5fcb51e7471ab05f15c486dbcc3c074c69ab0627c4845d56f",

	"b74aad34274fc6c328845899e1504400f82f632c0b79e5a72a91970a3a34f06f",
	"2e0ab2403e9cf3b909fc4f2c74715244c2a2d48f4b54ce158977bb1a5a44b58c",
	"fbb794ea32f07eece120ec7c516a136ae96f1bfe19a75393cf9c8476c7a4de2e",

	"dae2b85e5b6c6550cefb5b478cace23f7192d05d3bc4902f3fd7251f0e8fba6a",
	"a6d2fbd77623a7b52683e491d2a6d3e9c4c7bc78a96c95060faa12792c9f1a21",
	"54f67774b511a4301e18ceddfbc5d66690aaed090296a3f2bdd38fa59764e640",

	"138b69ce3c02c42384505d9e2e7c7ee9e2a45f8d825df084be48ced0b1b6fb4f",
	"c3849f4207ee75cba24a641a6ab9ebc7e613854bccbc385be35f27ea84123803",
	"9cb32db76b1073f0e6338288070f3090ceb1e59ac24a716840aaf7f408217a1e",

	"ac85820a0e9659376f241adbbab2da3fe797983096e7ee475a9ba2e49e74964d",
	"87c16396c42af29948cd1594f2396e59c45778ed36c739da0166c9cca2de2bd6",
	"bec4a7126f927e98807316f948fed476c2cd6ce131a8e8d1521e55cb1bf34f22",

	"51b4a3e99927cbc8bbb8584abafc0d1e915d2e5f0bf93617522593d5029ab935",
	"011916d017f7eeef2c95493d7b040f2fb5fa093fe0538294c854e74d861e5771",
	"38e5eff913b154c4b536de092f1734d59fdf0ae0f8e81af9fadd307b17259f59",

	"899c66330df8d48a7da39517f0babe7298fd18f79165f65218297258d2094879",
	"44c5603c739facc3a392c1277b39044e27c027eb2499e44796f38245e1f5affd",
	"4c6ec4dc4882256af0e362f849c3ebf5f79ff79bb79626d8516b301efbae247a",

	"e342ded1ed0f07137d7b2992e7fa4ea8e08bf45d0811c5e1d613e47938d5d21f",
	"d51a875780f9335ffb5a7aba312312388ba951bbe3980b8c79e390c0dcabb8c4",
	"dd72b3efd123911595d40909d5f27ed6b5b40e7677855f834366ad545f1fb34c",

	"9db6e81ed4fcf00e7b4a13509a976b25eb674b44ce6d47645cea7aa7f77b8e27",
	"03560cf6003a1a5bf34148c081a44ad345e6c3639c355b8ec63e396c77546fbc",
	"78f32da15f9ebe01b8d028eb3f4879204e1a141059d8d85b755947aa3b369073",

	"d8e11751381ad5cdccd6815ee6bdac4e386113dd24d2e17a48626a4280fcc94b",
	"a7b8933f86798c4bbba9589c100033ea2d0954802d1cc4d68f9282ea84bffb0e",
	"c961bcd739f530db91cbc83becc997c82e716dfe9ddf30c2946dd7eb57f62728",

	"6df18f24e7e180d23dcbac1f785023b8a11e55a2068438a745dc4301db1db301",
	"c843a3be8090195d566690801760a39626b12af2d398db96d24701197814cc7b",
	"46ed7c43da0dee7a80226e6e67f408f38a526c22e7593ad6e5933e4700af4954",

	"91d114745fa07c9ea408c091b0590ad2c921b83acba95270c948955c86619c55",
	"ff0f743577a47f9f73832ebf0e9b0d6dff89b208fc4b60e1ba18316af401b61d",
	"004397186d5503ed6eec0d6204c4462a12d590d2226d504c83add12bac59c414",

	"847610b28dfc4d72d52fdfc4e343e4508e65356319505868b643cb368775f07d",
	"ca761a610323e709e3014c5ec22272ce78eb4d5fafa0bd01b5d61aea7f3ffdc3",
	"ac114c9a50893d4b333cb05487bbf92b3c16228aacacbf51a41ea786b2f0ed25",

	"e1e1e70147d942f1ff27e837032a9bf34350209f4b562c6753e76039558ae728",
	"dc6b2d0b39ea749c36fed438e2e6f2d54f87f318a0ac80e8beb68c72f707945b",
	"0b1e250899cac9ea52e2837fbe1541a9a8d8ad5a29e7dca67cf1e9aa2b846f5e",

	NULL
};

static const char *KAT_DO255S_MC_POINT_MUL[] = {
	"a90da7dccc5a50916e7ce76e8292ddb9722a8bc210fb1a71258c50ea17e15722",

	"ae754a06df209494ece70bcffb24f0f4ddcd89c2674ddf9965d7ee14f69e4978",
	"fc34479bf9839114bc53daf081a155d68187d55809e684015a53f0f0cfdc9c0e",
	"c8c05caece52b2c29d6dc9f0bdd904ba42bcefe3f401098675a8ddb219ba3663",
	"88a4e1c958f1a24c016787c1caff2f6879b82380d5b02d92a29e2068bfa66328",
	"d50512d8b9580075fe99f77f78710900cc0bc21ec95cbe67e5de9e42e8c0f53c",
	"6da965ad8b472a62411bee2ebfc334868e3d6d4e7d82356f63cfbaaa86be3727",
	"85342d60cf733489500a66b5f130b5ed05e4fd35be61f177c4b8a30dddc8677f",
	"67d7674f7febe3e613f958010747a92d5724434f992b7a4446ef9c099b78b028",
	"eb0f2cb73d44fc21f8af9bb30502e44c8a07c3578feee67f0774d6541bf0ec13",
	"5ed693571fa186f9370e4eb4cea05ce4c28c7bba60d11f8c7ed0726d14a73f56",
	NULL
};

static const char *KAT_DO255E_MC_POINT_MUL[] = {
	"fae3f207045cdc8e1a73b743fbfd81b9d23bec1392603a524e4e53677ac8646d",

	"7c13b0849e3e5e26756e471cd2e8dd06c4d96e24c2edb547d39c9813579d0670",
	"bb140aaa8912e7379464635b220241b7eecd1cb21ac9ae63b9d6f2f2c076945a",
	"8a73932369250c2ced3c2c123558f150cf7b67320f18ea3e0a5271cb12dfca51",
	"a15c7e3bc340143bfa773b01530a034dfb0e9182d2bf5df120fbb8442e12ad7c",
	"4a2f85396ee4ef8b18c6fb5ac7ba7cdf4abb63e0084a675aa23ace64424c2f78",
	"11b2b72661412f1825b3c045906858345932542e04aaf916f5d0d65978883722",
	"6374912bf9223322bc6bab34b11440c779f1696f8d5031a4a09d51b6b4ce001d",
	"b5cd6da3d2465fa2ec920e09be4ff44857d80f62949ffa4bcdad91b60bfc8224",
	"f187c40e8f29bd7e825be662377d39f50e5d34e22a27ea910f5c8ac311ce1f71",
	"91bafd2379679082436e7c7ef2459f06e967225cea6ea37491b1ac901fab6f68",
	NULL
};

static void
test_do255s_mul(void)
{
	const char **str;
	uint8_t src[32], scalar[32], dst[32], tmp[32];
	do255s_point P;
	int i;
	shake_context rng;

	printf("Test do255s point mul: ");
	fflush(stdout);

	str = &KAT_DO255S_POINT_MUL[0];
	while (*str != NULL) {
		HEXTOBIN(src, *str ++);
		HEXTOBIN(scalar, *str ++);
		HEXTOBIN(dst, *str ++);
		if (!do255s_decode(&P, src)) {
			fprintf(stderr, "KAT decode\n");
			exit(EXIT_FAILURE);
		}
		do255s_mul(&P, &P, scalar);
		do255s_encode(tmp, &P);
		check_equals(dst, tmp, sizeof tmp, "KAT point mul");

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	shake_init(&rng, 128);
	shake_inject(&rng, "test do255s_mulgen", 18);
	shake_flip(&rng);
	for (i = 0; i < 1000; i ++) {
		shake_extract(&rng, scalar, 32);
		do255s_mulgen(&P, scalar);
		do255s_encode(dst, &P);
		do255s_mul(&P, &do255s_generator, scalar);
		do255s_encode(tmp, &P);
		check_equals(dst, tmp, sizeof tmp, "mul/mulgen");

		if (i % 100 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" ");
	fflush(stdout);

	str = &KAT_DO255S_MC_POINT_MUL[0];
	HEXTOBIN(src, *str ++);
	do255s_decode(&P, src);
	while (*str != NULL) {
		for (i = 0; i < 1000; i ++) {
			do255s_mul(&P, &P, src);
			do255s_encode(src, &P);
		}
		HEXTOBIN(dst, *str ++);
		check_equals(dst, src, sizeof src, "KAT MC point mul");
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_mul(void)
{
	const char **str;
	uint8_t src[32], scalar[32], dst[32], tmp[32];
	do255e_point P;
	int i;
	shake_context rng;

	printf("Test do255e point mul: ");
	fflush(stdout);

	str = &KAT_DO255E_POINT_MUL[0];
	while (*str != NULL) {
		HEXTOBIN(src, *str ++);
		HEXTOBIN(scalar, *str ++);
		HEXTOBIN(dst, *str ++);
		if (!do255e_decode(&P, src)) {
			fprintf(stderr, "KAT decode\n");
			exit(EXIT_FAILURE);
		}
		do255e_mul(&P, &P, scalar);
		do255e_encode(tmp, &P);
		check_equals(dst, tmp, sizeof tmp, "KAT point mul");

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	shake_init(&rng, 128);
	shake_inject(&rng, "test do255e_mulgen", 18);
	shake_flip(&rng);
	for (i = 0; i < 1000; i ++) {
		shake_extract(&rng, scalar, 32);
		do255e_mulgen(&P, scalar);
		do255e_encode(dst, &P);
		do255e_mul(&P, &do255e_generator, scalar);
		do255e_encode(tmp, &P);
		check_equals(dst, tmp, sizeof tmp, "mul/mulgen");

		if (i % 100 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" ");
	fflush(stdout);

	str = &KAT_DO255E_MC_POINT_MUL[0];
	HEXTOBIN(src, *str ++);
	do255e_decode(&P, src);
	while (*str != NULL) {
		for (i = 0; i < 1000; i ++) {
			do255e_mul(&P, &P, src);
			do255e_encode(src, &P);
		}
		HEXTOBIN(dst, *str ++);
		check_equals(dst, src, sizeof src, "KAT MC point mul");
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_verify_helper(void)
{
	shake_context rng;
	int i;

	printf("Test do255e verify helper: ");
	fflush(stdout);

	shake_init(&rng, 128);
	shake_inject(&rng, "test do255e_verify_helper", 25);
	shake_flip(&rng);

	for (i = 0; i < 1000; i ++) {
		uint8_t tmp[32], k0[32], k1[32];
		do255e_point P, T, U;

		shake_extract(&rng, tmp, 32);
		do255e_mulgen(&P, tmp);
		shake_extract(&rng, k0, 32);
		do255e_mulgen(&T, k0);
		shake_extract(&rng, k1, 32);
		do255e_mul(&U, &P, k1);
		do255e_add(&T, &T, &U);
		do255e_encode(tmp, &T);
		if (!do255e_verify_helper_vartime(k0, &P, k1, tmp)) {
			fprintf(stderr, "verify failed\n");
			exit(EXIT_FAILURE);
		}
		do255e_add(&T, &T, &do255e_generator);
		do255e_encode(tmp, &T);
		if (do255e_verify_helper_vartime(k0, &P, k1, tmp)) {
			fprintf(stderr, "verify should have failed\n");
			exit(EXIT_FAILURE);
		}

		if (i % 100 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255s_verify_helper(void)
{
	shake_context rng;
	int i;

	printf("Test do255s verify helper: ");
	fflush(stdout);

	shake_init(&rng, 128);
	shake_inject(&rng, "test do255s_verify_helper", 25);
	shake_flip(&rng);

	for (i = 0; i < 1000; i ++) {
		uint8_t tmp[32], k0[32], k1[32];
		do255s_point P, T, U;

		shake_extract(&rng, tmp, 32);
		do255s_mulgen(&P, tmp);
		shake_extract(&rng, k0, 32);
		do255s_mulgen(&T, k0);
		shake_extract(&rng, k1, 32);
		do255s_mul(&U, &P, k1);
		do255s_add(&T, &T, &U);
		do255s_encode(tmp, &T);
		if (!do255s_verify_helper_vartime(k0, &P, k1, tmp)) {
			fprintf(stderr, "verify failed\n");
			exit(EXIT_FAILURE);
		}
		do255s_add(&T, &T, &do255s_generator);
		do255s_encode(tmp, &T);
		if (do255s_verify_helper_vartime(k0, &P, k1, tmp)) {
			fprintf(stderr, "verify should have failed\n");
			exit(EXIT_FAILURE);
		}

		if (i % 100 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * KAT do255e keygen.
 * For i in 0..19, keygen with SHAKE256((byte)i) as source.
 */
static const char *const KAT_DO255E_KEYGEN[] = {
	"938b4583a72eb5382f3a2fa2ce57c3a4e5de0bbf30042ef0a86e36f4b8600d14",
	"b7d9c6c6adaf587b03120202c168831bc23b2a67f66d15aa0ce1382baaf14d52",

	"4a50b29655af442c8499c037f18ce8c278153d5b748baf796190856803d9773a",
	"ff2ac7224924892e240c731ec2e48309015d811568c9973d35ecdbc894a64704",

	"642e9115343df132c945255f0ff3d6969a9b5b7d7cb4adb2ca59716b5eacf811",
	"7f98f2bbf09789958c3613236c91b71762eb003778add1885490ee94c20aed47",

	"db4252337900d8ab7f609d170135d459a6798945d5a574d4a9db7b623c754d07",
	"5c9faa85d638a92d9e5c12afa46d49c2420efaab23190241530094b8efe72806",

	"87f1eb67b2f70780ff554843c8bd2980540cf15ed1d23499d013359759ef5927",
	"017536300dc562e33776b4ae9ee43a5e120667e292852e3aa978ab34fcc7820d",

	"41bc5a124a6e723eb004fb536a17dd1a577816bd83143f2397703f794ba3e130",
	"24d5334513be93d7eeac84f5694e9c20659107c17fe2055280e09d4c8819387a",

	"02c20375cc557dcb57eb02dc25f7f2390926071c3f4223cc321ba4cd0aec6502",
	"18dd6b94a95c55d29204751031bb915b6c6adc8878aecc305757b51e81dd983a",

	"493b59ef92575b18a6ad3ee1dbf2114c4aaa2571ae69543fa6a7a837532d0632",
	"f65a28c4afba44303bd8f1325c61d456e2e1f32b1525d0ccafa43dc36d204a27",

	"a662b8f9275568027b9bd157b6708b73ba8ca34a127e1c9f5d71d9ecb3afd423",
	"1cbe4b02bb3acbbc20764101c452061591b80d52eb89301424332c8a39ba9314",

	"288f7a3125bd12be279cb827b897c255267f97e4f7913713a80cb2630b646935",
	"59300e36d3d7d5a435c79b8dad6a974bffbeb3c9bc3ed54270e612dc10972451",

	"200f6e952677c538c9d831e3aeca6fc257a5b711d99cda07ff14042c4d8b833a",
	"9c9728c88778a9b6586d6f7049453160467829e16bb5ce4a7f63fac47e6f7f1f",

	"a56bc5fdf870fa3b77b648eb5ed73a6f0173659eededed0b7e6f6ebb886d392d",
	"26c4a638f7591e1a405a205594ce888992e24fbf5d9d8e14ea485d99c1908028",

	"47a53550f794a851c0dd80a0901cccbe981e40b1c4d1d69210f9f802e85fda1f",
	"289d5098634f96fe4e8ee706c580937946706db1b231458d4929cab8a482bf0a",

	"41ae0bf1dd0839c48045d6352642bb3119d92cff32063a1d3985f140fa2cfb2f",
	"c0a812552fc0516a4a3f131ed076fc95b4d98267683aeeed8d93655211cd5246",

	"4cf1f2bc5c4181cd0665fdae60d6e144a97e77a291188897de766170d6131503",
	"277d8ee6a02bfbfc276889c11ed2f5421713d7e61dd5e4dcd7bd080cbabc046e",

	"aabb07488ff9edd05d6a603b7791b60a16d45093608f1badc0c9cc9a9154f215",
	"da6a979fc42038aa88d27daf4df743322c4ab176c33ef1b310f5cd69f7f83b3c",

	"8d3e0b22c9824647109c1000aef645d43d83b327ef8d86caab540ef295edb33e",
	"2b2c16261d10532a98d7f88d68d5984b47e787842f158c9b46878c1a64e0bb55",

	"cd82320c3e9d68b16d40a3aefda68027950f0333bd845b66933d0bfabbec842f",
	"1e11506f4a8caef900bee6d2d9c80488ea8addedb1bfd663784a18a34615af04",

	"d3e72a448061f8c22a68661eef211387cd2851ebab1dc1cf4bff878a13686133",
	"d4bff10c04cfaa62793adb12318c9630361b4ebc62c17a4586585925f5570d27",

	"9ca015753a43f52834bab8182c96edbd895070a9998371ad3bab311da54c9716",
	"86ebdf40ca5df688a5534bd88066995f67502f1136818c609e1f92fc810d2f03",

	NULL
};

/*
 * KAT do255e ECDH.
 * Each group of five value is:
 *   private key
 *   public peer point (valid)
 *   secret from ECDH with valid peer point
 *   public peer point (invalid)
 *   secret from ECDH with invalid peer point
 */
static const char *const KAT_DO255E_ECDH[] = {
	"ef28a690f4129956e9309b080a4a9314a70f070bf23a4ac5a5c7ceccd79bf713",
	"1d80ce62a0a2b3916844eee076f0ede0f58c6b1569076c80d931f482af777023",
	"1958e2c1e90c16bd1ce032d65389acc03ca40e86d7cb7d2c3da4951d8bc1d053",
	"c72d892b575ba6bf2120d8d376db16a79c8bd835af80661bdcfd4b5149267907",
	"fac3236e3b3841f2f08274c58ecd704535a285ebbdd7767e1b45d33261083141",

	"c8e5bd9585757ca9e984c08b6d5f04bfeacca23ae70c08d6cf483566460fc22c",
	"b6544dcc56df2c151507957542bdbdd8cd89f905eeb7f586df21650bd3034928",
	"c00d4181e7f07693440b531c0c4c0321d11787155cb810955a07d7db6639dc6e",
	"efa30b7f745011a1c2179bdeb65bf676f6699ed7fac93753483f3a18a7b4e31c",
	"f8fbf93e12a1340a64e3540c9fd6868f9425c32e4d2dd8e05ed0befd643af586",

	"d65c50065aa45d6184712e758c95e7d73dca997dc5b670354e9be636bd31181a",
	"31bffbc837e9295be72eb55281c82a258663503b1846e8e65d492c6f78d33718",
	"4b220af3c0b0c53051fe3f392b7a750ca37bb193dbd8122922f8480cd10de520",
	"53b0e353877607a0c33d544da0cca0248fa932e9460e42734a25cfc404bd5364",
	"398c5dab9da224b7064be994198ad14b4b231fce7432fa92e301d9ceafa8ffd6",

	"52fc76c0f047480bfa05b3a5da31ab2c28828cb6b7465209b2635a9fd567b23e",
	"f6a806c9c087fd2b8af290f83446870447cdb1f9e702f8117da313ef96d97605",
	"76e190d84232eefcaf6140296b710c22d725176008f3aa7aadde7c54325bdc26",
	"8dde1bbcfc6a2be8665d07e992c387277638b86aab3f96cd7be78f6eb22da12e",
	"612d7131f7b769f918544cf187b95b30fb25366d4a3d03895e906319c171acbb",

	"4c0469f97b0245cfead8fbc104ca54d18a917010f50eee1a4fca31e6d0b44d11",
	"f4254ce377b6bfb57b6b78202d5cec1deb66b00ad7441c691ea5166156de7c7d",
	"7696af08c6d643d1e1473016a61d50edcdfa7d43cca6d8968a7b00b3a8f2e1ff",
	"d168edbebd43705c4d03cb3368d0e25aa333553da0637cf5ed7e6f74e6f23c0c",
	"b0dc53f8f55edd519f4106704e29c676df46f67d48039aa31debb30f0911a464",

	"61f3fa5e5c6aaa53815d2f4458bb55f9dc78dfa0b871bb507551538b7b7aa907",
	"80772b696144793e55dc706a007aff79c4df1718aafa8bc224ccb993a27d2a7a",
	"3a2d8c445b2c3905ab81524718801adb0a762e7c17a259427b367ddf92c3d4a0",
	"de029424e3c3f4fed8e42b36865a645895e2430373dfc5c309f55a130954a42a",
	"a6b8713d608f63f8afea100db16949071e5e257c4dbfa3ce9327285d12f2204d",

	"9754017a70bd45fbda023f0a315502d394a424e3191368228a6f2e8b62e15818",
	"8533475f7404f1fc9cb904c61768035293e3787308da15a601aaee3b60f8034e",
	"f260f1faba57ee8a8bac32a0595696eea9b9ab683c145f9274213d0717d4c108",
	"4f19f142b0fbd94563629227cd36f047ea7edbf973d540e778227966257c142b",
	"47560853d57536c2f910d26536d53e2a594432d02994c4a5ee4080273c0e7b21",

	"12604b60bb1bada85a247f47cbbf2a14647ab8828110353fee8180c7feb0f708",
	"222aaa94d01a6bbeca96eaa5af11b6e350e04a5590c7b2536606d6c363b4f41e",
	"77f5b1a5f608f96d952439feef5ee7547453e29f04499b123d3944d3e48439e1",
	"b79607b0cbeb197c8705e2999f255af712ebc583ce819f2bf0e40008f20eb23f",
	"32fe19e0ee4589f3d9fe0212c9a2ab3b8f61d1067867fb9cfb6b3caf457e12bf",

	"3cfd86fc7dbf092f33a0a3c92022684216c6d10576cd7b3f7f64dbf5d85a6726",
	"d4c1d69dd7868fc85a37716ef5621b6787a274180f27c9d6e7e3c08908b8d069",
	"274ca75b281d5fb6f3304fe716123dacf7222d7b9efe6453ac50e83ac145e9e0",
	"c6511b4cfbe214659a02642c1c0e7893ef73802b21ef39d10a0710c674a79e5d",
	"23f57d371414d13ff1a516fae7a9db71ea1633efd42c69d7a27f62c997b79803",

	"1db9a064623f90f985e3f66cf8ff46bfffdb64e2fb73c8dc31d62acd362df419",
	"d6ba5e02dff611438c8a82d01790f257dfa7dc0c576e8ddd27cb1e9075eae008",
	"d89bdfe6b78549e7abe8e94819ae2eb817cd57d3cf340a3ae962baf777745a66",
	"bb733ca5f060f156bc08ba9865adf28febc683f14146d91259327750bb49af62",
	"788ecb175be5a165bfef25a6fd296a059ccc432c354e55afbd92a261605c773a",

	"ca00d92e95336b0bff5db5cc624e6c58b0e259a8ec1b0f3e3df6b92b22751214",
	"cdc2f220bb4d0ee2839cda1cb2bd4ff3f8f3c0144e5bce7b040a543285016442",
	"5e6755dd5e3282729ab5255a9f244fb1716ee5f5ed86c3b0c65998fa36d6b770",
	"9ff0fc35d268faec8ac0facd3c8c07d35ed8e4015cc7a492e7035b1396d7d94b",
	"38eb314008e1514309fdd940880dbf14bab3dc6b45e6ee6ef72d5eea9cf1e288",

	"b56284ea283943c82b3326e7df3a7743b4cb5ffc8015465f3e5dfa14d1aa0627",
	"56ff5c7c341845844bac18776aad2eb36e76a2e740800ee5f9fa3e2189acd91c",
	"b9e5f0b180d56f6a92db793e0a3830cd9ea39f8cf8269f5c644ed47b22d61a3c",
	"d6f8a731ffd43a2da991fbf88b32cd7beeb77243d8e1d72dd966a17b21e5d657",
	"0620573ba4423ca1e6466d1dff02cccfd2895b89ff26a7426a1c903cb9970ae9",

	"ea394e65ffed2af6f374ee698cbed8fc13aab0af934a1b67b15ed7c6caf1a721",
	"e80856ed01be82403a378cf21f91eab1de8dea8b9fc2ed7da8916d2166add716",
	"3e98515fe823f348975cb5d3b58679d68efa0a3c063e5909b97b0ba4abdcf496",
	"373520a86245d80d6807cb1af60b912f1a4d3e277c7c05d60276484f76353c2d",
	"889ce022eebfd41a508e4acf21346698f35bfba35d5d284e00467d6610e586dc",

	"395446b0114c7c2f5fdf9c84b445738b4b6af5b421610a6234f151f58999aa16",
	"efb29c8aaa0701a1c71ec9077571bb6d97a21c8dbb580714d655c3aa50943b72",
	"5473acad9ee583c23912093e553db4a31b45a539333825bb81b878d6844b5a56",
	"0b49ed1ee6809c76e3b342fd54ee375e1a54cf297d179d590dc6927ce31f5c67",
	"caceca76026bd4f8553b261b2eb7e40f96785de03df5b3a263b152d7894b590d",

	"093b4502c12f43f1d68888211a083efc7e551146f7385fe0e0203c3a87b6f902",
	"2204f8dee58e75732662dc0fd8ae5c9aacd9c288e8520431ff34bfd229a0b856",
	"681eab3e0f73552bbcffe9e8acb22b8446e65447818be1c022a0e9d7f719e16b",
	"2461e08ebafc979aca3af695553fb1f510bc96bbcca1c409226ee642a301c112",
	"065e2afb921b257909e17cf7a9df7ab803a491d80315688aa94ef3b2d8e9fbaf",

	"296e9c430cf963e46fc2de1e48709a9f535e7e54184741598689b082f2fb712f",
	"e09b38d065ee00bebed7613ebc77392b54d72c1c59ac13975561bf00ae62be5c",
	"812047baacb81e4840bdf7c5fec3c0107fb2879eba441fda416736258e95d2bb",
	"42548ab4526e5176aff6631263881f57144329b132622b2ed0ecc49d1491f36f",
	"51d7cc6309ce07eeb509bef43a290364f2a285c0e8b99792a8711c40211c5a77",

	"d550a88023a96a010b7f0b7a93be2fa67bdb72cc2ef8a1addd9c27d1464fd211",
	"ab4e38805ce3870142e9327fb63fa51ea3ed9110f028bb225f4806fe34e3b165",
	"47a4065fbd3beb39477675449febeb0833736dd0585957019446b6b6f4d38003",
	"ffbf316c28bc263ea9e7ec7f19ab66ffdac6e751f4b9d5b59173c953c2f06f3d",
	"b89d3a6a00953f1586e1176254adc1dd098070b8e05eb54e3ccfcf9cfa17c2d9",

	"738b0034bcb054c33bde8b804a14c55e9035bd917cf09085133788835cbbb326",
	"d1c9fc36def70c49c4388050dbf568319a873bd86a5e17697064999015ef6e21",
	"28ae15ca72336e4fcd1bd5ad809abd7d4a7fa7c6a94356df2a199308a0c59746",
	"6c299b2954771e58e67cb03585e481c00a5a76bad24dcd34cfe988df8361d675",
	"4299fb5f1f0e1b9f3b02b0ae44eeb9e62275210257bbef944cb2aad78973236c",

	"303b84fec8756acdfc55d108800ed5c78795a79bfc718fc47cee55b1b2d1741d",
	"92a8956299d060f6bb3a1ae456d4e27a16054130460efefe3827a624a83e9a05",
	"3d7d64b38c30375ccc7b7252b432766e84b738e6ee40b48d722a94c47d85c246",
	"b3587506e4711f45e1f9de2355e1dd4ca403d892907ca3505ac842049556dc2e",
	"09ec0b7c6fb336f28d8b983dcfce94d0043e8090f1b74660c0b657a873d694b9",

	"95c863ec9ebdf963536500f4610964b96941d25e47f3098f5c6fd82704836703",
	"8129d1f4e3d8e63bcca497ba0fa234b6d915fb052dd6d75f7febaa1b0f31e80b",
	"64af14c411a4b660b3572382496a1b762fe5db9f049614aecfb80325d4ab0405",
	"2cb6327b08618354aaef4ed354c90875569277ee6ef8cc5e3eeae8e2091e1f22",
	"b8366704521d00c86b2cceab783b73135ab1bd80aa568aebfd3f24f073d80a75",

	NULL
};

/*
 * KAT do255e signature.
 * Each group of five value is:
 *   private key
 *   public key
 *   seed (can be empty!)
 *   data (SHA3-256 of "sample X" with X = 0, 1,...)
 *   signature
 */
static const char *const KAT_DO255E_SIGN[] = {
	"e7fefa881bd07b1c363b68a0f9e01ee04439f4aec941f402ae7be67993bcdf36",
	"4b383536c514fdd88fadd4699bb00f57cb16af174d086e65a291cb35957ee25f",
	"",
	"9eb7ccdd65cfb6e02f99210264c6424c911ab8890fcd93f064399d747b503ecf",
	"02df67f802428c4842a9337006f99bbe366e9d45330ebe59775fe233e24a1b471b450b81d63c76dfbf24b65dad69dcf19471a49346ae4e4239bc7846246bad29",

	"dd97542958b1929860ce00bbf16e70b96f585bdb40fe1c689c45089a25d2a33a",
	"2b2eaafa51faff7bb254388efe931a45efcccc9a18e0d2f22009ffac4a67dd33",
	"f6",
	"45515ea2aec5e64eab4ca100c26b030447c07b66a9fb6f532cba12e1f62ae627",
	"88bee51daef88d0122eefed288068e481d8600845232426c78833a2d3e2fb80d0d78997bedc3040667e65e37b6a73358a6b8be292c906dc924c0d87a1cf47c0c",

	"d621461ed0119fef78d82658d56db9a81e9de2a7dba43213941b150146cdd023",
	"82f078e1d7d4c27b87b44e94a5d65d10defa34b6c6187cdf27b6ba1cd0c8a416",
	"54e0",
	"61954339a2d5ff28e0d496af0716d757891010625c8859cf2ffde8fc75f3e74d",
	"aa362ed028a877b0e66c24ddfe04adff7168cafb6916e425f5997fe209de2e397eeed6838c1b2f37770612c2d816d1ee4b6dd72859131b92322ed6539403d035",

	"f1d2e620bdc765a7cca0da13c0181b15974d5ee3de8b66da9824b47706e8a90d",
	"245f4b6febec6d65a18dc721b92234b60f8cdd5cd6ec4bf379e2dbf1c4b3b960",
	"e2a81f",
	"e1bc41167127f1944e47d6c3e9775ea31919040ede6f1843baf1a039a363467e",
	"baa56a21d8d424da430e4b559faf3ce3d5c4a3500059ae91bd77b9c731b2a31657ca7fc36839b77045bedc022d17ebce315228dd7706225d3240cced47b0f21b",

	"7bf1ed60ed68b9ab68cd108f59ccf8c5244b3c67196ed682cdbdba199f2da826",
	"8d19cf9e1a466c45b25783295a5d7faa9c87aeebb66c3587b0a545c434c50062",
	"24e5f829",
	"68b268348810d58120705fe0ac3faaa11218c81f6560541891772c58305d3f6e",
	"f9dea4b04fdd1abcf5f363ac3baf8f9893c725ec7267257c7148fe6d9add5226d2ae40338d9b16826ce237b73f102271a5d3fc4c4d99b59824e4f90366ec1a38",

	"6a2c113c429adda15575b0df0e25423dcff7d999119e0259939fe13528976128",
	"242cc70b5bd1710e449864474cace3934da479971a462292ddfc86017533996c",
	"1d657c0787",
	"91cc12bdabea62f2fad09995b2177f3d01011d2b52fa51219d961e8801cf9b84",
	"afbf38a01391aabc062969fc90ae24844a05040804c254ab59c9e9ee4f27603a8702bfe8730c538ab1436e63b54dc73e1ea4826093a413b1940fd224de617507",

	"fdcb0a1f06aa24e5dfc01cd078e70509867d604275a1ac888cfccf45f4b4043e",
	"e84854f98bbed4bfc5cbb341ea31d74087818fb522828f0bc02300a9a3d84009",
	"5095bc20dffd",
	"0e298aef020cdd9cb9382d54f950d975b35f7807324d064c968cfb973f4d88a9",
	"5faef8cbc141c974bcb72a2174025730ba2cad3425e94302ab846c9b65ce935a6e8dfee1acfad9839506d6bbe7efc38b55e9bd8da0c591e9b206caaf4759f028",

	"e6f1ae77baf9b052eae31cae22ee3733fb976840f4b6ae4ad08188f00299f524",
	"2001c1321ef199a346ed03e826d93ea8d625a7aeb2e068abca08039409cd8c3b",
	"318a89d0101549",
	"bd6fa43ba1f3558c945ad97b86a16e46f18689d631832f36655808e0ae9625c4",
	"938102026e490a693d99f80c4385920006e6fb02a6fc751c3674fd739531e312aebf70726ee8454381d003310e7b3ea6010916283f7a54bf3ea8546558523709",

	"155e454cb96e895531ade93353db7b117fc99e6c4f5d3a7d66276e7f6230a42f",
	"1b758dfe2ac44b162db76c8763f0b942eb26dde011e6b69e5efdcbdbda5e9549",
	"b8000b24b6cad595",
	"15cad315726640f061a87262b030f8013d698829808ca191d29eab06e000d526",
	"c42e60dbaf1cc00342561dcbed70cc94bd80d37e261863721793682cab9e3502d94db97e792c2da4b0f3a386b8951b1f86dbda8371c1523dd67ca5c7ab5e101b",

	"4f73e0da62d285b6cb935de905cc86c5eb1172ce592b5a3f1b06329381db5e05",
	"42e38f43a66885ac8d5b55b096e1f2fbdac68c0cd734c95643cbbb1fc00a2e3c",
	"c37daa05a878983373",
	"a3a01a7f31c626a828d91a2ba34453fac4864a1c2418f22d135ac28158e12ab3",
	"75a0ae260f2fa0a64cd558a4adc7d3df092880201435d0181655453d3f6f715a57f9814be5265f9e087d47eef3b6cbe793cea5d090704161c5286df0cbec2e19",

	"ca209856bb248fccd018843b0b10e9fe1cb294a9213d78254a86a23ee5a79106",
	"7f8cb0ca92e9f8d0f1c9630833e6e28cfda0cb97d4728b1d28009abab545c234",
	"8b2bf5bec449e73da381",
	"3e2b0342ad12868adf498fc2a3c7bf168d115a3ca37fcb475541ead259f20c9a",
	"1699138979a0d5900446a16b31d7c3e30d438fbf482a6b026191f1fe11c37a7deb4033c7581705572a48b170d4e8c9ebbd3883f62f5bc5028bbc8b479e93501e",

	"9c5511323d7c143272ad5bbffb91a6af2f2be9a9d7a72643df3a7a790c944724",
	"6d76fe597859b316bec4f4c6a10403d84b8fd8797106ce8853a15e7a78318129",
	"4f3f552712d3c7ff7bb212",
	"8947bdb1b77d2117bc43d4933b89d835c28b0f5e228a2993f53711f0d6067499",
	"17330ef2c47365191905fa1ad9343b289340d62c44c94de0e1070dfc406fc335f12e800f4dc2600f18c61d9fa5cdb49f9af34a472d61720d346d402d9f2c6314",

	"4760c2b72a73b22cae763ce3903369fca15ab3fc3394aca4955e071213313820",
	"13d8e61e4a811bed22aec12fc4eb5176e05257c9e0ca62a75dedeb31b7d4832d",
	"867d02c2829dc1a05a553943",
	"ee4d9fd7516a121d065a71ba92bc70177ed7b12e1f62506d8f79127d790518c7",
	"a19560523ee01ad590ef3e6415713c3c0e7995cc2507de210e4ac4d20151a961a72c7ff988fe26ef7853780d6c87af552c22f61cfe0d083ce5d15a476484bd31",

	"6efd6cf3bc41af751a47a7a568c43fe05538c9bf4a4b56b92c87e384c369f701",
	"c790dfcf644d53bdb5c067df4a2d61e4d031d0422cdc6eae9e3a79572723a957",
	"9210040ba76b0d31aabbe1b79a",
	"1e713421d1a44b6ef1944e973d17764e44f580b91c53b9e11f4a1dc2783b435a",
	"6af33b5c8943e7e531a8fadb9a2a2da527f1219d6dcf9d56691dcccf505a4401f8e71e6710a5298005174c38e199782c7a26462219ab37fe0fb5919b1d5c8100",

	"9753da22261a95cdb29a35b2a53286ef3e091799fe88b80df7f081a05fbf8321",
	"cf8b7f0b66aad2ef6c53633c42f621406dd9c24ba97161a2efd451ed1422ae08",
	"4da18d848fa826de260ddc0c9cf6",
	"8584de70ef2a343366c334bc11b8d9e837ba63833dbb619127075bdb906ac403",
	"d9da40776f84ffbbdc0dd6dc765ffa3495e03502b5c628ee5a2f1d3ade7e2003cf43be376783a1fa6bf83644945f8036e08ab8c760fd78a456faa8794f77783a",

	"c87f825b563e86217e6eee2ed660982c7c9239c6cd90d011dc4a46fbba94e914",
	"bdd7e46cbba76e666be2a551564d449d3b3d10a480550c8a34c40ce696c82612",
	"338df4d100a8f105c7bc47ae7146cb",
	"018571598ecb3dc0b96a693a85de4927aba6e64ba6be14ffc9e9ff65896bc3dd",
	"79b67f4eef3383e68d53c63c76c185d87adec18f034ec476f416a12f5100a9050b059077f179902c0d8756e04e361f80a78b4596610f74d88fff7f7c58ccce0f",

	"d0cd7614980205ffdeb4cf87b81135d04a3975188ad1967b57f1a7b87d212a2a",
	"e7bedad9a67f17050b95d18978d1966cfc00300f39f4c9a706f7a1505e15ee34",
	"08495927e9336209338bf36dcf3aaa4e",
	"c02b87d23b05e56f0f7b96500368a41c23e83c5b0603642e5de35f4cda0d9b89",
	"d8effddf72c876db57beb1bc82b2e0feccd711bf5cc7384a92f2b9f2a987443744395abd62a55816333654cf1567a65d787da41444f80202d837dbce2c0d3805",

	"81db954da4ae7a5e0379a974094e64e6e23ff33e610de5e950d865f510adf325",
	"faf0d180629d98acbc28069304f932729689d72f72d4236b3fbbf189cf81141c",
	"fb52854932e67673d964ad8faad81bf1f9",
	"4739c6afd6c5b3f1f58efe807bb068940bdc6725057661660442af88f9050aa8",
	"6f2e0c0859922f76013ddd0830201723fb99cc58152edb121251e2a9c91c2b312bad2e2bbbe79ebb35b4b060843535b96cbe8a031a6eca6035e86706bebd1a16",

	"7f3365a6ac0fb684333b7d709b121369cd9a5cc6c191ee1df8cb40bda2635220",
	"9e07d0e40289c1c708f4d3f0bb3ee2fddafee94cf5bf3a441412edbc4a1e1c44",
	"da1eebfdebeb0b6e5c8999e3fbefd29d622a",
	"f8ec137d474d2fa4eebdd1f4acaa713743f0315f9480affad04af23ff7e4af4d",
	"cff66066d447a4d876590cfb0b02d1b553a7d942e6bbd65de3fdecf7ba6f2204861a5f4b112fc48cd0ff2f10ac987fbc2e2361c68da48f2cb3d4a05ec380561f",

	"3ea73a6a9e94de49026fa6c9385aed5a95b31643c0bc436e71ec15b55c7a6036",
	"78f4da98c85c33699114e047fcbb6d3cac9012bafe039e562356b1cae58bb53c",
	"e6efaa60b93f0c579dd6bd2f7bc6a60e4e9baa",
	"9403502de9faadba9fca6c0737f0f93058bcde07c7c0165714ace3084fb44086",
	"8dcb51309953408d891bd3a402dfb0dfced54130aabe5b490ad41614106d4c4aad6adb92ccc83a0fa95bdafd6084cf58b3d327525a80ca041a39e47bec759d01",

	NULL
};

static void
test_do255e_keygen(void)
{
	const char *const *s;
	int i;

	printf("Test do255e keygen: ");
	fflush(stdout);

	for (s = KAT_DO255E_KEYGEN, i = 0; *s != NULL; i ++) {
		do255e_private_key sk, ref_sk;
		do255e_public_key pk, ref_pk;
		shake_context sc;
		uint8_t x;

		HEXTOBIN(ref_sk.b, *s ++);
		HEXTOBIN(ref_pk.b, *s ++);
		x = (uint8_t)i;
		shake_init(&sc, 256);
		shake_inject(&sc, &x, 1);
		shake_flip(&sc);
		do255e_keygen(&sc, &sk, NULL);
		check_equals(sk.b, ref_sk.b, 32, "KAT keygen sk 1");
		memset(sk.b, 0, sizeof sk.b);
		shake_init(&sc, 256);
		shake_inject(&sc, &x, 1);
		shake_flip(&sc);
		do255e_keygen(&sc, &sk, &pk);
		check_equals(sk.b, ref_sk.b, 32, "KAT keygen sk 2");
		check_equals(pk.b, ref_pk.b, 32, "KAT keygen pk 1");
		memset(pk.b, 32, sizeof pk.b);
		do255e_make_public(&pk, &sk);
		check_equals(pk.b, ref_pk.b, 32, "KAT keygen pk 2");

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_ecdh(void)
{
	const char *const *s;

	printf("Test do255e ECDH: ");
	fflush(stdout);

	s = KAT_DO255E_ECDH;
	while (*s != NULL) {
		do255e_private_key sk;
		do255e_public_key pk1, pk2, pk_self;
		uint8_t sec1[32], sec2[32], tmp[32];

		HEXTOBIN(sk.b, *s ++);
		HEXTOBIN(pk1.b, *s ++);
		HEXTOBIN(sec1, *s ++);
		HEXTOBIN(pk2.b, *s ++);
		HEXTOBIN(sec2, *s ++);
		do255e_make_public(&pk_self, &sk);
		if (!do255e_key_exchange(tmp, 32, &sk, &pk_self, &pk1)) {
			fprintf(stderr, "ECDH failed\n");
			exit(EXIT_FAILURE);
		}
		check_equals(tmp, sec1, 32, "KAT ECDH 1");
		if (do255e_key_exchange(tmp, 32, &sk, &pk_self, &pk2)) {
			fprintf(stderr, "ECDH should have failed\n");
			exit(EXIT_FAILURE);
		}
		check_equals(tmp, sec2, 32, "KAT ECDH 2");

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255e_sign(void)
{
	const char *const *s;

	printf("Test do255e sign: ");
	fflush(stdout);

	s = KAT_DO255E_SIGN;
	while (*s != NULL) {
		do255e_private_key sk;
		do255e_public_key pk, pk_ref;
		do255e_signature sig, sig_ref;
		uint8_t seed[32], data[32];
		size_t seed_len;

		HEXTOBIN(sk.b, *s ++);
		HEXTOBIN(pk_ref.b, *s ++);
		seed_len = hextobin(seed, sizeof seed, *s ++);
		HEXTOBIN(data, *s ++);
		HEXTOBIN(sig_ref.b, *s ++);
		do255e_make_public(&pk, &sk);
		check_equals(pk.b, pk_ref.b, 32, "KAT sign pk");
		do255e_sign(&sig, &sk, &pk, DO255_OID_SHA3_256,
			data, 32, seed, seed_len);
		check_equals(sig.b, sig_ref.b, 64, "KAT sign sig");
		if (!do255e_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, data, 32))
		{
			fprintf(stderr, "KAT sign verify 1");
			exit(EXIT_FAILURE);
		}
		data[0] ^= 0x01;
		if (do255e_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, data, 32))
		{
			fprintf(stderr, "KAT sign verify 2");
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

/*
 * KAT do255s keygen.
 * For i in 0..19, keygen with SHAKE256((byte)i) as source.
 */
static const char *const KAT_DO255S_KEYGEN[] = {
	"f17dbcbef04a157b7e470b6563940017e5de0bbf30042ef0a86e36f4b8600d14",
	"07e4009d918a6c56f5a65896a8720ae26d4beeefcf67acd236886ca5475a1219",

	"0635a00de8e704b122b478bd1a0663a777153d5b748baf796190856803d9773a",
	"a6b4f6eff4a248044f7421d32c85ec7dec6677afc90df497beadd489f104f729",

	"c22008517d59517518530122a42f14099a9b5b7d7cb4adb2ca59716b5eacf811",
	"26b85f600b2f03b027024bb5156f16c6a39f72eb34c9133c3b50c14050218979",

	"db4252337900d8ab7f609d170135d459a6798945d5a574d4a9db7b623c754d07",
	"dc450fc99fe87dbec0f5f31edb85eb133ec9cbe6598481a09d2fb9bd92f48057",

	"87f1eb67b2f70780ff554843c8bd2980540cf15ed1d23499d013359759ef5927",
	"41492f6a1150193c44abae47ed6b164e27d808c10ac6d98c3d3ac6d87b7ca607",

	"fda04889dca632c34e1fb3d9939057ff557816bd83143f2397703f794ba3e130",
	"ef070d9bbff46988e6899511ebe3b532afb8caab3aa7a3bfe7280c23764a2559",

	"bea6f1eb5e8e3d50f605bb614f706d1e0826071c3f4223cc321ba4cd0aec6502",
	"0211f9839282ba0d86761aafb129e2c76da9312c492210aa96fe7d841b302c3e",

	"6312bea16eac7bdf93d5d2299aa8c9a248aa2571ae69543fa6a7a837532d0632",
	"8680740b543c20ce1dd3d2bc4950258cdad26b53a904fe531397e5e25e5aaf53",

	"c0391dac03aa88c968c365a0742643cab88ca34a127e1c9f5d71d9ecb3afd423",
	"e12f617a9656751a45db6f00746e90eb56e7b31933bf1afd3a4fa554ab31695f",

	"4266dfe30012338515c44c70764d7aac247f97e4f7913713a80cb2630b646935",
	"665bc6df189d53fc4d59d3b5e8cf45349ab7ae15f1ad1b0aa2705feef4741129",

	"7e01e5d06f93257b18e60da64307ad3457a5b711d99cda07ff14042c4d8b833a",
	"2e608d5fe618f81088d25e06d561f89ca90679c27673053d334a53dfee00817d",

	"6150b3748ba9bac015d100718850b5530073659eededed0b7e6f6ebb886d392d",
	"924e4fef4a4c467a01c0631ad3f3b2127e48ac4102eb806af3dee20b0ef2fd17",

	"a597ac8b40b108940feb5c6325590931981e40b1c4d1d69210f9f802e85fda1f",
	"27fceec7ac050d8353143392e2d63c38203fa1233af5ec37e8296639e0ba6770",

	"9fa0822c27259906d052b2f8ba7ef8a318d92cff32063a1d3985f140fa2cfb2f",
	"a23666c90f966c199ab0727a274abf39f4921d2f5df2d264da73a80856f4ef33",

	"66c8576f3896a194f48c91f71e8c999ba77e77a291188897de766170d6131503",
	"4675d07601fabed44c2835dd417efa263655c6fd0aca3fdfdf2dfc316240ef07",

	"aabb07488ff9edd05d6a603b7791b60a16d45093608f1badc0c9cc9a9154f215",
	"3ecbb71ee68ca84fbedd6167c48265aeee3df84dcf526b4042f52248649c0912",

	"a71570d4a4d7660efec3a4486cacfd2a3c83b327ef8d86caab540ef295edb33e",
	"05fbd65bf8cfc86403ef17f896c4396bb6aa3bcf9048c67330041a850ae7c119",

	"2b75a94787b9c8f3bc4d7f7192e3bd99940f0333bd845b66933d0bfabbec842f",
	"443294a320d54e81c54cd32d7f5b68eaee96331515887e19253cb897996c4260",

	"edbe8ff65bb6188a1890fa66add7caddcb2851ebab1dc1cf4bff878a13686133",
	"8f35d1a21a173ca61247229a6e16a98e77b37a4d52f7a0baa2e482ad3691b255",

	"fa928cb0835f556b83c794dbc0d22a30895070a9998371ad3bab311da54c9716",
	"a45d4fae9b1c3f0d553b2b181bd6f3f2362d59898b99b83fab2d0d5aa90d655d",

	NULL
};

/*
 * KAT do255s ECDH.
 * Each group of five value is:
 *   private key
 *   public peer point (valid)
 *   secret from ECDH with valid peer point
 *   public peer point (invalid)
 *   secret from ECDH with invalid peer point
 */
static const char *const KAT_DO255S_ECDH[] = {
	"0b60d9af0d907ad367c24916a2a5a7304d91627c062995a4247952fb90a67e07",
	"d642ebc9d2c25d04dd4f5a2b2fec9bb113c4a5e7cf18f133012a9ecc0c34805d",
	"496dc32968fa89024ccbbfa52a0fee974e5046dfdc61ea8b4c818e73e3bddf25",
	"76c6b5d3de8130ff877421680d399121e3940cb0de15f32b85ed660beeca5b05",
	"af8aa84ca1a92b5a493c09e10d015915659d53ee17eacaf30bcc6eef031c3a5c",

	"8f95a611095c4ef62bbff5655ef9709ec282624bcfcef0c6c397d2017e8c682b",
	"a80f03420c7a1f97f53e3ee1f7b250510870b47c05697aec5d0c6a78e5f38651",
	"5976a88675c86548925f463d02ca025199bae4b2b8eb543e60ba613f43650fd9",
	"bdbedc86551dac2f4e5aedd8dedd23f737bfecafc4d3ccef98c12434486acd21",
	"4edbf4dceef3fab819d65e748a8bac665d9a92e25d85678fd882186a7c5531ae",

	"5c1ddcb969ceff6d08376f37620e8db26d547988a9711bc0ceb81aecb287510e",
	"9be0a372230295c530ac8f1068b51798bb71e09db69a29e9bc6f24aca87c7313",
	"7727be9d3e2a49a12218501dd2f15f8133460e870d322f6793a5b037febd96ea",
	"23d40935d51eea778bbe72645278f0a738f33aee6b089a7a94c3426c519b893a",
	"db0b45c097684f8bbae4515a0b58aa4e1e42a5f0d883fedb758f6efb68e6581f",

	"58eb59334187c7b1177b19734e412115d8a83b87744d6f1c646f002d2fd56d2a",
	"3a52a1735586f1ae3f7bc97c3e04d342b861755874e92191191b03ae20ae9a3f",
	"02f70f575c05cc8f195c62f5c9371609f116c50aeddcfa037af678c3b032d7d4",
	"15483dea69a547c3a6ddbcb81e3374145344e64e38918ff40b3e0e8df5f78f0f",
	"55a7c3a1146e07ae19fab223e3d909672ac24e6c7949beadf6ee96493a7a6206",

	"6d206c6ef982d5de8b6142bf3d4e7a95dcb4e3db20a0cfd8fba72a5afe956531",
	"dd6e99b75e6b83422ce683f5a4081de22926767d3e34c21f51b4b88925256663",
	"7762e2dd5b5aed06ee17675931cc9540e9a917489e1ba05c82ca691e493f1f23",
	"a8b7403e0618b033ef62139deea274bbb7c7d455de8646bffd9f4a1ab3e3ba11",
	"389814e5e5c1479a8d7e7db177ee5b50029ab90ffe51e9f0c98685934c09a844",

	"e9cdb07c36fe154b43433e593f6f96cb2bd4352f66d530fb54b6401cdcb1942f",
	"15fbedede4fe2a5b48d74ed6bb44fbe85efc7da2648064d5e6e8c1c436de866e",
	"9279a6f3c4a6fcfe2be8de9ae27b7f4ef315cd360cad44e4def29b2e69aba62c",
	"a888375439a860c2f5f4b5e21366d7b27afd47c8e3c02f4577aa459c3db6a455",
	"3b9f382ebfbc5734f082867ebbacadea5e80e611b53400fce10312d657a73a74",

	"4a90bd7b42af4bddd61339010b58a87887990fd63883fb30825d15442f84de09",
	"b3d3629f8ebb06c6495c5e43da2dda4d70411e9b4dd2c5d547c1f71f4aa39b45",
	"eadc84acb8c474582062c9e028e35b06a027e19850dfbd6e5361dc02abb5f7a6",
	"5d54defb01469867de1586693abca766725a7614c0efaf3abd3971a365ee062f",
	"044ca8c19987d9e6147e70d2da59a55be892dc23e57e9b494cb01d4ede14ece5",

	"bedf5af3c0cacde96ec1774eddc4572bbd598fa53b9d8e91710260d67a4f782e",
	"a53dca0723c4fd22e8bc2ea0667a79658dad95c32b64f09b5fb5013e4e07391a",
	"4f0453d61d07dac237b44b6847408fc9757a8e6475b37d3275fdfcf740b34169",
	"5316bafd1104ea456a817c1a6967f2e3ad50257aad97ce4323090ab4809a7827",
	"5a2965d2e6cb6540c56d2e9a9e18a7f8b1f4ae904ca26b65cb228137c7c9ec16",

	"2a6fbfc17da763851cec7006abaee6a5890077891912801d6e7654cf4da0353b",
	"a7ee66e44141aec2ba0b01801db5785b686de1a1a1254f8ebc3de3c2505c730c",
	"64636087be0e692c209f302430febc61447436ebbd1b3bd92a1e89e971a07cc8",
	"8747b9b81a5152a075dfd12d6ccdb987cff6d069dd528fb895e5abdb9ca67526",
	"6568ce09e4f0390b060a2c8735d4137da14d46b39740b1ad503419b3e74a16fc",

	"978fe1674bb4875a12de35acdd4fa5bc703154f401e120c9c93f6a65f9015c28",
	"974a59c33edcf045159d8a8bcc4ea0dc8fad2ac8db0d2a6ebebcac8e503b7276",
	"dcab715ecccaf190a1d4bbca320908883c143df8c45055f32f83057b109938dc",
	"0d3c17ce3c7b8cb19d73a868e66892572688abd9ecd7c19b725efd3aa0138216",
	"119bb167339e0202cc240650e4e456b607cc19c54ce449e9d27047da7bde779f",

	"c888bb347d7a05609ba08dec3245db1af77dee4d854e07c0da9e20cee358473a",
	"faa6f8374e37d0d431f25b2b5249cfbad8b472b6aa931eac85ee3ac005326832",
	"ee942c37d828c77ad8cbd89effa1993628d336aca02b981f9df05b6e087537c8",
	"0b3aa099e0bca5abc5bd6736dcf6de294406e8e86f0af1e76db469fced569e26",
	"64cab2da182bf86382f4016fc902f5f9f79afffe20a0849a3597b86dc8eb13c2",

	"20074bd97f8f2d40ac4e1811bbf3405be9bc900cd7ac662b64f48d97f324bd1e",
	"8c1f8896cc56dfda73317b21403ba1ee38a9823230565fbe5706cb9d23a87f5d",
	"db25cdcd75c919ee5d10d1e9420061ebf4ea4e83a6dda8531c7ad63532981d41",
	"540910e3bb076075da39700f6b72cdc6c611a315ff0402a1471218a18e46991d",
	"79a4c3b6f27046ab7d7495dab1f1ee95b9da22e4660a1bb8555a9c7da6a9d760",

	"71fe7be0cefac677f315e53d9bbb5a317ee024ae8c82ea904aa7d5496dcd3a11",
	"f68bd5940810a1c3c85b50ac2b656eab983f0507bccedae3329adfc4876a7565",
	"ccf688fc2d007e78bc9f34ffb928839d6c60f6184cfc2307100c891439f4cd57",
	"865f2cf34ea2ee326d42b62d7868b9a12dbcb2b26858c333ea978df37b5d2870",
	"a33f09bef70a295702bfbfd992cea1fc31038667c2accc8b3873fd72ac143255",

	"c179c485ff5cb0537268906d66f36a25fb11b551662b799d761a4a9097f2fc32",
	"68bef03de98f175a15121a01234c6084882008eebab1c6cf95576724faa8e27c",
	"5ce14f7915ba47a14377c1bed66c79b4640cd0a243aa3b1748dd756a8fb68b55",
	"31596598561f475e5379c15129b1c9f833c831a376e4b639bb57ed77f7faa126",
	"6bcc12eee8e88b3ad3878e089bd07dd803256358ca71951cece61d65e62b3aa9",

	"6abe395062e48c5a4ff278c5aedfdca9a80ebe474ea5eb88849468802a58b83b",
	"b49af1e1ec15880fb7c5a0e5f97538a86130058f58758af5fb6d0272bf8f0237",
	"05145ccd3dbe28e76dfb7b48072c82a379e36fef1720d7343d257f488eb7a50e",
	"b355f9c90587debc900a7498491736785ec79c97de7b4b8ddcc3f9cfc8cc3f34",
	"4793c316b9ec1fc75fd9d83c9a050ace3b79940d2622208957c3ad8480f006e9",

	"6446de55796eb238e30fadff35977e2d62de45ce5a27c794017ed9790c869101",
	"2a0cdf633153be7cda181aa8e9967c92ce13212a7f13c337e25661fbc971c658",
	"82210570ba2a9783aa80ede04921e45c9c03c59076e8e36a28b735fcaf2dc3a0",
	"1097272b8bd51660c911f48e293198ef7b1dd83804d80f8728d83f1320b5400e",
	"01ad6b464607fe30fefb419c25709c386b30b58e79ba1075d7ba82f16d0baefc",

	"8ccc3f94c4c0f0fe6fee9b853a54907811a95f66696dc72843ebfc2ff8445536",
	"7dc8d90e4139704f345c6bd80d7e029af6e2c7208bf5d6a4b679a60da53ba87e",
	"0ac8906f35d3f43d973ebc85c51d2137f8e6c9271082d22f2f8183e0d97376ff",
	"154e473e1bbd2b3a4326602079aa69f066359ea2cb6a71455a8ff72da1f9f724",
	"c676c044aaf9dedb19b9b4113b41952234e9ba8d6d46e47295b4009a484b6b93",

	"50ddc9d34d0e9cbb63f877ca2724a29a02a5e5aa74eab98065a379c29f319d1c",
	"e5fa050131dfa6d36627441b248e3bce51553041cdb1292d4a6c1b08890c1313",
	"3fb891468269d14360b875a4f80071474e66cb8069f8086f15931ea69c943acc",
	"4fc6ef3f31cd4ab2cf4ba4f9ad9cbdf27e3c7819ff3d6556b1d3746dabfd651c",
	"4f6f90ec5f3634d888a5f0a395bb50924cde20baddb25eaf8c86461221912cb5",

	"8c707a5ca1995b331b83ef6fc2cf9d11106ed16da7b5283004cc76d4af708215",
	"2a92fc430bca0b9cc36f89dd9832e68a947723fc5839b117b8bfc47b12029e4b",
	"ad40abf4e50de087d601bc5ca8f62dbe8809f172fcdcbf074f1aa9a3ac86d73a",
	"008930651ab148f1aad3b9432a2258d15c3216bde9931ef1bfb1a06b29202007",
	"8bb1833a245779a2ff346bc231c6924d2b4c5a991336b5997de85c694153c518",

	"87553a36be57779d14f14f7a3e4112e1ceb051fe6a849f003f6a2bb9ac266e39",
	"7c2b52951827b02e740edf0ce42cc98a6c2d6422f1ccd2d7a873daa47533dd4b",
	"32213ceb8b8d81b101a004bf12ca17704a77dd7135e383769904425a10fd4b9b",
	"c6ea83ef596207d6493e07158da61d1e5ff5e323d00436e3d7c874fdd167572c",
	"a83506244231388ae7a981fbe8fc816d7ea350e3a5d5d52492ad03d0fbc911ac",

	NULL
};

/*
 * KAT do255s signature.
 * Each group of five value is:
 *   private key
 *   public key
 *   seed (can be empty!)
 *   data (SHA3-256 of "sample X" with X = 0, 1,...)
 *   signature
 */
static const char *const KAT_DO255S_SIGN[] = {
	"e7c1f12f0f7000d7f67a5a5c8d6a14700fba7a69d530c935c68105109ddce608",
	"e5a5789fe4447b6de68259cb8916a29d92b11bf817d51ba4423305c1ca22497d",
	"",
	"9eb7ccdd65cfb6e02f99210264c6424c911ab8890fcd93f064399d747b503ecf",
	"89d6cfd9a701d5a65b2f14570b340a21ba68d1d3c32e44a429b33d2db055910f0800d734bc10da80a0b2ab6e192d495eca1e71cc8ea53d8e70cbedf70a30643a",

	"115d54b74937d20f24816af4d1b35583f2c092285f59317b3b5c4ac4841ea113",
	"6bead75292150e0898e3404199f2d9490ac9ee2ff395d3a93458821bc29d2a1c",
	"b1",
	"45515ea2aec5e64eab4ca100c26b030447c07b66a9fb6f532cba12e1f62ae627",
	"59bb28b2be7735a883a2ce67b76352a9a993bd8dbf0739c4577b8ca6f6b7ae160dd09db54c27b3369405d6f5293628cb213c04fbc95550d7471c5465e4c1a239",

	"471329a46cd075adc0c4c36ac161a157e7900e7ce69457a58f9e2df4165cb401",
	"f306e28993bc5aa9ae3444bd59eb8ae3944c7b2b1393d53c182ab88566342210",
	"896d",
	"61954339a2d5ff28e0d496af0716d757891010625c8859cf2ffde8fc75f3e74d",
	"a71e987ee3a2e685376e2d91a8b169f6152b64ec59bc1a24ee070d9c5a06fd21db4173a6dd85818fe8d1f1c4136ef78e4e5112d8e4f8649bc1a1bda876df8d12",

	"cca6d9c557d0fec812da945f15258fd49288ac956ef6cfb46d5108ce14f2a604",
	"78e06a93e5d7e560f6a82695e7c5a3c74424f1ac0ea7cb9810f8d6aa7193ec7d",
	"356c39",
	"e1bc41167127f1944e47d6c3e9775ea31919040ede6f1843baf1a039a363467e",
	"57c144e06afe0e8434dd5619c6ad0edbc8117e42619267622a44df8efd8aea378d01662a382c49a6d6dafacb83762f6516f305150f1eb519740f758f2c604f12",

	"4f818f5503cac86c4988919dc6363e16646f677b35a8ffab5340d7144c301607",
	"7d528e8e6ba5fe7b942a31b881fe862febd67c53425a4eb1624bfadb70f1f430",
	"50c27940",
	"68b268348810d58120705fe0ac3faaa11218c81f6560541891772c58305d3f6e",
	"aa26f4b682f080f32cc265590c44cddede467d7454333dcf469687c629fb60668d6a8504e178e0492023996c05b4f8f320c0f2594f3fa76b6c21f4314c4c6425",

	"ec92eec04e4d608b7949b6673b223b0fada4b5cd0a45ed018d4d8341e7d73817",
	"ca8cfb6022395629da3be657808fd3985b34a5ecd20de61b023323b44e6eb711",
	"e47d5a9a3f",
	"91cc12bdabea62f2fad09995b2177f3d01011d2b52fa51219d961e8801cf9b84",
	"936599dbf96cdc7eadb145aa8e5c86c008d8810652cc91bfb26df3076b16e178328c28f89c4e4662be676d7994a02b96d9ae6c0f57e9702349c0c25955fe601e",

	"60f9d0bfbccf4b880b8c7961ef4c7f18c1424bfd38d68812ba288da371941d09",
	"903b6ce5c39c3a15f1f6a2782462cc39185192253fa1a8532e027e9a38fed13d",
	"e38710ceebd6",
	"0e298aef020cdd9cb9382d54f950d975b35f7807324d064c968cfb973f4d88a9",
	"28b268e0ce9ee196537ad9145ba51bb7911a21c383f698809b77f0a67fe355611559b223aaff4371d43d1c2f1eb15b4521b409062247a93b9e041089dfa21b1f",

	"38cab63ea60e7ccf213d04b3416c795e047f4897cc79968a3a28c2487c8dd51b",
	"00679269eebba60db4b8c2ef5dade2828d907ec582b8eefb13ea78b770402254",
	"e9ec0b55184b69",
	"bd6fa43ba1f3558c945ad97b86a16e46f18689d631832f36655808e0ae9625c4",
	"b250f5ee957d1fb1627e386947c8c6763797b3dfd818ec4601e169aff00db10ef88e829f7fc668c6f7d301b7c79f7d1feb45089ac4b4b7043185b2374a274815",

	"a0c852217d57f54af23b22c24663c08dfa39080a4ce479f15b29ca98aa1f7126",
	"7385259b9b2c90cb454e2a06ca0da0f7a2bd431a341d8df85592613f6519c06b",
	"89f04c834dc56e08",
	"15cad315726640f061a87262b030f8013d698829808ca191d29eab06e000d526",
	"3f699ff1c85c8355e95250a08cfa319ae38c5631a060878fd3feaae1a16c4d3a02d5f9c610306ecfc6be5692c01ae31dc915df2d864343d507093b0c396c3e29",

	"5f2aee1f839af4437708be00dc9fdd4b464f78cce503a2a46b7e49c7f8801e31",
	"cb6a952d21f8dac3e73719c5a2e25a5c0d25d296baa914a91017d4e2fd607543",
	"6ffcf67ae37c48bd44",
	"a3a01a7f31c626a828d91a2ba34453fac4864a1c2418f22d135ac28158e12ab3",
	"fc4e13e24ca421d18cd2abcac996c5583da63dab668862870e8338c984797d3f3f1e8ffabab56fe7b845d8f635d25f410105286244c288ab63ea73c2d75ba32e",

	"5180783172ddc9a101ea14d267606e7b086b9f199d285b9d963ad74e9c6a6129",
	"316e6c618476ff7f51315d9d4d7f2ec63099818b5d41fe64f109ad57365fc86b",
	"b5ec13448a6074acbe72",
	"3e2b0342ad12868adf498fc2a3c7bf168d115a3ca37fcb475541ead259f20c9a",
	"3fda8fcb7d32f98ecc7f7e1ed218bca9005774260e9db9ec7f26c0b9bef7dd0a1924bdca462df14986e740e8fb6fb498efb55765dc8de83eac6fb8cd9c8bed02",

	"6bfbb1a80603bb69eb67e30ec92a9643009223ae066317cd580763c74a5fee06",
	"e5625e553805421578ba2ee5399a5e6665ab343a67d0915cd4f789d22ee07652",
	"af57e5cd20151b84cb9b50",
	"8947bdb1b77d2117bc43d4933b89d835c28b0f5e228a2993f53711f0d6067499",
	"d7c5d21aa9fcd5632c5c5e79accdb1002e9930615f9d1f07c5f443f82ebed07a052060c8aa6b367bf9a4fda0b2cf2b488372ff57b9b4c6c48ade2d9dfbd9d024",

	"cebef0642b7dddf73345d784b85018b8c1f26b4fb9f1bfbb717bfc153e82c60f",
	"bda06c4aa5091c2bc2b7d81bf6962935f1e16fa4858bfe0bdd36ceaca8553200",
	"df6d2a11e996f8622c06c698",
	"ee4d9fd7516a121d065a71ba92bc70177ed7b12e1f62506d8f79127d790518c7",
	"ecc43c40ab39fbdaac1217d8407b540fcd41af0179ac06b595afaf1602ee1f564d9501f96e000ee249728df187dc8c3143e30f4d2967be37ebde283889f54939",

	"d163603878563d7a36f840d22c7a850363a6e8b094e471f60f7bfdc75da18e0a",
	"0721d84e87361db69031bf970e8bd5f3ae9d86c751b59417ee7139230f88da43",
	"fd6326d2da4c6683622d212e37",
	"1e713421d1a44b6ef1944e973d17764e44f580b91c53b9e11f4a1dc2783b435a",
	"0262368df64a673dcbdf5ba3791a1068e941e482d7d5aa06bdfa75847a4b5e2ed4b005878774c8884822bc8cd9f19d9f55a6aafb2f2b44899dd1c80e3e57bc0d",

	"93ae6d9673f5c48e26b37047d79c95bb1fe9aa1b496b5aaa1eb49eca880f2339",
	"7151af8a2e5cdaf205eb64172e88b23c257c73b0c1826a09ef5f02af37955118",
	"4f5b169cd6ae70f1490291708205",
	"8584de70ef2a343366c334bc11b8d9e837ba63833dbb619127075bdb906ac403",
	"41c78cf12f0c1fe16ab9d180d086d2faee58607fc891ab8132e641cefc9ea25e4e340928d90259419a08a9ebfa789ceff91bd86caf49d3875e83573a631fe226",

	"3b70f2cc6d753938d23ecdfab9621db34d1b606ddaafef479bdc0cf5fda83a18",
	"e73f4c27bfe24e2fda57ad76aab76554dc9f410589e78c93f4d152b5c4438807",
	"1ca32967015de8ec8276e1408ff3ce",
	"018571598ecb3dc0b96a693a85de4927aba6e64ba6be14ffc9e9ff65896bc3dd",
	"2e523cd538a1500265144b5d8c5d75e614045d398ffbaa44a9a65cee8f9b5717aa1e6e7ae16ec82682eddb547a4a2cb02feaafdecfa07d102b9802358c412f2d",

	"8becb40aff59d9d93289bb56b1957675794211bc3427882d1c54c7f3563d6c2d",
	"e29619933aaa702f0b733626eb7b7d3ae54047095d2448cf9a0281c1eb20114b",
	"38493f0d9d3f93d250a46ccfac0ad4f5",
	"c02b87d23b05e56f0f7b96500368a41c23e83c5b0603642e5de35f4cda0d9b89",
	"d31254216e4fd07da82d0c5ab5cc904dd475f44177e73d38e96551de6ec54c17c8b5a7239f05abea895e5831a58ad7c16246bdd6f40d6cba0880e18420e6612b",

	"25f0dfdeee9e6091cf0d2dc21c511c9a83cc341b62580e046aaaac8635aff502",
	"5b2244cfe0fdbc2dd9ca51c788b2ecf633c9877536e2ce2a9c444b58469a7b53",
	"a7e6bbd4302d191956d3e23a247b46b0ba",
	"4739c6afd6c5b3f1f58efe807bb068940bdc6725057661660442af88f9050aa8",
	"195773730c0136df3c4544437af5a8f10ec973e597990d2319b46eae58710a02c6c219c2afc30aaf96d769dbaf2b5157d4664299ea6e5b1533047d2b22e40911",

	"bd3c1a07a2961e85dd0ccca39045567a1f1f4c4e4f41e1190c8d5712be9a5b15",
	"b2faec833d71736c0fc8d28c5261bd4b61f29679d726fc17df3d68599e14435e",
	"c15641937d84e33443da70b5b5458a767097",
	"f8ec137d474d2fa4eebdd1f4acaa713743f0315f9480affad04af23ff7e4af4d",
	"a289f0ebeae6106d6e9d33721814c640a4cdaff4f856cdf582ba41819b958049a1eff69cdcbe022dfaa6889924687d5ed72e94982d0d2213f6a8d8e76d04931e",

	"56beec9f18c40a560e280e84d25480a68b6cabb1771f81828e3eac64d3d7530a",
	"723f50876ceccbd10f706ec389e34ecddc54552a1d6b5b2307879143dca1a365",
	"ace83bb731fdd197c1ac260a0b243057ca5f3c",
	"9403502de9faadba9fca6c0737f0f93058bcde07c7c0165714ace3084fb44086",
	"ebd4f2c77aa76448041f9c4e405fd4f1aef2fcc9a150a2257c60dd9dd68761657e1f4e2b7e331bede740598c35247401b251412261d4144832eaf2be18e1061b",

	NULL
};

static void
test_do255s_keygen(void)
{
	const char *const *s;
	int i;

	printf("Test do255s keygen: ");
	fflush(stdout);

	for (s = KAT_DO255S_KEYGEN, i = 0; *s != NULL; i ++) {
		do255s_private_key sk, ref_sk;
		do255s_public_key pk, ref_pk;
		shake_context sc;
		uint8_t x;

		HEXTOBIN(ref_sk.b, *s ++);
		HEXTOBIN(ref_pk.b, *s ++);
		x = (uint8_t)i;
		shake_init(&sc, 256);
		shake_inject(&sc, &x, 1);
		shake_flip(&sc);
		do255s_keygen(&sc, &sk, NULL);
		check_equals(sk.b, ref_sk.b, 32, "KAT keygen sk 1");
		memset(sk.b, 0, sizeof sk.b);
		shake_init(&sc, 256);
		shake_inject(&sc, &x, 1);
		shake_flip(&sc);
		do255s_keygen(&sc, &sk, &pk);
		check_equals(sk.b, ref_sk.b, 32, "KAT keygen sk 2");
		check_equals(pk.b, ref_pk.b, 32, "KAT keygen pk 1");
		memset(pk.b, 32, sizeof pk.b);
		do255s_make_public(&pk, &sk);
		check_equals(pk.b, ref_pk.b, 32, "KAT keygen pk 2");

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255s_ecdh(void)
{
	const char *const *s;

	printf("Test do255s ECDH: ");
	fflush(stdout);

	s = KAT_DO255S_ECDH;
	while (*s != NULL) {
		do255s_private_key sk;
		do255s_public_key pk1, pk2, pk_self;
		uint8_t sec1[32], sec2[32], tmp[32];

		HEXTOBIN(sk.b, *s ++);
		HEXTOBIN(pk1.b, *s ++);
		HEXTOBIN(sec1, *s ++);
		HEXTOBIN(pk2.b, *s ++);
		HEXTOBIN(sec2, *s ++);
		do255s_make_public(&pk_self, &sk);
		if (!do255s_key_exchange(tmp, 32, &sk, &pk_self, &pk1)) {
			fprintf(stderr, "ECDH failed\n");
			exit(EXIT_FAILURE);
		}
		check_equals(tmp, sec1, 32, "KAT ECDH 1");
		if (do255s_key_exchange(tmp, 32, &sk, &pk_self, &pk2)) {
			fprintf(stderr, "ECDH should have failed\n");
			exit(EXIT_FAILURE);
		}
		check_equals(tmp, sec2, 32, "KAT ECDH 2");

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

static void
test_do255s_sign(void)
{
	const char *const *s;

	printf("Test do255s sign: ");
	fflush(stdout);

	s = KAT_DO255S_SIGN;
	while (*s != NULL) {
		do255s_private_key sk;
		do255s_public_key pk, pk_ref;
		do255s_signature sig, sig_ref;
		uint8_t seed[32], data[32];
		size_t seed_len;

		HEXTOBIN(sk.b, *s ++);
		HEXTOBIN(pk_ref.b, *s ++);
		seed_len = hextobin(seed, sizeof seed, *s ++);
		HEXTOBIN(data, *s ++);
		HEXTOBIN(sig_ref.b, *s ++);
		do255s_make_public(&pk, &sk);
		check_equals(pk.b, pk_ref.b, 32, "KAT sign pk");
		do255s_sign(&sig, &sk, &pk, DO255_OID_SHA3_256,
			data, 32, seed, seed_len);
		check_equals(sig.b, sig_ref.b, 64, "KAT sign sig");
		if (!do255s_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, data, 32))
		{
			fprintf(stderr, "KAT sign verify 1");
			exit(EXIT_FAILURE);
		}
		data[0] ^= 0x01;
		if (do255s_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, data, 32))
		{
			fprintf(stderr, "KAT sign verify 2");
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}
	printf(" done.\n");
	fflush(stdout);
}

#if DO_BENCH86

static inline uint64_t
core_cycles(void)
{
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}

static int
cmp_u64(const void *p1, const void *p2)
{
	uint64_t v1, v2;

	v1 = *(const uint64_t *)p1;
	v2 = *(const uint64_t *)p2;
	if (v1 < v2) {
		return -1;
	} else if (v1 == v2) {
		return 0;
	} else {
		return 1;
	}
}

static void
speed_do255s_decode(void)
{
	size_t u;
	uint64_t tt[1000];
	do255s_point P;
	uint8_t buf[32];

	do255s_encode(buf, &do255s_generator);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255s_decode(&P, buf);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s decode:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_decode(void)
{
	size_t u;
	uint64_t tt[1000];
	do255e_point P;
	uint8_t buf[32];

	do255e_encode(buf, &do255e_generator);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255e_decode(&P, buf);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e decode:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_encode(void)
{
	size_t u;
	uint64_t tt[1000];
	do255s_point P;
	uint8_t buf[32];

	P = do255s_generator;
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255s_encode(buf, &P);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s encode:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_encode(void)
{
	size_t u;
	uint64_t tt[1000];
	do255e_point P;
	uint8_t buf[32];

	P = do255e_generator;
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255e_encode(buf, &P);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e encode:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_mul(void)
{
	size_t u;
	uint64_t tt[1000];
	do255s_point P;
	uint8_t scalar[32];

	P = do255s_generator;
	memset(scalar, 'T', sizeof scalar);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255s_mul(&P, &P, scalar);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s mul:            %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_mulgen(void)
{
	size_t u;
	uint64_t tt[1000];
	do255s_point P;
	uint8_t scalar[32];

	memset(scalar, 'T', sizeof scalar);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255s_mulgen(&P, scalar);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s mulgen:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_mul(void)
{
	size_t u;
	uint64_t tt[1000];
	do255e_point P;
	uint8_t scalar[32];

	P = do255e_generator;
	memset(scalar, 'T', sizeof scalar);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255e_mul(&P, &P, scalar);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e mul:            %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_mulgen(void)
{
	size_t u;
	uint64_t tt[1000];
	do255e_point P;
	uint8_t scalar[32];

	memset(scalar, 'T', sizeof scalar);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		begin = core_cycles();
		do255e_mulgen(&P, scalar);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e mulgen:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_verify_helper(void)
{
	size_t u;
	uint64_t tt[1000];
	uint8_t k0[32], k1[32], P_enc[32];
	shake_context rng;
	do255s_point P;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255s verify_helper", 26);
	shake_flip(&rng);
	shake_extract(&rng, k0, 32);
	do255s_mulgen(&P, k0);
	do255s_encode(P_enc, &P);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		shake_extract(&rng, k0, 32);
		shake_extract(&rng, k1, 32);
		begin = core_cycles();
		do255s_verify_helper_vartime(k0, &P, k1, P_enc);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s verify_helper:  %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_verify_helper(void)
{
	size_t u;
	uint64_t tt[1000];
	uint8_t k0[32], k1[32], P_enc[32];
	shake_context rng;
	do255e_point P;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255e verify_helper", 26);
	shake_flip(&rng);
	shake_extract(&rng, k0, 32);
	do255e_mulgen(&P, k0);
	do255e_encode(P_enc, &P);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;

		shake_extract(&rng, k0, 32);
		shake_extract(&rng, k1, 32);
		begin = core_cycles();
		do255e_verify_helper_vartime(k0, &P, k1, P_enc);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e verify_helper:  %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_keygen(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255e keygen", 19);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char seed[32];
		shake_context sc;
		do255e_private_key sk;
		do255e_public_key pk;

		shake_extract(&rng, seed, sizeof seed);
		shake_init(&sc, 256);
		shake_inject(&sc, seed, sizeof seed);
		shake_flip(&sc);
		shake_extract(&sc, seed, 8);
		begin = core_cycles();
		do255e_keygen(&sc, &sk, &pk);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e keygen:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_ecdh(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255e ecdh", 17);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char sec[32];
		do255e_private_key sk1, sk2;
		do255e_public_key pk1, pk2;

		do255e_keygen(&rng, &sk1, &pk1);
		do255e_keygen(&rng, &sk2, &pk2);
		begin = core_cycles();
		do255e_key_exchange(sec, sizeof sec, &sk1, &pk1, &pk2);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e key exchange:   %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_sign(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255e sign", 17);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char hv[32];
		do255e_private_key sk;
		do255e_public_key pk;
		do255e_signature sig;

		do255e_keygen(&rng, &sk, &pk);
		shake_extract(&rng, hv, sizeof hv);
		begin = core_cycles();
		do255e_sign(&sig, &sk, &pk, DO255_OID_SHA3_256, hv, sizeof hv,
			NULL, 0);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e sign:           %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255e_verify(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255e verify", 19);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char hv[32];
		do255e_private_key sk;
		do255e_public_key pk;
		do255e_signature sig;

		do255e_keygen(&rng, &sk, &pk);
		shake_extract(&rng, hv, sizeof hv);
		do255e_sign(&sig, &sk, &pk, DO255_OID_SHA3_256, hv, sizeof hv,
			NULL, 0);
		begin = core_cycles();
		do255e_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, hv, sizeof hv);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255e verify:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_keygen(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255s keygen", 19);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char seed[32];
		shake_context sc;
		do255s_private_key sk;
		do255s_public_key pk;

		shake_extract(&rng, seed, sizeof seed);
		shake_init(&sc, 256);
		shake_inject(&sc, seed, sizeof seed);
		shake_flip(&sc);
		shake_extract(&sc, seed, 8);
		begin = core_cycles();
		do255s_keygen(&sc, &sk, &pk);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s keygen:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_ecdh(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255s ecdh", 17);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char sec[32];
		do255s_private_key sk1, sk2;
		do255s_public_key pk1, pk2;

		do255s_keygen(&rng, &sk1, &pk1);
		do255s_keygen(&rng, &sk2, &pk2);
		begin = core_cycles();
		do255s_key_exchange(sec, sizeof sec, &sk1, &pk1, &pk2);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s key exchange:   %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_sign(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255s sign", 17);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char hv[32];
		do255s_private_key sk;
		do255s_public_key pk;
		do255s_signature sig;

		do255s_keygen(&rng, &sk, &pk);
		shake_extract(&rng, hv, sizeof hv);
		begin = core_cycles();
		do255s_sign(&sig, &sk, &pk, DO255_OID_SHA3_256, hv, sizeof hv,
			NULL, 0);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s sign:           %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

static void
speed_do255s_verify(void)
{
	size_t u;
	uint64_t tt[1000];
	shake_context rng;

	shake_init(&rng, 128);
	shake_inject(&rng, "speed do255s verify", 19);
	shake_flip(&rng);
	for (u = 0; u < 2000; u ++) {
		uint64_t begin, end;
		unsigned char hv[32];
		do255s_private_key sk;
		do255s_public_key pk;
		do255s_signature sig;

		do255s_keygen(&rng, &sk, &pk);
		shake_extract(&rng, hv, sizeof hv);
		do255s_sign(&sig, &sk, &pk, DO255_OID_SHA3_256, hv, sizeof hv,
			NULL, 0);
		begin = core_cycles();
		do255s_verify_vartime(&sig, &pk,
			DO255_OID_SHA3_256, hv, sizeof hv);
		end = core_cycles();
		if (u >= 1000) {
			tt[u - 1000] = end - begin;
		}
	}
	qsort(tt, (sizeof tt) / sizeof tt[0], sizeof tt[0], &cmp_u64);
	printf("do255s verify:         %9lu (%lu .. %lu)\n",
		(unsigned long)tt[500],
		(unsigned long)tt[100],
		(unsigned long)tt[900]);
	fflush(stdout);
}

#endif

int
main(void)
{
	test_do255e_scalar();
	test_do255s_scalar();
	test_do255e_decode();
	test_do255s_decode();
	test_do255e_map_to_curve();
	test_do255s_map_to_curve();
	test_do255e_add();
	test_do255s_add();
	test_do255e_mul();
	test_do255s_mul();
	test_do255e_verify_helper();
	test_do255s_verify_helper();
	test_do255e_keygen();
	test_do255e_ecdh();
	test_do255e_sign();
	test_do255s_keygen();
	test_do255s_ecdh();
	test_do255s_sign();
#if DO_BENCH86
	speed_do255e_decode();
	speed_do255s_decode();
	speed_do255e_encode();
	speed_do255s_encode();
	speed_do255e_mul();
	speed_do255s_mul();
	speed_do255e_mulgen();
	speed_do255s_mulgen();
	speed_do255e_verify_helper();
	speed_do255s_verify_helper();
	printf("\n");
	speed_do255e_keygen();
	speed_do255s_keygen();
	speed_do255e_ecdh();
	speed_do255s_ecdh();
	speed_do255e_sign();
	speed_do255s_sign();
	speed_do255e_verify();
	speed_do255s_verify();
#endif
	return 0;
}
