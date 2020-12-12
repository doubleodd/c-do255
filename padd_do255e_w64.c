/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations, including gf_sqrt() or gf_issquare()
 *
 * This file is for all implementations of do255e that use 64-bit limbs.
 * It defines:
 *  - do255e_neutral
 *  - do255e_generator
 *  - do255e_decode()
 *  - do255e_encode()
 *  - do255e_is_neutral()
 *  - do255e_eq()
 *  - do255e_add()
 *  - do255e_double()
 *  - do255e_double_x()
 */

/* a */
static const gf CURVE_A = {
	0, 0, 0, 0
};

/* 4*b */
static const gf CURVE_4B = {
	0xFFFFFFFFFFFFB71D,
	0xFFFFFFFFFFFFFFFF,
	0xFFFFFFFFFFFFFFFF,
	0x7FFFFFFFFFFFFFFF
};

/* see do255.h */
const do255e_point do255e_generator = {
	{ { 2, 0, 0, 0 } },
	{ { 1, 0, 0, 0 } },
	{ { 1, 0, 0, 0 } }
};

/*
 * Custom structure for a point in affine coordinates:
 *  - if X = 0, then this is the neutral (W is ignored);
 *  - if X != 0, then coordinate Z is implicitly equal to 1.
 */
typedef struct {
	do255_int256 X, W;
} do255e_point_affine;

/*
 * We get do255e_neutral, do255e_decode(), do255e_encode(),
 * do255e_is_neutral() and do255e_eq() from pcore_w64.c.
 */
#include "pcore_w64.c"

/* see do255.h */
void
do255e_add(do255e_point *P3,
	const do255e_point *P1, const do255e_point *P2)
{
	gf t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, X3, W3, Z3;
	uint64_t fz1, fz2;

	/*
	 * Test whether P1 and/or P2 is neutral.
	 */
	fz1 = gf_iszero(&P1->Z.w64);
	fz2 = gf_iszero(&P2->Z.w64);

	/* t1 <- Z1^2 */
	gf_sqr_inline(&t1, &P1->Z.w64);

	/* t2 <- Z2^2 */
	gf_sqr_inline(&t2, &P2->Z.w64);

	/* t3 <- ((Z1 + Z2)^2 - t1 - t2)/2 */
	gf_add(&t3, &P1->Z.w64, &P2->Z.w64);
	gf_sqr_inline(&t3, &t3);
	gf_sub2(&t3, &t3, &t1, &t2);
	gf_half(&t3, &t3);

	/* t4 <- t3^2 */
	gf_sqr_inline(&t4, &t3);

	/* t5 <- W1*W2 */
	gf_mul_inline(&t5, &P1->W.w64, &P2->W.w64);

	/* t6 <- X1*X2 */
	gf_mul_inline(&t6, &P1->X.w64, &P2->X.w64);

	/* t7 <- (W1 + Z1)*(W2 + Z2) - t3 - t5 */
	gf_add(&t7, &P1->W.w64, &P1->Z.w64);
	gf_add(&t8, &P2->W.w64, &P2->Z.w64);
	gf_mul_inline(&t7, &t7, &t8);
	gf_sub2(&t7, &t7, &t3, &t5);

	/* t8 <- (X1 + t1)*(X2 + t2) - t4 - t6 */
	gf_add(&t8, &P1->X.w64, &t1);
	gf_add(&t9, &P2->X.w64, &t2);
	gf_mul_inline(&t8, &t8, &t9);
	gf_sub2(&t8, &t8, &t4, &t6);

	/* Z3 <- (t6 - b*t4)*t7
	   Also, replace t4 with -b*t4 */
	gf_mul2(&t4, &t4);
	gf_add(&t9, &t6, &t4);
	gf_mul_inline(&Z3, &t9, &t7);

	/* t9 <- t7^4 */
	gf_sqr_inline(&t9, &t7);
	gf_sqr_inline(&t9, &t9);

	/* X3 <- b*t6*t9 */
	gf_mul_inline(&X3, &t6, &t9);
	gf_mul2(&X3, &X3);
	gf_neg(&X3, &X3);

	/* t10 <- (t5 + a*t3)*(t6 + b*t4)
	   a = 0
	   -b*t4 was already computed (in t4)
	   We overwrite t5 and t6, which we won't need anymore */
	gf_sub(&t6, &t6, &t4);
	gf_mul_inline(&t10, &t5, &t6);

	/* W3 <- -t10 - 2*b*t3*t8
	   b = 1/2, hence 2*b = 1.
	   We overwrite t8. */
	gf_mul_inline(&t8, &t3, &t8);
	gf_mul4(&t8, &t8);
	gf_sub(&W3, &t8, &t10);

	/*
	 * If P1 is neutral, replace P3 with P2.
	 * If P2 is neutral, replace P3 with P1.
	 */
	gf_sel3(&P3->X.w64, &P2->X.w64, &P1->X.w64, &X3, fz1, fz2);
	gf_sel3(&P3->W.w64, &P2->W.w64, &P1->W.w64, &W3, fz1, fz2);
	gf_sel3(&P3->Z.w64, &P2->Z.w64, &P1->Z.w64, &Z3, fz1, fz2);
}

/*
 * Point addition, with the second point in affine coordinates.
 */
UNUSED
static void
do255e_add_mixed(do255e_point *P3,
	const do255e_point *P1, const do255e_point_affine *P2)
{
	gf t1, t5, t6, t7, t8, t9, t10, X3, W3, Z3;
	uint64_t fz1, fz2;

	/*
	 * Test whether P1 and/or P2 is neutral.
	 */
	fz1 = gf_iszero(&P1->Z.w64);
	fz2 = gf_iszero(&P2->X.w64);

	/* t1 <- Z1^2 */
	gf_sqr_inline(&t1, &P1->Z.w64);

	/* t2 = 1 */
	/* t3 = Z1 */
	/* t4 = t1 */

	/* t5 <- W1*W2 */
	gf_mul_inline(&t5, &P1->W.w64, &P2->W.w64);

	/* t6 <- X1*X2 */
	gf_mul_inline(&t6, &P1->X.w64, &P2->X.w64);

	/* t7 <- W1 + W2*Z1 */
	gf_mul_inline(&t7, &P1->Z.w64, &P2->W.w64);
	gf_add(&t7, &t7, &P1->W.w64);

	/* t8 <- X1 + X2*t1 */
	gf_mul_inline(&t8, &t1, &P2->X.w64);
	gf_add(&t8, &t8, &P1->X.w64);

	/* Z3 <- (t6 - b*t1)*t7
	   Also, replace t1 with -b*t1 */
	gf_mul2(&t1, &t1);
	gf_add(&t9, &t6, &t1);
	gf_mul_inline(&Z3, &t9, &t7);

	/* t9 <- t7^4 */
	gf_sqr_inline(&t9, &t7);
	gf_sqr_inline(&t9, &t9);

	/* X3 <- b*t6*t9 */
	gf_mul_inline(&X3, &t6, &t9);
	gf_mul2(&X3, &X3);
	gf_neg(&X3, &X3);

	/* t10 <- (t5 + a*t3)*(t6 + b*t1)
	   a = 0
	   -b*t1 was already computed (in t1)
	   We overwrite t5 and t6, which we won't need anymore */
	gf_sub(&t6, &t6, &t1);
	gf_mul_inline(&t10, &t5, &t6);

	/* W3 <- -t10 - 2*b*t3*t8
	   b = -2, hence 2*b = -4.
	   We overwrite t8. */
	gf_mul_inline(&t8, &P1->Z.w64, &t8);
	gf_mul4(&t8, &t8);
	gf_sub(&W3, &t8, &t10);

	/*
	 * If P1 is neutral, replace P3 with P2.
	 * If P2 is neutral, replace P3 with P1.
	 * If both are neutral, then we want to use P1 as source, whose
	 * Z coordinate is then zero; this allows us to assume that
	 * P2.Z = 1 here.
	 */
	gf_sel3(&P3->X.w64, &P1->X.w64, &P2->X.w64, &X3, fz2, fz1);
	gf_sel3(&P3->W.w64, &P1->W.w64, &P2->W.w64, &W3, fz2, fz1);
	gf_sel3(&P3->Z.w64, &P1->Z.w64, &GF_ONE, &Z3, fz2, fz1);
}

/* see do255.h */
void
do255e_double(do255e_point *P3, const do255e_point *P1)
{
	gf tX, tW, tZ, t4, t5;

	/*
	 * X' = W^4         (in tX)
	 * W' = W^2 - 2*X   (in tW)
	 * Z' = W*Z         (in tZ)
	 */
	gf_sqr_inline(&tX, &P1->W.w64);
	gf_mul2(&tW, &P1->X.w64);
	gf_mul_inline(&tZ, &P1->W.w64, &P1->Z.w64);
	gf_sub(&tW, &tX, &tW);
	gf_sqr_inline(&tX, &tX);

	/*
	 * X3 = 4*b*Z'^4 = -8*Z'^4 (because b = -2)
	 * W3 = X' - (1/2)*W'^2
	 * Z3 = W'*Z'
	 *
	 * We could compute Z3 = (1/2)*((W'+Z')^2 - W'^2 - Z'^2) but
	 * on a Coffe Lake CPU this seems to be slower than a mul.
	 */
	gf_sqr_inline(&t4, &tZ);
	gf_sqr_inline(&t5, &tW);
	gf_mul_inline(&P3->Z.w64, &tW, &tZ);
	gf_sqr_inline(&t4, &t4);
	gf_half(&t5, &t5);
	gf_sub(&P3->W.w64, &tX, &t5);
	gf_mul8(&t4, &t4);
	gf_neg(&P3->X.w64, &t4);
}

/* see do255.h */
void
do255e_double_x(do255e_point *P3, const do255e_point *P1, unsigned n)
{
	gf tX, tW, tZ, t1, t2;

	if (n == 0) {
		*P3 = *P1;
		return;
	}
	tX = P1->X.w64;
	tW = P1->W.w64;
	tZ = P1->Z.w64;
	while (n -- > 1) {
		/*
		 * Here we use 2M+4S per iteration. We could change that
		 * to 1M+5S by noticing that:
		 *   Z'' = 2*W*W'*Z
		 *       = ((W+W')^2 - W^2 - W'^2)*Z
		 * but on a big x86 (Coffee Lake core) this happens to be
		 * somewhat slower, because of the extra subtractions
		 * and/or less favourable dependencies.
		 */

		/*
		 * X' = W^4
		 * W' = W^2 - 2*X
		 * Z' = W*Z
		 */
		gf_sqr_inline(&t1, &tW);
		gf_mul2(&t2, &tX);
		gf_mul_inline(&tZ, &tW, &tZ);
		gf_sub(&tW, &t1, &t2);
		gf_sqr_inline(&tX, &t1);

		/*
		 * X'' = W'^4
		 * W'' = W'^2 - 2*X'
		 * Z'' = 2*W'*Z'
		 */
		gf_sqr_inline(&t1, &tW);
		gf_mul2(&t2, &tX);
		gf_mul_inline(&tZ, &tW, &tZ);
		gf_sub(&tW, &t1, &t2);
		gf_sqr_inline(&tX, &t1);
		gf_mul2(&tZ, &tZ);
	}

	/*
	 * X' = W^4
	 * W' = W^2 - 2*X
	 * Z' = W*Z
	 */
	gf_sqr_inline(&t1, &tW);
	gf_mul2(&t2, &tX);
	gf_mul_inline(&tZ, &tW, &tZ);
	gf_sub(&tW, &t1, &t2);
	gf_sqr_inline(&tX, &t1);

	/*
	 * X3 = 4*b*Z'^4 = -8*Z'^4 (because b = -2)
	 * W3 = X' - (1/2)*W'^2
	 * Z3 = W'*Z'
	 *
	 * We could compute Z3 = (1/2)*((W'+Z')^2 - W'^2 - Z'^2) but
	 * on a Coffe Lake CPU this seems to be slower than a mul.
	 */
	gf_sqr_inline(&t1, &tZ);
	gf_sqr_inline(&t2, &tW);
	gf_mul_inline(&P3->Z.w64, &tW, &tZ);
	gf_sqr_inline(&t1, &t1);
	gf_half(&t2, &t2);
	gf_sub(&P3->W.w64, &tX, &t2);
	gf_mul8(&t1, &t1);
	gf_neg(&P3->X.w64, &t1);
}

/*
 * Custom structures for points in (x,u) coordinates (fractional and
 * affine).
 */
typedef struct {
	do255_int256 X, Z, U, T;
} do255e_point_xu;
typedef struct {
	do255_int256 X, U;
} do255e_point_affine_xu;

static const do255e_point_affine_xu window_G_xu[] = {
	/* 1 */
	{
		{ { 0x0000000000000002, 0x0000000000000000,
		    0x0000000000000000, 0x0000000000000000 } },
		{ { 0x0000000000000001, 0x0000000000000000,
		    0x0000000000000000, 0x0000000000000000 } }
	},
	/* 2 */
	{
		{ { 0xE38E38E38E38AAE3, 0x8E38E38E38E38E38,
		    0x38E38E38E38E38E3, 0x638E38E38E38E38E } },
		{ { 0xB6DB6DB6DB6D97A3, 0xDB6DB6DB6DB6DB6D,
		    0x6DB6DB6DB6DB6DB6, 0x36DB6DB6DB6DB6DB } }
	},
	/* 3 */
	{
		{ { 0x0000000000000152, 0x0000000000000000,
		    0x0000000000000000, 0x0000000000000000 } },
		{ { 0xC2F21347C4043E79, 0x6B1CEBA6066D4156,
		    0xAB617909A3E20224, 0x12358E75D30336A0 } }
	},
	/* 4 */
	{
		{ { 0xEA5E1BA07D1A36D1, 0xD0AF77073CCEA916,
		    0xC68A748FF2D92037, 0x43554439A8DE7571 } },
		{ { 0x65A29F71130DB4AD, 0x9F71130DFA47C8BB,
		    0x130DFA47C8BB65A2, 0x7A47C8BB65A29F71 } }
	},
	/* 5 */
	{
		{ { 0x9B7D88CD74D7D3CA, 0x0E31B461193896AC,
		    0x93464506E97D44DB, 0x0ABC8AC61DCD9949 } },
		{ { 0x1F2B6B08DA5B43EE, 0xE40F8B8BC44A0C63,
		    0x5866F1F8B35FB70C, 0x185034D250F768D7 } }
	},
	/* 6 */
	{
		{ { 0x19348B5724DF12CF, 0xAB7572F66232945A,
		    0xFFAA89042A63BF4B, 0x3024ED85633EF2E2 } },
		{ { 0x0BD0C5F1F91D6B18, 0xBB4A410D263610A7,
		    0xA1AB0B9D98F35F00, 0x4FA6D8B6AFDDC92B } }
	},
	/* 7 */
	{
		{ { 0x33E381251F43C3D5, 0x49BF0E2E71C6FE8F,
		    0x6AF69CC116BAEF18, 0x36199FBAD0C9585E } },
		{ { 0x7EB52414159EF4EA, 0xB885C9D1EB4CC9E1,
		    0x350914B3EE64BF7F, 0x6DD8CDFA520AED5A } }
	},
	/* 8 */
	{
		{ { 0xDDBCA65ECDCBDFAE, 0xD19C744A47F38FBD,
		    0x15099FB55653455E, 0x1FD0DC2A961EC01F } },
		{ { 0x99DA8C93EB513A8B, 0x0706B8B95DEDFC87,
		    0xC54D8F471F778CE9, 0x4766315BFA2E63E5 } }
	},
	/* 9 */
	{
		{ { 0x8DF7F2C9CE799A9C, 0x7AC7F7C3AFDD04F9,
		    0x915FE4A27D833740, 0x1ED67871986F29BA } },
		{ { 0xA84A27A9D0A08E61, 0x27E9084D132CCAC1,
		    0x498C7D8B01F68C40, 0x6957FDFF940E4159 } }
	},
	/* 10 */
	{
		{ { 0xE944946059A9387C, 0x32BB0464A75287B8,
		    0x122E571F46C8845D, 0x0D05AC0126E0A481 } },
		{ { 0x3AA366BBB889903E, 0x55838146CC140A37,
		    0x4AA37581A9B6AD5E, 0x7B37113C916F803C } }
	},
	/* 11 */
	{
		{ { 0x6562064C442E3709, 0xD013EB4D114A7267,
		    0x166892C716D5320A, 0x2824BCCA3B493396 } },
		{ { 0xA9A8911D864E7F82, 0x65CF6B9CAB741725,
		    0x8C133221E772B327, 0x158521078CD1F209 } }
	},
	/* 12 */
	{
		{ { 0x4AE1AE2876FFE733, 0x55A43A11F9D28845,
		    0xBAACD8A3E4990483, 0x37B39256440F5C21 } },
		{ { 0xE4C1C725087640AA, 0xFB902D6A3EF5D5E0,
		    0x53EF35932E1297EB, 0x67E65CF7E1787343 } }
	},
	/* 13 */
	{
		{ { 0x90DF642868789634, 0x267A28B9CB72C6CA,
		    0x27BE4B2B937625B5, 0x62003971A89B844F } },
		{ { 0x90F8839881061965, 0x67D0394FF2BFCB98,
		    0x913200FCCD1396D8, 0x17F96D76306A3580 } }
	},
	/* 14 */
	{
		{ { 0xC76258B805A821DA, 0xDC3C29F024B5765F,
		    0xB646CAAD30897EED, 0x46F594DEF8D35CB9 } },
		{ { 0x05B49673D2AC4172, 0xA016A6890D77E4E6,
		    0x7C6DAA970635E1C0, 0x42C8034547A6A04A } }
	},
	/* 15 */
	{
		{ { 0x496055CDC3DFB745, 0x0673F992F547B770,
		    0x8EFAE8F99B3E5BB3, 0x33C76A13E12C07DD } },
		{ { 0x7FFA4AF719120727, 0x705D12571BF74984,
		    0x4AD1FA649FAE1F07, 0x2F4CA2B6265D7456 } }
	},
	/* 16 */
	{
		{ { 0xE1DAC53644F243F9, 0xE1154ECDDFFD59A1,
		    0x1731585F8B8C6649, 0x1BFB93F1365D2CAF } },
		{ { 0x2D808316E1227049, 0x15064C9132683177,
		    0x706D8A1F41E90ED8, 0x251A19311A6DB76E } }
	}
};
static const do255e_point_affine_xu window_G65_xu[] = {
	/* 1 */
	{
		{ { 0xABEF504D87FDEB41, 0x3A2D867D250A6B59,
		    0x5906F21AAB3E16A0, 0x58327D52081B8A67 } },
		{ { 0xAD5F8FBFC596FA71, 0x415893549DE223FF,
		    0x395D2181E50A4384, 0x1B313D36A8A7626E } }
	},
	/* 2 */
	{
		{ { 0x15356A7D39689FC6, 0xF52A5DE2A0967CE6,
		    0xAD2738B8D02C8707, 0x75B28A23988AA077 } },
		{ { 0x2B7C2F71889E3F33, 0x4CA4C049A5E65CF4,
		    0x5CD27D909976BFE7, 0x0BE56F359985D602 } }
	},
	/* 3 */
	{
		{ { 0xDC3DE207745F399F, 0xD7C5D084CA47E2C5,
		    0xDFA47415EF678E9B, 0x5F4DD2E4049C2FE8 } },
		{ { 0x4A504A5DED61CB7F, 0xAF41508342D7801D,
		    0x0519A68AAB4295EB, 0x098D3AB90B09C2B4 } }
	},
	/* 4 */
	{
		{ { 0x977017D6C9BF2E76, 0x17CE586894D3EC82,
		    0xA07F1106A288282D, 0x1334610CA14E17A8 } },
		{ { 0xD98BB32E81B4D89F, 0x5A0FCD3AE1067F08,
		    0x3B845EF3AFAD3191, 0x3540EB32AA6E2C23 } }
	},
	/* 5 */
	{
		{ { 0x592F8CC3E5028371, 0x45BC1E36D117EB2E,
		    0x531CFBF930BEEB8F, 0x2829EA75A68A6C5E } },
		{ { 0xDB2EB15DD80EB7AB, 0x695C5863633AC0BB,
		    0x9DD36D925FB45810, 0x4B36B3BE374E74CB } }
	},
	/* 6 */
	{
		{ { 0x94156556D9FB736E, 0x3FCA76C281E2244C,
		    0x7A11BA953BA95F28, 0x2E6BCF6A3DD91BBF } },
		{ { 0x2473D6DB1425C001, 0x763067186E58108F,
		    0x0C0776C9F6ED32DE, 0x143CF8A090ACB085 } }
	},
	/* 7 */
	{
		{ { 0xE24BA37DA32EFB19, 0xB319B89EBE6E4852,
		    0xB99C333F17938861, 0x053979BEC7A8BA46 } },
		{ { 0xE95E218127D26209, 0x33088969E8FB612B,
		    0xBAC06821F2BF788C, 0x7EE7B8C1D61C83DD } }
	},
	/* 8 */
	{
		{ { 0x6B70CE277F74345E, 0x5FB603F5B5D2FF08,
		    0x9A28F63978F82349, 0x1E7BF7585E894EE6 } },
		{ { 0x4C033B2301944CB8, 0x3260DA08A5D337CB,
		    0x34DFEFB2A6C39FFC, 0x235872299B5B142E } }
	},
	/* 9 */
	{
		{ { 0xEDC150B81CDB24C5, 0xCB26D5421BE70FF6,
		    0xF19BB5BF01AA63A4, 0x316C4EA24603AF0D } },
		{ { 0x6742B67C6086DF92, 0x53193A718D75331F,
		    0xADD356D1CA14E352, 0x07DF0CC5D7C9E88D } }
	},
	/* 10 */
	{
		{ { 0xE36D90729B572DFF, 0xC4F6284A7B9548E7,
		    0x0FDCACBCFEDA82E0, 0x40A3C68531C84199 } },
		{ { 0x7982C54051C4DF08, 0x476A1D01EBBFE2BB,
		    0x3F6064CC7E287D25, 0x79E82E420C17C457 } }
	},
	/* 11 */
	{
		{ { 0xC10298FEAF5791A3, 0x36AD705BB82FDC4C,
		    0xA8764B80FAC163BB, 0x1388EA98CCAF5744 } },
		{ { 0x00FF128613403E58, 0x854CBABF2D02F041,
		    0xB97069DD65593890, 0x3CC8BC3168B5A376 } }
	},
	/* 12 */
	{
		{ { 0x230E07753775607D, 0xDA9A2205BF1C2AE5,
		    0x975954E5C2DD10BB, 0x28652C57B5420892 } },
		{ { 0x3688782C6588D284, 0xA50B5985163B51DB,
		    0x23447D42F7311FA5, 0x51CBB2B65B751D79 } }
	},
	/* 13 */
	{
		{ { 0x4CBEDE49F7C4233B, 0x4EC8BC06D529F82D,
		    0x75FB83AA2630BCB9, 0x7F3F2E391D0902B6 } },
		{ { 0x6AF277E67A9FAA68, 0xFDEE555C9009703F,
		    0xAA36C89F8D5B502A, 0x052051EDAFC87131 } }
	},
	/* 14 */
	{
		{ { 0x2C9FA67ECEC074F3, 0x077855CE1E217BAC,
		    0xAE612A7D92BCDFF8, 0x7B5805F9F2BC9375 } },
		{ { 0x86A494FC08EB56DD, 0xF694DF6ED62219A0,
		    0x214282B022098B86, 0x7182E92FA95620EB } }
	},
	/* 15 */
	{
		{ { 0xD48702C2FB7A4ACE, 0x9009655147735647,
		    0xD98938B066322FCE, 0x08B7C64434681715 } },
		{ { 0xD25E4076E71417FD, 0xAE53B9A6D617E037,
		    0xEDF6131ABB08D620, 0x406E88577206C6C4 } }
	},
	/* 16 */
	{
		{ { 0xD86B084C27BD12CA, 0x772D4452961D7D94,
		    0x109ACCAAED9EC998, 0x575C5B51EF95AAE4 } },
		{ { 0x7D37BBB6D0781C73, 0x1ADF4CFE82DAF6AB,
		    0x6CB90D1DDFBA601C, 0x7B916A823C090D35 } }
	}
};
static const do255e_point_affine_xu window_G130_xu[] = {
	/* 1 */
	{
		{ { 0xC796B2402F4B5577, 0xB848DE143824131E,
		    0x35479390B17AF31F, 0x0EB9DD4F3DFA2F98 } },
		{ { 0x6A6F8767CAF58BCB, 0xEF9A2E5ACD520CD9,
		    0x2B998E19EE40437C, 0x1E3A7692F3E02AB1 } }
	},
	/* 2 */
	{
		{ { 0x9B55C0CC7A0BD2F8, 0x73D4BA09593578E5,
		    0x89A107F02B82E751, 0x23C4CFE49D8C1DB2 } },
		{ { 0x97643EA7C28259E8, 0x64C33BBEA0416456,
		    0xAC5EBA85AFFFBCEB, 0x1DE0359F8936CEEA } }
	},
	/* 3 */
	{
		{ { 0x29BC60641F6F9439, 0xA4470A0D45E0AB6F,
		    0x3CF550B4CBA43D48, 0x098DBD0C76126C5C } },
		{ { 0x4A8F6858185A63C2, 0xC29D0227105E6338,
		    0x4FA122A357313D72, 0x62BAF3EF3C842009 } }
	},
	/* 4 */
	{
		{ { 0xE5AFFE37730A6B2D, 0x4AD5E971E27ED8D1,
		    0xCD5F1FF510C3D246, 0x05D71BC1B236B80B } },
		{ { 0x17625F88FFC35397, 0xEF697901DA783099,
		    0xDCCCD3E1EF458EDC, 0x37EDC360CFBBEDE2 } }
	},
	/* 5 */
	{
		{ { 0x15EE9577D26B0A12, 0x796D2F22498E2395,
		    0xEEBF418FA5DC2FFF, 0x6EFC71DFB61B892C } },
		{ { 0x055D6CCDE9EC95DF, 0xEB86F24ADFCADBF5,
		    0x6DDAC4AB9F7AE17C, 0x282E209409ADD692 } }
	},
	/* 6 */
	{
		{ { 0x1D2065C52CC05E06, 0xE136A646AD8CCE04,
		    0x9EEE9C6F51A1AD2C, 0x29BD9CB831800380 } },
		{ { 0xBA9DA77EDCAD7528, 0xB6FA24A8892C93FF,
		    0xF406AD1FFA1ACD6B, 0x7E4F46D5131E195C } }
	},
	/* 7 */
	{
		{ { 0x3764EFB34F88E075, 0x932413347C635C0D,
		    0xC9DA4905A98FA573, 0x3B23E4AB52D89A2F } },
		{ { 0x9E343B6FEBEF1413, 0xC1293137B7F95D67,
		    0x6FB5672D05CA0B5B, 0x372239482AEC172D } }
	},
	/* 8 */
	{
		{ { 0xC2D4B896B99ECE2D, 0x12CFC73515B85FA6,
		    0x68977BE4A000B16A, 0x1A7EDDB0A3FD9E15 } },
		{ { 0x21232B19EA446499, 0x0F110F3BEDDC580B,
		    0x7BFEE5167F928B59, 0x4E1CF6CCC48289F4 } }
	},
	/* 9 */
	{
		{ { 0x66DD5CBDC8683D52, 0x0E57F934D5717796,
		    0x06E7AFEFE940A0DE, 0x5E68DFB004FFD454 } },
		{ { 0xFFFA0532A28FF940, 0x06134380E791A9D0,
		    0x24120A75934E696A, 0x6F671B022C83BC57 } }
	},
	/* 10 */
	{
		{ { 0x382654D53F6219C8, 0x5B068A05F61C4345,
		    0xC106F1243C0212B6, 0x7C5492368C243A54 } },
		{ { 0x0EE1B5003B071E37, 0x3EC320448D80BCA7,
		    0x108D5A1C5EC7534A, 0x64ECAAC0EBA1A511 } }
	},
	/* 11 */
	{
		{ { 0xBFB4842BF8CB93CF, 0xF994337DE38CC098,
		    0xB9BB68C291517A22, 0x56B3E07FE03B3A5A } },
		{ { 0x6CC3477AD4AB70A1, 0x0BC92225712BA4D4,
		    0xD51C733D4564D8F4, 0x029CAFA5599B372A } }
	},
	/* 12 */
	{
		{ { 0xA0A8890960C40BEC, 0xCE8ECA1A58BBBBC9,
		    0x076AB31BA9BBF6DE, 0x0E1B23A02743D405 } },
		{ { 0xCB9E1D5CB7854025, 0x960067A5C064493E,
		    0x41F420B5A91ED717, 0x6E81F9FA7E661D8D } }
	},
	/* 13 */
	{
		{ { 0x56380CDCFD534861, 0xBF23AAE1C53EC995,
		    0xEA53DCFF88E21C0D, 0x44019A922EE00C9D } },
		{ { 0x1A722773E489DDB8, 0xF94983893D4AABD6,
		    0x59F3D4C5BB3DFDCC, 0x653EE371D2801E6A } }
	},
	/* 14 */
	{
		{ { 0x1AEE789343714BC2, 0x5027559B5A242FA0,
		    0xD35C383123797CC4, 0x7E1BA15EC794EA34 } },
		{ { 0x460E164D09693F50, 0x96FABF7744D22EC2,
		    0x216A1928595E868E, 0x50E1BEE9AC402680 } }
	},
	/* 15 */
	{
		{ { 0x0C6BCDDDE880E3FE, 0x9D52DF1A80CA1D01,
		    0x6D28A55FD9814394, 0x07CD186F1749AC3B } },
		{ { 0x442C914C64C6EE61, 0x5486463AAA3D41CD,
		    0x2323BA05744FB271, 0x5CE94782B63D2983 } }
	},
	/* 16 */
	{
		{ { 0xF633589461CE1D8F, 0x43656EE4CD988663,
		    0x601D9906BC58B752, 0x585EFD1D9E157C20 } },
		{ { 0x20DDE2D9560BD063, 0x68337F979386B815,
		    0x9CAE33A6B5F9B94C, 0x0F2ED8418B17674E } }
	}
};
static const do255e_point_affine_xu window_G195_xu[] = {
	/* 1 */
	{
		{ { 0xC9D74325CBE23870, 0x3C8F799792EB5981,
		    0xB4D3EBF20E4ACC9D, 0x5D505019F9CD4639 } },
		{ { 0xC6C4EA52864B8022, 0x3FACF03027F2FE05,
		    0x5A78F8FDAFE0F2B2, 0x7A2058682117A352 } }
	},
	/* 2 */
	{
		{ { 0xB5FFEA9757F8DCE6, 0x992EAADD0950F49E,
		    0xC30EE566764DB296, 0x77BEE9FA736A26DC } },
		{ { 0xE7A9C455FDCCE69F, 0xB043C24E23C52866,
		    0xCBC1DD8A3179B032, 0x597FE7EC4E366C38 } }
	},
	/* 3 */
	{
		{ { 0xF4762F3FF4FB8115, 0x60D22515C8F29371,
		    0xE64A746FA6B9C81A, 0x107BB7D6EED2E10E } },
		{ { 0xEBF1C192782AD7E7, 0xDAC867CE0228990B,
		    0xB0C0AEB839C2A9BB, 0x5D529C2B2E3222F2 } }
	},
	/* 4 */
	{
		{ { 0xA0C19926FC389D68, 0xD93A4C18F4C2CCCD,
		    0x09EB7E080DF8E02C, 0x5CD950BF71F691C5 } },
		{ { 0xE5AD1FDA050CE7CF, 0xD5179DCB398FE9E2,
		    0x880F0F9CA2B23DE9, 0x73E9DA1D7C583AB6 } }
	},
	/* 5 */
	{
		{ { 0x97A1BBD1E5AA8841, 0x09B24BDC0DC8FFB6,
		    0xA57657C8DCE4DE79, 0x228ECFC1B5307822 } },
		{ { 0x81533FD0CED00FE0, 0x2B41B323457375B0,
		    0x3428954D0B0B6412, 0x3FB05C6B656FCDE7 } }
	},
	/* 6 */
	{
		{ { 0xF03AEB994C6B5021, 0xB0156B2F0414CE7E,
		    0x64B75C8B8346FCE3, 0x4F54E4B6B9C3FD25 } },
		{ { 0xABECE8DEAA4DEFF3, 0xA6B25F5370FF8BED,
		    0xC70C1F018B95875D, 0x1EEE7F4380019FB8 } }
	},
	/* 7 */
	{
		{ { 0x9854F37D39978A1E, 0xA8401C2863D1E85D,
		    0x021EDD635FDF6914, 0x317884D08D246053 } },
		{ { 0xC2513D53CE5A6CD2, 0x8AF5B5BD4C9ADB58,
		    0xDF748C1856292D78, 0x1C54D437C147EB47 } }
	},
	/* 8 */
	{
		{ { 0x09F1B77F26BC9F8E, 0x219F33E838C90E61,
		    0x320CDCEF213ECF00, 0x67131909DEA4A881 } },
		{ { 0x961441BFA4853698, 0x76396E28425429D3,
		    0x23187D9C49399AB4, 0x47EC89341C754A72 } }
	},
	/* 9 */
	{
		{ { 0x1BE2A51DB2720321, 0x645B8A7DD6376B36,
		    0x686F85695A8133B5, 0x07E6F607B6D91397 } },
		{ { 0x790DDCB656DA5EBB, 0x899CA30F8A6C1157,
		    0xB055E943E160FF52, 0x0C4E4B67E97A3F02 } }
	},
	/* 10 */
	{
		{ { 0x5940A4043221C587, 0x5FA59E201799740F,
		    0x73197D2BCF61AC84, 0x303417596A2CF352 } },
		{ { 0x7CBAD6A204D1F6B9, 0xEB725D50999A4399,
		    0x7C80807D104B0670, 0x44B8942C5DF07889 } }
	},
	/* 11 */
	{
		{ { 0x97C6A80A64F51D16, 0xEAC3C45ABBB9912B,
		    0x72FB0626EAAD16C4, 0x5DA3C2E5773227DE } },
		{ { 0x30FE96716E7DF796, 0xB6A214969A98317F,
		    0xEB5423DB7543C3F6, 0x7DD81F0AD475BB65 } }
	},
	/* 12 */
	{
		{ { 0x0FE33479735C7A13, 0xA6F8DF8C4AC20F15,
		    0xD19D08B9A74DDF08, 0x7D923FBD990B2A82 } },
		{ { 0x041881DC4BB15593, 0x0620628690F070A8,
		    0xF6647FFFA6239BFD, 0x60406418E0A8E484 } }
	},
	/* 13 */
	{
		{ { 0x9D7AE8C91CE6917B, 0xB5DB6A7F7DEA58CC,
		    0x9B1C372C39EEA275, 0x58742A848CC4EB79 } },
		{ { 0x02A5C58CE250DF40, 0x80A9A62960313C1E,
		    0xEDC81102A9A286D9, 0x03C8B0610E5DE932 } }
	},
	/* 14 */
	{
		{ { 0xDD44633AC9F7D040, 0x2632460F3A278E77,
		    0x092C7F41C75E959E, 0x00D4C6E23ED1D27A } },
		{ { 0x2359C265188EE74D, 0x3498D65BDB7611FC,
		    0x97D9FFD6286A6BD1, 0x63F775224E36165F } }
	},
	/* 15 */
	{
		{ { 0x515E61CB24DAF0E2, 0x606D332C7B076CF7,
		    0x4792A1E865D47BC0, 0x5A839AAC1191AB83 } },
		{ { 0x14BD04EFD34FC573, 0x58F89E267B42BEA3,
		    0xDD4D47FA083CC9BF, 0x5C69FCC38DA29629 } }
	},
	/* 16 */
	{
		{ { 0xF04AA6935B6381F5, 0xA05C72D30E2E31F2,
		    0x5414DCBFD3642F59, 0x75A8FE2604285EC3 } },
		{ { 0x805E028481B3D10D, 0x3E6A069F39FCEFCB,
		    0xDA636B907FED771B, 0x162581D9B675A4E1 } }
	}
};

/*
 * Doubling in fractional (x,u) coordinates.
 */
UNUSED
static void
do255e_double_xu(do255e_point_xu *P3, const do255e_point_xu *P1)
{
	gf tX, tW, tZ, t1, t2;

	/*
	 * First half-doubling, combined with conversion from fractional (x,u)
	 * to Jacobian (x,w); output in E[r](-2*a,a^2-4*b).
	 *   X' = Z^2*T^4
	 *   W' = (2*X + a*Z)*U^2 - Z*T^2
	 *   Z' = Z*U*T
	 * Note that a = 0 for curve do255e.
	 * Cost: 4M+2S
	 */
	gf_sqr_inline(&tW, &P1->U.w64);              /* tW <- U^2 */
	gf_mul_inline(&t1, &P1->Z.w64, &P1->T.w64);  /* t1 <- Z*T */
	gf_mul2(&tW, &tW);                           /* tW <- 2*U^2 */
	gf_mul_inline(&t2, &t1, &P1->T.w64);         /* t2 <- Z*T^2 */
	gf_mul_inline(&tW, &tW, &P1->X.w64);         /* tW <- 2*X*U^2 */
	gf_sqr_inline(&tX, &t2);                     /* tX <- Z^2*T^4 */
	gf_sub(&tW, &t2, &tW);                       /* tW <- Z*T^2 - 2*X*U^2 */
	gf_mul_inline(&tZ, &t1, &P1->U.w64);         /* tZ <- Z*U*T */

	/*
	 * Second half-doubling, combined with conversion back to
	 * fractional (x,u) coordinates.
	 *   X' = 4*b*Z^2
	 *   Z' = W^2
	 *   U' = 2*W*Z
	 *   T' = 2*X - 2*a*Z^2 - W^2
	 * Note that a = 0 and b = -2 for curve do255e.
	 * Cost: 3S  (with 2*W*Z = (W+Z)^2 - W^2 - Z^2)
	 */
	gf_sqr_inline(&t2, &tW);             /* t2 <- W^2 */
	gf_add(&t1, &tW, &tZ);               /* t1 <- W + Z */
	gf_sqr_inline(&tZ, &tZ);             /* tZ <- Z^2 */
	gf_sqr_inline(&t1, &t1);             /* t1 <- (W + Z)^2 */
	gf_mul2(&tX, &tX);                   /* tX <- 2*X */
	gf_sub2(&P3->U.w64, &t1, &t2, &tZ);  /* U3 <- 2*W*Z */
	gf_mul8(&tZ, &tZ);                   /* tZ <- 8*Z^2 */
	gf_neg(&P3->X.w64, &tZ);             /* X3 <- -8*Z^2 */
	P3->Z.w64 = t2;                      /* Z3 <- W^2 */
	gf_sub(&P3->T.w64, &tX, &t2);        /* T3 <- 2*X - W^2 */
}

/*
 * Repeated doublings in fractional (x,u) coordinates.
 */
UNUSED
static void
do255e_double_x_xu(do255e_point_xu *P3, const do255e_point_xu *P1, unsigned n)
{
	gf tX, tW, tZ, t1, t2;

	if (n == 0) {
		*P3 = *P1;
		return;
	}

	/*
	 * First half-doubling, combined with conversion from fractional (x,u)
	 * to Jacobian (x,w); output in E[r](-2*a,a^2-4*b).
	 *   X' = Z^2*T^4
	 *   W' = Z*T^2 - (2*X + a*Z)*U^2
	 *   Z' = Z*U*T
	 * Note that a = 0 for curve do255e.
	 * Cost: 4M+2S
	 */
	gf_sqr_inline(&tW, &P1->U.w64);              /* tW <- U^2 */
	gf_mul_inline(&t1, &P1->Z.w64, &P1->T.w64);  /* t1 <- Z*T */
	gf_mul2(&tW, &tW);                           /* tW <- 2*U^2 */
	gf_mul_inline(&t2, &t1, &P1->T.w64);         /* t2 <- Z*T^2 */
	gf_mul_inline(&tW, &tW, &P1->X.w64);         /* tW <- 2*X*U^2 */
	gf_sqr_inline(&tX, &t2);                     /* tX <- Z^2*T^4 */
	gf_sub(&tW, &t2, &tW);                       /* tW <- Z*T^2 - 2*X*U^2 */
	gf_mul_inline(&tZ, &t1, &P1->U.w64);         /* tZ <- Z*U*T */

	/*
	 * For n-1 doublings, apply psi_1/2() then psi_1().
	 */
	while (n -- > 1) {
		/*
		 * Here we use 2M+4S per iteration. We could change that
		 * to 1M+5S by noticing that:
		 *   Z'' = 2*W*W'*Z
		 *       = ((W+W')^2 - W^2 - W'^2)*Z
		 * but on a big x86 (Coffee Lake core) this happens to be
		 * somewhat slower, because of the extra subtractions
		 * and/or less favourable dependencies.
		 */

		/*
		 * X' = W^4
		 * W' = W^2 - 2*X
		 * Z' = 2*W*Z
		 */
		gf_sqr_inline(&t1, &tW);
		gf_mul2(&t2, &tX);
		gf_mul_inline(&tZ, &tW, &tZ);
		gf_sub(&tW, &t1, &t2);
		gf_sqr_inline(&tX, &t1);
		gf_mul2(&tZ, &tZ);

		/*
		 * X'' = W'^4
		 * W'' = W'^2 - 2*X'
		 * Z'' = W'*Z'
		 */
		gf_sqr_inline(&t1, &tW);
		gf_mul2(&t2, &tX);
		gf_mul_inline(&tZ, &tW, &tZ);
		gf_sub(&tW, &t1, &t2);
		gf_sqr_inline(&tX, &t1);
	}

	/*
	 * Final half-doubling, combined with conversion back to
	 * fractional (x,u) coordinates.
	 *   X' = 4*b*Z^2
	 *   Z' = W^2
	 *   U' = 2*W*Z
	 *   T' = 2*X - 2*a*Z^2 - W^2
	 * Note that a = 0 and b = -2 for curve do255e.
	 * Cost: 3S  (with 2*W*Z = (W+Z)^2 - W^2 - Z^2)
	 */
	gf_sqr_inline(&t2, &tW);             /* t2 <- W^2 */
	gf_add(&t1, &tW, &tZ);               /* t1 <- W + Z */
	gf_sqr_inline(&tZ, &tZ);             /* tZ <- Z^2 */
	gf_sqr_inline(&t1, &t1);             /* t1 <- (W + Z)^2 */
	gf_mul2(&tX, &tX);                   /* tX <- 2*X */
	gf_sub2(&P3->U.w64, &t1, &t2, &tZ);  /* U3 <- 2*W*Z */
	gf_mul8(&tZ, &tZ);                   /* tZ <- 8*Z^2 */
	gf_neg(&P3->X.w64, &tZ);             /* X3 <- -8*Z^2 */
	P3->Z.w64 = t2;                      /* Z3 <- W^2 */
	gf_sub(&P3->T.w64, &tX, &t2);        /* T3 <- 2*X - W^2 */
}

/*
 * Mixed addition of a point in fractional (x,u) coordinates with a
 * point in affine coordinates.
 */
UNUSED
static void
do255e_add_mixed_xu(do255e_point_xu *P3,
	const do255e_point_xu *P1, const do255e_point_affine_xu *P2)
{
	gf t1, t3, t5, t6, t7, t8, t9, t10;

	/* t1 <- X1*X2 */
	gf_mul_inline(&t1, &P1->X.w64, &P2->X.w64);

	/* t2 <- Z1*Z2 = Z1 */

	/* t3 <- U1*U2 */
	gf_mul_inline(&t3, &P1->U.w64, &P2->U.w64);

	/* t4 <- T1*T2 = T1 */

	/* t5 <- X1*Z2 + X2*Z1 = X1 + X2*Z1 */
	gf_mul_inline(&t5, &P1->Z.w64, &P2->X.w64);
	gf_add(&t5, &t5, &P1->X.w64);

	/* t6 <- U1*T2 + U2*T1 = U1 + U2*T1 */
	gf_mul_inline(&t6, &P1->T.w64, &P2->U.w64);
	gf_add(&t6, &t6, &P1->U.w64);

	/* t7 <- t1 + b*t2  (with b = -2 and t2 = Z1) */
	gf_sub2(&t7, &t1, &P1->Z.w64, &P1->Z.w64);

	/* t8 <- t4*t7  (with t4 = T1) */
	gf_mul_inline(&t8, &P1->T.w64, &t7);

	/* t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2) */
	gf_mul_inline(&t9, &t3, &t5);
	gf_mul4(&t9, &t9);
	gf_neg(&t9, &t9);

	/* t10 <- (t4 + alpha*t3)*(t5 + t7)  (with t4 = T1 and alpha = 2) */
	gf_add(&t5, &t5, &t7);
	gf_mul2(&t3, &t3);
	gf_add(&t3, &t3, &P1->T.w64);
	gf_mul_inline(&t10, &t3, &t5);

	/* U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*Z1) */
	gf_neg(&t6, &t6);
	gf_add(&t1, &t1, &P1->Z.w64);
	gf_add(&t1, &t1, &P1->Z.w64);
	gf_mul_inline(&P3->U.w64, &t6, &t1);

	/* Z3 <- t8 - t9 */
	gf_sub(&P3->Z.w64, &t8, &t9);

	/* T3 <- t8 + t9 */
	gf_add(&P3->T.w64, &t8, &t9);

	/* X3 <- b*(t10 - t8 + beta*t9)  (with b = -2 and beta = 1/2)
	         = 2*(t8 - t10) - t9 */
	gf_sub(&t8, &t8, &t10);
	gf_mul2(&t8, &t8);
	gf_sub(&P3->X.w64, &t8, &t9);
}

/*
 * Point addition in fractional (x,u) coordinates.
 */
UNUSED
static void
do255e_add_xu(do255e_point_xu *P3,
	const do255e_point_xu *P1, const do255e_point_xu *P2)
{
	gf t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

	/* t1 <- X1*X2 */
	gf_mul_inline(&t1, &P1->X.w64, &P2->X.w64);

	/* t2 <- Z1*Z2 */
	gf_mul_inline(&t2, &P1->Z.w64, &P2->Z.w64);

	/* t3 <- U1*U2 */
	gf_mul_inline(&t3, &P1->U.w64, &P2->U.w64);

	/* t4 <- T1*T2 */
	gf_mul_inline(&t4, &P1->T.w64, &P2->T.w64);

	/* t5 <- (X1 + Z1)*(X2 + Z2) - t1 - t2 */
	gf_add(&t5, &P1->X.w64, &P1->Z.w64);
	gf_add(&t8, &P2->X.w64, &P2->Z.w64);
	gf_mul_inline(&t5, &t5, &t8);
	gf_sub2(&t5, &t5, &t1, &t2);

	/* t6 <- (U1 + T1)*(U2 + T2) - t3 - t4 */
	gf_add(&t6, &P1->U.w64, &P1->T.w64);
	gf_add(&t9, &P2->U.w64, &P2->T.w64);
	gf_mul_inline(&t6, &t6, &t9);
	gf_sub2(&t6, &t6, &t3, &t4);

	/* t7 <- t1 + b*t2  (with b = -2) */
	gf_sub2(&t7, &t1, &t2, &t2);

	/* t8 <- t4*t7 */
	gf_mul_inline(&t8, &t4, &t7);

	/* t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2) */
	gf_mul_inline(&t9, &t3, &t5);
	gf_mul4(&t9, &t9);
	gf_neg(&t9, &t9);

	/* t10 <- (t4 + alpha*t3)*(t5 + t7)  (with alpha = 2) */
	gf_add(&t5, &t5, &t7);
	gf_mul2(&t3, &t3);
	gf_add(&t3, &t3, &t4);
	gf_mul_inline(&t10, &t3, &t5);

	/* U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*t2) */
	gf_neg(&t6, &t6);
	gf_mul2(&t2, &t2);
	gf_add(&t1, &t1, &t2);
	gf_mul_inline(&P3->U.w64, &t6, &t1);

	/* Z3 <- t8 - t9 */
	gf_sub(&P3->Z.w64, &t8, &t9);

	/* T3 <- t8 + t9 */
	gf_add(&P3->T.w64, &t8, &t9);

	/* X3 <- b*(t10 - t8 + beta*t9)  (with b = -2 and beta = 1/2)
	         = 2*(t8 - t10) - t9 */
	gf_sub(&t8, &t8, &t10);
	gf_mul2(&t8, &t8);
	gf_sub(&P3->X.w64, &t8, &t9);
}
