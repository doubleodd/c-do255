/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations, including gf_sqrt() or gf_issquare()
 *
 * This file is for all implementations of do255e that use 32-bit limbs.
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
static const gf CURVE_A = { {
	0, 0, 0, 0, 0, 0, 0, 0
} };

/* 4*b */
static const gf CURVE_4B = { {
	0xFFFFB71D, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF
} };

/* see do255.h */
const do255e_point do255e_generator = {
	{ .w32 = { { 2, 0, 0, 0, 0, 0, 0, 0 } } },
	{ .w32 = { { 1, 0, 0, 0, 0, 0, 0, 0 } } },
	{ .w32 = { { 1, 0, 0, 0, 0, 0, 0, 0 } } }
};

/*
 * Custom structure for a point in affine coordinates:
 *  - if X = 0, then this is the neutral (W is ignored);
 *  - if X != 0, then coordinate Z is implicitly equal to 1.
 */
typedef struct {
	do255_int256 X, W;
} do255e_point_affine;

/* Precomputed windows for the generator. */
static const do255e_point_affine window_G[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x00000002, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} },
		{ .w32 = {{ 0x00000001, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x8E38AAE3, 0xE38E38E3, 0x38E38E38, 0x8E38E38E,
		            0xE38E38E3, 0x38E38E38, 0x8E38E38E, 0x638E38E3 }} },
		{ .w32 = {{ 0x55554932, 0x55555555, 0x55555555, 0x55555555,
		            0x55555555, 0x55555555, 0x55555555, 0x15555555 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x00000152, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} },
		{ .w32 = {{ 0x3B139548, 0xB13B13B1, 0x13B13B13, 0x3B13B13B,
		            0xB13B13B1, 0x13B13B13, 0x3B13B13B, 0x313B13B1 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x7D1A36D1, 0xEA5E1BA0, 0x3CCEA916, 0xD0AF7707,
		            0xF2D92037, 0xC68A748F, 0xA8DE7571, 0x43554439 }} },
		{ .w32 = {{ 0xBC930B68, 0xA736E93C, 0x614EDC4A, 0x93CBC932,
		            0xC4AA736E, 0x932614ED, 0x36E93CBC, 0x2EDC4AA7 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x74D7D3CA, 0x9B7D88CD, 0x193896AC, 0x0E31B461,
		            0xE97D44DB, 0x93464506, 0x1DCD9949, 0x0ABC8AC6 }} },
		{ .w32 = {{ 0x6A3F5D77, 0xFB43801A, 0x9087618B, 0xEDA4968D,
		            0x0D73F849, 0x10D79E3C, 0xFBBDE55B, 0x609CA158 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x24DF12CF, 0x19348B57, 0x6232945A, 0xAB7572F6,
		            0x2A63BF4B, 0xFFAA8904, 0x633EF2E2, 0x3024ED85 }} },
		{ .w32 = {{ 0x28B85FAB, 0x262D3C74, 0x998A5A0A, 0x9DF21226,
		            0xE92A27C8, 0xDFD986E1, 0x2417F9F7, 0x18BF7AAC }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x1F43C3D5, 0x33E38125, 0x71C6FE8F, 0x49BF0E2E,
		            0x16BAEF18, 0x6AF69CC1, 0xD0C9585E, 0x36199FBA }} },
		{ .w32 = {{ 0x63844C67, 0xC5EB6106, 0xD13F51E5, 0x6C2EEAF2,
		            0xAAE3769B, 0xCB5992FC, 0x39CF7FD9, 0x7FE818CE }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xCDCBDFAE, 0xDDBCA65E, 0x47F38FBD, 0xD19C744A,
		            0x5653455E, 0x15099FB5, 0x961EC01F, 0x1FD0DC2A }} },
		{ .w32 = {{ 0x4678FCF8, 0x2636BC87, 0xB0F9793B, 0xDE73D7EE,
		            0xA7C16795, 0xC255284A, 0x83360C10, 0x288B55D8 }} }
	}
};
static const do255e_point_affine window_G64[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x0B89842A, 0xF288EF16, 0x78434147, 0x29A63B0E,
		            0x560CEB19, 0xE0CF3EC6, 0xD8D61F33, 0x01EB989A }} },
		{ .w32 = {{ 0x5A47EAA4, 0x8B86DEF1, 0x1750050F, 0xAFB67C34,
		            0xAE7A2FF3, 0x8F287F28, 0xC765F5EC, 0x036A33B4 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x87FDEB41, 0xABEF504D, 0x250A6B59, 0x3A2D867D,
		            0xAB3E16A0, 0x5906F21A, 0x081B8A67, 0x58327D52 }} },
		{ .w32 = {{ 0xBE4BCA79, 0x901C9C88, 0x2CC89A6B, 0xABFA8A1E,
		            0x82BC00EB, 0x685DB6BD, 0x126CF2D2, 0x678BABB7 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x6C66AD1E, 0x82670D06, 0x92F2D505, 0x73870F4A,
		            0x54269292, 0x86539F60, 0xCCE5F35C, 0x729C7248 }} },
		{ .w32 = {{ 0x2CAEE5F4, 0x0B13F655, 0x8CE6DF73, 0xCE86AFAA,
		            0xB6734D10, 0xBB436375, 0xB12AD294, 0x6B5E41D5 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x39689FC6, 0x15356A7D, 0xA0967CE6, 0xF52A5DE2,
		            0xD02C8707, 0xAD2738B8, 0x988AA077, 0x75B28A23 }} },
		{ .w32 = {{ 0x33B9389F, 0x87A93451, 0xB96C58BC, 0xB6C29B4F,
		            0x1D21802A, 0xF9F51FD9, 0xBE98C7C1, 0x254BD9E2 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x35736AF3, 0x912E4A40, 0x82FDA79B, 0x07140C0F,
		            0xC841F59F, 0xD2B0BE5C, 0xBAAB26EE, 0x29EF71A2 }} },
		{ .w32 = {{ 0xF8C60628, 0xFAFD5722, 0x7FAEBFC2, 0x141F54D7,
		            0x6FCCD1F9, 0x27A15253, 0xB9EBBC2E, 0x7CBD0296 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x745F399F, 0xDC3DE207, 0xCA47E2C5, 0xD7C5D084,
		            0xEF678E9B, 0xDFA47415, 0x049C2FE8, 0x5F4DD2E4 }} },
		{ .w32 = {{ 0x594577E3, 0xA10D33EA, 0xF0329B3F, 0xBE0DEDD7,
		            0xCCC4B58D, 0xD5067AAE, 0xE8E4962D, 0x2C413D98 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xD13300C2, 0xBD95F2D0, 0xCC6D2EB2, 0xB45A06E8,
		            0x9C33B4F5, 0x0FE48D17, 0x66F67471, 0x799D6EDE }} },
		{ .w32 = {{ 0xF547D795, 0x3E3AD5AF, 0x514C3F3D, 0x5A751E7D,
		            0x86AE7788, 0x9932746E, 0x60B1F263, 0x7EA491EC }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xC9BF2E76, 0x977017D6, 0x94D3EC82, 0x17CE5868,
		            0xA288282D, 0xA07F1106, 0xA14E17A8, 0x1334610C }} },
		{ .w32 = {{ 0xE790CC6B, 0x7163483F, 0xB99D871A, 0xE36106D7,
		            0x728AE136, 0x112CBBAB, 0xB08F9F64, 0x310817C6 }} }
	}
};
static const do255e_point_affine window_G128[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x412DC1CF, 0xF4C0EDA0, 0x42CAE105, 0x922CFBCE,
		            0xA9AB696A, 0x3DA00EE7, 0x2D1685D3, 0x4901F275 }} },
		{ .w32 = {{ 0x960CE525, 0x73F8A80C, 0xE05AF970, 0xCAA5C24A,
		            0xF397CDC2, 0x9E41072A, 0x9C380AF6, 0x12429062 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x24642BCE, 0x56C2951A, 0xE3F95859, 0x237D2D2F,
		            0x3923ED4E, 0xECDDD0D2, 0xB4E9ACF6, 0x54B0FC59 }} },
		{ .w32 = {{ 0xB06E0E19, 0xD111D48B, 0x94B35170, 0x3077E018,
		            0xD3791EBC, 0xF3052C06, 0x92F22EBD, 0x7EC3EF93 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x0D4BC8C6, 0xAE76901F, 0x650411FC, 0x639CF3CA,
		            0x7A6D29DD, 0xBD65C636, 0x11B00C5F, 0x7A9E514A }} },
		{ .w32 = {{ 0x452DB961, 0xAB6AA839, 0x8E437067, 0x30516EEB,
		            0x93D6924C, 0x8325BF70, 0x65D85744, 0x1B3303F4 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x2F4B5577, 0xC796B240, 0x3824131E, 0xB848DE14,
		            0xB17AF31F, 0x35479390, 0x3DFA2F98, 0x0EB9DD4F }} },
		{ .w32 = {{ 0xA98FE0D9, 0xF028FEB8, 0x62B98251, 0x5847F0A6,
		            0x8C18D93D, 0xCA08C438, 0x31E9F31D, 0x7C1F86F7 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x90F5B3EF, 0xD0188F1C, 0xCDAA6603, 0x6EB97D5C,
		            0x128A881A, 0x3A02A52F, 0xD48913C2, 0x24B1E144 }} },
		{ .w32 = {{ 0xE494709A, 0xD00E7396, 0x1916D0CA, 0xAF46F90C,
		            0x5400DD81, 0x1FC2A062, 0x23DF0998, 0x737B0452 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0xAA4F7318, 0x3A8904BA, 0xD54F3120, 0x20706347,
		            0x810ECB00, 0x018E979E, 0x5DD1DD32, 0x2DC835C3 }} },
		{ .w32 = {{ 0x6FB31533, 0x20CFB1A1, 0xE639371C, 0x26359A57,
		            0x52D3CB49, 0x1DF94CD1, 0x5DB32796, 0x1C597178 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xAB722BF0, 0xE00A7990, 0x28802BA6, 0x6AECDBE1,
		            0xC1B63688, 0x29CBAEC1, 0x08514327, 0x0A35E900 }} },
		{ .w32 = {{ 0x8C4767E5, 0x67D86103, 0x53B2D34F, 0xF8650F1E,
		            0x3AB8E979, 0x89E74C47, 0x9D470823, 0x1489C34C }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0x7A0BD2F8, 0x9B55C0CC, 0x593578E5, 0x73D4BA09,
		            0x2B82E751, 0x89A107F0, 0x9D8C1DB2, 0x23C4CFE4 }} },
		{ .w32 = {{ 0xD5475CF5, 0x8E395956, 0xB08F567F, 0x5416F288,
		            0x0F99D290, 0x77EACEDD, 0x6B0407F8, 0x7AB28933 }} }
	}
};
static const do255e_point_affine window_G192[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x49A387B1, 0x299651B1, 0x0EA619DD, 0xDBBBF3CA,
		            0x0447FC6D, 0x32B9793D, 0x48FD74BD, 0x5B0539CA }} },
		{ .w32 = {{ 0x464B007F, 0x7A9A2C16, 0x62B189C8, 0xE508F756,
		            0x67ADAB60, 0x03392299, 0xBF5D806A, 0x21E5BFE4 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x3C884959, 0xF64B79D9, 0xAC0E0B1E, 0xB251DCFE,
		            0xC8983EDD, 0x58A9A996, 0x4BE960EE, 0x6336494A }} },
		{ .w32 = {{ 0xFF86D6E3, 0x0C96BE28, 0xCA1E25A2, 0x0A09592D,
		            0x464938F7, 0x695E57BE, 0x63961EB4, 0x4F61E6B8 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x47BDD5E2, 0x778CBF55, 0x9E30D01E, 0x3A27EC8F,
		            0xEFF3D017, 0x8C5EBE33, 0x462227E7, 0x0E09365D }} },
		{ .w32 = {{ 0x3177B04C, 0x48564E6A, 0x7B3A521E, 0x6F6E360D,
		            0x82DF23BB, 0x65633314, 0x112B0B14, 0x09D6D0BC }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x62013ABB, 0xEB1402F8, 0xCB9A3CAF, 0x726CE071,
		            0x7E5C377D, 0x05A0C0E3, 0xDE4F58DE, 0x2F430D8B }} },
		{ .w32 = {{ 0xC861FA49, 0xA579EB87, 0x98CEB220, 0xD493015D,
		            0x3717BF69, 0xB3BA5F65, 0x55B15EF8, 0x144BDE48 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x22FFC155, 0xE8ACF6B3, 0x395D1744, 0x96FCF4C8,
		            0x5FA28430, 0xA8172E14, 0xAF670D78, 0x51120313 }} },
		{ .w32 = {{ 0x02639F0C, 0x74D74596, 0xCC7BEE01, 0x038FB618,
		            0x15084D21, 0xC9738E8F, 0xF3AA3A38, 0x6C0461C3 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x89911FA5, 0x0FF4B381, 0xE4696279, 0x323F3CAC,
		            0xB754FE64, 0x8B16307F, 0x8E2F729F, 0x12C68BEF }} },
		{ .w32 = {{ 0x0367C664, 0x8591E1D7, 0xD984C4C0, 0x7784D4A1,
		            0x4FB3AF7F, 0xFC0B24DD, 0x61CE91B8, 0x1B8C613A }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xB3D218FF, 0xAD3D5EA0, 0xD9E27F66, 0x6FC9D7A8,
		            0x360755AA, 0xCDA645EA, 0xA66796A7, 0x629AD6AF }} },
		{ .w32 = {{ 0xCFC1C383, 0x21A115AE, 0xF1928F9C, 0xA88FF348,
		            0xC7454E17, 0x7F075335, 0xA4900326, 0x4914ADF6 }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xCBE23870, 0xC9D74325, 0x92EB5981, 0x3C8F7997,
		            0x0E4ACC9D, 0xB4D3EBF2, 0xF9CD4639, 0x5D505019 }} },
		{ .w32 = {{ 0x5E47E695, 0x8FA137B8, 0xC53F8D60, 0x15A6AA40,
		            0x4D13094E, 0xCE4451C4, 0x48AE6218, 0x18744EBC }} }
	}
};
static const do255e_point_affine window_odd_G128[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x412DC1CF, 0xF4C0EDA0, 0x42CAE105, 0x922CFBCE,
		            0xA9AB696A, 0x3DA00EE7, 0x2D1685D3, 0x4901F275 }} },
		{ .w32 = {{ 0x960CE525, 0x73F8A80C, 0xE05AF970, 0xCAA5C24A,
		            0xF397CDC2, 0x9E41072A, 0x9C380AF6, 0x12429062 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x0D4BC8C6, 0xAE76901F, 0x650411FC, 0x639CF3CA,
		            0x7A6D29DD, 0xBD65C636, 0x11B00C5F, 0x7A9E514A }} },
		{ .w32 = {{ 0x452DB961, 0xAB6AA839, 0x8E437067, 0x30516EEB,
		            0x93D6924C, 0x8325BF70, 0x65D85744, 0x1B3303F4 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x90F5B3EF, 0xD0188F1C, 0xCDAA6603, 0x6EB97D5C,
		            0x128A881A, 0x3A02A52F, 0xD48913C2, 0x24B1E144 }} },
		{ .w32 = {{ 0xE494709A, 0xD00E7396, 0x1916D0CA, 0xAF46F90C,
		            0x5400DD81, 0x1FC2A062, 0x23DF0998, 0x737B0452 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xAB722BF0, 0xE00A7990, 0x28802BA6, 0x6AECDBE1,
		            0xC1B63688, 0x29CBAEC1, 0x08514327, 0x0A35E900 }} },
		{ .w32 = {{ 0x8C4767E5, 0x67D86103, 0x53B2D34F, 0xF8650F1E,
		            0x3AB8E979, 0x89E74C47, 0x9D470823, 0x1489C34C }} }
	}
};

/*
 * We get do255e_neutral, do255e_decode(), do255e_encode(),
 * do255e_is_neutral() and do255e_eq() from pcore_w32.c.
 */
#include "pcore_w32.c"

/* 2*r */
static const uint32_t R2[] = {
	0xE9B08A4A, 0x3EA5915C, 0xA80F18A6, 0x3A19261E,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF
};

/* 3*r */
static const uint32_t R3[] = {
	0x5E88CF6F, 0x5DF85A0B, 0xFC16A4F9, 0xD725B92D,
	0xFFFFFFFE, 0xFFFFFFFF, 0xFFFFFFFF, 0xBFFFFFFF
};

/* see do255.h */
void
do255e_add(do255e_point *P3,
	const do255e_point *P1, const do255e_point *P2)
{
	gf t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, X3, W3, Z3;
	uint32_t fz1, fz2;

	/*
	 * Test whether P1 and/or P2 is neutral.
	 */
	fz1 = gf_iszero(&P1->Z.w32);
	fz2 = gf_iszero(&P2->Z.w32);

	/* t1 <- Z1^2 */
	gf_sqr_inline(&t1, &P1->Z.w32);

	/* t2 <- Z2^2 */
	gf_sqr_inline(&t2, &P2->Z.w32);

	/* t3 <- ((Z1 + Z2)^2 - t1 - t2)/2 */
	gf_add(&t3, &P1->Z.w32, &P2->Z.w32);
	gf_sqr_inline(&t3, &t3);
	gf_sub2(&t3, &t3, &t1, &t2);
	gf_half(&t3, &t3);

	/* t4 <- t3^2 */
	gf_sqr_inline(&t4, &t3);

	/* t5 <- W1*W2 */
	gf_mul_inline(&t5, &P1->W.w32, &P2->W.w32);

	/* t6 <- X1*X2 */
	gf_mul_inline(&t6, &P1->X.w32, &P2->X.w32);

	/* t7 <- (W1 + Z1)*(W2 + Z2) - t3 - t5 */
	gf_add(&t7, &P1->W.w32, &P1->Z.w32);
	gf_add(&t8, &P2->W.w32, &P2->Z.w32);
	gf_mul_inline(&t7, &t7, &t8);
	gf_sub2(&t7, &t7, &t3, &t5);

	/* t8 <- (X1 + t1)*(X2 + t2) - t4 - t6 */
	gf_add(&t8, &P1->X.w32, &t1);
	gf_add(&t9, &P2->X.w32, &t2);
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
	gf_sel3(&P3->X.w32, &P2->X.w32, &P1->X.w32, &X3, fz1, fz2);
	gf_sel3(&P3->W.w32, &P2->W.w32, &P1->W.w32, &W3, fz1, fz2);
	gf_sel3(&P3->Z.w32, &P2->Z.w32, &P1->Z.w32, &Z3, fz1, fz2);
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
	uint32_t fz1, fz2;

	/*
	 * Test whether P1 and/or P2 is neutral.
	 */
	fz1 = gf_iszero(&P1->Z.w32);
	fz2 = gf_iszero(&P2->X.w32);

	/* t1 <- Z1^2 */
	gf_sqr_inline(&t1, &P1->Z.w32);

	/* t2 = 1 */
	/* t3 = Z1 */
	/* t4 = t1 */

	/* t5 <- W1*W2 */
	gf_mul_inline(&t5, &P1->W.w32, &P2->W.w32);

	/* t6 <- X1*X2 */
	gf_mul_inline(&t6, &P1->X.w32, &P2->X.w32);

	/* t7 <- W1 + W2*Z1 */
	gf_mul_inline(&t7, &P1->Z.w32, &P2->W.w32);
	gf_add(&t7, &t7, &P1->W.w32);

	/* t8 <- X1 + X2*t1 */
	gf_mul_inline(&t8, &t1, &P2->X.w32);
	gf_add(&t8, &t8, &P1->X.w32);

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
	gf_mul_inline(&t8, &P1->Z.w32, &t8);
	gf_mul4(&t8, &t8);
	gf_sub(&W3, &t8, &t10);

	/*
	 * If P1 is neutral, replace P3 with P2.
	 * If P2 is neutral, replace P3 with P1.
	 * If both are neutral, then we want to use P1 as source, whose
	 * Z coordinate is then zero; this allows us to assume that
	 * P2.Z = 1 here.
	 */
	gf_sel3(&P3->X.w32, &P1->X.w32, &P2->X.w32, &X3, fz2, fz1);
	gf_sel3(&P3->W.w32, &P1->W.w32, &P2->W.w32, &W3, fz2, fz1);
	gf_sel3(&P3->Z.w32, &P1->Z.w32, &GF_ONE, &Z3, fz2, fz1);
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
	gf_sqr_inline(&tX, &P1->W.w32);
	gf_mul2(&tW, &P1->X.w32);
	gf_mul_inline(&tZ, &P1->W.w32, &P1->Z.w32);
	gf_sub(&tW, &tX, &tW);
	gf_sqr_inline(&tX, &tX);

	/*
	 * X3 = 4*b*Z'^4 = -8*Z'^4 (because b = -2)
	 * W3 = X' - (1/2)*W'^2
	 * Z3 = W'*Z'
	 *
	 * We can also compute Z3 = (1/2)*((W'+Z')^2 - W'^2 - Z'^2)
	 * (not done here).
	 */
	gf_sqr_inline(&t4, &tZ);
	gf_sqr_inline(&t5, &tW);
	gf_mul_inline(&P3->Z.w32, &tW, &tZ);
	gf_sqr_inline(&t4, &t4);
	gf_half(&t5, &t5);
	gf_sub(&P3->W.w32, &tX, &t5);
	gf_mul8(&t4, &t4);
	gf_neg(&P3->X.w32, &t4);
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
	tX = P1->X.w32;
	tW = P1->W.w32;
	tZ = P1->Z.w32;
	while (n -- > 1) {
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
	 * We can also compute Z3 = (1/2)*((W'+Z')^2 - W'^2 - Z'^2)
	 * (not done here).
	 */
	gf_sqr_inline(&t1, &tZ);
	gf_sqr_inline(&t2, &tW);
	gf_mul_inline(&P3->Z.w32, &tW, &tZ);
	gf_sqr_inline(&t1, &t1);
	gf_half(&t2, &t2);
	gf_sub(&P3->W.w32, &tX, &t2);
	gf_mul8(&t1, &t1);
	gf_neg(&P3->X.w32, &t1);
}
