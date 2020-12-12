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

/*
 * We get do255e_neutral, do255e_decode(), do255e_encode(),
 * do255e_is_neutral() and do255e_eq() from pcore_w64.c.
 */
#include "pcore_w32.c"

/* do255e_add() is implemented in assembly. */
/* do255e_double() is implemented in assembly. */
/* do255e_double_x() is implemented in assembly. */

/*
 * Point addition, with the second point in affine coordinates.
 * (implemented in assembly)
 */
void
do255e_add_mixed(do255e_point *P3,
	const do255e_point *P1, const do255e_point_affine *P2);

/*
 * Custom structures for points in (x.u) coordinates (fractional and
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
		{ .w32 = {{ 0x00000002, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} },
		{ .w32 = {{ 0x00000001, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x8E38AAE3, 0xE38E38E3, 0x38E38E38, 0x8E38E38E,
		            0xE38E38E3, 0x38E38E38, 0x8E38E38E, 0x638E38E3 }} },
		{ .w32 = {{ 0xDB6D97A3, 0xB6DB6DB6, 0x6DB6DB6D, 0xDB6DB6DB,
		            0xB6DB6DB6, 0x6DB6DB6D, 0xDB6DB6DB, 0x36DB6DB6 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x00000152, 0x00000000, 0x00000000, 0x00000000,
		            0x00000000, 0x00000000, 0x00000000, 0x00000000 }} },
		{ .w32 = {{ 0xC4043E79, 0xC2F21347, 0x066D4156, 0x6B1CEBA6,
		            0xA3E20224, 0xAB617909, 0xD30336A0, 0x12358E75 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x7D1A36D1, 0xEA5E1BA0, 0x3CCEA916, 0xD0AF7707,
		            0xF2D92037, 0xC68A748F, 0xA8DE7571, 0x43554439 }} },
		{ .w32 = {{ 0x130DB4AD, 0x65A29F71, 0xFA47C8BB, 0x9F71130D,
		            0xC8BB65A2, 0x130DFA47, 0x65A29F71, 0x7A47C8BB }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x74D7D3CA, 0x9B7D88CD, 0x193896AC, 0x0E31B461,
		            0xE97D44DB, 0x93464506, 0x1DCD9949, 0x0ABC8AC6 }} },
		{ .w32 = {{ 0xDA5B43EE, 0x1F2B6B08, 0xC44A0C63, 0xE40F8B8B,
		            0xB35FB70C, 0x5866F1F8, 0x50F768D7, 0x185034D2 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x24DF12CF, 0x19348B57, 0x6232945A, 0xAB7572F6,
		            0x2A63BF4B, 0xFFAA8904, 0x633EF2E2, 0x3024ED85 }} },
		{ .w32 = {{ 0xF91D6B18, 0x0BD0C5F1, 0x263610A7, 0xBB4A410D,
		            0x98F35F00, 0xA1AB0B9D, 0xAFDDC92B, 0x4FA6D8B6 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x1F43C3D5, 0x33E38125, 0x71C6FE8F, 0x49BF0E2E,
		            0x16BAEF18, 0x6AF69CC1, 0xD0C9585E, 0x36199FBA }} },
		{ .w32 = {{ 0x159EF4EA, 0x7EB52414, 0xEB4CC9E1, 0xB885C9D1,
		            0xEE64BF7F, 0x350914B3, 0x520AED5A, 0x6DD8CDFA }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xCDCBDFAE, 0xDDBCA65E, 0x47F38FBD, 0xD19C744A,
		            0x5653455E, 0x15099FB5, 0x961EC01F, 0x1FD0DC2A }} },
		{ .w32 = {{ 0xEB513A8B, 0x99DA8C93, 0x5DEDFC87, 0x0706B8B9,
		            0x1F778CE9, 0xC54D8F47, 0xFA2E63E5, 0x4766315B }} }
	}
};
static const do255e_point_affine_xu window_G64_xu[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x0B89842A, 0xF288EF16, 0x78434147, 0x29A63B0E,
		            0x560CEB19, 0xE0CF3EC6, 0xD8D61F33, 0x01EB989A }} },
		{ .w32 = {{ 0x039C1F18, 0x2F13AC07, 0xA2C065C1, 0x2C672A38,
		            0xAC39AC80, 0x8C35ED12, 0x4F4720D7, 0x704582D8 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x87FDEB41, 0xABEF504D, 0x250A6B59, 0x3A2D867D,
		            0xAB3E16A0, 0x5906F21A, 0x081B8A67, 0x58327D52 }} },
		{ .w32 = {{ 0xC596FA71, 0xAD5F8FBF, 0x9DE223FF, 0x41589354,
		            0xE50A4384, 0x395D2181, 0xA8A7626E, 0x1B313D36 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x6C66AD1E, 0x82670D06, 0x92F2D505, 0x73870F4A,
		            0x54269292, 0x86539F60, 0xCCE5F35C, 0x729C7248 }} },
		{ .w32 = {{ 0xA9DF1CB3, 0xA138B8CB, 0xBD96F5C9, 0x32EB5467,
		            0xDDA8C801, 0xCE4D1C30, 0xF046C173, 0x4A7E597C }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x39689FC6, 0x15356A7D, 0xA0967CE6, 0xF52A5DE2,
		            0xD02C8707, 0xAD2738B8, 0x988AA077, 0x75B28A23 }} },
		{ .w32 = {{ 0x889E3F33, 0x2B7C2F71, 0xA5E65CF4, 0x4CA4C049,
		            0x9976BFE7, 0x5CD27D90, 0x9985D602, 0x0BE56F35 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x35736AF3, 0x912E4A40, 0x82FDA79B, 0x07140C0F,
		            0xC841F59F, 0xD2B0BE5C, 0xBAAB26EE, 0x29EF71A2 }} },
		{ .w32 = {{ 0xDD093EB1, 0x4E99731A, 0x73F0A9B7, 0xE579B5FA,
		            0x77360995, 0xF98F264D, 0xC0EC0878, 0x1763183C }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x745F399F, 0xDC3DE207, 0xCA47E2C5, 0xD7C5D084,
		            0xEF678E9B, 0xDFA47415, 0x049C2FE8, 0x5F4DD2E4 }} },
		{ .w32 = {{ 0xED61CB7F, 0x4A504A5D, 0x42D7801D, 0xAF415083,
		            0xAB4295EB, 0x0519A68A, 0x0B09C2B4, 0x098D3AB9 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xD13300C2, 0xBD95F2D0, 0xCC6D2EB2, 0xB45A06E8,
		            0x9C33B4F5, 0x0FE48D17, 0x66F67471, 0x799D6EDE }} },
		{ .w32 = {{ 0x78219368, 0x660115E6, 0xC85936C8, 0x7FAE38D5,
		            0x8F33CAF0, 0x8EAD56CC, 0x410FF469, 0x51D5808D }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xC9BF2E76, 0x977017D6, 0x94D3EC82, 0x17CE5868,
		            0xA288282D, 0xA07F1106, 0xA14E17A8, 0x1334610C }} },
		{ .w32 = {{ 0x81B4D89F, 0xD98BB32E, 0xE1067F08, 0x5A0FCD3A,
		            0xAFAD3191, 0x3B845EF3, 0xAA6E2C23, 0x3540EB32 }} }
	}
};
static const do255e_point_affine_xu window_G128_xu[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x412DC1CF, 0xF4C0EDA0, 0x42CAE105, 0x922CFBCE,
		            0xA9AB696A, 0x3DA00EE7, 0x2D1685D3, 0x4901F275 }} },
		{ .w32 = {{ 0xBDA5CF10, 0x62D452F8, 0xF6EE3087, 0x50C25327,
		            0xF5A910DF, 0x9A9FE56C, 0xCCC0E440, 0x02C7AD54 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x24642BCE, 0x56C2951A, 0xE3F95859, 0x237D2D2F,
		            0x3923ED4E, 0xECDDD0D2, 0xB4E9ACF6, 0x54B0FC59 }} },
		{ .w32 = {{ 0x48BCEEAB, 0xDBB77AE4, 0x7C28884E, 0x49786096,
		            0xA8B0E7C8, 0xE77D5DAD, 0x1D343E4E, 0x12EBB188 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x0D4BC8C6, 0xAE76901F, 0x650411FC, 0x639CF3CA,
		            0x7A6D29DD, 0xBD65C636, 0x11B00C5F, 0x7A9E514A }} },
		{ .w32 = {{ 0xF7639FFE, 0xF238B8D6, 0xEF6F1127, 0xC258F59D,
		            0x870B4AD5, 0x61E2FFEB, 0xA762BCCA, 0x58F51A2E }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x2F4B5577, 0xC796B240, 0x3824131E, 0xB848DE14,
		            0xB17AF31F, 0x35479390, 0x3DFA2F98, 0x0EB9DD4F }} },
		{ .w32 = {{ 0xCAF58BCB, 0x6A6F8767, 0xCD520CD9, 0xEF9A2E5A,
		            0xEE40437C, 0x2B998E19, 0xF3E02AB1, 0x1E3A7692 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x90F5B3EF, 0xD0188F1C, 0xCDAA6603, 0x6EB97D5C,
		            0x128A881A, 0x3A02A52F, 0xD48913C2, 0x24B1E144 }} },
		{ .w32 = {{ 0x35F16389, 0x4F01C519, 0xBC0E42F2, 0x1E248A4C,
		            0x8D7A98FD, 0x6A1E8ECF, 0x47025D2C, 0x221D2D6B }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0xAA4F7318, 0x3A8904BA, 0xD54F3120, 0x20706347,
		            0x810ECB00, 0x018E979E, 0x5DD1DD32, 0x2DC835C3 }} },
		{ .w32 = {{ 0x481251FF, 0x3316AC55, 0xD5F9B4CE, 0xC3485D13,
		            0xA2893D9C, 0x4C001EB1, 0x3EA38958, 0x4E7BD8D1 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xAB722BF0, 0xE00A7990, 0x28802BA6, 0x6AECDBE1,
		            0xC1B63688, 0x29CBAEC1, 0x08514327, 0x0A35E900 }} },
		{ .w32 = {{ 0xEB88BDFF, 0xFEA35DFF, 0x46D4FCA2, 0x4D9CA40E,
		            0x7664057F, 0x4E24F2C2, 0x80DD251B, 0x5843DA12 }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0x7A0BD2F8, 0x9B55C0CC, 0x593578E5, 0x73D4BA09,
		            0x2B82E751, 0x89A107F0, 0x9D8C1DB2, 0x23C4CFE4 }} },
		{ .w32 = {{ 0xC28259E8, 0x97643EA7, 0xA0416456, 0x64C33BBE,
		            0xAFFFBCEB, 0xAC5EBA85, 0x8936CEEA, 0x1DE0359F }} }
	}
};
static const do255e_point_affine_xu window_G192_xu[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x49A387B1, 0x299651B1, 0x0EA619DD, 0xDBBBF3CA,
		            0x0447FC6D, 0x32B9793D, 0x48FD74BD, 0x5B0539CA }} },
		{ .w32 = {{ 0xA3207925, 0xC2CF1BD2, 0xAF2D7854, 0x0DF0A41F,
		            0x792ECF9F, 0xA9DEC81F, 0x328C94CB, 0x553D899B }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x3C884959, 0xF64B79D9, 0xAC0E0B1E, 0xB251DCFE,
		            0xC8983EDD, 0x58A9A996, 0x4BE960EE, 0x6336494A }} },
		{ .w32 = {{ 0xCCCD2ECD, 0x2FE6DF2E, 0xD495A757, 0x88E041F9,
		            0x711E4DC4, 0x867B0E5B, 0x917E839F, 0x29A11311 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x47BDD5E2, 0x778CBF55, 0x9E30D01E, 0x3A27EC8F,
		            0xEFF3D017, 0x8C5EBE33, 0x462227E7, 0x0E09365D }} },
		{ .w32 = {{ 0xC447B151, 0xA706C8C0, 0xDAEE410A, 0xEBBE7804,
		            0xD6BBE4DD, 0x48919032, 0x4CD5F6C3, 0x4BD05ACD }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x62013ABB, 0xEB1402F8, 0xCB9A3CAF, 0x726CE071,
		            0x7E5C377D, 0x05A0C0E3, 0xDE4F58DE, 0x2F430D8B }} },
		{ .w32 = {{ 0x73E39618, 0xA5168209, 0x137CAEAD, 0x69B80498,
		            0x24D25DF1, 0xBEBE79BE, 0xBA617D3F, 0x73449587 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x22FFC155, 0xE8ACF6B3, 0x395D1744, 0x96FCF4C8,
		            0x5FA28430, 0xA8172E14, 0xAF670D78, 0x51120313 }} },
		{ .w32 = {{ 0x9A56CFE6, 0x8718291B, 0x46914057, 0xC1B24382,
		            0x41443C33, 0xF6034589, 0x4B78B1BD, 0x63CE6A92 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x89911FA5, 0x0FF4B381, 0xE4696279, 0x323F3CAC,
		            0xB754FE64, 0x8B16307F, 0x8E2F729F, 0x12C68BEF }} },
		{ .w32 = {{ 0xDA3331BE, 0xFC68B3CE, 0x998694A9, 0xBB412817,
		            0x0AD015B5, 0xF7FD91E0, 0xAC50C135, 0x6E777DBB }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0xB3D218FF, 0xAD3D5EA0, 0xD9E27F66, 0x6FC9D7A8,
		            0x360755AA, 0xCDA645EA, 0xA66796A7, 0x629AD6AF }} },
		{ .w32 = {{ 0x054248BD, 0x48DE2BE6, 0x0C4C74F1, 0x777196A0,
		            0xD5A4F9BD, 0xC874A7B2, 0xAE285A2D, 0x58F14D66 }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xCBE23870, 0xC9D74325, 0x92EB5981, 0x3C8F7997,
		            0x0E4ACC9D, 0xB4D3EBF2, 0xF9CD4639, 0x5D505019 }} },
		{ .w32 = {{ 0x864B8022, 0xC6C4EA52, 0x27F2FE05, 0x3FACF030,
		            0xAFE0F2B2, 0x5A78F8FD, 0x2117A352, 0x7A205868 }} }
	}
};

/*
 * Repeated doublings in fractional (x,u) coordinates.
 * (implemented in assembly)
 */
void do255e_double_x_xu(do255e_point_xu *P3,
	const do255e_point_xu *P1, unsigned n);

/*
 * Mixed addition in (x,u) coordinates.
 * (implemented in assembly)
 */
void do255e_add_mixed_xu(do255e_point_xu *P3,
	const do255e_point_xu *P1, const do255e_point_affine_xu *P2);
