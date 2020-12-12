/*
 * This file is meant to be included, not compiled by itself.
 * Caller must have included/defined the following prior to inclusion:
 *
 *  - included "do255.h"
 *  - defined gf and operations, including gf_sqrt() or gf_issquare()
 *
 * This file is for all implementations of do255s that use 32-bit limbs.
 * It defines:
 *  - do255s_neutral
 *  - do255s_generator
 *  - do255s_decode()
 *  - do255s_encode()
 *  - do255s_is_neutral()
 *  - do255s_eq()
 *  - do255s_add()
 *  - do255s_double()
 *  - do255s_double_x()
 */

/* a */
static const gf CURVE_A = { {
	0xFFFFF08A, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,
	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x7FFFFFFF
} };

/* 4*b */
static const gf CURVE_4B = { {
	2, 0, 0, 0, 0, 0, 0, 0
} };

/* see do255.h */
const do255s_point do255s_generator = {
	{ .w32 = { {
		0x33B156B1, 0x4803AC7D, 0x5840B591, 0x3EF83226,
		0xCB010B9D, 0x213759EC, 0x1783FB6D, 0x39BD7265
	} } },
	{ .w32 = { {
		0xAAAAA584, 0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA,
		0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA, 0x2AAAAAAA
	} } },
	{ .w32 = { { 1, 0, 0, 0, 0, 0, 0, 0 } } }
};

/*
 * Custom structure for a point in affine coordinates:
 *  - if X = 0, then this is the neutral (W is ignored);
 *  - if X != 0, then coordinate Z is implicitly equal to 1.
 */
typedef struct {
	do255_int256 X, W;
} do255s_point_affine;

/* Precomputed windows for the generator. */
static const do255s_point_affine window_G[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x33B156B1, 0x4803AC7D, 0x5840B591, 0x3EF83226,
		            0xCB010B9D, 0x213759EC, 0x1783FB6D, 0x39BD7265 }} },
		{ .w32 = {{ 0xAAAAA584, 0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA,
		            0xAAAAAAAA, 0xAAAAAAAA, 0xAAAAAAAA, 0x2AAAAAAA }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0x318C5721, 0x8C6318C6, 0x6318C631, 0x18C6318C,
		            0xC6318C63, 0x318C6318, 0x8C6318C6, 0x6318C631 }} },
		{ .w32 = {{ 0x5C3E5011, 0x7311544E, 0x2F5E574B, 0x3667EAC4,
		            0xA4176AA4, 0x6DC59919, 0x3F062A69, 0x701BB820 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x5AF542C8, 0xC92EA6A4, 0xEB2EBB62, 0xACA2B1F0,
		            0x772D49BF, 0x16AA9AB4, 0x68993FDF, 0x586E82D4 }} },
		{ .w32 = {{ 0x07647DD6, 0xB0E2616D, 0xD9B70D8F, 0x79835A85,
		            0x733292C6, 0x39545AAB, 0x3EE92301, 0x27A19C07 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x0F24B43F, 0xB21E437C, 0x2A191529, 0xE352D46A,
		            0x8E691D57, 0x5105F27F, 0x85ABD1BC, 0x203E7421 }} },
		{ .w32 = {{ 0xAFC4CC85, 0x0AA6E498, 0x51FC11F0, 0x34CABA83,
		            0x2BD1700A, 0x09F86D5B, 0x0924022E, 0x537BA3D4 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x8FAAFD54, 0xD4C60438, 0xCA27D2A0, 0x478144CA,
		            0xA3381C38, 0x4B4B7554, 0xA8C94117, 0x78CC3E71 }} },
		{ .w32 = {{ 0xB48175C6, 0x9819FEFF, 0xA19AC3DD, 0xEB222706,
		            0xBE5D5845, 0x410723C5, 0xF4EA9551, 0x3C18C25D }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x5D8AEC75, 0x0F8CCA9A, 0xF00CFF35, 0x2337B2C5,
		            0x60F3A520, 0x4C38E257, 0xB38EC25A, 0x5B8B19C1 }} },
		{ .w32 = {{ 0x77A6F073, 0x5E273AB2, 0x5019A510, 0xDFCD0E54,
		            0x7DE61215, 0xB4D344CB, 0xE6CE27C4, 0x20CB4C0B }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x9916AE56, 0xE4574479, 0x444F94C0, 0xAB719188,
		            0xF4727AB0, 0x8BCF3B53, 0x4119EE4B, 0x239A8F01 }} },
		{ .w32 = {{ 0x646D9E07, 0x6D72E617, 0x9844C50F, 0xE08B7A2B,
		            0x18FDB104, 0x164E30CF, 0xF1A5F7AF, 0x44C18E58 }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xD4089F3F, 0x486D0418, 0x7E21EC49, 0x415C3064,
		            0x6C060489, 0xDB2F693C, 0x5E1CF582, 0x475F06F5 }} },
		{ .w32 = {{ 0x88DF33A3, 0x716647F9, 0xA7CEE964, 0xFA6FD869,
		            0xF220CA25, 0xDB2F9858, 0x03825F9C, 0x7F45D8FB }} }
	}
};
static const do255s_point_affine window_G64[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x60E01FEA, 0xE13252EB, 0xD231ECB5, 0x9197370C,
		            0xABBECE1A, 0xF31CFB1E, 0x31F15CB8, 0x11BFABD6 }} },
		{ .w32 = {{ 0xDD7C174B, 0xB28D51CF, 0x098E1C01, 0xDE06E300,
		            0xD9B8B95C, 0x6452EFF8, 0xCD23E562, 0x752F0B51 }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0xF7B2D34D, 0xE3AD260D, 0x64E44AEC, 0xEF4C98E5,
		            0x25234018, 0x41B749AF, 0x5D77E307, 0x0722EF99 }} },
		{ .w32 = {{ 0xBB760285, 0x02784B47, 0x675F77EC, 0x674B8473,
		            0xD6CEEED4, 0x5F71BB6F, 0x58D994F1, 0x34D98636 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x33514CC7, 0xE8DD6EC9, 0x1DB42CF4, 0x0202CF4A,
		            0x80156AEC, 0x3F149B87, 0x8499AA74, 0x0C60AC54 }} },
		{ .w32 = {{ 0x5AF8352D, 0x30FE0F15, 0x7B8C5B86, 0x34869785,
		            0x97D45FFC, 0xACF44E15, 0x72745EA9, 0x3716511D }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0xEECB8F77, 0x2ADA4E59, 0x354A3FA3, 0xD577FCEA,
		            0xC41ABDB8, 0xE3E9AA4A, 0x0D3E7390, 0x0FA4D0D0 }} },
		{ .w32 = {{ 0x6CF3E045, 0x903B4E22, 0xC979F4C6, 0x8F6A1B46,
		            0x19C37B18, 0x751E1AD6, 0xB6C5C58B, 0x2BF70C62 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0xA82DE754, 0x263A3655, 0xEB125279, 0x5D06767F,
		            0x0D36949E, 0xBF69A9F9, 0x5330C479, 0x016525E0 }} },
		{ .w32 = {{ 0x1F8CB310, 0x47F9C71D, 0x84526B7E, 0xC09602D6,
		            0xD3D0B456, 0x7B98E8A5, 0x7DACED6A, 0x1740308C }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x44A76584, 0x83437B3B, 0x5EA37019, 0x85337A5B,
		            0xBACBCA6E, 0x0AF40A3C, 0x7ECE3AE0, 0x6129FB01 }} },
		{ .w32 = {{ 0x8AFDFB16, 0xE2AF05AD, 0xDA602A72, 0x00D46F68,
		            0x7B2FCA11, 0x06730D41, 0xAE0247FC, 0x4A9C28F6 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x217CF53D, 0x231E2747, 0xC981B66C, 0x55A000CE,
		            0xB0324CF0, 0xDEB03F43, 0xB9D86B7A, 0x6856DF69 }} },
		{ .w32 = {{ 0x314FD2F0, 0x42E5774F, 0x32E45BFE, 0xDB2CD756,
		            0xA21D4CDF, 0x339FDC58, 0xA88DAA47, 0x4F727EDE }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0xA8BA33DF, 0x417B5BE5, 0xC729E5BE, 0x1146B966,
		            0xC9837F83, 0xF19DEF75, 0x0C99B2FC, 0x4B20562C }} },
		{ .w32 = {{ 0x8B0DF192, 0x7C1122BA, 0x0975002A, 0x3BD15C55,
		            0x5A7A5DB3, 0x2186EFEC, 0xEB3A7184, 0x0BFBF827 }} }
	}
};
static const do255s_point_affine window_G128[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x7635D9CF, 0x4D5B37E4, 0xB313CF11, 0xAE7E7A89,
		            0xA5963E8B, 0xA277DD52, 0xEAAC2050, 0x6B925003 }} },
		{ .w32 = {{ 0xCC1D108C, 0x5DE53387, 0xC4102FDF, 0x1B92AB3F,
		            0xFFC65BC6, 0x3B38F632, 0xBDEA8465, 0x01AAC5EE }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0xD6481616, 0x8CF76BBE, 0xD3896825, 0xA5DECE3F,
		            0x947BA3BA, 0x92DC7116, 0xC2BF4CAD, 0x615D397B }} },
		{ .w32 = {{ 0xBBF61353, 0x5B833529, 0x8CB4847B, 0x60A31BDE,
		            0xC5400E3B, 0x7127A393, 0x95631BD2, 0x052888DC }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x777B8F2F, 0xFB21E7ED, 0x40239F37, 0x1EE37518,
		            0xEC00783A, 0xBE214D43, 0xE49E0AF7, 0x6D6E7DAD }} },
		{ .w32 = {{ 0xDF608FBC, 0x8373BB18, 0xC95742CE, 0xCD100440,
		            0x1E917572, 0xC3F95DCD, 0x9CDCECBD, 0x25C18129 }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x313E150A, 0xD9AB2AB4, 0x923D8F48, 0xAF3680D8,
		            0x1CE8DFAF, 0x8FE12E04, 0xAEB11CAB, 0x254561A8 }} },
		{ .w32 = {{ 0x9BC36E11, 0x9FA5CE6B, 0xC4C6B597, 0x365A4759,
		            0x218EE6B1, 0xC381D92E, 0xFFFB0E9F, 0x61878956 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x2AF21E15, 0x8CAB4C3A, 0x9E9632CD, 0xC5AE64DF,
		            0x0FE2775D, 0x6A931D2B, 0x9C0049FA, 0x0F36BEC8 }} },
		{ .w32 = {{ 0xEB3CF496, 0x993D9FD0, 0xC0BB21FB, 0xC53D312D,
		            0x15D40A91, 0xAE76833F, 0x2725408B, 0x780E0C23 }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x62C08A1B, 0x9741CAC3, 0x6C3BAE68, 0x2E9F8753,
		            0x7CBB757C, 0xD3108E6E, 0xB1298012, 0x5880D4DE }} },
		{ .w32 = {{ 0x0BFA72C1, 0x09531A9B, 0x571D5B75, 0x4B256BDC,
		            0x0D75DA62, 0x20AADE43, 0x1721E73B, 0x621B11AF }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x9B2990E5, 0xF1800926, 0x279C58A8, 0x61572BF5,
		            0x5E0A0DD4, 0x1E021898, 0x252D00F2, 0x7286BEC0 }} },
		{ .w32 = {{ 0x1C2DC846, 0xBCEA20E4, 0xBEFC73C5, 0xA59EA984,
		            0x2EEC48A2, 0x283FF87E, 0xF9D49A0F, 0x4EC4844F }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0x7F6C8053, 0x9C35A563, 0xD6021602, 0x5706387F,
		            0x0C1697F6, 0x0AF3C147, 0xA7420D24, 0x7847E5A5 }} },
		{ .w32 = {{ 0x693C216C, 0x247FA839, 0x227D9686, 0x9F1A1C6F,
		            0xE0BCB5DA, 0xC9754112, 0xA79828C4, 0x5184C659 }} }
	}
};
static const do255s_point_affine window_G192[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x0A435303, 0x25E4F86C, 0x2C92D6F9, 0xFB1350AF,
		            0xA9209A23, 0xBB42CF3A, 0x3FE7D419, 0x5CA75997 }} },
		{ .w32 = {{ 0xDB95C9B1, 0xDB63AE7E, 0x624E03F1, 0xADFC025C,
		            0x752FDD08, 0x08185CED, 0xFF3B79C9, 0x5DB3D29A }} }
	},
	/* 2 */
	{
		{ .w32 = {{ 0xAEDC3F8C, 0x9A53A899, 0x898348F7, 0xAF074A24,
		            0x90D94DED, 0x90B71D2F, 0xF702F9AA, 0x4A7E9271 }} },
		{ .w32 = {{ 0xA64D45B7, 0x6E055713, 0x14F85099, 0x8052CB7C,
		            0x86425CD2, 0xB63BB7F8, 0x6BF10114, 0x35451641 }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x37DC3DE4, 0xB67F7ABB, 0x929F1AD1, 0x18FA3503,
		            0x569E43DA, 0x925D2711, 0x09D1C53B, 0x7A1B7478 }} },
		{ .w32 = {{ 0xB121D266, 0x6A6BAD93, 0xE45B328B, 0x8D5B43DA,
		            0x37AB40DE, 0xA0D5F082, 0xEE791C15, 0x58CDFB5C }} }
	},
	/* 4 */
	{
		{ .w32 = {{ 0x80BC554A, 0xB202F9EA, 0x1B4F5DBC, 0x582E9CC2,
		            0xEAB2F0CE, 0x35C48881, 0xDC57D3A5, 0x7B6750B1 }} },
		{ .w32 = {{ 0x3158C2D5, 0x50AA49C4, 0x7B85FDBF, 0x38FE415E,
		            0x7A7DD9BE, 0x0E7883AF, 0x5D0B9459, 0x69B8E78D }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x51B600E4, 0xDDF898D1, 0xCB1749C1, 0xDEBA13EF,
		            0xD36CEA4A, 0x958A72DF, 0x26B50975, 0x16D05067 }} },
		{ .w32 = {{ 0xA984B080, 0x12052CAB, 0xA527F662, 0xF6F6C0A3,
		            0x4002659A, 0x633EF292, 0x1EA30B6B, 0x7E3429FF }} }
	},
	/* 6 */
	{
		{ .w32 = {{ 0x5E2EBD62, 0xD9514371, 0x89448F48, 0xF4D671FD,
		            0x4C7281E9, 0x6CD2F819, 0x14CDE8FA, 0x125707C4 }} },
		{ .w32 = {{ 0x2FB9E4A9, 0x1D8A92C5, 0x3791F474, 0x65F9A9D4,
		            0xD2DF1771, 0x9D2F8902, 0xFABEBDB3, 0x5DE484BC }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x0C18BF9C, 0xE4046DEE, 0xA4840D07, 0x1890CEAF,
		            0x6FB8A6E6, 0xA32DD6DA, 0xFA22AABD, 0x66A9ACFF }} },
		{ .w32 = {{ 0x43412540, 0x365B9775, 0xA1160301, 0x56064313,
		            0x69147780, 0xA6B97610, 0x45655A99, 0x44F29938 }} }
	},
	/* 8 */
	{
		{ .w32 = {{ 0x8B607110, 0xA5B332AE, 0x8C5C2A7B, 0xC6A1F588,
		            0x5F176F78, 0x4686BCEF, 0x0771051D, 0x2EEA2C71 }} },
		{ .w32 = {{ 0x6AE34331, 0x6D0E9F71, 0x7EBB560C, 0xA2D563E3,
		            0x13FC17A6, 0x02C47EC2, 0x5FEC969B, 0x6A77A305 }} }
	}
};
static const do255s_point_affine window_odd_G128[] = {
	/* 1 */
	{
		{ .w32 = {{ 0x7635D9CF, 0x4D5B37E4, 0xB313CF11, 0xAE7E7A89,
		            0xA5963E8B, 0xA277DD52, 0xEAAC2050, 0x6B925003 }} },
		{ .w32 = {{ 0xCC1D108C, 0x5DE53387, 0xC4102FDF, 0x1B92AB3F,
		            0xFFC65BC6, 0x3B38F632, 0xBDEA8465, 0x01AAC5EE }} }
	},
	/* 3 */
	{
		{ .w32 = {{ 0x777B8F2F, 0xFB21E7ED, 0x40239F37, 0x1EE37518,
		            0xEC00783A, 0xBE214D43, 0xE49E0AF7, 0x6D6E7DAD }} },
		{ .w32 = {{ 0xDF608FBC, 0x8373BB18, 0xC95742CE, 0xCD100440,
		            0x1E917572, 0xC3F95DCD, 0x9CDCECBD, 0x25C18129 }} }
	},
	/* 5 */
	{
		{ .w32 = {{ 0x2AF21E15, 0x8CAB4C3A, 0x9E9632CD, 0xC5AE64DF,
		            0x0FE2775D, 0x6A931D2B, 0x9C0049FA, 0x0F36BEC8 }} },
		{ .w32 = {{ 0xEB3CF496, 0x993D9FD0, 0xC0BB21FB, 0xC53D312D,
		            0x15D40A91, 0xAE76833F, 0x2725408B, 0x780E0C23 }} }
	},
	/* 7 */
	{
		{ .w32 = {{ 0x9B2990E5, 0xF1800926, 0x279C58A8, 0x61572BF5,
		            0x5E0A0DD4, 0x1E021898, 0x252D00F2, 0x7286BEC0 }} },
		{ .w32 = {{ 0x1C2DC846, 0xBCEA20E4, 0xBEFC73C5, 0xA59EA984,
		            0x2EEC48A2, 0x283FF87E, 0xF9D49A0F, 0x4EC4844F }} }
	}
};

/*
 * We get do255s_neutral, do255s_decode(), do255s_encode(),
 * do255s_is_neutral() and do255s_eq() from pcore_w64.c.
 */
#include "pcore_w32.c"

/* r */
static const uint32_t R1[] = {
	0x396152C7, 0xDCF2AC65, 0x912B7F03, 0x2ACF567A,
	0x00000000, 0x00000000, 0x00000000, 0x40000000
};

/* 2*r */
static const uint32_t R2[] = {
	0x72C2A58E, 0xB9E558CA, 0x2256FE07, 0x559EACF5,
	0x00000000, 0x00000000, 0x00000000, 0x80000000
};

/* see do255.h */
void
do255s_add(do255s_point *P3,
	const do255s_point *P1, const do255s_point *P2)
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
	   Also, replace t4 with b*t4 */
	gf_half(&t4, &t4);
	gf_sub(&t9, &t6, &t4);
	gf_mul_inline(&Z3, &t9, &t7);

	/* t9 <- t7^4 */
	gf_sqr_inline(&t9, &t7);
	gf_sqr_inline(&t9, &t9);

	/* X3 <- b*t6*t9 */
	gf_mul_inline(&X3, &t6, &t9);
	gf_half(&X3, &X3);

	/* t10 <- (t5 + a*t3)*(t6 + b*t4)
	   a = -1
	   b*t4 was already computed (in t4)
	   We overwrite t5 and t6, which we won't need anymore */
	gf_sub(&t5, &t5, &t3);
	gf_add(&t6, &t6, &t4);
	gf_mul_inline(&t10, &t5, &t6);

	/* W3 <- -t10 - 2*b*t3*t8
	   b = 1/2, hence 2*b = 1.
	   We overwrite t8. */
	gf_mul_inline(&t8, &t3, &t8);
	gf_sub2(&W3, &GF_ZERO, &t10, &t8);

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
do255s_add_mixed(do255s_point *P3,
	const do255s_point *P1, const do255s_point_affine *P2)
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
	   Also, replace t1 with b*t1 */
	gf_half(&t1, &t1);
	gf_sub(&t9, &t6, &t1);
	gf_mul_inline(&Z3, &t9, &t7);

	/* t9 <- t7^4 */
	gf_sqr_inline(&t9, &t7);
	gf_sqr_inline(&t9, &t9);

	/* X3 <- b*t6*t9 */
	gf_mul_inline(&X3, &t6, &t9);
	gf_half(&X3, &X3);

	/* t10 <- (t5 + a*t3)*(t6 + b*t1)
	   a = -1
	   b*t1 was already computed (in t1)
	   We overwrite t5 and t6, which we won't need anymore */
	gf_sub(&t5, &t5, &P1->Z.w32);
	gf_add(&t6, &t6, &t1);
	gf_mul_inline(&t10, &t5, &t6);

	/* W3 <- -t10 - 2*b*t3*t8
	   b = 1/2, hence 2*b = 1.
	   We overwrite t8. */
	gf_mul_inline(&t8, &P1->Z.w32, &t8);
	gf_sub2(&W3, &GF_ZERO, &t10, &t8);

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

static inline void
do255s_double_inline(do255s_point *P3, const do255s_point *P1)
{
	gf t1, t2, t3, t4, X3;

	/* t1 <- W*Z */
	gf_mul_inline(&t1, &P1->W.w32, &P1->Z.w32);

	/* t2 <- t1^2 */
	gf_sqr_inline(&t2, &t1);

	/* X' <- 8*t2^2 */
	gf_sqr_inline(&t4, &t2);
	gf_mul8(&X3, &t4);

	/* t3 <- (W + Z)^2 - 2*t1 */
	gf_add(&t3, &P1->W.w32, &P1->Z.w32);
	gf_sqr_inline(&t3, &t3);
	gf_sub2(&t3, &t3, &t1, &t1);

	/* W' <- 2*t2 - t3^2
	   We reuse t2 as scratch */
	gf_sqr_inline(&t4, &t3);
	gf_mul2(&t2, &t2);
	gf_sub(&P3->W.w32, &t2, &t4);

	/* Z' <- 2*t1*(2*X - t3)
	   We reuse t2 and t3 as scratch */
	gf_mul2(&t2, &P1->X.w32);
	gf_sub(&t2, &t2, &t3);
	gf_mul_inline(&t3, &t2, &t1);
	gf_mul2(&P3->Z.w32, &t3);

	/* Writing to P3->X was delayed because we still needed P1->X,
	   which might be the same value. */
	P3->X.w32 = X3;
}

/* see do255.h */
void
do255s_double(do255s_point *P3, const do255s_point *P1)
{
	do255s_double_inline(P3, P1);
}

/* see do255.h */
void
do255s_double_x(do255s_point *P3, const do255s_point *P1, unsigned n)
{
	*P3 = *P1;
	while (n -- > 0) {
		do255s_double_inline(P3, P3);
	}
}
