@ =======================================================================
@ Assembly module for curve do255e on ARM Cortex M4.
@ =======================================================================

	.syntax	unified
	.cpu	cortex-m4
	.file	"asm_do255e_cm4.s"
	.text

#define CN(name)   do255e_ ## name
#define MQ         18651
#define INVT510    0x4AF83D33, 0x750356D7, 0x1CA0B754, 0x8ECD5CE6, 0x003BB5AB, 0xCFB32C8A, 0x63A3EB02, 0x713CD715

@ =======================================================================

#include "asm_gf_cm4.S"
#include "asm_padd_do255e_cm4.S"
#include "asm_icore_cm4.S"
#include "asm_scalar_do255e_cm4.S"
#include "asm_pmul_cm4.S"
#include "asm_pmul_do255e_cm4.S"
