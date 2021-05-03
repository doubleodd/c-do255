@ =======================================================================
@ Assembly module for curve do255s on ARM Cortex M4.
@ =======================================================================

	.syntax	unified
	.cpu	cortex-m4
	.file	"asm_do255s_cm4.s"
	.text

#define CN(name)   do255s_ ## name
#define MQ         3957
#define INVT510    0x0035E7B3, 0xF1F837B1, 0xC7DBEC43, 0xB32AF51D, 0x1ADD2F9B, 0x6503C8E8, 0xCCD520D9, 0x4711614B

@ =======================================================================

#include "asm_gf_cm4.S"
#include "asm_padd_do255s_cm4.S"
#include "asm_icore_cm4.S"
#include "asm_scalar_do255s_cm4.S"
#include "asm_pmul_cm4.S"
#include "asm_pmul_do255s_cm4.S"
#include "asm_lagrange_cm4.S"
#include "asm_pvrfy_do255s_cm4.S"
