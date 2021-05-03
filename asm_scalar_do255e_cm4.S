@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It implements modular reduction for scalars (curve do255e).
@ The following macro must have been defined:
@   CN(name)   makes a symbol name suitable for exporting
@ =======================================================================

@ =======================================================================
@ void modr_reduce256_partial_inner(i256 *d, const i256 *a, uint32_t ah)
@
@ ABI: all registers are consumed.
@
@ Cost: 40
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce256_partial_inner, %function
modr_reduce256_partial_inner:
	@ Load the source value completely.
	ACC_LOAD  r1

	@ Extract the top two bits of the input, and combine with 'ah'.
	lsls	r2, #2
	orr	r3, r2, r12, lsr #30
	bic	r12, r12, #0xC0000000

	@ Multiply the extra bits (in r2) with R0, and add to the current
	@ value.
	ldr	r1, const_modr_reduce256_partial_R0
	eors	r2, r2
	umaal	r4, r2, r1, r3
	ldr	r1, const_modr_reduce256_partial_R0 + 4
	umaal	r5, r2, r1, r3
	ldr	r1, const_modr_reduce256_partial_R0 + 8
	umaal	r6, r2, r1, r3
	ldr	r1, const_modr_reduce256_partial_R0 + 12
	umaal	r7, r2, r1, r3
	adds	r8, r2
	adcs	r10, #0
	adcs	r11, #0
	adc	r12, #0

	@ Write the result.
	ACC_STORE  r0

	bx	lr
	.align	2
const_modr_reduce256_partial_R0:
	.long	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
	.size	modr_reduce256_partial_inner, .-modr_reduce256_partial_inner

@ =======================================================================
@ void modr_reduce256_finish_inner(i256 *d, const i256 *a)
@
@ ABI: all registers are consumed.
@
@ Cost: 54
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce256_finish_inner, %function
modr_reduce256_finish_inner:
	@ Load the whole value.
	ACC_LOAD  r1

	@ Subtract R = 2^254 - R0 (i.e. add R0 but subtract 2^254).
	.align	2
	ldrd	r1, r2, const_modr_reduce256_finish_R0
	adds	r4, r1
	adcs	r5, r2
	.align	2
	ldrd	r1, r2, const_modr_reduce256_finish_R0 + 8
	adcs	r6, r1
	adcs	r7, r2
	adcs	r8, #0
	adcs	r10, #0
	adcs	r11, #0
	adc	r12, #0xC0000000

	@ If the value is negative, then add back R.
	lsr	r3, r12, #31
	ldr	r1, const_modr_reduce256_finish_minusR0
	eors	r2, r2
	umaal	r4, r2, r3, r1
	ldr	r1, const_modr_reduce256_finish_minusR0 + 4
	umaal	r5, r2, r3, r1
	ldr	r1, const_modr_reduce256_finish_minusR0 + 8
	umaal	r6, r2, r3, r1
	ldr	r1, const_modr_reduce256_finish_minusR0 + 12
	umaal	r7, r2, r3, r1
	rsbs	r1, r3, #0
	umaal	r8, r2, r3, r1
	umaal	r10, r2, r3, r1
	umaal	r11, r2, r3, r1
	lsrs	r1, #2
	umaal	r12, r2, r3, r1

	@ Store the result.
	ACC_STORE  r0

	bx	lr
	.align	2
const_modr_reduce256_finish_R0:
	.long	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
const_modr_reduce256_finish_minusR0:
	.long	0x74D84525, 0x1F52C8AE, 0x54078C53, 0x9D0C930F
	.size	modr_reduce256_finish_inner, .-modr_reduce256_finish_inner

@ =======================================================================
@ void modr_reduce384_partial_inner(i256 *d, const i384 *a)
@
@ ABI: all registers are consumed.
@
@ Cost: 72 + cost(CURVE_mul128x128) + cost(modr_reduce256_partial_inner)
@       = 176
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce384_partial_inner, %function
modr_reduce384_partial_inner:
	push	{ r0, r1, lr }
	sub	sp, #36

	@ Multiply the high third (ah) by R0, into stack (0)
	mov	r0, sp
	adds	r1, #32
	adr	r2, const_modr_reduce384_partial_R0
	bl	CN(mul128x128)

	@ Load ah*R0 and multiply it by 4 (extra bits in r14).
	ACC_LOAD  sp
	movs	r2, #3
	lsr	r14, r4, #30
	lsls	r4, r4, #2
	umaal	r5, r14, r2, r5
	umaal	r6, r14, r2, r6
	umaal	r7, r14, r2, r7
	umaal	r8, r14, r2, r8
	umaal	r10, r14, r2, r10
	umaal	r11, r14, r2, r11
	umaal	r12, r14, r2, r12

	@ Add 4*r0*ah to (a mod 2^256) into stack (0) (extra bits in r2).
	ldr	r0, [sp, #40]
	ldm	r0!, { r1, r2, r3 }
	adds	r4, r1
	adcs	r5, r2
	adcs	r6, r3
	ldm	r0!, { r1, r2, r3 }
	adcs	r7, r1
	adcs	r8, r2
	adcs	r10, r3
	ldm	r0!, { r1, r2 }
	adcs	r11, r1
	adcs	r12, r2
	adcs	r2, r14, #0
	ACC_STORE  sp

	@ Perform partial reduction, write into the output buffer.
	ldr	r0, [sp, #36]
	mov	r1, sp
	bl	modr_reduce256_partial_inner

	add	sp, #44
	pop	{ pc }
	.align	2
const_modr_reduce384_partial_R0:
	.long	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
	.size	modr_reduce384_partial_inner, .-modr_reduce384_partial_inner

@ =======================================================================
@ void scalar_decode_inner(uint32_t *d, const void *src)
@
@ Source scalar may be unaligned. Output is fully reduced. Input and
@ output MUST NOT overlap.
@
@ ABI: all registers are consumed.
@
@ Cost: 46 + cost(modr_reduce256_partial_inner)
@       + cost(modr_reduce256_finish_inner)
@       = 140  (if source is aligned)
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	scalar_decode_inner, %function
scalar_decode_inner:
	push	{ r0, lr }
	ACC_LOAD_UNALIGNED  r1
	ACC_STORE  r0
	ldr	r0, [sp]
	movs	r1, r0
	movs	r2, #0
	bl	modr_reduce256_partial_inner
	ldr	r0, [sp]
	movs	r1, r0
	bl	modr_reduce256_finish_inner
	pop	{ r0, pc }
	.size	scalar_decode_inner, .-scalar_decode_inner

@ =======================================================================
@ int CURVE_scalar_is_reduced(const void *a)
@
@ This function conforms to the public ABI.
@
@ Cost: 52  (if source is aligned)
@ =======================================================================

	.align	1
	.global	CN(scalar_is_reduced)
	.thumb
	.thumb_func
	.type	CN(scalar_is_reduced), %function
CN(scalar_is_reduced):
	push	{ r4, r5, r6, r7, r8, r10, r11 }
	ACC_LOAD_UNALIGNED  r0
	adr	r0, const_scalar_is_reduced_r
	ldm	r0, { r0, r1, r2, r3 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	mvn	r0, #0
	mvn	r1, #0xC0000000
	sbcs	r8, r0
	sbcs	r10, r0
	sbcs	r11, r0
	sbcs	r12, r1
	sbcs	r0, r0
	rsbs	r0, r0, #0
	pop	{ r4, r5, r6, r7, r8, r10, r11 }
	bx	lr
	.align	2
const_scalar_is_reduced_r:
	.long	0x74D84525, 0x1F52C8AE, 0x54078C53, 0x9D0C930F
	.size	CN(scalar_is_reduced), .-CN(scalar_is_reduced)

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_modr_reduce256_partial(i256 *d, const i256 *a, uint32_t ah)
@ =======================================================================

	.align	1
	.global	CN(modr_reduce256_partial)
	.thumb
	.thumb_func
	.type	CN(modr_reduce256_partial), %function
CN(modr_reduce256_partial):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	modr_reduce256_partial_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(modr_reduce256_partial), .-CN(modr_reduce256_partial)

@ =======================================================================
@ void CURVE_modr_reduce256_finish(i256 *d, const i256 *a)
@ =======================================================================

	.align	1
	.global	CN(modr_reduce256_finish)
	.thumb
	.thumb_func
	.type	CN(modr_reduce256_finish), %function
CN(modr_reduce256_finish):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	modr_reduce256_finish_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(modr_reduce256_finish), .-CN(modr_reduce256_finish)

@ =======================================================================
@ void CURVE_modr_reduce384_partial(i256 *d, const i384 *a)
@ =======================================================================

	.align	1
	.global	CN(modr_reduce384_partial)
	.thumb
	.thumb_func
	.type	CN(modr_reduce384_partial), %function
CN(modr_reduce384_partial):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	modr_reduce384_partial_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(modr_reduce384_partial), .-CN(modr_reduce384_partial)
