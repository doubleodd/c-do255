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
@ ABI: registers r0-r7 are consumed; high registers are untouched.
@
@ Cost: 95
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce256_partial_inner, %function
modr_reduce256_partial_inner:
	@ Extract the top two bits of the input, and combine with 'ah'.
	@ Also, keep a copy of the top word into a stack slot (truncated
	@ to 30 bits).
	ldr	r4, [r1, #28]
	lsrs	r3, r4, #30
	lsls	r2, r2, #2
	orrs	r2, r3
	lsls	r4, r4, #2
	lsrs	r4, r4, #2
	push	{ r1, r4 }

	@ Multiply the extra bits (in r2) with R0, into r4:r5:r6:r7:r1
	adr	r1, const_modr_reduce256_partial_R0
	ldm	r1!, { r4, r5, r6, r7 }

.macro MMah  rx, carry_in
	lsrs	r3, \rx, #16
	uxth	\rx, \rx
	muls	\rx, r2
	muls	r3, r2
	.if (\carry_in) != 0
	adds	\rx, r1
	.endif
	lsrs	r1, \rx, #16
	adds	r3, r1
	uxth	\rx, \rx
	lsls	r1, r3, #16
	orrs	\rx, r1
	lsrs	r1, r3, #16
.endm

	MMah  r4, 0
	MMah  r5, 1
	MMah  r6, 1
	MMah  r7, 1

	@ Add the result (5 words) to the input (truncated to 254 bits).
	pop	{ r2 }
	ldm	r2!, { r3 }
	adds	r3, r4
	stm	r0!, { r3 }
	ldm	r2!, { r3, r4 }
	adcs	r3, r5
	adcs	r4, r6
	stm	r0!, { r3, r4 }
	ldm	r2, { r2, r3, r4, r5 }
	adcs	r2, r7
	adcs	r3, r1
	eors	r1, r1
	adcs	r4, r1
	adcs	r5, r1
	pop	{ r6 }    @ a[7] truncated to 30 bits
	adcs	r6, r1
	stm	r0!, { r2, r3, r4, r5, r6 }

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
@ Cost: 94
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce256_finish_inner, %function
modr_reduce256_finish_inner:
	@ Subtract R = 2^254 - R0 (i.e. add R0 but subtract 2^254),
	@ into r8:r10:r11:r12:r4:r5:r6:r7. Since the source value is
	@ less than 2*R < 2^255, the result can fit in 256 bits with
	@ a sign bit.
	adr	r2, const_modr_reduce256_finish_R0
	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	adds	r4, r6
	adcs	r5, r7
	mov	r8, r4
	mov	r10, r5
	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	adcs	r4, r6
	adcs	r5, r7
	mov	r11, r4
	mov	r12, r5
	ldm	r1!, { r4, r5, r6, r7 }
	eors	r2, r2
	adcs	r4, r2
	adcs	r5, r2
	adcs	r6, r2
	adcs	r7, r2
	movs	r2, #1
	lsls	r2, r2, #30
	subs	r7, r2

	@ If the result is negative, then we must replace the computed
	@ subtraction with the source value. Otherwise, we keep it.
	asrs	r3, r7, #31
	bics	r4, r3
	bics	r5, r3
	bics	r6, r3
	bics	r7, r3
	subs	r1, #16
	ldm	r1!, { r2 }
	ands	r2, r3
	orrs	r4, r2
	ldm	r1!, { r2 }
	ands	r2, r3
	orrs	r5, r2
	ldm	r1!, { r2 }
	ands	r2, r3
	orrs	r6, r2
	ldm	r1!, { r2 }
	ands	r2, r3
	orrs	r7, r2
	adds	r0, #16
	stm	r0!, { r4, r5, r6, r7 }
	subs	r0, #32
	subs	r1, #32
	mov	r4, r8
	mov	r5, r10
	bics	r4, r3
	bics	r5, r3
	ldm	r1!, { r6, r7 }
	ands	r6, r3
	ands	r7, r3
	orrs	r4, r6
	orrs	r5, r7
	stm	r0!, { r4, r5 }
	mov	r4, r11
	mov	r5, r12
	bics	r4, r3
	bics	r5, r3
	ldm	r1!, { r6, r7 }
	ands	r6, r3
	ands	r7, r3
	orrs	r4, r6
	orrs	r5, r7
	stm	r0!, { r4, r5 }

	bx	lr
	.align	2
const_modr_reduce256_finish_R0:
	.long	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
	.size	modr_reduce256_finish_inner, .-modr_reduce256_finish_inner

@ =======================================================================
@ void modr_reduce384_partial_inner(i256 *d, const i384 *a)
@
@ ABI: all registers are consumed.
@
@ Cost: 105 + cost(iXX_mul_inner[128x128->256])
@       + cost(modr_reduce256_partial_inner)
@       = 887
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	modr_reduce384_partial_inner, %function
modr_reduce384_partial_inner:
	push	{ r0, r1, lr }
	sub	sp, #36

	@ Multiply the high third (ah) by R0, into stack (0)
	adds	r1, #32
	movs	r2, r1
	movs	r3, #4
	adr	r4, const_modr_reduce384_partial_R0
	movs	r5, #4
	mov	r0, sp
	movs	r1, #8
	bl	iXX_mul_inner

	@ Multiply ah*R0 by 4 (into r2:r3:r4:r5:r8:r10:r11:r12, extra
	@ bits in r7).
	add	r0, sp, #12
	ldm	r0!, { r2, r3, r4, r5, r6 }
	lsrs	r7, r2, #30
	lsrs	r1, r3, #30
	lsls	r3, r3, #2
	orrs	r3, r7
	mov	r8, r3
	lsrs	r7, r4, #30
	lsls	r4, r4, #2
	orrs	r4, r1
	mov	r10, r4
	lsrs	r1, r5, #30
	lsls	r5, r5, #2
	orrs	r5, r7
	mov	r11, r5
	lsrs	r7, r6, #30
	lsls	r6, r6, #2
	orrs	r6, r1
	mov	r12, r6
	lsls	r5, r2, #2
	mov	r0, sp
	ldm	r0!, { r2, r3, r4 }
	lsrs	r6, r4, #30
	orrs	r5, r6
	lsls	r4, r4, #2
	lsrs	r6, r3, #30
	orrs	r4, r6
	lsls	r3, r3, #2
	lsrs	r6, r2, #30
	orrs	r3, r6
	lsls	r2, r2, #2

	@ Add 4*r0*ah to (a mod 2^256) into stack (0) (extra bits in r2).
	ldr	r0, [sp, #40]
	ldm	r0!, { r1, r6 }
	adds	r1, r2
	adcs	r6, r3
	str	r1, [sp]
	add	r1, sp, #4
	stm	r1!, { r6 }
	ldm	r0!, { r2, r3, r6 }
	adcs	r2, r4
	adcs	r3, r5
	mov	r4, r8
	adcs	r6, r4
	stm	r1!, { r2, r3, r6 }
	ldm	r0!, { r2, r3, r4 }
	mov	r5, r10
	adcs	r2, r5
	mov	r5, r11
	adcs	r3, r5
	mov	r5, r12
	adcs	r4, r5
	stm	r1!, { r2, r3, r4 }
	eors	r2, r2
	adcs	r2, r7

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
@ Cost: 25 + cost(modr_reduce256_partial_inner)
@       + cost(modr_reduce256_finish_inner) + cost(memcpy[32])
@       = 214 + cost(memcpy[32])
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	scalar_decode_inner, %function
scalar_decode_inner:
	push	{ r0, lr }
	movs	r2, #32
	bl	memcpy
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
@ Cost: 56 + cost(memcpy[32])
@ =======================================================================

	.align	1
	.global	CN(scalar_is_reduced)
	.thumb
	.thumb_func
	.type	CN(scalar_is_reduced), %function
CN(scalar_is_reduced):
	push	{ r4, r5, r6, r7, lr }
	sub	sp, #36

	movs	r1, r0
	mov	r0, sp
	movs	r2, #32
	bl	memcpy
	mov	r0, sp
	adr	r1, const_scalar_is_reduced_r
	ldm	r0!, { r2, r3, r4 }
	ldm	r1!, { r5, r6, r7 }
	subs	r2, r5
	sbcs	r3, r6
	sbcs	r4, r7
	ldm	r0!, { r2, r3, r4 }
	ldm	r1!, { r5, r6, r7 }
	sbcs	r2, r5
	sbcs	r3, r6
	sbcs	r4, r7
	ldm	r0!, { r2, r3 }
	ldm	r1!, { r5, r6 }
	sbcs	r2, r5
	sbcs	r3, r6
	sbcs	r0, r0
	rsbs	r0, r0, #0

	add	sp, #36
	pop	{ r4, r5, r6, r7, pc }
	.align	2
const_scalar_is_reduced_r:
	.long	0x74D84525, 0x1F52C8AE, 0x54078C53, 0x9D0C930F
	.long	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x3FFFFFFF
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
	push	{ r4, r5, r6, r7, lr }
	bl	modr_reduce256_partial_inner
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(modr_reduce256_partial), .-CN(modr_reduce256_partial)

@ =======================================================================
@ void CURVE_modr_reduce256_finish(i256 *d, const i256 *a, uint32_t ah)
@ =======================================================================

	.align	1
	.global	CN(modr_reduce256_finish)
	.thumb
	.thumb_func
	.type	CN(modr_reduce256_finish), %function
CN(modr_reduce256_finish):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	modr_reduce256_finish_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(modr_reduce256_finish), .-CN(modr_reduce256_finish)

@ =======================================================================
@ void CURVE_modr_reduce384_partial(i256 *d, const i256 *a, uint32_t ah)
@ =======================================================================

	.align	1
	.global	CN(modr_reduce384_partial)
	.thumb
	.thumb_func
	.type	CN(modr_reduce384_partial), %function
CN(modr_reduce384_partial):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	modr_reduce384_partial_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(modr_reduce384_partial), .-CN(modr_reduce384_partial)
