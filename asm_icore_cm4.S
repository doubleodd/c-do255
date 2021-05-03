@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It implements basic integer operations for scalars.
@ The following macro must have been defined:
@   CN(name)   makes a symbol name suitable for exporting
@ =======================================================================

@ =======================================================================
@ void mul128x128trunc_inner(const u128 *a, const i128 *b)
@
@ ABI: all registers consumed. Result is returned in r0:r1:r2:r3.
@
@ Cost: 23
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	mul128x128trunc_inner, %function
mul128x128trunc_inner:
	ldm	r0, { r8, r10, r11, r12 }
	ldm	r1!, { r4, r5, r6, r7 }

	@ a0*b3 + a3*b0 -> r3
	mul	r3, r7, r8
	umlal	r3, r2, r4, r12

	@ r7 and r12 are free.

	@ a0*b0
	umull	r0, r1, r4, r8

	@ a0*b1 and a1*b0
	umull	r7, r2, r4, r10
	umaal	r1, r7, r5, r8

	@ carries: r2, r7

	@ a0*b2, a1*b1 and a2*b0
	umull	r12, r4, r4, r11
	umaal	r2, r12, r5, r10
	umaal	r2, r7, r6, r8

	@ carries: r3, r4, r7, r12

	@ a1*b2 and a2*b1
	umaal	r3, r4, r5, r11
	umaal	r3, r7, r6, r10
	add	r3, r12

	bx	lr
	.size	mul128x128trunc_inner, .-mul128x128trunc_inner

@ =======================================================================
@ void CURVE_mul128x128trunc(i128 *d, const i128 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 32 + cost(mul128x128trunc_inner) = 55
@ =======================================================================

	.align	1
	.global	CN(mul128x128trunc)
	.thumb
	.thumb_func
	.type	CN(mul128x128trunc), %function
CN(mul128x128trunc):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	movs	r0, r2
	bl	mul128x128trunc_inner
	pop	{ r4 }
	stm	r4!, { r0, r1, r2, r3 }
	pop	{ r4, r5, r6, r7, r8, r10, r11, lr }
	.size	CN(mul128x128trunc), .-CN(mul128x128trunc)

@ =======================================================================
@ void CURVE_mul128x128(i256 *d, const i128 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 64
@ =======================================================================

	.align	1
	.global	CN(mul128x128)
	.thumb
	.thumb_func
	.type	CN(mul128x128), %function
CN(mul128x128):
	push	{ r4, r5, r6, r7, r8, r10, r11 }

	@ Read b0..b3.
	ldm	r2, { r8, r10, r11, r12 }

	@ Process a0.
	ldm	r1!, { r2 }
	umull	r3, r4, r2, r8
	eors	r5, r5
	umull	r6, r7, r5, r5
	umlal	r4, r5, r2, r10
	umlal	r5, r6, r2, r11
	umlal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a1.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a2.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a3.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3, r4, r5, r6, r7 }

	pop	{ r4, r5, r6, r7, r8, r10, r11 }
	bx	lr
	.size	CN(mul128x128), .-CN(mul128x128)

@ =======================================================================
@ void CURVE_mul256x128(i384 *d, const i256 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 100
@ =======================================================================

	.align	1
	.global	CN(mul256x128)
	.thumb
	.thumb_func
	.type	CN(mul256x128), %function
CN(mul256x128):
	push	{ r4, r5, r6, r7, r8, r10, r11 }

	@ Read b0..b3.
	ldm	r2, { r8, r10, r11, r12 }

	@ Process a0.
	ldm	r1!, { r2 }
	umull	r3, r4, r2, r8
	eors	r5, r5
	umull	r6, r7, r5, r5
	umlal	r4, r5, r2, r10
	umlal	r5, r6, r2, r11
	umlal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a1.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a2.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a3.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a4.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a5.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a6.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3 }

	@ Process a7.
	ldm	r1!, { r2 }
	eors	r3, r3
	umaal	r3, r4, r2, r8
	umaal	r4, r5, r2, r10
	umaal	r5, r6, r2, r11
	umaal	r6, r7, r2, r12
	stm	r0!, { r3, r4, r5, r6, r7 }

	pop	{ r4, r5, r6, r7, r8, r10, r11 }
	bx	lr
	.size	CN(mul256x128), .-CN(mul256x128)

@ =======================================================================
@ void CURVE_mul256x256(i512 *d, const i256 *a, const i256 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 82 + 2*cost(CURVE_mul256x128) = 282
@ =======================================================================

	.align	1
	.global	CN(mul256x256)
	.thumb
	.thumb_func
	.type	CN(mul256x256), %function
CN(mul256x256):
	@ Since this function is not time-critical, we implement it
	@ by calling mul256x128() twice, which saves code space.
	push	{ r4, r5, r6, lr }
	sub	sp, #48
	movs	r4, r0
	movs	r5, r1
	movs	r6, r2
	bl	CN(mul256x128)
	mov	r0, sp
	movs	r1, r5
	add	r2, r6, #16
	bl	CN(mul256x128)

	adds	r4, #16
	pop	{ r0, r1, r2 }
	ldm	r4, { r3, r5, r6 }
	adds	r0, r3
	adcs	r1, r5
	adcs	r2, r6
	stm	r4!, { r0, r1, r2 }
	pop	{ r0, r1, r2 }
	ldm	r4, { r3, r5, r6 }
	adcs	r0, r3
	adcs	r1, r5
	adcs	r2, r6
	stm	r4!, { r0, r1, r2 }
	pop	{ r0, r1, r2 }
	ldm	r4, { r3, r5 }
	eors	r6, r6
	adcs	r0, r3
	adcs	r1, r5
	adcs	r2, r6
	stm	r4!, { r0, r1, r2 }
	pop	{ r0, r1, r2 }
	adcs	r0, r6
	adcs	r1, r6
	adcs	r2, r6
	stm	r4!, { r0, r1, r2 }

	pop	{ r4, r5, r6, pc }
	.size	CN(mul256x256), .-CN(mul256x256)

@ =======================================================================
@ uint32_t CURVE_i256_add(i256 *d, const i256 *a, const i256 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 56
@ =======================================================================

	.align	1
	.global	CN(i256_add)
	.thumb
	.thumb_func
	.type	CN(i256_add), %function
CN(i256_add):
	push	{ r4, r5, r6, r7 }

	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	adds	r4, r6
	adcs	r5, r7
	stm	r0!, { r4, r5 }
	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	adcs	r4, r6
	adcs	r5, r7
	stm	r0!, { r4, r5 }
	ldm	r1, { r1, r3, r4, r5 }
	ldm	r2!, { r6, r7 }
	adcs	r1, r6
	adcs	r3, r7
	ldm	r2!, { r6, r7 }
	adcs	r4, r6
	adcs	r5, r7
	stm	r0!, { r1, r3, r4, r5 }
	sbcs	r0, r0
	adds	r0, #1

	pop	{ r4, r5, r6, r7 }
	bx	lr
	.size	CN(i256_add), .-CN(i256_add)

@ =======================================================================
@ uint32_t CURVE_i256_sub(i256 *d, const i256 *a, const i256 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 56
@ =======================================================================

	.align	1
	.global	CN(i256_sub)
	.thumb
	.thumb_func
	.type	CN(i256_sub), %function
CN(i256_sub):
	push	{ r4, r5, r6, r7 }

	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	subs	r4, r6
	sbcs	r5, r7
	stm	r0!, { r4, r5 }
	ldm	r1!, { r4, r5 }
	ldm	r2!, { r6, r7 }
	sbcs	r4, r6
	sbcs	r5, r7
	stm	r0!, { r4, r5 }
	ldm	r1, { r1, r3, r4, r5 }
	ldm	r2!, { r6, r7 }
	sbcs	r1, r6
	sbcs	r3, r7
	ldm	r2!, { r6, r7 }
	sbcs	r4, r6
	sbcs	r5, r7
	stm	r0!, { r1, r3, r4, r5 }
	sbcs	r0, r0
	rsbs	r0, r0, #0

	pop	{ r4, r5, r6, r7 }
	bx	lr
	.size	CN(i256_sub), .-CN(i256_sub)
