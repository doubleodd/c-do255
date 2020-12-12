@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It implements basic integer operations for scalars.
@ The following macro must have been defined:
@   CN(name)   makes a symbol name suitable for exporting
@ =======================================================================

@ =======================================================================
@ iXX_mul_inner(uint32_t *d, size_t d_len,
@               const uint32_t *a, size_t a_len,
@               const uint32_t *b, size_t b_len)
@
@ This function computes the product of 'a' (a_len limbs) with 'b'
@ (b_len limbs); the result has size d_len limbs. d_len MUST NOT be
@ larger than a_len + b_len. d_len MUST NOT be lower than a_len or b_len.
@ a_len and b_len MUST NOT be zero.
@
@ Output MUST NOT overlap with either input.
@
@ ABI:
@  - input parameters are expected in r0..r5
@  - all registers are consumed
@
@ Cost: If d_len == a_len + b_len:  34*a_len*b_len + 27*a_len + 5*b_len + 15
@   128x128->256    687
@   256x128->384   1339
@   128x256->384   1251
@   256x256->512   2447
@ If truncating, then cost is lower:
@   128x128->128    461
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	iXX_mul_inner, %function
iXX_mul_inner:
	@ Clear destination buffer.
	lsls	r7, r1, #2
	eors	r6, r6
.LiXX_mul_inner_loop0:
	subs	r7, #4
	str	r6, [r0, r7]
	bne	.LiXX_mul_inner_loop0

	@ Compute into the stack the limit difference for b_len.
	subs	r1, r3
	lsls	r1, r1, #2
	subs	r1, #4
	push	{ r1, lr }

	@ Adjust a_len and b_len to count bytes instead of limbs.
	lsls	r3, r3, #2
	lsls	r5, r5, #2

.LiXX_mul_inner_loop1:
	@ Input conditions:
	@   r0   pointer to output (adjusted)
	@   r3   remaining limbs in a (non-zero)
	@   r4   pointer to second operand (adjusted)

	@ Save number of remaining outer iterations (x4) in r10.
	mov	r10, r3

	@ Save number of iterations of inner loop (x4) into r14.
	@ This will be used to adjust pointers afterwards.
	mov	r14, r5

	@ At this point:
	@   r0    running pointer to output
	@   r2    current pointer to first operand
	@   r4    running pointer to second operand
	ldm	r2!, { r6 }
	mov	r8, r6
	mov	r12, r2
	movs	r7, #0

.LiXX_mul_inner_loop2:
	@ Current contents:
	@   r0   running pointer to output
	@   r4   running pointer to second operand
	@   r5   number of remaining inner iterations (x4)
	@   r7   carry word
	@   r8   multiplier for this loop
	mov	r11, r5
	ldm	r4!, { r1 }
	mov	r3, r8
	MUL32x32  r1, r2, r3, r5, r6
	adds	r1, r7
	eors	r7, r7
	adcs	r2, r7
	ldr	r3, [r0]
	adds	r1, r3
	adcs	r7, r2
	stm	r0!, { r1 }
	mov	r5, r11
	subs	r5, #4
	bne	.LiXX_mul_inner_loop2

	@ Recover the pointer to the first operand.
	mov	r2, r12

	@ Get in r5 the number of inner iterations we just did (x4).
	@ We must rewind r4 by that many bytes. We also rewind r0 by
	@ r5 bytes.
	mov	r5, r14
	subs	r4, r5
	subs	r0, r5

	@ If we have reached the end of the output buffer, then we
	@ must decrement the size of b for the next iteration. Otherwise,
	@ we must store the carry word. We also loop until we get to the
	@ end of a.
	mov	r3, r10
	ldr	r6, [sp]
	adds	r6, r3
	cmp	r6, r5
	bcc	.LiXX_mul_inner_cont2
	str	r7, [r0, r5]
	adds	r0, #4
	subs	r3, #4
	bne	.LiXX_mul_inner_loop1
	b	.LiXX_mul_inner_exit
.LiXX_mul_inner_cont2:
	subs	r5, #4
	adds	r0, #4
	subs	r3, #4
	bne	.LiXX_mul_inner_loop1

.LiXX_mul_inner_exit:
	pop	{ r0, pc }
	.size	iXX_mul_inner, .-iXX_mul_inner

@ =======================================================================
@ void CURVE_mul128x128trunc(i128 *d, const i128 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 53 + cost(iXX_mul_inner[128x128->128]) = 514
@ =======================================================================

	.align	1
	.global	CN(mul128x128trunc)
	.thumb
	.thumb_func
	.type	CN(mul128x128trunc), %function
CN(mul128x128trunc):
	@ Since output has the same type of the inputs, it may be
	@ at the same address; we thus need an intermediate stack
	@ buffer.
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r0, r4, r5, r6 }
	sub	sp, #20

	movs	r4, r2
	movs	r2, r1
	mov	r0, sp
	movs	r1, #4
	movs	r3, #4
	movs	r5, #4
	bl	iXX_mul_inner
	ldr	r0, [sp, #20]
	mov	r1, sp
	ldm	r1!, { r4, r5, r6, r7 }
	stm	r0!, { r4, r5, r6, r7 }

	add	sp, #24
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(mul128x128trunc), .-CN(mul128x128trunc)

@ =======================================================================
@ void CURVE_mul128x128(i256 *d, const i128 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 36 + cost(iXX_mul_inner[128x128->256]) = 723
@ =======================================================================

	.align	1
	.global	CN(mul128x128)
	.thumb
	.thumb_func
	.type	CN(mul128x128), %function
CN(mul128x128):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	movs	r4, r2
	movs	r2, r1
	movs	r1, #8
	movs	r3, #4
	movs	r5, #4
	bl	iXX_mul_inner

	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(mul128x128), .-CN(mul128x128)

@ =======================================================================
@ void CURVE_mul256x128(i384 *d, const i256 *a, const i128 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 36 + cost(iXX_mul_inner[256x128->384]) = 1375
@ =======================================================================

	.align	1
	.global	CN(mul256x128)
	.thumb
	.thumb_func
	.type	CN(mul256x128), %function
CN(mul256x128):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	movs	r4, r2
	movs	r2, r1
	movs	r1, #12
	movs	r3, #8
	movs	r5, #4
	bl	iXX_mul_inner

	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(mul256x128), .-CN(mul256x128)

@ =======================================================================
@ void CURVE_mul256x256(i512 *d, const i256 *a, const i256 *b)
@
@ This function conforms to the public ABI.
@
@ Cost: 36 + cost(iXX_mul_inner[256x256->512]) = 2483
@ =======================================================================

	.align	1
	.global	CN(mul256x256)
	.thumb
	.thumb_func
	.type	CN(mul256x256), %function
CN(mul256x256):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	movs	r4, r2
	movs	r2, r1
	movs	r1, #16
	movs	r3, #8
	movs	r5, #8
	bl	iXX_mul_inner

	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
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
