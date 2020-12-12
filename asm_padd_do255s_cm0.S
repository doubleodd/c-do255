@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines functions for point addition and doubling in curve do255s.
@ =======================================================================

@ =======================================================================
@ void point_double_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Cost: 102 + 2*cost(gf_mul_inner) + 4*cost(gf_sqr_inner)
@       + 4*cost(gf_mul2_inner) + 3*cost(gf_sub_inner)
@       + cost(gf_add_inner)
@       = 7533
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_double_inner, %function
point_double_inner:
	push	{ r0, r1, lr }
	sub	sp, #100

	@ Stack offsets:
	@   t1 = 0
	@   t2 = 32
	@   t3 = 64
	@ destination pointer:  [sp, #100]
	@ source pointer:       [sp, #104]

	@ t1 <- W*Z
	mov	r0, sp
	ldr	r1, [sp, #104]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_mul_inner

	@ t2 <- t1^2
	add	r0, sp, #32
	mov	r1, sp
	bl	gf_sqr_inner

	@ t3 <- (W + Z)^2 - 2*t1
	@ (also replace t1 with 2*t1)
	mov	r0, sp
	mov	r1, sp
	bl	gf_mul2_inner
	add	r0, sp, #64
	ldr	r1, [sp, #104]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_add_inner
	add	r0, sp, #64
	add	r1, sp, #64
	bl	gf_sqr_inner
	add	r0, sp, #64
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_sub_inner

	@ Z' <- 2*t1*(2*X - t3)
	@ The formal value '2*t1' is already in t1.
	ldr	r0, [sp, #100]
	adds	r0, #64
	ldr	r1, [sp, #104]
	bl	gf_mul2_inner
	ldr	r0, [sp, #100]
	adds	r0, #64
	movs	r1, r0
	add	r2, sp, #64
	bl	gf_sub_inner
	ldr	r0, [sp, #100]
	adds	r0, #64
	movs	r1, r0
	mov	r2, sp
	bl	gf_mul_inner

	@ W' <- 2*t2 - t3^2
	@ We replace t2 with 2*t2.
	@ We replace t3 with t3^2.
	add	r0, sp, #32
	add	r1, sp, #32
	bl	gf_mul2_inner
	add	r0, sp, #64
	add	r1, sp, #64
	bl	gf_sqr_inner
	ldr	r0, [sp, #100]
	adds	r0, #32
	add	r1, sp, #32
	add	r2, sp, #64
	bl	gf_sub_inner

	@ X' <- 8*t2^2
	@ The formal value '2*t2' is already in t2.
	add	r0, sp, #32
	add	r1, sp, #32
	bl	gf_sqr_inner
	ldr	r0, [sp, #100]
	add	r1, sp, #32
	bl	gf_mul2_inner

	add	sp, #108
	pop	{ pc }
	.size	point_double_inner, .-point_double_inner

@ =======================================================================
@ void point_add_inner(CURVE_point *P3,
@                      const CURVE_point *P1, const CURVE_point *P2)
@
@ Cost: 260 + 8*cost(gf_mul_inner) + 6*cost(gf_sqr_inner)
@       + 2*cost(gf_iszero_inner) + 7*cost(gf_add_inner)
@       + 3*cost(gf_half_inner) + 8*cost(gf_sub_inner)
@       + cost(gf_neg_inner) + 3*cost(gf_sel3_inner)
@       = 19759
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_add_inner, %function
point_add_inner:
	push	{ r0, r1, r2, lr }
	sub	sp, #328

	@ Stack offsets:
	@   t1    0
	@   t2    32
	@   t3    64
	@   t4    96
	@   t5    128
	@   t6    160
	@   t7    192
	@   t8    0     (shared with t1)
	@   t9    32    (shared with t2)
	@   t10   192   (shared with t7)
	@   X3    224   (X3, W3 and Z3 must be consecutive in RAM)
	@   W3    256
	@   Z3    288
	@   fz1   320   (1 word)
	@   fz2   324   (1 word)
	@   &P3   328   (1 word)
	@   &P1   332   (1 word)
	@   &P2   336   (1 word)

	@ Test whether P1 and/or P2 is neutral; save the result.
	mov	r8, r2
	movs	r0, r1
	bl	gf_iszero_inner
	str	r0, [sp, #320]
	mov	r0, r8
	bl	gf_iszero_inner
	str	r0, [sp, #324]

	@ t1 <- Z1^2
	mov	r0, sp
	ldr	r1, [sp, #332]
	adds	r1, #64
	bl	gf_sqr_inner

	@ t2 <- Z2^2
	add	r0, sp, #32
	ldr	r1, [sp, #336]
	adds	r1, #64
	bl	gf_sqr_inner

	@ t3 <- ((Z1 + Z2)^2 - t1 - t2)/2
	add	r0, sp, #64
	ldr	r1, [sp, #332]
	adds	r1, #64
	ldr	r2, [sp, #336]
	adds	r2, #64
	bl	gf_add_inner
	add	r0, sp, #64
	add	r1, sp, #64
	bl	gf_sqr_inner
	add	r0, sp, #64
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_sub_inner
	add	r0, sp, #64
	add	r1, sp, #64
	add	r2, sp, #32
	bl	gf_sub_inner
	add	r0, sp, #64
	add	r1, sp, #64
	bl	gf_half_inner

	@ t4 <- t3^2
	add	r0, sp, #96
	add	r1, sp, #64
	bl	gf_sqr_inner

	@ t5 <- W1*W2
	add	r0, sp, #128
	ldr	r1, [sp, #332]
	adds	r1, #32
	ldr	r2, [sp, #336]
	adds	r2, #32
	bl	gf_mul_inner

	@ t6 <- X1*X2
	add	r0, sp, #160
	ldr	r1, [sp, #332]
	ldr	r2, [sp, #336]
	bl	gf_mul_inner

	@ t7 <- (W1 + Z1)*(W2 + Z2) - t3 - t5
	add	r0, sp, #192
	ldr	r1, [sp, #332]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_add_inner
	add	r0, sp, #224    @ We use X3 (on the stack) as temporary
	ldr	r1, [sp, #336]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_add_inner
	add	r0, sp, #192
	add	r1, sp, #192
	add	r2, sp, #224
	bl	gf_mul_inner
	add	r0, sp, #192
	add	r1, sp, #192
	add	r2, sp, #64
	bl	gf_sub_inner
	add	r0, sp, #192
	add	r1, sp, #192
	add	r2, sp, #128
	bl	gf_sub_inner

	@ t8 <- (X1 + t1)*(X2 + t2) - t4 - t6
	mov	r0, sp
	ldr	r1, [sp, #332]
	mov	r2, sp
	bl	gf_add_inner
	add	r0, sp, #32
	ldr	r1, [sp, #336]
	add	r2, sp, #32
	bl	gf_add_inner
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #32
	bl	gf_mul_inner
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #96
	bl	gf_sub_inner
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #160
	bl	gf_sub_inner

	@ Z3 <- (t6 - b*t4)*t7
	@ Also, replace t4 with b*t4.
	add	r0, sp, #96
	add	r1, sp, #96
	bl	gf_half_inner
	add	r0, sp, #32
	add	r1, sp, #160
	add	r2, sp, #96
	bl	gf_sub_inner
	add	r0, sp, #288
	add	r1, sp, #32
	add	r2, sp, #192
	bl	gf_mul_inner

	@ t9 <- t7^4
	add	r0, sp, #32
	add	r1, sp, #192
	bl	gf_sqr_inner
	add	r0, sp, #32
	add	r1, sp, #32
	bl	gf_sqr_inner

	@ X3 <- b*t6*t9
	add	r0, sp, #224
	add	r1, sp, #160
	add	r2, sp, #32
	bl	gf_mul_inner
	add	r0, sp, #224
	add	r1, sp, #224
	bl	gf_half_inner

	@ t10 <- (t5 + a*t3)*(t6 + b*t4)
	@ We have a = -1. b*t4 is already computed (in t4).
	@ We use t5 and t6 as scratch.
	add	r0, sp, #128
	add	r1, sp, #128
	add	r2, sp, #64
	bl	gf_sub_inner
	add	r0, sp, #160
	add	r1, sp, #160
	add	r2, sp, #96
	bl	gf_add_inner
	add	r0, sp, #192
	add	r1, sp, #128
	add	r2, sp, #160
	bl	gf_mul_inner

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = 1/2, hence 2*b = 1.
	@ We use t8 as scratch.
	mov	r0, sp
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_mul_inner
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #192
	bl	gf_add_inner
	add	r0, sp, #256
	mov	r1, sp
	bl	gf_neg_inner

	@ If P1 is neutral, we replace P3 with P2.
	@ If P2 is neutral, we replace P3 with P1.
	ldr	r0, [sp, #328]
	ldr	r1, [sp, #332]
	ldr	r2, [sp, #336]
	add	r3, sp, #224
	ldr	r4, [sp, #324]
	ldr	r5, [sp, #320]
	mov	r8, r4
	mov	r10, r5
	bl	gf_sel3_inner
	mov	r4, r8
	mov	r5, r10
	bl	gf_sel3_inner
	mov	r4, r8
	mov	r5, r10
	bl	gf_sel3_inner

	add	sp, #340
	pop	{ pc }
	.size	point_add_inner, .-point_add_inner

@ =======================================================================
@ void point_add_mixed_inner(CURVE_point *P3,
@                            const CURVE_point *P1,
@                            const CURVE_point_affine *P2)
@
@ Cost: 185 + 8*cost(gf_mul_inner) + 3*cost(gf_sqr_inner)
@       + 2*cost(gf_iszero_inner) + 4*cost(gf_add_inner)
@       + 2*cost(gf_half_inner) + 2*cost(gf_sub_inner)
@       + cost(gf_neg_inner) + 3*cost(gf_sel3_inner)
@       = 15980
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_add_mixed_inner, %function
point_add_mixed_inner:
	push	{ r0, r1, r2, lr }
	sub	sp, #264

	@ Stack offsets:
	@   t1    0
	@   t5    32
	@   t6    64
	@   t7    96
	@   t8    128
	@   t9    96    (shared with t7)
	@   t10   0     (shared with t1)
	@   X3    160   (X3, W3 and Z3 must be consecutive in RAM)
	@   W3    192
	@   Z3    224
	@   fz1   256   (1 word)
	@   fz2   260   (1 word)
	@   &P3   264   (1 word)
	@   &P1   268   (1 word)
	@   &P2   272   (1 word)

	@ Test whether P1 and/or P2 is neutral; save the result.
	mov	r8, r2
	movs	r0, r1
	bl	gf_iszero_inner
	str	r0, [sp, #256]
	mov	r0, r8
	bl	gf_iszero_inner
	str	r0, [sp, #260]

	@ t1 <- Z1^2
	mov	r0, sp
	ldr	r1, [sp, #268]
	adds	r1, #64
	bl	gf_sqr_inner

	@ t5 <- W1*W2
	add	r0, sp, #32
	ldr	r1, [sp, #268]
	adds	r1, #32
	ldr	r2, [sp, #272]
	adds	r2, #32
	bl	gf_mul_inner

	@ t6 <- X1*X2
	add	r0, sp, #64
	ldr	r1, [sp, #268]
	ldr	r2, [sp, #272]
	bl	gf_mul_inner

	@ t7 <- W1 + W2*Z1
	add	r0, sp, #96
	ldr	r1, [sp, #272]
	adds	r1, #32
	ldr	r2, [sp, #268]
	adds	r2, #64
	bl	gf_mul_inner
	add	r0, sp, #96
	ldr	r1, [sp, #268]
	adds	r1, #32
	add	r2, sp, #96
	bl	gf_add_inner

	@ t8 <- X1 + X2*t1
	add	r0, sp, #128
	ldr	r1, [sp, #272]
	mov	r2, sp
	bl	gf_mul_inner
	add	r0, sp, #128
	ldr	r1, [sp, #268]
	add	r2, sp, #128
	bl	gf_add_inner

	@ Z3 <- (t6 - b*t1)*t7
	@ Also, replace t1 with b*t1.
	mov	r0, sp
	mov	r1, sp
	bl	gf_half_inner
	add	r0, sp, #224
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_sub_inner
	add	r0, sp, #224
	add	r1, sp, #224
	add	r2, sp, #96
	bl	gf_mul_inner

	@ t9 <- t7^4
	add	r0, sp, #96
	add	r1, sp, #96
	bl	gf_sqr_inner
	add	r0, sp, #96
	add	r1, sp, #96
	bl	gf_sqr_inner

	@ X3 <- b*t6*t9
	add	r0, sp, #160
	add	r1, sp, #64
	add	r2, sp, #96
	bl	gf_mul_inner
	add	r0, sp, #160
	add	r1, sp, #160
	bl	gf_half_inner

	@ t10 <- (t5 + a*t3)*(t6 + b*t1)
	@ We have a = -1 and t3 = Z1. b*t1 is already computed (in t1).
	@ We use t5 and t6 as scratch.
	add	r0, sp, #32
	add	r1, sp, #32
	ldr	r2, [sp, #268]
	adds	r2, #64
	bl	gf_sub_inner
	add	r0, sp, #64
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_add_inner
	mov	r0, sp
	add	r1, sp, #32
	add	r2, sp, #64
	bl	gf_mul_inner

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = 1/2, hence 2*b = 1.
	@ We use t8 as scratch.
	add	r0, sp, #128
	ldr	r1, [sp, #268]
	adds	r1, #64
	add	r2, sp, #128
	bl	gf_mul_inner
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #128
	bl	gf_add_inner
	add	r0, sp, #192
	mov	r1, sp
	bl	gf_neg_inner

	@ If P1 is neutral, we replace P3 with P2.
	@ If P2 is neutral, we replace P3 with P1.
	ldr	r0, [sp, #264]
	ldr	r1, [sp, #268]
	ldr	r2, [sp, #272]
	add	r3, sp, #160
	ldr	r4, [sp, #260]
	ldr	r5, [sp, #256]
	mov	r8, r4
	mov	r10, r5
	bl	gf_sel3_inner
	mov	r4, r8
	mov	r5, r10
	bl	gf_sel3_inner
	adr	r2, const_point_add_mixed_gfone
	mov	r4, r8
	mov	r5, r10
	bl	gf_sel3_inner

	add	sp, #276
	pop	{ pc }
	.align	2
const_point_add_mixed_gfone:
	.long	1, 0, 0, 0, 0, 0, 0, 0
	.size	point_add_mixed_inner, .-point_add_mixed_inner

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_double(CURVE_point *P3, const CURVE_point *P1)
@ =======================================================================

	.align	1
	.global	CN(double)
	.thumb
	.thumb_func
	.type	CN(double), %function
CN(double):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	point_double_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(double), .-CN(double)

@ =======================================================================
@ void CURVE_double_x(CURVE_point *P3, const CURVE_point *P1, unsigned n)
@ =======================================================================

	.align	1
	.global	CN(double_x)
	.thumb
	.thumb_func
	.type	CN(double_x), %function
CN(double_x):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	@ If n == 0, no doubling; just copy the point (if necessary).
	cmp	r2, #0
	bne	.Lpoint_double_x_step2
	cmp	r0, r1
	beq	.Lpoint_double_x_exit
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	b	.Lpoint_double_x_exit

	@ Perform requested doublings. First one still has r1 pointing
	@ to P1; afterwards, we use destination point (P3).
.Lpoint_double_x_step2:
	push	{ r0, r2 }
	bl	point_double_inner
.Lpoint_double_x_loop:
	pop	{ r0, r2 }
	subs	r2, #1
	beq	.Lpoint_double_x_exit
	push	{ r0, r2 }
	movs	r1, r0
	bl	point_double_inner
	b	.Lpoint_double_x_loop

.Lpoint_double_x_exit:
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(double), .-CN(double)

@ =======================================================================
@ void CURVE_add(CURVE_point *P3,
@                const CURVE_point *P1, const CURVE_point *P2)
@ =======================================================================

	.align	1
	.global	CN(add)
	.thumb
	.thumb_func
	.type	CN(add), %function
CN(add):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	point_add_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(add), .-CN(add)

@ =======================================================================
@ void CURVE_add_mixed(CURVE_point *P3,
@                      const CURVE_point *P1, const CURVE_point_affine *P2)
@ =======================================================================

	.align	1
	.global	CN(add_mixed)
	.thumb
	.thumb_func
	.type	CN(add_mixed), %function
CN(add_mixed):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	point_add_mixed_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(add_mixed), .-CN(add_mixed)
