@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines functions for point addition and doubling in curve do255s.
@ =======================================================================

@ =======================================================================
@ void point_double_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Cost: 240 + cost(gf_mul_to_acc_inner) + cost(gf_mul_to_acc_inner_altentry)
@       + 4*cost(gf_sqr_acc_inner) + cost(gf_add_acc_inner)
@       + 2*cost(gf_sub_acc_inner) + cost(gf_rsub_acc_inner)
@       = 1186
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
	adds	r1, #32
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	ACC_STORE  sp

	@ t2 <- t1^2
	bl	gf_sqr_acc_inner
	add	r0, sp, #32
	ACC_STORE  r0

	@ replace t1 with 2*t1
	ACC_LOAD  sp
	ACC_MULCONST  2, r2, r3
	ACC_STORE  sp

	@ t3 <- (W + Z)^2 - 2*t1
	@ (2*t1 is already computed)
	ldr	r0, [sp, #104]
	adds	r0, #32
	ACC_LOAD_UPDATE  r0
	bl	gf_add_acc_inner
	bl	gf_sqr_acc_inner
	mov	r0, sp
	bl	gf_sub_acc_inner
	add	r0, sp, #64
	ACC_STORE  r0

	@ Z' <- 2*t1*(2*X - t3)
	@ The formal value '2*t1' is already in t1.
	ldr	r0, [sp, #104]
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #64
	bl	gf_sub_acc_inner
	ldr	r1, [sp, #100]
	adds	r1, #80
	stm	r1, { r8, r10, r11, r12 }
	mov	r2, sp
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #100]
	adds	r0, #64
	ACC_STORE  r0

	@ X' <- 8*t2^2
	@ We also store 2*t2 in t1.
	add	r0, sp, #32
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	ACC_STORE  sp
	bl	gf_sqr_acc_inner
	ACC_MULCONST  2, r2, r3
	ldr	r0, [sp, #100]
	ACC_STORE  r0

	@ W' <- 2*t2 - t3^2
	add	r0, sp, #64
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	mov	r0, sp
	bl	gf_rsub_acc_inner
	ldr	r0, [sp, #100]
	adds	r0, #32
	ACC_STORE  r0

	add	sp, #108
	pop	{ pc }
	.size	point_double_inner, .-point_double_inner

@ =======================================================================
@ void point_add_inner(CURVE_point *P3,
@                      const CURVE_point *P1, const CURVE_point *P2)
@
@ Cost: 476 + 2*cost(gf_iszero_acc_inner) + 5*cost(gf_sqr_acc_inner)
@       + 4*cost(gf_mul_to_acc_inner) + 5*cost(gf_mul_to_acc_inner_altentry)
@       + 6*cost(gf_add_acc_inner) + 5*cost(gf_sub_acc_inner)
@       + cost(gf_rsub_acc_inner) + 2*cost(gf_half_acc_inner)
@       + cost(gf_neg_acc_inner) + 3*cost(gf_sel3_inner)
@       = 3310
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_add_inner, %function
point_add_inner:
	@ We use 9M+5S formulas, computing t3 = Z1*Z2 with a multiplication
	@ instead of ((Z1+Z2)^2 - Z1^2 - Z2^2)/2, because the extra cost
	@ of a multiplication (compared to a squaring) is more than
	@ compensated by avoiding two additions, one subtraction and one
	@ halving.

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

	@ t1 <- Z1^2
	@ Also test whether P1 is neutral (saved on stack).
	adds	r1, #64
	ACC_LOAD  r1
	bl	gf_iszero_acc_inner
	str	r0, [sp, #320]
	bl	gf_sqr_acc_inner
	ACC_STORE  sp

	@ t2 <- Z2^2
	@ Also test whether P2 is neutral (saved on stack).
	ldr	r2, [sp, #336]
	adds	r2, #64
	ACC_LOAD  r2
	bl	gf_iszero_acc_inner
	str	r0, [sp, #324]
	bl	gf_sqr_acc_inner
	add	r2, sp, #32
	ACC_STORE  r2

	@ t3 <- Z1*Z2
	ldrd	r1, r2, [sp, #332]
	adds	r1, #64
	adds	r2, #64
	bl	gf_mul_to_acc_inner
	add	r0, sp, #64
	ACC_STORE  r0

	@ t4 <- t3^2
	bl	gf_sqr_acc_inner
	add	r0, sp, #96
	ACC_STORE  r0

	@ t5 <- W1*W2
	ldrd	r1, r2, [sp, #332]
	adds	r1, #32
	adds	r2, #32
	bl	gf_mul_to_acc_inner
	add	r0, sp, #128
	ACC_STORE  r0

	@ t6 <- X1*X2
	ldrd	r1, r2, [sp, #332]
	bl	gf_mul_to_acc_inner
	add	r0, sp, #160
	ACC_STORE  r0

	@ t7 <- (W1 + Z1)*(W2 + Z2) - t3 - t5
	@ We use X3 (stack) as temporary.
	ldr	r0, [sp, #332]
	adds	r0, #32
	ACC_LOAD  r0
	adds	r0, #32
	bl	gf_add_acc_inner
	add	r0, sp, #192
	ACC_STORE  r0
	ldr	r0, [sp, #336]
	adds	r0, #32
	ACC_LOAD  r0
	adds	r0, #32
	bl	gf_add_acc_inner
	add	r1, sp, #240
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #192
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #64
	bl	gf_sub_acc_inner
	add	r0, sp, #128
	bl	gf_sub_acc_inner
	add	r0, sp, #192
	ACC_STORE  r0

	@ t8 <- (X1 + t1)*(X2 + t2) - t4 - t6
	@ t1 and t2 won't be used afterwards; t8 is an alias on t1,
	@ and we can use t2 as temporary.
	ACC_LOAD  sp
	ldr	r0, [sp, #332]
	bl	gf_add_acc_inner
	ACC_STORE  sp
	add	r0, sp, #32
	ACC_LOAD  r0
	ldr	r0, [sp, #336]
	bl	gf_add_acc_inner
	add	r1, sp, #48
	stm	r1, { r8, r10, r11, r12 }
	mov	r2, sp
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #96
	bl	gf_sub_acc_inner
	add	r0, sp, #160
	bl	gf_sub_acc_inner
	ACC_STORE  sp

	@ Z3 <- (t6 - b*t4)*t7
	@ Also, replace t4 with b*t4.
	add	r0, sp, #96
	ACC_LOAD  r0
	bl	gf_half_acc_inner
	add	r0, sp, #96
	ACC_STORE  r0
	add	r0, sp, #160
	bl	gf_rsub_acc_inner
	add	r1, sp, #304
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #192
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #288
	ACC_STORE  r0

	@ t9 <- t7^4  (in acc, not stored)
	add	r0, sp, #192
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner

	@ X3 <- b*t6*t9
	add	r1, sp, #48
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #160
	bl	gf_mul_to_acc_inner_altentry
	bl	gf_half_acc_inner
	add	r0, sp, #224
	ACC_STORE  r0

	@ t10 <- (t5 + a*t3)*(t6 + b*t4)
	@ We have a = -1. b*t4 is already computed (in t4).
	@ t5 is consumed.
	add	r0, sp, #128
	ACC_LOAD  r0
	add	r0, sp, #64
	bl	gf_sub_acc_inner
	add	r0, sp, #128
	ACC_STORE  r0
	add	r0, sp, #160
	ACC_LOAD  r0
	add	r0, sp, #96
	bl	gf_add_acc_inner
	add	r1, sp, #208
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #128
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #192
	ACC_STORE  r0

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = 1/2.
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_mul_to_acc_inner
	add	r0, sp, #192
	bl	gf_add_acc_inner
	bl	gf_neg_acc_inner
	add	r0, sp, #256
	ACC_STORE  r0

	@ If P1 is neutral, we replace P3 with P2.
	@ If P2 is neutral, we replace P3 with P1.
	ldrd	r0, r1, [sp, #328]
	ldr	r2, [sp, #336]
	add	r3, sp, #224
	ldrd	r5, r4, [sp, #320]
	bl	gf_sel3_inner
	ldrd	r0, r1, [sp, #328]
	adds	r0, #32
	adds	r1, #32
	ldr	r2, [sp, #336]
	adds	r2, #32
	add	r3, sp, #256
	ldrd	r5, r4, [sp, #320]
	bl	gf_sel3_inner
	ldrd	r0, r1, [sp, #328]
	adds	r0, #64
	adds	r1, #64
	ldr	r2, [sp, #336]
	adds	r2, #64
	add	r3, sp, #288
	ldrd	r5, r4, [sp, #320]
	bl	gf_sel3_inner

	add	sp, #340
	pop	{ pc }
	.size	point_add_inner, .-point_add_inner

@ =======================================================================
@ void point_add_mixed_inner(CURVE_point *P3,
@                            const CURVE_point *P1,
@                            const CURVE_point_affine *P2)
@
@ Cost: 319 + 2*cost(gf_iszero_acc_inner) + 3*cost(gf_sqr_acc_inner)
@       + 5*cost(gf_mul_to_acc_inner) + 3*cost(gf_mul_to_acc_inner_altentry)
@       + 4*cost(gf_add_acc_inner) + cost(gf_sub_acc_inner)
@       + cost(gf_rsub_acc_inner) + 2*cost(gf_half_acc_inner)
@       + cost(gf_neg_acc_inner) + 3*cost(gf_sel3_inner)
@       = 2576
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

	@ Test whether P2 is neutral; save the result.
	ACC_LOAD  r2
	bl	gf_iszero_acc_inner
	str	r0, [sp, #260]

	@ t1 <- Z1^2
	@ Also test whether P1 is neutral (saved on stack).
	ldr	r0, [sp, #268]
	adds	r0, #64
	ACC_LOAD  r0
	bl	gf_iszero_acc_inner
	str	r0, [sp, #256]
	bl	gf_sqr_acc_inner
	ACC_STORE  sp

	@ t6 <- X1*X2
	ldrd	r1, r2, [sp, #268]
	bl	gf_mul_to_acc_inner
	add	r0, sp, #64
	ACC_STORE  r0

	@ t7 <- W1 + W2*Z1
	ldrd	r1, r2, [sp, #268]
	adds	r1, #64
	adds	r2, #32
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #268]
	adds	r0, #32
	bl	gf_add_acc_inner
	add	r0, sp, #96
	ACC_STORE  r0

	@ t8 <- X1 + X2*t1
	ldr	r1, [sp, #272]
	mov	r2, sp
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #268]
	bl	gf_add_acc_inner
	add	r0, sp, #128
	ACC_STORE  r0

	@ Z3 <- (t6 - b*t1)*t7
	@ Also, replace t1 with b*t1.
	ACC_LOAD  sp
	bl	gf_half_acc_inner
	ACC_STORE  sp
	add	r0, sp, #64
	bl	gf_rsub_acc_inner
	add	r1, sp, #240
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #96
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #224
	ACC_STORE  r0

	@ t9 <- t7^4  (into acc only)
	add	r0, sp, #96
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner

	@ X3 <- b*t6*t9
	add	r1, sp, #176
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #64
	bl	gf_mul_to_acc_inner_altentry
	bl	gf_half_acc_inner
	add	r0, sp, #160
	ACC_STORE  r0

	@ t5 <- W1*W2  (into acc only)
	ldrd	r1, r2, [sp, #268]
	adds	r1, #32
	adds	r2, #32
	bl	gf_mul_to_acc_inner

	@ t10 <- (t5 + a*t3)*(t6 + b*t1)  (into acc)
	@ We have a = -1 and t3 = Z1. b*t1 is already computed (in t1).
	@ We use t5 as scratch. t1 is consumed.
	ldr	r0, [sp, #268]
	adds	r0, #64
	bl	gf_sub_acc_inner
	add	r0, sp, #32
	ACC_STORE  r0
	ACC_LOAD  sp
	add	r0, sp, #64
	bl	gf_add_acc_inner
	add	r1, sp, #16
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #32
	bl	gf_mul_to_acc_inner_altentry
	ACC_STORE  sp

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = 1/2, hence 2*b = 1.
	@ We use t8 as scratch.
	ldr	r1, [sp, #268]
	adds	r1, #64
	add	r2, sp, #128
	bl	gf_mul_to_acc_inner
	mov	r0, sp
	bl	gf_add_acc_inner
	bl	gf_neg_acc_inner
	add	r0, sp, #192
	ACC_STORE  r0

	@ If P1 is neutral, we replace P3 with P2.
	@ If P2 is neutral, we replace P3 with P1.
	ldrd	r0, r1, [sp, #264]
	ldr	r2, [sp, #272]
	add	r3, sp, #160
	ldrd	r5, r4, [sp, #256]
	bl	gf_sel3_inner
	ldrd	r0, r1, [sp, #264]
	adds	r0, #32
	adds	r1, #32
	ldr	r2, [sp, #272]
	adds	r2, #32
	add	r3, sp, #192
	ldrd	r5, r4, [sp, #256]
	bl	gf_sel3_inner
	ldrd	r0, r1, [sp, #264]
	adds	r0, #64
	adds	r1, #64
	adr	r2, const_point_add_mixed_gfone
	add	r3, sp, #224
	ldrd	r5, r4, [sp, #256]
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	point_double_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }

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
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	point_add_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	point_add_mixed_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(add_mixed), .-CN(add_mixed)
