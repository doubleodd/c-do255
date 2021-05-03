@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines functions for point addition and doubling in curve do255e.
@ =======================================================================

@ =======================================================================
@ void psi_x2_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Compute P3 = iso(psi(psi(P1)))
@ psi() is the isogeny from group E(a,b) into group E_QR(-2*a,a^2-4*b).
@ iso() is the isomorphism from E(4*a,16*b) into E(a,b)
@ This function is _almost_ a double, except that its output is in
@ E_QR(a,b) instead of our usual group E_nonQR(a,b).
@
@ Since a = 0, and b does not intervene in the computations, this
@ function also work for an almost-double from E(-2*a,a^2-4*b)
@ to E(-2*a,a^2-4*b).
@
@ Cost: 192 + 4*cost(gf_sqr_acc_inner) + cost(gf_mul_to_acc_inner)
@       + cost(gf_mul_to_acc_inner_altentry) + cost(gf_sub_acc_inner)
@       + cost(gf_rsub_acc_inner)
@       = 1073
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	psi_x2_inner, %function
psi_x2_inner:
	@ Formulas:
	@   X' = W^4
	@   W' = W^2 - 2*X
	@   Z' = W*Z
	@   X'' = W'^4
	@   W'' = W'^2 - 2*X'
	@   Z'' = 2*W'*Z'
	@ On the ARM Cortex M4, the relative costs of multiplications,
	@ squarings and additions/subtractions make it so that it is not
	@ worthwhile to compute Z'' as:
	@   Z'' = 2*Z*((W + W')^2 - W^2 - W'^2)
	@ Thus, we use 2M+4S formulas instead of 1M+5S.

	push	{ r0, r1, lr }
	sub	sp, #4

	@ Compute W^2 and save it to stack(64)
	adds	r1, #32
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	ACC_PUSH

	@ 2*X' = 2*W^4 -> stack(32)
	bl	gf_sqr_acc_inner
	ACC_MULCONST  2, r2, r3
	ACC_PUSH

	@ W' = W^2 - 2*X -> acc, stack(64)
	ldr	r0, [sp, #72]
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #32
	bl	gf_rsub_acc_inner
	add	r0, sp, #32
	ACC_STORE  r0

	@ W'^2 -> acc, stack(0)
	bl	gf_sqr_acc_inner
	ACC_PUSH

	@ X'' = W'^4 -> output
	bl	gf_sqr_acc_inner
	ldr	r0, [sp, #100]
	ACC_STORE  r0

	@ Z'' = 2*W*W'*Z -> output
	ldr	r1, [sp, #104]
	adds	r1, #32
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	@ We need to write down the upper accumulator words into the
	@ output buffer, for the alternate entry to the multiplication
	@ routine.
	ldr	r1, [sp, #100]
	adds	r1, #80
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #64
	bl	gf_mul_to_acc_inner_altentry
	ACC_MULCONST  2, r2, r3
	ldr	r0, [sp, #100]
	adds	r0, #64
	ACC_STORE  r0

	@ W'' = W'^2 - 2*X' -> output
	ACC_LOAD  sp
	add	r0, sp, #32
	bl	gf_sub_acc_inner
	ldr	r0, [sp, #100]
	adds	r0, #32
	ACC_STORE  r0

	add	sp, #108
	pop	{ pc }
	.size	psi_x2_inner, .-psi_x2_inner

@ =======================================================================
@ void point_double_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Cost: 207 + 5*cost(gf_sqr_acc_inner) + 2*cost(gf_mul_to_acc_inner)
@       + 2*cost(gf_rsub_acc_inner) + cost(gf_neg_acc_inner)
@       = 1226
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_double_inner, %function
point_double_inner:
	push	{ r0, r1, lr }
	sub	sp, #36

	@ Formulas:
	@   X' = W^4
	@   W' = W^2 - 2*X
	@   Z' = W*Z
	@   X'' = 16*b*Z'^4
	@   W'' = 2*X' - W'^2
	@   Z'' = 2*W'*Z'
	@ Total cost is 2M+5S. We could do it in 1M+6S by noticing that
	@ we internally compute Z'^2 and W'^2, so we could compute:
	@   Z'' = 2*W'*Z' = (W'+Z')^2 - Z'^2 - W'^2
	@ However, add+add+sub = 95 cycles while replacing a mul by
	@ a square saves only 57 cycles.

	@ Load W and compute W^2 (into stack(32)) and X' = W^4 (acc).
	adds	r1, #32
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	ACC_PUSH
	bl	gf_sqr_acc_inner

	@ 2*X' -> stack(0)
	ACC_MULCONST  2, r2, r3
	ACC_PUSH

	@ Load X and compute W' = W^2 - 2*X (acc and stack(32)).
	ldr	r0, [sp, #104]
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #32
	bl	gf_rsub_acc_inner
	add	r0, sp, #32
	ACC_STORE  r0

	@ W'' = 2*X' - W'^2 -> stack(0)
	bl	gf_sqr_acc_inner
	mov	r0, sp
	bl	gf_rsub_acc_inner
	ACC_STORE  sp

	@ 2*Z' = 2*W*Z -> acc and stack(64)
	ldr	r1, [sp, #104]
	adds	r1, #32
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #64
	ACC_STORE  r0

	@ X'' = -32*Z'^4 -> output.
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	ACC_MULCONST  2, r2, r3
	bl	gf_neg_acc_inner
	ldr	r0, [sp, #100]
	ACC_STORE  r0

	@ W'' -> output.
	ACC_LOAD  sp
	adds	r0, #32
	ACC_STORE  r0

	@ Z'' = 2*Z'*W'  (2*Z' is in stack(64), and W' is in stack(32))
	add	r1, sp, #64
	add	r2, sp, #32
	bl	gf_mul_to_acc_inner

	@ Z'' -> output.
	ldr	r0, [sp, #100]
	adds	r0, #64
	ACC_STORE  r0

	add	sp, #108
	pop	{ pc }
	.size	point_double_inner, .-point_double_inner

@ =======================================================================
@ void CURVE_double_x(CURVE_point *P3, const CURVE_point *P1, unsigned n)
@
@ This function is global.
@
@ First n-1 doublings use psi_x2_inner(), last one uses point_double_inner().
@ =======================================================================

	.align	1
	.global	CN(double_x)
	.thumb
	.thumb_func
	.type	CN(double_x), %function
CN(double_x):
	@ If n <= 1 (unusual case), jump to specific code (at the end).
	@ Also set r2 to n-1.
	subs	r2, #1
	bls	.Ldouble_x_single

	@ Make first n-1 doublings with psi_x2_inner().
	push	{ r0, r2, r4, r5, r6, r7, r8, r10, r11, lr }
.Ldouble_x_loop:
	bl	psi_x2_inner
	ldrd	r0, r2, [sp]
	movs	r1, r0
	subs	r2, #1
	beq	.Ldouble_x_last
	str	r2, [sp, #4]
	b	.Ldouble_x_loop

.Ldouble_x_last:
	@ Last doubling uses the dedicated function.
	bl	point_double_inner
	add	sp, #8
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }

.Ldouble_x_single:
	@ If a single doubling is requested, use the dedicated function
	@ (tail call).
	bne	.Ldouble_x_none
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	point_double_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }

.Ldouble_x_none:
	@ If there is no doubling, we may still have to copy the source
	@ to the destination (if addresses don't match).
	cmp	r0, r1
	beq	.Ldouble_x_none_exit
	push	{ r4, r5, r6, r7, r8 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r12 }
	pop	{ r4, r5, r6, r7, r8 }
.Ldouble_x_none_exit:
	bx	lr
	.size	CN(double_x), .-CN(double_x)

@ =======================================================================
@ void CURVE_double_x_xu(CURVE_point_xu *P3,
@                        const CURVE_point_xu *P1, unsigned n)
@
@ This function is global.
@
@ Cost: 314 + 3*cost(gf_mul_to_acc_inner)
@       + 2*cost(gf_mul_to_acc_inner_altentry) + 4*cost(gf_sqr_acc_inner)
@       + cost(gf_add_acc_inner) + cost(gf_sub_acc_inner)
@       + cost(gf_rsub_acc_inner) + cost(gf_neg_acc_inner)
@       + (n-1)*(14 + cost(psi_x2_inner))
@       = 670 + n*1088
@ =======================================================================

	.align	1
	.global	CN(double_x_xu)
	.thumb
	.thumb_func
	.type	CN(double_x_xu), %function
CN(double_x_xu):
	push	{ r0, r1, r2, r4, r5, r6, r7, r8, r10, r11, lr }
	sub	sp, #36

	@ Stack layout:
	@     0   t1
	@    36   pointer to P3
	@    40   pointer to P1
	@    44   number of doubles to compute

	@ If n == 0, jump to specific code that copies the input to the
	@ output (at the end).
	cmp	r2, #0
	beq	.Lpoint_double_x_xu_copy

	@ First half-doubling, with conversion to Jacobian (x,w).
	@   X' = Z^2*T^4
	@   W' = Z*T^2 - (2*X + a*Z)*U^2
	@   Z' = Z*U*T
	@ (a = 0 on curve do255e).
	@ Cost: 4M+2S
	@
	@ We store X', W' and Z' in the three first slots of P3. We take
	@ care to make operations in an order that still works if P1 and
	@ P3 are the same point.

	@ 2*X*U^2 -> stack(0)
	adds	r1, #64
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #16
	stm	r1, { r8, r10, r11, r12 }
	ldr	r2, [sp, #40]
	bl	gf_mul_to_acc_inner_altentry
	ACC_MULCONST  2, r2, r3
	ACC_STORE  sp

	@ Z*T -> acc, out:X (may overwrite in:X)
	ldr	r1, [sp, #40]
	adds	r1, #32
	add	r2, r1, #64
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #36]
	ACC_STORE  r1

	@ Z' = Z*U*T -> out:U
	adds	r1, #16
	ldr	r2, [sp, #40]
	adds	r2, #64
	bl	gf_mul_to_acc_inner_altentry
	ldr	r1, [sp, #36]
	add	r2, r1, #64
	ACC_STORE  r2

	@ Z*T^2 -> acc, out:X
	ldr	r2, [sp, #40]
	adds	r2, #96
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #36]
	ACC_STORE  r1

	@ W' = Z*T^2 - 2*X*U^2 -> out:Z
	mov	r0, sp
	bl	gf_sub_acc_inner
	ldr	r0, [sp, #36]
	add	r1, r0, #32
	ACC_STORE  r1

	@ X' = Z^2*T^4 -> out:X
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	ldr	r0, [sp, #36]
	ACC_STORE  r0

	@ Do n-1 inner doublings (with psi(psi)). We temporarily
	@ release our stack buffer to save space.
	add	sp, #32
.Lpoint_double_x_xu_loop:
	ldr	r2, [sp, #12]
	subs	r2, #1
	beq	.Lpoint_double_x_xu_finish
	str	r2, [sp, #12]
	ldr	r0, [sp, #4]
	movs	r1, r0
	bl	psi_x2_inner
	b	.Lpoint_double_x_xu_loop

.Lpoint_double_x_xu_finish:
	@ Final half doubling, with conversion back to fractional (x,u).
	@   X' = 4*b*Z^2
	@   Z' = W^2
	@   U' = 2*W*Z
	@   T' = 2*X - 2*a*Z^2 - W^2
	@ (a = 0 and b = -2 for curve do255e).
	@
	@ Nominal cost is 3S (using 2*W*Z = (W+Z)^2 - W^2 - Z^2) but on
	@ the ARM Cortex M4, the relative costs of multiplications and
	@ squarings are such that it's cheaper to do 2*W*Z directly.

	@ Reallocate the stack buffer.
	sub	sp, #32

	@ U' = 2*W*Z -> stack(0)
	ldr	r1, [sp, #36]
	adds	r1, #32
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	ACC_MULCONST  2, r2, r3
	ACC_STORE  sp

	@ Z' = W^2 -> acc, out:Z  (overwrites in:W)
	ldr	r0, [sp, #36]
	adds	r0, #32
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	ldr	r0, [sp, #36]
	adds	r0, #32
	ACC_STORE  r0

	@ T' = 2*X - W^2 -> out:T
	ldr	r0, [sp, #36]
	bl	gf_rsub_acc_inner
	ldr	r0, [sp, #36]
	bl	gf_add_acc_inner
	ldr	r0, [sp, #36]
	adds	r0, #96
	ACC_STORE  r0

	@ X' = -8*Z^2 -> out:X  (overwrites in:X)
	ldr	r0, [sp, #36]
	adds	r0, #64
	ACC_LOAD  r0
	bl	gf_sqr_acc_inner
	ACC_MULCONST  8, r2, r3
	bl	gf_neg_acc_inner
	ldr	r0, [sp, #36]
	ACC_STORE  r0

	@ U' -> out:U  (overwrites in:Z)
	ACC_LOAD  sp
	adds	r0, #64
	ACC_STORE  r0

	@ Done, go to exit.
	b	.Lpoint_double_x_xu_exit

.Lpoint_double_x_xu_copy:
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	stm	r0!, { r2, r3, r4, r5, r6, r7, r8, r10 }
	@ Fall through to exit.

.Lpoint_double_x_xu_exit:
	add	sp, #48
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(double_x_xu), .-CN(double_x_xu)

@ =======================================================================
@ void point_add_inner(CURVE_point *P3,
@                      const CURVE_point *P1, const CURVE_point *P2)
@
@ Cost: 484 + 2*cost(gf_iszero_acc_inner) + 5*cost(gf_sqr_acc_inner)
@       + 4*cost(gf_mul_to_acc_inner) + 5*cost(gf_mul_to_acc_inner_altentry)
@       + 5*cost(gf_add_acc_inner) + 6*cost(gf_sub_acc_inner)
@       + cost(gf_neg_acc_inner) + 3*cost(gf_sel3_inner)
@       = 3267
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
	@ Also, replace t4 with -b*t4.
	add	r0, sp, #96
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	ACC_STORE  r0
	add	r0, sp, #160
	bl	gf_add_acc_inner
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
	ACC_MULCONST  2, r2, r3
	bl	gf_neg_acc_inner
	add	r0, sp, #224
	ACC_STORE  r0

	@ t10 <- (t5 + a*t3)*(t6 + b*t4)
	@ We have a = 0. -b*t4 is already computed (in t4).
	@ We use t6 as scratch.
	add	r0, sp, #160
	ACC_LOAD  r0
	add	r0, sp, #96
	bl	gf_sub_acc_inner
	add	r1, sp, #208
	stm	r1, { r8, r10, r11, r12 }
	add	r2, sp, #128
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #192
	ACC_STORE  r0

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = -2, hence 2*b = -4.
	@ We use t8 as scratch.
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_mul_to_acc_inner
	ACC_MULCONST  4, r2, r3
	add	r0, sp, #192
	bl	gf_sub_acc_inner
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
@ void point_add_mixed_xu_inner(CURVE_point_xu *P3,
@                            const CURVE_point_xu *P1,
@                            const CURVE_point_affine_xu *P2)
@
@ Cost: 366 + 5*cost(gf_mul_to_acc_inner)
@       + 3*cost(gf_mul_to_acc_inner_altentry) + 7*cost(gf_add_acc_inner)
@       + 2*cost(gf_sub_acc_inner) + cost(gf_rsub_acc_inner)
@       + cost(gf_neg_acc_inner)
@       = 2074
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_add_mixed_xu_inner, %function
point_add_mixed_xu_inner:
	push	{ r0, r1, r2, lr }
	sub	sp, #128

	@ Stack offsets:
	@   v1    0
	@   v3    32
	@   v5    64
	@   v6    alias on X3 (not on stack)
	@   v7    0     (shared with v1)
	@   v8    96
	@   v9    64    (shared with v5)
	@   v10   32    (shared with v3)
	@ parameters:
	@   [sp, #128]   pointer to P3 (X3:Z3:U3:T3)
	@   [sp, #132]   pointer to P1 (X1:Z1:U1:T1)
	@   [sp, #136]   pointer to P2 (X2:U2)

	@ t1 <- X1*X2  (into v1)
	bl	gf_mul_to_acc_inner
	ACC_STORE  sp

	@ t2 = Z1*Z2 = Z1

	@ t3 <- U1*U2  (into v3)
	ldrd	r1, r2, [sp, #132]
	adds	r1, #64
	adds	r2, #32
	bl	gf_mul_to_acc_inner
	add	r0, sp, #32
	ACC_STORE  r0

	@ t4 = T1*T2 = T1

	@ t5 <- X1*Z2 + X2*Z1 = X1 + X2*Z1  (into v5)
	ldrd	r1, r2, [sp, #132]
	adds	r1, #32
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #132]
	bl	gf_add_acc_inner
	add	r0, sp, #64
	ACC_STORE  r0

	@ t6 <- U1*T2 + U2*T1 = U1 + U2*T1
	@ U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*Z1)
	@ We use X3 as scratch (X1 won't be needed any more)
	@ We also want 2*Z1 into v8 on output

	@ *** X3 <- -t6 = -(U1 + U2*T1)
	ldrd	r1, r2, [sp, #132]
	adds	r1, #96
	adds	r2, #32
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #132]
	adds	r0, #64
	bl	gf_add_acc_inner
	bl	gf_neg_acc_inner
	ldr	r0, [sp, #128]
	ACC_STORE  r0
	@ *** U3 <- -t6*(t1 + 2*Z1)
	@ with copy of 2*Z1 into v8
	ldr	r0, [sp, #132]
	adds	r0, #32
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #96
	ACC_STORE  r0
	mov	r0, sp
	bl	gf_add_acc_inner
	ldr	r2, [sp, #128]
	add	r1, r2, #80
	stm	r1, { r8, r10, r11, r12 }
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #128]
	adds	r0, #64
	ACC_STORE  r0

	@ t7 <- t1 + b*t2  (with b = -2 and t2 = Z1)  (into v7)
	@ 2*Z1 is still in v8
	mov	r0, sp
	ACC_LOAD  sp
	add	r0, sp, #96
	bl	gf_sub_acc_inner
	ACC_STORE  sp

	@ t8 <- t4*t7  (with t4 = T1)  (into v8)
	add	r1, sp, #16     @ t7 already fully stored there
	ldr	r2, [sp, #132]
	adds	r2, #96
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #96
	ACC_STORE  r0

	@ t5 + t7  (into v7)
	ACC_LOAD  sp
	add	r0, sp, #64
	bl	gf_add_acc_inner
	ACC_STORE  sp

	@ t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2)
	@ We actually store -t9 into v9.
	add	r1, sp, #32
	add	r2, sp, #64
	bl	gf_mul_to_acc_inner
	ACC_MULCONST  4, r2, r3
	add	r0, sp, #64
	ACC_STORE  r0

	@ t10 <- (t4 + alpha*t3)*(t5 + t7)  (t4 = T1, alpha = 2)  (into acc)
	add	r0, sp, #32
	ACC_LOAD  r0
	ACC_MULCONST  2, r2, r3
	ldr	r0, [sp, #132]
	adds	r0, #96
	bl	gf_add_acc_inner
	add	r1, sp, #48
	stm	r1, { r8, r10, r11, r12 }
	mov	r2, sp
	bl	gf_mul_to_acc_inner_altentry

	@ X3 <- b*(t10 - t8 + beta*t9) = 2*(t8 - t10) - t9  (v9 contains -t9)
	add	r0, sp, #96
	bl	gf_rsub_acc_inner
	ACC_MULCONST  2, r2, r3
	add	r0, sp, #64
	bl	gf_add_acc_inner
	ldr	r0, [sp, #128]
	ACC_STORE  r0

	@ Z3 <- t8 - t9   (v9 contains -t9)
	add	r0, sp, #96
	ACC_LOAD  r0
	add	r0, sp, #64
	bl	gf_add_acc_inner
	ldr	r0, [sp, #128]
	adds	r0, #32
	ACC_STORE  r0

	@ T3 <- t8 + t9   (v9 contains -t9)
	add	r0, sp, #96
	ACC_LOAD  r0
	add	r0, sp, #64
	bl	gf_sub_acc_inner
	ldr	r0, [sp, #128]
	adds	r0, #96
	ACC_STORE  r0

	add	sp, #140
	pop	{ pc }
	.size	point_add_mixed_xu_inner, .-point_add_mixed_xu_inner

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
@ void CURVE_add_mixed_xu(CURVE_point_xu *P3, const CURVE_point_xu *P1,
@                         const CURVE_point_affine_xu *P2)
@ =======================================================================

	.align	1
	.global	CN(add_mixed_xu)
	.thumb
	.thumb_func
	.type	CN(add_mixed_xu), %function
CN(add_mixed_xu):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	point_add_mixed_xu_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(add_mixed_xu), .-CN(add_mixed_xu)
