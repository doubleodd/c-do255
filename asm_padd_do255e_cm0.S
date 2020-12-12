@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines functions for point addition and doubling in curve do255e.
@ =======================================================================

@ =======================================================================
@ void psi_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Compute P3 = psi(P1)
@ psi() is the isogeny from group E(a,b) into group E_QR(-2*a,a^2-4*b).
@
@ Cost: 46 + cost(gf_mul_inner) + 2*cost(gf_sqr_inner)
@       + cost(gf_mul2_inner) + cost(gf_sub_inner)
@       = 3633
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	psi_inner, %function
psi_inner:
	push	{ r0, r1, lr }
	sub	sp, #36

	@ Z' = W*Z
	adds	r0, #64
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_mul_inner

	@ W' = W^2 - 2*X
	@ (also store W^2 into a stack buffer)
	mov	r0, sp
	ldr	r1, [sp, #40]
	adds	r1, #32
	bl	gf_sqr_inner
	ldr	r0, [sp, #36]
	adds	r0, #32
	ldr	r1, [sp, #40]
	bl	gf_mul2_inner
	ldr	r0, [sp, #36]
	adds	r0, #32
	mov	r1, sp
	movs	r2, r0
	bl	gf_sub_inner

	@ X' = W^4
	ldr	r0, [sp, #36]
	mov	r1, sp
	bl	gf_sqr_inner

	add	sp, #44
	pop	{ pc }
	.size	psi_inner, .-psi_inner

@ =======================================================================
@ void phi_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Compute P3 = phi(P1)
@ phi() is the isogeny from group E(a,b) into group E_nonQR(-a/2,(a^2/4)-b).
@
@ Cost: 98 + 4*cost(gf_sqr_inner) + cost(gf_add_inner)
@       + 3*cost(gf_sub_inner) + 2*cost(gf_half_inner)
@       + cost(gf_mul8_inner) + cost(gf_neg_inner)
@       = 4642
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	phi_inner, %function
phi_inner:
	push	{ r0, r1, lr }
	sub	sp, #36

	@ Compute Z^2 into stack buffer.
	mov	r0, sp
	adds	r1, #64
	bl	gf_sqr_inner

	@ Compute (W+Z)^2 (into Z').
	ldr	r0, [sp, #36]
	adds	r0, #64
	ldr	r1, [sp, #40]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_add_inner
	ldr	r0, [sp, #36]
	adds	r0, #64
	movs	r1, r0
	bl	gf_sqr_inner

	@ Compute W^2 (into W').
	ldr	r0, [sp, #36]
	adds	r0, #32
	ldr	r1, [sp, #40]
	adds	r1, #32
	bl	gf_sqr_inner

	@ Z' = W*Z = ((W+Z)^2 - W^2 - Z^2)/2
	@ (W+Z)^2 is in Z', W^2 in W', and Z^2 in the stack buffer.
	ldr	r0, [sp, #36]
	adds	r0, #64
	movs	r1, r0
	movs	r2, r0
	subs	r2, #32
	bl	gf_sub_inner
	ldr	r0, [sp, #36]
	adds	r0, #64
	movs	r1, r0
	mov	r2, sp
	bl	gf_sub_inner
	ldr	r0, [sp, #36]
	adds	r0, #64
	movs	r1, r0
	bl	gf_half_inner

	@ W' = X - (1/2)*W^2
	ldr	r0, [sp, #36]
	adds	r0, #32
	movs	r1, r0
	bl	gf_half_inner
	ldr	r0, [sp, #36]
	adds	r0, #32
	ldr	r1, [sp, #40]
	movs	r2, r0
	bl	gf_sub_inner

	@ X' = 4*b*Z^4 = -8*Z^4  (because b = -2 on curve do255e)
	mov	r0, sp
	mov	r1, sp
	bl	gf_sqr_inner
	mov	r0, sp
	mov	r1, sp
	bl	gf_mul8_inner
	ldr	r0, [sp, #36]
	mov	r1, sp
	bl	gf_neg_inner

	add	sp, #44
	pop	{ pc }
	.size	phi_inner, .-phi_inner

@ =======================================================================
@ void psi_x2_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Compute P3 = iso(psi(psi(P1)))
@ psi() is the isogeny from group E(a,b) into group E_QR(-2*a,a^2-4*b).
@ iso() is the isomorphism from E(4*a,16*b) into E(a,b)
@ This function is _almost_ a double, except that its output is in
@ E_QR(a,b) instead of our usual group E_nonQR(a,b).
@
@ Cost: 93 + cost(gf_mul_inner) + 5*cost(gf_sqr_inner)
@       + 2*cost(gf_mul2_inner) + 4*cost(gf_sub_inner)
@       + cost(gf_add_inner)
@       = 7015
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	psi_x2_inner, %function
psi_x2_inner:
	push	{ r0, r1, lr }
	sub	sp, #68

	@ W^2 (into stack (0))
	mov	r0, sp
	adds	r1, #32
	bl	gf_sqr_inner

	@ W' = W^2 - 2*X (into stack (32))
	add	r0, sp, #32
	ldr	r1, [sp, #72]
	bl	gf_mul2_inner
	add	r0, sp, #32
	mov	r1, sp
	add	r2, sp, #32
	bl	gf_sub_inner

	@ W'^2 (into out[X])
	ldr	r0, [sp, #68]
	add	r1, sp, #32
	bl	gf_sqr_inner

	@ (W+W')^2 (into stack (32))
	add	r0, sp, #32
	ldr	r1, [sp, #72]
	adds	r1, #32
	add	r2, sp, #32
	bl	gf_add_inner
	add	r0, sp, #32
	add	r1, sp, #32
	bl	gf_sqr_inner

	@ Z'' = 2*W*W'*Z = ((W+W')^2 - W^2 - W'^2)*Z (into out[Z])
	add	r0, sp, #32
	add	r1, sp, #32
	mov	r2, sp
	bl	gf_sub_inner
	add	r0, sp, #32
	add	r1, sp, #32
	ldr	r2, [sp, #68]
	bl	gf_sub_inner
	ldr	r0, [sp, #68]
	adds	r0, #64
	ldr	r1, [sp, #72]
	adds	r1, #64
	add	r2, sp, #32
	bl	gf_mul_inner

	@ 2*X' = 2*W^4 (into stack (0))
	mov	r0, sp
	mov	r1, sp
	bl	gf_sqr_inner
	mov	r0, sp
	mov	r1, sp
	bl	gf_mul2_inner

	@ W'' = W'^2 - 2*X' (into out[W])
	ldr	r1, [sp, #68]
	mov	r0, r1
	adds	r0, #32
	mov	r2, sp
	bl	gf_sub_inner

	@ X'' = W'^4 (into out[X])
	ldr	r0, [sp, #68]
	mov	r1, r0
	bl	gf_sqr_inner

	add	sp, #76
	pop	{ pc }
	.size	psi_x2_inner, .-psi_x2_inner

@ =======================================================================
@ void point_double_inner(CURVE_point *P3, const CURVE_point *P1)
@
@ Cost: 18 + cost(psi_inner) + cost(phi_inner)
@       = 8293
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	point_double_inner, %function
point_double_inner:
	push	{ r0, lr }

	@ (P+N)*(P+N) = phi(psi(P))
	bl	psi_inner
	ldr	r0, [sp]
	movs	r1, r0
	bl	phi_inner

	pop	{ r0, pc }
	.size	point_double_inner, .-point_double_inner

@ =======================================================================
@ void CURVE_double_x(CURVE_point *P3, const CURVE_point *P1, unsigned n)
@
@ This function is global.
@
@ First n-1 doublings use two psi() (and one iso()), for a cost per
@ iteration of 7028 cycles. Last doubling uses phi(psi()), which has
@ cost 8293 cycles.
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

	@ If n <= 1, jump to specific code (at the end).
	cmp	r2, #2
	bcc	.Lpoint_double_x_single

	@ Generic case. We first do n-1 doublings into E_QR(a,b)
	@ using iso(psi(psi(P))) = P+P (on the curve).
	subs	r2, #1
	push	{ r0, r2 }
.Lpoint_double_x_loop:
	bl	psi_x2_inner
	ldr	r0, [sp]
	movs	r1, r0
	ldr	r2, [sp, #4]
	subs	r2, #1
	str	r2, [sp, #4]
	bne	.Lpoint_double_x_loop

	@ Final doubling, with phi(psi(P))
	bl	psi_inner
	ldr	r0, [sp]
	movs	r1, r0
	bl	phi_inner
	add	sp, #8
	b	.Lpoint_double_x_exit

.Lpoint_double_x_single:
	@ If n == 1, this is a simple doubling.
	cmp	r2, #0
	beq	.Lpoint_double_x_none

	push	{ r0, r1 }
	bl	psi_inner
	ldr	r0, [sp]
	movs	r1, r0
	bl	phi_inner
	add	sp, #8
	b	.Lpoint_double_x_exit

.Lpoint_double_x_none:
	@ If n == 0, there is no doubling; we simply copy the input to
	@ the output if the pointers are distinct.
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
	@ Fall through to exit.

.Lpoint_double_x_exit:
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(double_x), .-CN(double_x)

@ =======================================================================
@ void CURVE_double_x_xu(CURVE_point_xu *P3,
@                        const CURVE_point_xu *P1, unsigned n)
@
@ This function is global.
@
@ Cost: 177 + 5*cost(gf_sqr_inner) + 4*cost(gf_mul_inner)
@       + 2*cost(gf_mul2_inner) + cost(gf_mul8_inner) + cost(gf_neg_inner)
@       + cost(gf_add_inner) + 4*cost(gf_sub_inner)
@       + (n-1)*(14 + cost(psi_x2_inner))
@       = 4609 + n*7015
@ =======================================================================

	.align	1
	.global	CN(double_x_xu)
	.thumb
	.thumb_func
	.type	CN(double_x_xu), %function
CN(double_x_xu):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r0, r1, r2, r4, r5, r6 }
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

	@ t1 <- 2*X*U^2
	mov	r0, sp
	adds	r1, #64
	bl	gf_sqr_inner
	mov	r0, sp
	mov	r1, sp
	ldr	r2, [sp, #40]
	bl	gf_mul_inner
	mov	r0, sp
	mov	r1, sp
	bl	gf_mul2_inner

	@ X' <- Z*T  (may overwrite X)
	ldr	r0, [sp, #36]
	ldr	r1, [sp, #40]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #64
	bl	gf_mul_inner

	@ Z' <- X'*U = Z*U*T  (may overwrite U)
	ldr	r1, [sp, #36]
	movs	r0, r1
	adds	r0, #64
	ldr	r2, [sp, #40]
	adds	r2, #64
	bl	gf_mul_inner

	@ X' <- X'*T = Z*T^2
	ldr	r0, [sp, #36]
	movs	r1, r0
	ldr	r2, [sp, #40]
	adds	r2, #96
	bl	gf_mul_inner

	@ W' <- X' - t1 = Z*T^2 - 2*X*U^2  (may overwrite Z)
	ldr	r1, [sp, #36]
	movs	r0, r1
	adds	r0, #32
	mov	r2, sp
	bl	gf_sub_inner

	@ X' <- X'^2 = Z^2*T^4
	ldr	r0, [sp, #36]
	movs	r1, r0
	bl	gf_sqr_inner

	@ Do n-1 inner doublings. Each has cost 1M+5S. We temporarily
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
	@ Cost: 3S

	@ Reallocate the stack buffer.
	sub	sp, #32

	@ T' <- 2*X
	ldr	r1, [sp, #36]
	movs	r0, r1
	adds	r0, #96
	bl	gf_mul2_inner

	@ t1 <- (W + Z)^2
	mov	r0, sp
	ldr	r1, [sp, #36]
	adds	r1, #32
	movs	r2, r1
	adds	r2, #32
	bl	gf_add_inner
	mov	r0, sp
	mov	r1, sp
	bl	gf_sqr_inner

	@ Z' <- W^2  (overwrites W)
	ldr	r0, [sp, #36]
	adds	r0, #32
	movs	r1, r0
	bl	gf_sqr_inner

	@ T' <- T' - Z' = 2*X - W^2
	ldr	r2, [sp, #36]
	adds	r2, #32
	movs	r0, r2
	adds	r0, #64
	movs	r1, r0
	bl	gf_sub_inner

	@ X' <- Z^2  (overwrites X)
	ldr	r0, [sp, #36]
	movs	r1, r0
	adds	r1, #64
	bl	gf_sqr_inner

	@ U' <- t1 - Z' - X' = 2*W*Z  (overwrites Z)
	ldr	r0, [sp, #36]
	movs	r2, r0
	adds	r0, #64
	mov	r1, sp
	adds	r2, #32
	bl	gf_sub_inner
	ldr	r0, [sp, #36]
	movs	r2, r0
	adds	r0, #64
	movs	r1, r0
	bl	gf_sub_inner

	@ X' <- -8*X' = 4*b*Z^2
	ldr	r0, [sp, #36]
	movs	r1, r0
	bl	gf_mul8_inner
	ldr	r0, [sp, #36]
	movs	r1, r0
	bl	gf_neg_inner

	@ Done, go to exit.
	b	.Lpoint_double_x_xu_exit

.Lpoint_double_x_xu_copy:
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3 }
	stm	r0!, { r2, r3 }
	@ Fall through to exit.

.Lpoint_double_x_xu_exit:
	add	sp, #48
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(double_x_xu), .-CN(double_x_xu)

@ =======================================================================
@ void point_add_inner(CURVE_point *P3,
@                      const CURVE_point *P1, const CURVE_point *P2)
@
@ Cost: 259 + 8*cost(gf_mul_inner) + 6*cost(gf_sqr_inner)
@       + 2*cost(gf_iszero_inner) + 6*cost(gf_add_inner)
@       + 8*cost(gf_sub_inner) + cost(gf_half_inner)
@       + 2*cost(gf_mul2_inner) + cost(gf_mul4_inner)
@       + cost(gf_neg_inner) + 3*cost(gf_sel3_inner)
@       = 19731
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
	@ Also, replace t4 with -b*t4.
	add	r0, sp, #96
	add	r1, sp, #96
	bl	gf_mul2_inner
	add	r0, sp, #32
	add	r1, sp, #160
	add	r2, sp, #96
	bl	gf_add_inner
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
	bl	gf_mul2_inner
	add	r0, sp, #224
	add	r1, sp, #224
	bl	gf_neg_inner

	@ t10 <- (t5 + a*t3)*(t6 + b*t4)
	@ We have a = 0. -b*t4 is already computed (in t4).
	@ We use t6 as scratch.
	add	r0, sp, #160
	add	r1, sp, #160
	add	r2, sp, #96
	bl	gf_sub_inner
	add	r0, sp, #192
	add	r1, sp, #128
	add	r2, sp, #160
	bl	gf_mul_inner

	@ W3 <- -t10 - 2*b*t3*t8
	@ We have b = -2, hence 2*b = -4.
	@ We use t8 as scratch.
	mov	r0, sp
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_mul_inner
	mov	r0, sp
	mov	r1, sp
	bl	gf_mul4_inner
	add	r0, sp, #256
	mov	r1, sp
	add	r2, sp, #192
	bl	gf_sub_inner

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
@ void point_add_mixed_xu_inner(CURVE_point_xu *P3,
@                            const CURVE_point_xu *P1,
@                            const CURVE_point_affine_xu *P2)
@
@ Cost: 174 + 8*cost(gf_mul_inner)
@       + 7*cost(gf_add_inner) + 3*cost(gf_mul2_inner)
@       + cost(gf_neg_inner) + 3*cost(gf_sub_inner)
@       + cost(gf_mul4_inner)
@       = 12903
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
	@   v7    0     (shared with t1)
	@   v8    96
	@   v9    64    (shared with t5)
	@   v10   32    (shared with t3)
	@ parameters:
	@   [sp, #128]   pointer to P3 (X3:Z3:U3:T3)
	@   [sp, #132]   pointer to P1 (X1:Z1:U1:T1)
	@   [sp, #136]   pointer to P2 (X2:U2)

	@ t1 <- X1*X2  (into v1)
	mov	r0, sp
	bl	gf_mul_inner

	@ t2 = Z1*Z2 = Z1

	@ t3 <- U1*U2  (into v3)
	add	r0, sp, #32
	ldr	r1, [sp, #132]
	adds	r1, #64
	ldr	r2, [sp, #136]
	adds	r2, #32
	bl	gf_mul_inner

	@ t4 = T1*T2 = T1

	@ t5 <- X1*Z2 + X2*Z1 = X1 + X2*Z1  (into v5)
	add	r0, sp, #64
	ldr	r1, [sp, #136]
	ldr	r2, [sp, #132]
	adds	r2, #32
	bl	gf_mul_inner
	add	r0, sp, #64
	add	r1, sp, #64
	ldr	r2, [sp, #132]
	bl	gf_add_inner

	@ t6 <- U1*T2 + U2*T1 = U1 + U2*T1
	@ U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*Z1)
	@ We use X3 as scratch (X1 won't be needed any more)
	@ We also want 2*Z1 into v8 on output
	@ *** v8 <- U2*T1
	add	r0, sp, #96
	ldr	r1, [sp, #132]
	adds	r1, #96
	ldr	r2, [sp, #136]
	adds	r2, #32
	bl	gf_mul_inner
	@ *** U3 <- U1 + v8 = t6 = U1 + U2*T1
	ldr	r0, [sp, #128]
	adds	r0, #64
	ldr	r1, [sp, #132]
	adds	r1, #64
	add	r2, sp, #96
	bl	gf_add_inner
	@ *** v8 <- 2*Z1
	add	r0, sp, #96
	ldr	r1, [sp, #132]
	adds	r1, #32
	bl	gf_mul2_inner
	@ *** X3 <- v1 + v8 = t1 + 2*Z1
	ldr	r0, [sp, #128]
	mov	r1, sp
	add	r2, sp, #96
	bl	gf_add_inner
	@ *** U3 <- U3*X3 = t6*(t1 + 2*Z1)
	ldr	r0, [sp, #128]
	movs	r1, r0
	adds	r0, #64
	movs	r2, r0
	bl	gf_mul_inner
	@ *** U3 <- -U3 = -t6*(t1 + 2*Z1)
	ldr	r0, [sp, #128]
	adds	r0, #64
	movs	r1, r0
	bl	gf_neg_inner

	@ t7 <- t1 + b*t2  (with b = -2 and t2 = Z1)  (into v7)
	@ 2*Z1 is still in v8
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #96
	bl	gf_sub_inner

	@ t8 <- t4*t7  (with t4 = T1)  (into v8)
	add	r0, sp, #96
	mov	r1, sp
	ldr	r2, [sp, #132]
	adds	r2, #96
	bl	gf_mul_inner

	@ t5 + t7  (into v7)
	mov	r0, sp
	add	r1, sp, #64
	mov	r2, sp
	bl	gf_add_inner

	@ t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2)
	@ We actually store -t9 into v9.
	add	r0, sp, #64
	add	r1, sp, #32
	add	r2, sp, #64
	bl	gf_mul_inner
	add	r0, sp, #64
	add	r1, sp, #64
	bl	gf_mul4_inner

	@ t10 <- (t4 + alpha*t3)*(t5 + t7)  (t4 = T1, alpha = 2)  (into v10)
	add	r0, sp, #32
	add	r1, sp, #32
	bl	gf_mul2_inner
	add	r0, sp, #32
	add	r1, sp, #32
	ldr	r2, [sp, #132]
	adds	r2, #96
	bl	gf_add_inner
	add	r0, sp, #32
	add	r1, sp, #32
	mov	r2, sp
	bl	gf_mul_inner

	@ Z3 <- t8 - t9   (v9 contains -t9)
	ldr	r0, [sp, #128]
	adds	r0, #32
	add	r1, sp, #96
	add	r2, sp, #64
	bl	gf_add_inner

	@ T3 <- t8 + t9   (v9 contains -t9)
	ldr	r0, [sp, #128]
	adds	r0, #96
	add	r1, sp, #96
	add	r2, sp, #64
	bl	gf_sub_inner

	@ X3 <- b*(t10 - t8 + beta*t9) = 2*(t8 - t10) - t9  (v9 contains -t9)
	add	r0, sp, #96
	add	r1, sp, #96
	add	r2, sp, #32
	bl	gf_sub_inner
	add	r0, sp, #96
	add	r1, sp, #96
	bl	gf_mul2_inner
	ldr	r0, [sp, #128]
	add	r1, sp, #96
	add	r2, sp, #64
	bl	gf_add_inner

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
@ void CURVE_add_mixed_xu(CURVE_point_xu *P3, const CURVE_point_xu *P1,
@                         const CURVE_point_affine_xu *P2)
@ =======================================================================

	.align	1
	.global	CN(add_mixed_xu)
	.thumb
	.thumb_func
	.type	CN(add_mixed_xu), %function
CN(add_mixed_xu):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	point_add_mixed_xu_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(add_mixed_xu), .-CN(add_mixed_xu)
