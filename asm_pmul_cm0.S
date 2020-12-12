@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines support functions for performing point multiplication (called
@ from pmul_arm.c).
@ =======================================================================

@ =======================================================================

@ Window lookup step. Field element at address r0 is masked with r4
@ (which has value 0 or -1), then added limb-wise to the accumulator
@ (in r5:r6:r7:r8:r10:r11:r12:r14).
@ Input:
@   r0        pointer to field element (updated)
@   r4        mask (0 or -1) (unmodified)
@   r5..r14   accumulator (updated)
@ r1, r2 and r3 are scratch.
@ On output, r0 has been incremented by 32.
@ Cost: 28
.macro WINDOW_LOOKUP_STEP
	ldm	r0!, { r2, r3 }
	ands	r2, r4
	ands	r3, r4
	adds	r5, r2
	adds	r6, r3
	ldm	r0!, { r2, r3 }
	ands	r2, r4
	ands	r3, r4
	adds	r7, r2
	add	r8, r3
	ldm	r0!, { r2, r3 }
	ands	r2, r4
	ands	r3, r4
	add	r10, r2
	add	r11, r3
	ldm	r0!, { r2, r3 }
	ands	r2, r4
	ands	r3, r4
	add	r12, r2
	add	r14, r3
.endm

@ Window lookup for one coordinate.
@ Input:
@   r0        pointer to the coordinate in the first window element (updated)
@   r2        index - 1 (updated)
@ Output:
@   r5..r14   coordinate value (or zero)
@ On output:
@   r0 has been incremented by 480
@   r2 has been decremented by 8
@ r1, r3 and r4 are scratch.
@ Cost: 63 + 8*cost(WINDOW_LOOKUP_STEP) = 287
.macro WINDOW_LOOKUP_COORDINATE
	eors	r5, r5
	eors	r6, r6
	eors	r7, r7
	mov	r8, r5
	mov	r10, r5
	mov	r11, r5
	mov	r12, r5
	mov	r14, r5

	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	pop	{ r2 }
	subs	r2, #1
	push	{ r2 }
	sbcs	r4, r4
	WINDOW_LOOKUP_STEP
	pop	{ r2 }
.endm

@ =======================================================================
@ void window_lookup_8_affine_inner(CURVE_point_affine *P,
@                                   const CURVE_point_affine *win,
@                                   size_t index)
@
@ ABI: all registers are consumed.
@
@ Cost: 47 + 2*cost(WINDOW_LOOKUP_COORDINATE) = 621
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_lookup_8_affine_inner, %function
window_lookup_8_affine_inner:
	push	{ r0, lr }

	@ Subtracting 1 yields no borrow except when going from 0 to -1.

	@ Lookup X coordinate. The accumulator is set to zero. The value
	@ corresponding to the index parameter is added to the accumulator;
	@ if none matches (i.e. index was zero), then the accumulator
	@ remains at zero, which is correct.
	subs	r2, #1
	movs	r0, r1
	WINDOW_LOOKUP_COORDINATE
	ldr	r1, [sp]
	stm	r1!, { r5, r6, r7 }
	mov	r3, r8
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	stm	r1!, { r3, r4, r5, r6, r7 }

	@ Lookup W coordinate.
	subs	r0, #224
	subs	r0, #224
	adds	r2, #8
	WINDOW_LOOKUP_COORDINATE
	pop	{ r1 }
	adds	r1, #32
	stm	r1!, { r5, r6, r7 }
	mov	r3, r8
	mov	r4, r10
	mov	r5, r11
	mov	r6, r12
	mov	r7, r14
	stm	r1!, { r3, r4, r5, r6, r7 }

	pop	{ pc }
	.size	window_lookup_8_affine_inner, .-window_lookup_8_affine_inner

@ =======================================================================
@ void invert_multiple_inner(gf *dx, const gf *sx, size_t num)
@
@ This function inverts num+1 field elements. Most input and output is
@ done on the stack:
@   sp              element x[1] on input, receives 1/x[1]
@   sp+32           element x[2] on input, receives 1/x[2]
@   sp+64           element x[3] on input, receives 1/x[3]
@   ...
@   sp+32*(num-1)   element x[num] on input, receives 1/x[num]
@   r1 (sx)         element x[0] on input
@   r0 (dx)         receives 1/x[0]
@ The function assumes that at least 32*(num+1) bytes are available on
@ the stack beyond x[num].
@
@ num MUST be non-zero. The caller is responsible for making sure that
@ large enough stack buffers have been allocated.
@
@ It is acceptable that dx == sx, or that dx points to one of the
@ temporary stack slots except the first two (i.e. we can have
@ dx == sp + 32*num + 64, but addresses sp to sp + 32*num + 63 are
@ forbidden).
@
@ ABI: all registers are consumed.
@
@ Cost: 9 + num * 78 + cost(gf_inv_inner) + (3*num) * cost(gf_mul_inner)
@       = 54802 + 4470*num
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_inv_multiple_inner, %function
gf_inv_multiple_inner:
	@ Save parameters.
	lsls	r3, r2, #5
	add	r3, sp, r3
	push	{ r0, r1, r2, r3, r4, lr }

	@ Stack layout:
	@    0   dx (pointer to output for 1/x[0])
	@    4   sx (pointer to input for 1/x[0])
	@    8   num
	@   12   pointer to temporary stack buffer
	@   16   temporary slot (loop counter)
	@   20   saved linked register
	@   24   input/output buffer (32 * num)

	@ Step 1:
	@ Compute m[i] for i = 1..num:
	@    m[1] = x[0] * x[1]
	@    m[i] = m[i-1] * x[i] for i = 2..num
	@ Values m[] go to the stack buffer (m[i] goes to offset 32*i)
	@ We also formally define m[0] = x[0] (implicit).

	@ m[1]
	ldr	r0, [sp, #12]
	adds	r0, #32
	add	r2, sp, #24
	bl	gf_mul_inner

	@ m[i] for i = 2..num
	movs	r0, #1
.Lgf_inv_multiple_loop1:
	ldr	r1, [sp, #8]
	cmp	r0, r1
	bcs	.Lgf_inv_multiple_cont1
	str	r0, [sp, #16]
	@ r0 contains i-1
	lsls	r3, r0, #5
	ldr	r1, [sp, #12]
	adds	r1, r3
	add	r2, sp, #24
	adds	r2, r3
	movs	r0, r1
	adds	r0, #32
	bl	gf_mul_inner
	ldr	r0, [sp, #16]
	adds	r0, #1
	b	.Lgf_inv_multiple_loop1

	@ Invert m[num].
.Lgf_inv_multiple_cont1:
	ldr	r3, [sp, #8]
	lsls	r3, r3, #5
	ldr	r0, [sp, #12]
	adds	r0, r3
	movs	r1, r0
	bl	gf_inv_inner

	@ For i = num down to 1, we have:
	@   1/m[i-1] = x[i]*(1/m[i])
	@   1/x[i] = m[i-1]*(1/m[i])
	@ and 1/x[0] = 1/m[0].
	@ Loop below handles i = num down to 2.
	ldr	r0, [sp, #8]
.Lgf_inv_multiple_loop2:
	subs	r0, #1
	beq	.Lgf_inv_multiple_cont2
	str	r0, [sp, #16]
	@ r0 contains i-1 (from num-1 down to 1)
	@ Compute 1/x[i] = m[i-1]*(1/m[i]) (into temporary slot)
	lsls	r3, r0, #5
	ldr	r0, [sp, #12]
	movs	r1, r0
	adds	r1, r3
	movs	r2, r1
	adds	r2, #32
	bl	gf_mul_inner
	@ Compute 1/m[i-1] = x[i]*(1/m[i]) (into slot of m[i-1])
	ldr	r0, [sp, #16]
	lsls	r3, r0, #5
	ldr	r0, [sp, #12]
	adds	r0, r3
	add	r1, sp, #24
	adds	r1, r3
	movs	r2, r0
	adds	r2, #32
	bl	gf_mul_inner
	@ Copy the computed 1/x[i] into its destination.
	ldr	r0, [sp, #16]
	lsls	r7, r0, #5
	add	r7, sp, r7
	adds	r7, #24
	ldr	r1, [sp, #12]
	ldm	r1!, { r3, r4, r5, r6 }
	stm	r7!, { r3, r4, r5, r6 }
	ldm	r1!, { r3, r4, r5, r6 }
	stm	r7!, { r3, r4, r5, r6 }
	@ Loop. Current counter is already in r0.
	b	.Lgf_inv_multiple_loop2

	@ When we reach this point, 1/x[i] has been computed and written
	@ for i = 2..num. We have 1/m[1] in its stack slot. Since
	@ m[1] = x[0]*x[1], we can invert x[0] and x[1] with one
	@ multiplication each.
.Lgf_inv_multiple_cont2:
	@ 1/x[1] = x[0]*(1/m[1]) (into stack temporary)
	ldr	r0, [sp, #12]
	movs	r1, r0
	adds	r1, #32
	ldr	r2, [sp, #4]
	bl	gf_mul_inner
	@ 1/x[0] = x[1]*(1/m[1])
	ldr	r0, [sp]
	ldr	r1, [sp, #12]
	adds	r1, #32
	add	r2, sp, #24
	bl	gf_mul_inner
	@ Copy 1/x[1] into its destination.
	add	r0, sp, #24
	ldr	r1, [sp, #12]
	ldm	r1!, { r4, r5, r6, r7 }
	stm	r0!, { r4, r5, r6, r7 }
	ldm	r1!, { r4, r5, r6, r7 }
	stm	r0!, { r4, r5, r6, r7 }

	@ Exit.
	pop	{ r0, r1, r2, r3, r4, pc }
	.size	gf_inv_multiple_inner, .-gf_inv_multiple_inner

@ =======================================================================
@ void window_fill_8_inner(CURVE_point_affine *win,
@                          const CURVE_point *P)
@
@ This routine computes 8 window points in Jacobian coordinates.
@ The Z coordinate for points 2*P to 8*P goes into a stack buffer.
@ This routine assumes that the stack buffer has already been allocated,
@ and fills it (and the window) as follows:
@           window (X, W)   stack (Z)
@   1*P        0     32         -   (warning: not actually copied)
@   2*P       64     96         0
@   3*P      128    160        32
@   4*P      192    224        64
@   5*P      256    288        96
@   6*P      320    352       128
@   7*P      384    416       160
@   8*P      448    480       192
@
@ WARNING: coordinates of 1*P are NOT copied. The caller should use the
@ source point instead.
@
@ Parameters are supposed to be on the stack AND in registers:
@   offset   register   value
@     480               free slot (used by this routine)
@     484               pointer to window (unmodified)
@     488       r1      pointer to P (unmodified)
@     492               free slot (used by this routine)
@
@ ABI: all registers are consumed.
@
@ Cost: 637 + 4*cost(point_double_inner) + 3*cost(point_add_inner)
@       = 93002 (do255e) or 90046 (do255s)
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_fill_8_inner, %function
window_fill_8_inner:
	@ Save return address in one of the free slots.
	mov	r0, lr
	str	r0, [sp, #492]

	@ 2*P -> stack (224), then copy X and W to window (64), Z to stack (0)
	@ r1 points to source point P.
	add	r0, sp, #224
	bl	point_double_inner
	ldr	r0, [sp, #484]
	adds	r0, #64
	add	r1, sp, #224
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5 }
	mov	r0, sp
	stm	r0!, { r6, r7 }
	ldm	r1!, { r2, r3, r4, r5, r6, r7 }
	stm	r0!, { r2, r3, r4, r5, r6, r7 }

	@ For i = 2 to 4, compute (2*i-1)*P and (2*i)*P, with X and W
	@ copied to the window, and Z to the stack. The complete (2*i)*P
	@ is left in the stack (224) for the next iteration.

	@ We store 2*i-2 in r7 (value from 2 to 6).
	movs	r7, #2
.Lwindow_fill_8_inner_loop:
	@ Odd index: add P to previous value. Result goes to window
	@ buffer, and Z is copied out of it.
	str	r7, [sp, #480]
	ldr	r0, [sp, #484]
	lsls	r3, r7, #6
	adds	r0, r3
	ldr	r1, [sp, #488]
	add	r2, sp, #224
	bl	point_add_inner
	ldr	r7, [sp, #480]
	lsls	r0, r7, #5
	subs	r0, #32
	add	r0, sp, r0
	lsls	r3, r7, #6
	ldr	r1, [sp, #484]
	adds	r1, r3
	adds	r1, #64
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }

	@ Even index: copy point i*P from window + stack into stack (224),
	@ then double, and copy resulting X and W to window, Z to stack.
	add	r0, sp, #224
	lsls	r3, r7, #5
	ldr	r1, [sp, #484]
	adds	r1, r3
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	lsls	r1, r7, #4
	add	r1, sp, r1
	subs	r1, #32
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	add	r0, sp, #224
	add	r1, sp, #224
	bl	point_double_inner
	ldr	r7, [sp, #480]
	lsls	r3, r7, #6
	ldr	r0, [sp, #484]
	adds	r0, r3
	adds	r0, #64
	add	r1, sp, #224
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	lsls	r0, r7, #5
	add	r0, sp, r0
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }

	@ Loop control.
	adds	r7, #2
	cmp	r7, #8
	bne	.Lwindow_fill_8_inner_loop

	@ Exit function.
	ldr	r0, [sp, #492]
	bx	r0
	.size	window_fill_8_inner, .-window_fill_8_inner

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_window_lookup_8_affine(CURVE_point_affine *P,
@                                   const CURVE_point_affine *win,
@                                   size_t index)
@ =======================================================================

	.align	1
	.global	CN(window_lookup_8_affine)
	.thumb
	.thumb_func
	.type	CN(window_lookup_8_affine), %function
CN(window_lookup_8_affine):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	window_lookup_8_affine_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(window_lookup_8_affine), .-CN(window_lookup_8_affine)
