@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines support functions for performing point multiplication (called
@ from pmul_arm.c).
@ =======================================================================

@ =======================================================================

@ Window lookup, inner step.
@   r0    current pointer (updated: +32)
@   r14   index (updated: -1)
@ Value is set in the accumulator. r1, r2 and r3 are consumed.
@ Cost: 22
.macro WINDOW_LOOKUP_STEP
	subs	r14, #1
	sbcs	r1, r1
	uadd8	r1, r1, r1
	ldm	r0!, { r1, r2, r3 }
	sel	r4, r1, r4
	sel	r5, r2, r5
	sel	r6, r3, r6
	ldm	r0!, { r1, r2, r3 }
	sel	r7, r1, r7
	sel	r8, r2, r8
	sel	r10, r3, r10
	ldm	r0!, { r1, r2 }
	sel	r11, r1, r11
	sel	r12, r2, r12
.endm

@ Window lookup, for one coordinate.
@   r0    pointer to first table entry (updated: +480)
@   r14   index-1 (update: -8)
@ Value is set in the accumulator (0 if index == 0).
@ r1, r2 and r3 are consumed.
@ Cost: 12 + 8*cost(WINDOW_LOOKUP_STEP) = 188
.macro WINDOW_LOOKUP_COORDINATE
	eors	r4, r4
	eors	r5, r5
	umull	r6, r7, r4, r4
	umull	r8, r10, r4, r4
	umull	r11, r12, r4, r4
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
	adds	r0, #32
	WINDOW_LOOKUP_STEP
.endm

@ =======================================================================
@ void window_lookup_8_affine_inner(CURVE_point_affine *P,
@                                   const CURVE_point_affine *win,
@                                   size_t index)
@
@ ABI: all registers are consumed.
@
@ Cost: 35 + 2*cost(WINDOW_LOOKUP_COORDINATE) = 411
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_lookup_8_affine_inner, %function
window_lookup_8_affine_inner:
	push	{ r0, lr }

	@ We set r14 to index-1. Each lookup decrements it by 1; this
	@ triggers a borrow when going from 0 to -1.
	sub	r14, r2, #1

	@ Lookup first coordinate.
	movs	r0, r1
	WINDOW_LOOKUP_COORDINATE
	ldr	r2, [sp]
	ACC_STORE  r2

	@ Lookup second coordinate. Window pointer and index must be
	@ adjusted to cancel the side effects of the first lookup, and
	@ to start at the second coordinate.
	subw	r0, r0, #448
	add	r14, r14, #8
	WINDOW_LOOKUP_COORDINATE
	ldr	r2, [sp]
	adds	r2, #32
	ACC_STORE  r2

	pop	{ r0, pc }
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
@ Cost: 6 + num * 98 + cost(gf_inv_inner)
@       + (2*num+1) * cost(gf_mul_to_acc_inner)
@       + (num-1) * cost(gf_mul_to_acc_inner_altentry)
@       = 28458 + 615*num
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_inv_multiple_inner, %function
gf_inv_multiple_inner:
	@ Save parameters.
	add	r3, sp, r2, lsl #5
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
	add	r2, sp, #24
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #12]
	adds	r1, #32
	ACC_STORE  r1

	@ m[i] for i = 2..num
	movs	r0, #1
.Lgf_inv_multiple_loop1:
	@ r0 contains i-1.
	@ r1 points to where we stored m[i-1].
	@ m[i-1] is also in the accumulator.
	ldr	r3, [sp, #8]
	cmp	r0, r3
	bcs	.Lgf_inv_multiple_cont1
	str	r0, [sp, #16]
	adds	r1, #16
	add	r2, sp, #24
	add	r2, r2, r0, lsl #5
	bl	gf_mul_to_acc_inner_altentry
	ldrd	r1, r0, [sp, #12]
	adds	r0, #1
	add	r1, r1, r0, lsl #5
	ACC_STORE  r1
	b	.Lgf_inv_multiple_loop1

.Lgf_inv_multiple_cont1:
	@ Invert m[num] (in place); r1 points to m[num].
	movs	r0, r1
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
	ldr	r1, [sp, #12]
	add	r1, r1, r0, lsl #5
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	ldrd	r1, r3, [sp, #12]
	ACC_STORE  r1
	@ Compute 1/m[i-1] = x[i]*(1/m[i]) (into slot of m[i-1])
	@ We have i-1 in r3.
	add	r1, sp, #24
	add	r1, r1, r3, lsl #5
	ldr	r2, [sp, #12]
	add	r2, r2, r3, lsl #5
	adds	r2, #32
	bl	gf_mul_to_acc_inner
	ldrd	r0, r1, [sp, #12]
	add	r0, r0, r1, lsl #5
	ACC_STORE  r0
	@ Copy the computed 1/x[i] into its destination.
	ldrd	r1, r0, [sp, #12]
	add	r2, sp, #24
	add	r2, r2, r0, lsl #5
	ACC_LOAD  r1
	ACC_STORE  r2
	@ Loop. Current counter is already in r0.
	b	.Lgf_inv_multiple_loop2

	@ When we reach this point, 1/x[i] has been computed and written
	@ for i = 2..num. We have 1/m[1] in its stack slot. Since
	@ m[1] = x[0]*x[1], we can invert x[0] and x[1] with one
	@ multiplication each.
.Lgf_inv_multiple_cont2:
	@ 1/x[1] = x[0]*(1/m[1]) (into stack temporary)
	ldr	r1, [sp, #12]
	adds	r1, #32
	ldr	r2, [sp, #4]
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #12]
	ACC_STORE  r1
	@ 1/x[0] = x[1]*(1/m[1])
	adds	r1, #32
	add	r2, sp, #24
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp]
	ACC_STORE  r0
	@ Copy 1/x[1] into its destination.
	add	r0, sp, #24
	ldr	r1, [sp, #12]
	ACC_LOAD  r1
	ACC_STORE  r0

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
@ Cost: 617 + 4*cost(point_double_inner) + 3*cost(point_add_inner)
@       = 15342 (do255e) or 15291 (do255s)
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_fill_8_inner, %function
window_fill_8_inner:
	@ Save return address in one of the free slots.
	str	lr, [sp, #492]

	@ 2*P -> stack (224), then copy X and W to window (64), Z to stack (0)
	@ r1 points to source point P.
	add	r0, sp, #224
	bl	point_double_inner
	ldr	r0, [sp, #484]
	adds	r0, #64
	add	r1, sp, #224
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE  r0
	ACC_LOAD  r1
	ACC_STORE  sp

	@ For i = 2 to 4, compute (2*i-1)*P and (2*i)*P, with X and W
	@ copied to the window, and Z to the stack. The complete (2*i)*P
	@ is left in the stack (224) for the next iteration.

	@ We store 2*i-2 in r3 (value from 2 to 6).
	movs	r3, #2
.Lwindow_fill_8_inner_loop:
	@ Odd index: add P to previous value. Result goes to window
	@ buffer, and Z is copied out of it.
	str	r3, [sp, #480]
	ldrd	r0, r1, [sp, #484]
	add	r0, r0, r3, lsl #6
	add	r2, sp, #224
	bl	point_add_inner
	ldrd	r3, r1, [sp, #480]
	add	r0, sp, r3, lsl #5
	subs	r0, #32
	add	r1, r1, r3, lsl #6
	adds	r1, #64
	ACC_LOAD  r1
	ACC_STORE  r0
	ACC_LOAD  r1
	ACC_STORE  r0

	@ Even index: copy point i*P from window + stack into stack (224),
	@ then double, and copy resulting X and W to window, Z to stack.
	add	r0, sp, #224
	ldr	r1, [sp, #484]
	add	r1, r1, r3, lsl #5
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD  r1
	ACC_STORE_UPDATE  r0
	add	r1, sp, r3, lsl #4
	subs	r1, #32
	ACC_LOAD  r1
	ACC_STORE  r0
	add	r0, sp, #224
	add	r1, sp, #224
	bl	point_double_inner
	ldrd	r3, r0, [sp, #480]
	add	r0, r0, r3, lsl #6
	adds	r0, #64
	add	r1, sp, #224
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE  r0
	add	r0, sp, r3, lsl #5
	ACC_LOAD  r1
	ACC_STORE  r0

	@ Loop control.
	adds	r3, #2
	cmp	r3, #8
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	window_lookup_8_affine_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(window_lookup_8_affine), .-CN(window_lookup_8_affine)
