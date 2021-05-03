@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines support functions for performing signature verification in
@ curve do255s (called from pvrfy_do255s_arm.c).
@ =======================================================================

@ =======================================================================
@ void window_fill_8odd_x2_affine_inner(CURVE_point_affine *win1,
@                                       CURVE_point_affine *win2,
@                                       const CURVE_point *P1,
@                                       const CURVE_point *P2)
@
@ ABI: all registers are consumed.
@
@ Cost: 882 + 2*cost(point_double_inner) + 6*cost(point_add_inner)
@       + cost(gf_inv_multiple_inner[num=7])
@       + 8*cost(gf_sqr_acc_inner) + 8*cost(gf_mul_to_acc_inner)
@       + 8*cost(gf_mul_to_acc_inner_altentry)
@       = 59549
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_fill_8odd_x2_affine_inner, %function
window_fill_8odd_x2_affine_inner:
	@ Stack layout:
	@    0   first buffer for Z coordinates (3*P1, 5*P1, 7*P1)
	@   96   second buffer for Z coordinates (1*P2, 3*P2, 5*P2, 7*P2)
	@  224   temporary (256 bytes)
	@  480   temporary (4 bytes)
	@  484   pointer to first window (win1)
	@  488   pointer to second window (win2)
	@  492   pointer to P1
	@  496   pointer to P2
	push	{ r0, r1, r2, r3, lr }
	sub	sp, #484

	@ Step 1: compute window points in Jacobian coordinates. We store
	@ the Z coordinates in a stack buffer.
	@ Destination of (X:W:Z) for each point:
	@            win1 (X, W)    stack (Z)
	@   1*P1       0     32         -
	@   3*P1      64     96         0
	@   5*P1     128    160        32
	@   7*P1     192    224        64

	@ 2*P1 -> stack (0)
	mov	r0, sp
	movs	r1, r2
	bl	point_double_inner

	@ 3*P1 -> win1 (64)
	ldr	r0, [sp, #484]
	adds	r0, #64
	mov	r1, sp
	ldr	r2, [sp, #492]
	bl	point_add_inner

	@ 5*P1 -> stack (96)
	add	r0, sp, #96
	ldr	r1, [sp, #484]
	adds	r1, #64
	mov	r2, sp
	bl	point_add_inner

	@ 7*P -> stack (0)
	mov	r0, sp
	mov	r1, sp
	add	r2, sp, #96
	bl	point_add_inner

	@ At this point, we have:
	@   window (win1):
	@      0   -
	@     32   -
	@     64   X(3*P1)   (correct)
	@     96   W(3*P1)   (correct)
	@    128   Z(3*P1)
	@    160   -
	@    192   -
	@    224   -
	@
	@   stack:
	@      0   X(7*P1)
	@     32   W(7*P1)
	@     64   Z(7*P1)   (correct)
	@     96   X(5*P1)
	@    128   W(5*P1)
	@    160   Z(5*P1)

	@ Copy X(7*P1) and W(7*P1) to win1 (192)
	ldr	r0, [sp, #484]
	adds	r0, #192
	mov	r1, sp
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy Z(3*P1) to stack (0)
	movs	r1, r0
	subs	r1, #128
	mov	r0, sp
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy Z(5*P1) to stack (32)
	add	r1, sp, #160
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy X(5*P1) and W(5*P1) to win1 (128)
	ldr	r0, [sp, #484]
	adds	r0, #128
	add	r1, sp, #96
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ We now similarly compute 3*P2, 5*P2 and 7*P2, with Z coordinates
	@ copied to the stack, starting at offset 128. We reserve a stack
	@ slot (96) for a copy of the Z coordinate of P2.
	@            win2 (X, W)    stack (Z)
	@   1*P2       0     32        96
	@   3*P2      64     96       128
	@   5*P2     128    160       160
	@   7*P2     192    224       192

	@ 2*P2 -> stack (128)
	add	r0, sp, #128
	ldr	r1, [sp, #496]
	bl	point_double_inner

	@ 3*P2 -> win2 (64)
	ldr	r0, [sp, #488]
	adds	r0, #64
	add	r1, sp, #128
	ldr	r2, [sp, #496]
	bl	point_add_inner

	@ 5*P2 -> stack (224)
	add	r0, sp, #224
	ldr	r1, [sp, #488]
	adds	r1, #64
	add	r2, sp, #128
	bl	point_add_inner

	@ 7*P2 -> stack (128)
	add	r0, sp, #128
	add	r1, sp, #128
	add	r2, sp, #224
	bl	point_add_inner

	@ At this point, we have:
	@   window (win2):
	@      0   -
	@     32   -
	@     64   X(3*P2)   (correct)
	@     96   W(3*P2)   (correct)
	@    128   Z(3*P2)
	@    160   -
	@    192   -
	@    224   -
	@
	@   stack:
	@    128   X(7*P2)
	@    160   W(7*P2)
	@    192   Z(7*P2)   (correct)
	@    224   X(5*P2)
	@    256   W(5*P2)
	@    288   Z(5*P2)

	@ Copy X(7*P2) and W(7*P2) to win2 (192)
	ldr	r0, [sp, #488]
	adds	r0, #192
	add	r1, sp, #128
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy Z(3*P2) to stack (128)
	movs	r1, r0
	subs	r1, #128
	add	r0, sp, #128
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy Z(5*P2) to stack (160)
	add	r1, sp, #288
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy X(5*P2) and W(5*P2) to win2 (128)
	ldr	r0, [sp, #488]
	adds	r0, #128
	add	r1, sp, #224
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Copy Z(P2) to stack (96)
	add	r0, sp, #96
	ldr	r1, [sp, #496]
	adds	r1, #64
	ACC_LOAD_UPDATE  r1
	ACC_STORE_UPDATE  r0

	@ Step 2: We now have the X and W of 3*P1, 5*P1 and 7*P1 in the
	@ first window (win1) and the X and W of 3*P2, 5*P2 and 7*P2 in
	@ the second window (win2). The stack contains the Z coordinates
	@ of 3*P1, 5*P1, 7*P1, 1*P2, 3*P2, 5*P2 and 7*P2, in that order.
	@ We invert all the Z coordinates, including the Z coordinate of
	@ the first source point (P1).
	@ The inverse of the Z of P1 goes to the stack (288).
	add	r0, sp, #288
	ldr	r1, [sp, #492]
	adds	r1, #64
	movs	r2, #7
	bl	gf_inv_multiple_inner

	@ Step 3: apply inverted Z coordinates to all points in the
	@ window to convert them to affine coordinates. If a source
	@ point was the neutral, then we get (0,0) for all corresponding
	@ affine coordinates, which is correct.

	@ Convert 1*P1 to affine.
	ldr	r1, [sp, #492]
	adds	r1, #32
	add	r2, sp, #288
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #484]
	adds	r0, #32
	ACC_STORE  r0
	add	r1, sp, #288
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #304
	stm	r1, { r8, r10, r11, r12 }
	ldr	r2, [sp, #492]
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #484]
	ACC_STORE  r0

	@ Convert 3*P1 to 7*P1 to affine.
	movs	r3, #0
.Lwindow_fill_8odd_x2_affine_inner_loop1:
	str	r3, [sp, #480]
	ldr	r1, [sp, #484]
	adds	r1, #96
	add	r1, r1, r3, lsl #6
	add	r2, sp, r3, lsl #5
	bl	gf_mul_to_acc_inner
	ldrd	r3, r0, [sp, #480]
	adds	r0, #96
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0
	add	r1, sp, r3, lsl #5
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	ldrd	r3, r2, [sp, #480]
	add	r1, sp, r3, lsl #5
	adds	r1, #16
	stm	r1, { r8, r10, r11, r12 }
	adds	r2, #64
	add	r2, r2, r3, lsl #6
	bl	gf_mul_to_acc_inner_altentry
	ldrd	r3, r0, [sp, #480]
	adds	r0, #64
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0
	adds	r3, #1
	cmp	r3, #3
	bne	.Lwindow_fill_8odd_x2_affine_inner_loop1

	@ Convert 1*P2 to affine.
	ldr	r1, [sp, #496]
	adds	r1, #32
	add	r2, sp, #96
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #488]
	adds	r0, #32
	ACC_STORE  r0
	add	r1, sp, #96
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #112
	stm	r1, { r8, r10, r11, r12 }
	ldr	r2, [sp, #496]
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #488]
	ACC_STORE  r0

	@ Convert 3*P2 to 7*P2 to affine.
	movs	r3, #0
.Lwindow_fill_8odd_x2_affine_inner_loop2:
	str	r3, [sp, #480]
	ldr	r1, [sp, #488]
	adds	r1, #96
	add	r1, r1, r3, lsl #6
	add	r2, sp, r3, lsl #5
	adds	r2, #128
	bl	gf_mul_to_acc_inner
	ldr	r3, [sp, #480]
	ldr	r0, [sp, #488]
	adds	r0, #96
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0
	add	r1, sp, r3, lsl #5
	adds	r1, #128
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	ldr	r3, [sp, #480]
	ldr	r2, [sp, #488]
	add	r1, sp, r3, lsl #5
	adds	r1, #144
	stm	r1, { r8, r10, r11, r12 }
	adds	r2, #64
	add	r2, r2, r3, lsl #6
	bl	gf_mul_to_acc_inner_altentry
	ldr	r3, [sp, #480]
	ldr	r0, [sp, #488]
	adds	r0, #64
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0
	adds	r3, #1
	cmp	r3, #3
	bne	.Lwindow_fill_8odd_x2_affine_inner_loop2

	@ All done, exit.
	add	sp, #500
	pop	{ pc }
	.size	window_fill_8odd_x2_affine_inner, .-window_fill_8odd_x2_affine_inner

@ =======================================================================
@ void neg_i128(i128 *x);
@
@ ABI: This function conforms to the public ABI (registers r4 to r11 are
@ preserved).
@
@ Cost: 24
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.global	CN(neg_i128)
	.type	CN(neg_i128), %function
CN(neg_i128):
	push	{ r4, r5 }
	ldm	r0!, { r1, r2, r3, r4 }
	eors	r5, r5
	subs	r1, r5, r1
	sbcs	r2, r5, r2
	sbcs	r3, r5, r3
	sbcs	r4, r5, r4
	subs	r0, #16
	stm	r0!, { r1, r2, r3, r4 }
	pop	{ r4, r5 }
	bx	lr
	.size	CN(neg_i128), .-CN(neg_i128)

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_window_fill_8odd_x2_affine(CURVE_point_affine *win1,
@                                       CURVE_point_affine *win2,
@                                       const CURVE_point *P1,
@                                       const CURVE_point *P2)
@ =======================================================================

	.align	1
	.global	CN(window_fill_8odd_x2_affine)
	.thumb
	.thumb_func
	.type	CN(window_fill_8odd_x2_affine), %function
CN(window_fill_8odd_x2_affine):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	window_fill_8odd_x2_affine_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(window_fill_8odd_x2_affine), .-CN(window_fill_8odd_x2_affine)
