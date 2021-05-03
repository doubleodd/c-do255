@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines support functions for performing point multiplication (called
@ from pmul_do255s_arm.c).
@ =======================================================================

@ =======================================================================
@ void window_fill_8_affine_inner(CURVE_point_affine *win,
@                                 const CURVE_point *P)
@
@ ABI: all registers are consumed.
@
@ Cost: 567 + cost(window_fill_8_inner)
@       + cost(gf_inv_multiple_inner[num=7])
@       + 8*cost(gf_sqr_acc_inner) + 8*cost(gf_mul_to_acc_inner)
@       + 8*cost(gf_mul_to_acc_inner_altentry)
@       = 52301
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_fill_8_affine_inner, %function
window_fill_8_affine_inner:
	@ Stack layout:
	@    0   buffer for Z coordinates (2*P ... 8*P) (224 bytes)
	@  224   temporary (256 bytes)
	@  480   temporary (4 bytes)
	@  484   pointer to window
	@  488   pointer to P
	@  492   temporary (4 bytes)
	@  496   temporary (4 bytes)
	push	{ r0, r1, r2, r3, lr }
	sub	sp, #484

	@ Step 1: compute window points in Jacobian coordinates.
	bl	window_fill_8_inner

	@ Step 2: We now have the X and W of 2*P to 8*P in the window,
	@ and the Z of 2*P to 8*P in the stack. We invert all the Z
	@ coordinates, including the Z coordinate of the source (P).
	@ Inverse of the Z of P goes to stack (320).
	add	r0, sp, #320
	ldr	r1, [sp, #488]
	adds	r1, #64
	movs	r2, #7
	bl	gf_inv_multiple_inner

	@ Step 3: apply inverted Z coordinates to all points in the
	@ window to convert them to affine coordinates. If the source
	@ point was the neutral, then we get (0,0) for all affine
	@ coordinates, which is correct.

	@ Convert 1*P to affine.
	ldr	r1, [sp, #488]
	adds	r1, #32
	add	r2, sp, #320
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #484]
	adds	r0, #32
	ACC_STORE  r0
	add	r1, sp, #320
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #336
	stm	r1, { r8, r10, r11, r12 }
	ldr	r2, [sp, #488]
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #484]
	ACC_STORE  r0

	@ Convert 2*P to 8*P to affine.
	movs	r3, #0
.Lwindow_fill_8_affine_inner_loop2:
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
	ldrd	r3, r1, [sp, #480]
	adds	r1, #64
	add	r1, r1, r3, lsl #6
	ACC_STORE  r1
	adds	r3, #1
	cmp	r3, #7
	bne	.Lwindow_fill_8_affine_inner_loop2

	@ All done, exit.
	add	sp, #500
	pop	{ pc }
	.size	window_fill_8_affine_inner, .-window_fill_8_affine_inner

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_window_fill_8_affine(CURVE_point_affine *win,
@                                 const CURVE_point *P)
@ =======================================================================

	.align	1
	.global	CN(window_fill_8_affine)
	.thumb
	.thumb_func
	.type	CN(window_fill_8_affine), %function
CN(window_fill_8_affine):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	window_fill_8_affine_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(window_fill_8_affine), .-CN(window_fill_8_affine)
