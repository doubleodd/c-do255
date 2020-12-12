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
@ Cost: 326 + cost(window_fill_8_inner)
@       + cost(gf_inv_multiple_inner[num=7])
@       + 8*cost(gf_sqr_inner) + 16*cost(gf_mul_inner)
@       = 207864
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
	ldr	r0, [sp, #484]
	adds	r0, #32
	ldr	r1, [sp, #488]
	adds	r1, #32
	add	r2, sp, #320
	bl	gf_mul_inner
	add	r0, sp, #320
	add	r1, sp, #320
	bl	gf_sqr_inner
	ldr	r0, [sp, #484]
	ldr	r1, [sp, #488]
	add	r2, sp, #320
	bl	gf_mul_inner

	@ Convert 2*P to 8*P to affine.
	movs	r7, #0
.Lwindow_fill_8_affine_inner_loop2:
	str	r7, [sp, #480]
	ldr	r0, [sp, #484]
	lsls	r3, r7, #6
	adds	r0, r3
	adds	r0, #96
	movs	r1, r0
	lsls	r2, r7, #5
	add	r2, sp, r2
	bl	gf_mul_inner
	ldr	r7, [sp, #480]
	lsls	r0, r7, #5
	add	r0, sp, r0
	movs	r1, r0
	bl	gf_sqr_inner
	ldr	r7, [sp, #480]
	ldr	r0, [sp, #484]
	lsls	r3, r7, #6
	adds	r0, r3
	adds	r0, #64
	movs	r1, r0
	lsls	r2, r7, #5
	add	r2, sp, r2
	bl	gf_mul_inner
	ldr	r7, [sp, #480]
	adds	r7, #1
	cmp	r7, #7
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
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }
	bl	window_fill_8_affine_inner
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }
	.size	CN(window_fill_8_affine), .-CN(window_fill_8_affine)
