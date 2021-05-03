@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines support functions for performing point multiplication (called
@ from pmul_arm.c).
@ =======================================================================

@ =======================================================================
@ uint32_t recode4_small_inner(uint8_t *sd, const i128 *s)
@
@ ABI: all registers consumed except r12. Returned value in r0.
@
@ Cost: 459
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	recode4_small_inner, %function
recode4_small_inner:
	@ Compute x = abs(s); remember sign (in r8).
	ldm	r1!, { r4, r5, r6, r7 }
	asrs	r3, r7, #31
	eors	r4, r3
	eors	r5, r3
	eors	r6, r3
	eors	r7, r3
	subs	r4, r3
	sbcs	r5, r3
	sbcs	r6, r3
	sbcs	r7, r3
	lsr	r8, r3, #31

	@ Carry is in r3; it is 0 or 1. Initial value is zero.
	eors	r3, r3

	@ Outer loop counter in r10.
	movw	r10, #4

.Lrecode4_small_inner_loop1:
	@ Inner loop counter in r11.
	movw	r11, #8

.Lrecode4_small_inner_loop2:
	@ Get next four bits, and add carry.
	ubfx	r2, r4, #0, #4
	lsrs	r4, r4, #4
	adds	r2, r3
	@ If result is greater than 8, subtract it from 16.
	rsbs	r3, r2, #8
	lsls	r1, r3, #1
	lsrs	r3, r3, #31
	umlal	r2, r1, r1, r3
	@ Set top bit of the output byte if there was a carry, and
	@ write the byte.
	orr	r2, r2, r3, lsl #7
	strb	r2, [r0], #1

	@ Inner loop control.
	subs	r11, #1
	bne	.Lrecode4_small_inner_loop2

	@ Outer loop control. Arrange for next source word to be in r4.
	movs	r4, r5
	movs	r5, r6
	movs	r6, r7
	subs	r10, #1
	bne	.Lrecode4_small_inner_loop1

	@ Return original sign of s.
	mov	r0, r8
	bx	lr
	.size	recode4_small_inner, .-recode4_small_inner

@ =======================================================================
@ void mul_divr_rounded_inner(i128 *y, const i256 *k, const i128 *e)
@
@ ABI: all registers consumed.
@
@ Cost: 150 + cost(CURVE_mul256x128) + cost(CURVE_mul128x128) = 314
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	mul_divr_rounded_inner, %function
mul_divr_rounded_inner:
	push	{ r0, lr }
	sub	sp, #64

	@ Stack layout:
	@   0   z   48 bytes
	@  32   t   32 bytes (overlaps top part of z)
	@  64   y   destination pointer

	@ z <- k*e
	mov	r0, sp
	bl	CN(mul256x128)

	@ z <- z + (r-1)/2
	ldm	sp, { r0, r1, r2, r3 }
	adr	r4, const_mul_divr_rounded_hr_lo
	ldm	r4, { r4, r5, r6, r7 }
	adds	r0, r4
	adcs	r1, r5
	adcs	r2, r6
	adcs	r3, r7
	mov	r4, sp
	stm	r4!, { r0, r1, r2, r3 }
	ldm	r4!, { r0, r1, r2, r3 }
	adr	r4, const_mul_divr_rounded_hr_hi
	ldm	r4, { r4, r5, r6, r7 }
	adcs	r0, r4
	adcs	r1, r5
	adcs	r2, r6
	adcs	r3, r7
	add	r4, sp, #16
	stm	r4!, { r0, r1, r2, r3 }
	ldm	r4!, { r0, r1, r2, r3 }
	eors	r7, r7
	adcs	r0, r7
	adcs	r1, r7
	adcs	r2, r7
	adcs	r3, r7
	add	r4, sp, #32
	stm	r4!, { r0, r1, r2, r3 }

	@ y <- floor(z / 2^254) + 1
	add	r1, sp, #28
	ldm	r1, { r1, r2, r3, r4, r5 }
	eors	r6, r6
	movs	r7, #3
	umaal	r1, r6, r7, r1
	adds	r6, #1
	umaal	r2, r6, r7, r2
	umaal	r3, r6, r7, r3
	umaal	r4, r6, r7, r4
	umaal	r5, r6, r7, r5
	ldr	r0, [sp, #64]
	stm	r0!, { r2, r3, r4, r5 }

	@ Compute y*r0 into stack (32). This implicitly overwrites the
	@ top part of z.
	add	r0, sp, #32
	ldr	r1, [sp, #64]
	adr	r2, const_mul_divr_rounded_r0
	bl	CN(mul128x128)

	@ Add y*r0 (currently in stack (32)) to z mod 2^254. We are
	@ only interested in the top limb of the result.
	mov	r0, sp
	add	r1, sp, #32
	ldm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r6, r7, r8, r10 }
	adds	r2, r6
	adcs	r3, r7
	adcs	r4, r8
	adcs	r5, r10
	ldm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r6, r7, r8, r10 }
	bic	r5, r5, #0xC0000000
	adcs	r2, r6
	adcs	r3, r7
	adcs	r4, r8
	adcs	r5, r10

	@ The high limb (in r5) has value less than 2^31. If it is
	@ lower than 2^30, then y must be decremented.
	ldr	r0, [sp, #64]
	ldm	r0!, { r3, r4, r6, r7 }
	lsrs	r5, r5, #30
	subs	r5, #1
	adds	r3, r5
	adcs	r4, r5
	adcs	r6, r5
	adcs	r7, r5
	subs	r0, #16
	stm	r0!, { r3, r4, r6, r7 }

	add	sp, #68
	pop	{ pc }
	.align	2
const_mul_divr_rounded_hr_lo:
	.long	0x3A6C2292, 0x8FA96457, 0xAA03C629, 0xCE864987
const_mul_divr_rounded_hr_hi:
	.long	0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x1FFFFFFF
const_mul_divr_rounded_r0:
	.long	0x8B27BADB, 0xE0AD3751, 0xABF873AC, 0x62F36CF0
	.size	mul_divr_rounded_inner, .-mul_divr_rounded_inner

@ =======================================================================
@ void split_inner()
@
@ WARNING: This routine has a special ABI:
@  - It expects a 68-byte buffer on the stack.
@  - The source value k (32 bytes) should already be loaded in the stack
@    buffer (starting at [sp]). That value MUST be fully reduced modulo r.
@  - The computed values are produced in the stack buffer at addresses
@    sp (for k0) and sp+16 (for k1); they replace the source scalar k.
@  - All registers are consumed.
@
@ ABI: all registers consumed.
@
@ Cost: 88 + 2*cost(mul_divr_rounded_inner) + 4*cost(mul128x128trunc_inner)
@       = 808
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	split_inner, %function
split_inner:
	str	lr, [sp, #64]

	@ Stack layout
	@   0   k    32 bytes
	@  32   c    16 bytes
	@  48   d    16 bytes
	@  64   lr   4 bytes (saved return address)

	@ c = round(k*v / r)
	@ d = round(k*u / r)
	add	r0, sp, #32
	mov	r1, sp
	adr	r2, const_split_v
	bl	mul_divr_rounded_inner
	add	r0, sp, #48
	mov	r1, sp
	adr	r2, const_split_u
	bl	mul_divr_rounded_inner

	@ k0 = k - d*u - c*v  (computed modulo 2^128)
	add	r0, sp, #48
	adr	r1, const_split_u
	bl	mul128x128trunc_inner
	pop	{ r4, r5, r6, r7 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	push	{ r4, r5, r6, r7 }
	add	r0, sp, #32
	adr	r1, const_split_v
	bl	mul128x128trunc_inner
	pop	{ r4, r5, r6, r7 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	push	{ r4, r5, r6, r7 }

	@ k1 = d*v - c*u  (computed modulo 2^128)

	add	r0, sp, #48
	adr	r1, const_split_v
	bl	mul128x128trunc_inner
	add	r4, sp, #16
	stm	r4!, { r0, r1, r2, r3 }
	add	r0, sp, #32
	adr	r1, const_split_u
	bl	mul128x128trunc_inner
	add	r4, sp, #16
	ldm	r4, { r4, r5, r6, r7 }
	subs	r4, r0
	sbcs	r5, r1
	sbcs	r6, r2
	sbcs	r7, r3
	add	r0, sp, #16
	stm	r0!, { r4, r5, r6, r7 }

	ldr	r0, [sp, #64]
	bx	r0
	.align	2
const_split_u:
	.long	0xC93F6111, 0x2ACCF9DE, 0x53C2C6E6, 0x1A509F7A
const_split_v:
	.long	0x5466F77E, 0x0B7A3130, 0xFFBB3A93, 0x7D440C6A
	.size	split_inner, .-split_inner

@ =======================================================================
@ void window_fill_8_to_xu_affine_inner(CURVE_point_affine_xu *win,
@                                       const CURVE_point *P)
@
@ This function fills the window with points 1*P, 2*P... 8*P in
@ affine (x,u) coordinates. The source point is in Jacobian (x,w)
@ coordinates.
@
@ ABI: all registers are consumed.
@
@ Cost: 1309 + cost(window_fill_8_inner) + cost(gf_iszero_acc_inner)
@       + 8*cost(gf_sqr_acc_inner) + 31*cost(gf_mul_to_acc_inner)
@       + 9*cost(gf_mul_to_acc_inner_altentry)
@       + cost(gf_inv_multiple_inner[num=7])
@       = 57296
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	window_fill_8_to_xu_affine_inner, %function
window_fill_8_to_xu_affine_inner:
	@ Stack layout:
	@    0   buffer for Z coordinates (2*P ... 8*P) (224 bytes)
	@  224   temporary (256 bytes)
	@  480   temporary (4 bytes)
	@  484   pointer to window
	@  488   pointer to P
	@  492   temporary (4 bytes)
	@  496   temporary (4 bytes)
	@ This is the layout expected by window_fill_8_inner, but we
	@ need an extra 32-byte buffer, which we will use later on.
	@ That buffer is at offset 492.
	push	{ lr }
	sub	sp, #40
	push	{ r0, r1 }
	sub	sp, #484

	@ Step 1: compute window points in Jacobian coordinates.
	bl	window_fill_8_inner

	@ Step 2: convert points to projective (x,u) coordinates
	@   (X:W:Z) -> (X':Z':U':T') with Z' = T'
	@   X' = X*W
	@   U' = Z^3
	@   Z' = W*Z^2
	@ For the source point P, we write Z' in the extra buffer
	@ (at offset 496). Another buffer (tmp) is used at offset 256
	@ for scratch values.
	@
	@ We must take care to handle the neutral point properly: if
	@ the source point is N, then the formulas above set all
	@ three coordinates to 0, and we need to set Z' to a non-zero
	@ value. If P = N, then all points are set to N. If P != N,
	@ then none of 2*P ... 8*P will be the neutral either.

	@ Store P == N into stack slot at offset 492.
	ldr	r0, [sp, #488]
	ACC_LOAD  r0
	bl	gf_iszero_acc_inner
	str	r0, [sp, #492]

	@ Point P itself.
	@ *** X' <- X*W
	ldr	r1, [sp, #488]
	add	r2, r1, #32
	adds	r1, #16
	bl	gf_mul_to_acc_inner_altentry
	ldr	r0, [sp, #484]
	ACC_STORE  r0
	@ *** tmp <- Z^2
	ldr	r1, [sp, #488]
	adds	r1, #64
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #256
	ACC_STORE  r1
	@ *** Z' <- W*tmp = W*Z^2
	adds	r1, #16
	ldr	r2, [sp, #488]
	adds	r2, #32
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #496
	ACC_STORE  r0
	@ *** U' <- tmp*Z = Z^3
	add	r1, sp, #256
	ldr	r2, [sp, #488]
	adds	r2, #64
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #484]
	adds	r0, #32
	ACC_STORE  r0
	@ Handle P == N situation.
	ldrd	r0, r1, [sp, #492]
	orrs	r0, r1
	str	r0, [sp, #496]

	@ Points 2*P ... 8*P
	movs	r7, #0
.Lwindow_fill_8_to_xu_affine_inner_loop1:
	str	r7, [sp, #480]

	@ We convert point (r7+2)*P:
	@   X   win + 64*r7 + 64  (also X')
	@   W   win + 64*r7 + 96  (also U')
	@   Z   sp + 32*r7        (also Z')
	@ We now need two temporary buffers (t1 and t2) at offsets 256 and 288.

	@ *** X' <- X*W  (overwrites X)
	ldr	r1, [sp, #484]
	adds	r1, #64
	add	r1, r1, r7, lsl #6
	add	r2, r1, #32
	bl	gf_mul_to_acc_inner
	ldrd	r3, r0, [sp, #480]
	adds	r0, #64
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0

	@ *** t1 <- Z^2
	add	r1, sp, r3, lsl #5
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	add	r1, sp, #256
	ACC_STORE  r1

	@ *** t2 <- W*t1 = W*Z^2
	adds	r1, #16
	ldr	r2, [sp, #484]
	adds	r2, #96
	ldr	r3, [sp, #480]
	add	r2, r2, r3, lsl #6
	bl	gf_mul_to_acc_inner_altentry
	add	r0, sp, #288
	ACC_STORE  r0

	@ *** U' <- Z*t1 = Z^3  (overwrites W)
	ldr	r3, [sp, #480]
	add	r1, sp, r3, lsl #5
	add	r2, sp, #256
	bl	gf_mul_to_acc_inner
	ldrd	r3, r0, [sp, #480]
	adds	r0, #96
	add	r0, r0, r3, lsl #6
	ACC_STORE  r0

	@ *** Z' <- t2  (overwrites Z)
	@ Also handle P == N situation.
	add	r0, sp, #288
	ACC_LOAD  r0
	ldr	r3, [sp, #492]
	orrs	r4, r3
	ldr	r3, [sp, #480]
	add	r0, sp, r3, lsl #5
	ACC_STORE  r0

	@ Loop control.
	ldr	r7, [sp, #480]
	adds	r7, #1
	cmp	r7, #7
	bne	.Lwindow_fill_8_to_xu_affine_inner_loop1

	@ Step 3: We now have the X and U of 1*P to 8*P in the window,
	@ and the Z=T of 1*P to 8*P in the stack (Z of 1*P is in the
	@ extra buffer at offset 496). We invert all these Z coordinates.
	@ Inverse of the Z of P goes to stack (320).
	add	r0, sp, #320
	add	r1, sp, #496
	movs	r2, #7
	bl	gf_inv_multiple_inner

	@ Step 4: apply inverted Z coordinates to all points in the
	@ window to convert them to affine coordinates.

	@ Convert 1*P to affine.
	ldr	r1, [sp, #484]
	add	r2, sp, #320
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #484]
	ACC_STORE_UPDATE  r1
	add	r2, sp, #320
	bl	gf_mul_to_acc_inner
	ldr	r1, [sp, #484]
	adds	r1, #32
	ACC_STORE  r1

	@ Convert 2*P to 8*P to affine.
	movs	r3, #0
.Lwindow_fill_8_to_xu_affine_inner_loop2:
	str	r3, [sp, #480]
	ldr	r1, [sp, #484]
	adds	r1, #64
	add	r1, r1, r3, lsl #6
	add	r2, sp, r3, lsl #5
	bl	gf_mul_to_acc_inner
	ldrd	r3, r1, [sp, #480]
	adds	r1, #64
	add	r1, r1, r3, lsl #6
	ACC_STORE_UPDATE  r1
	add	r2, sp, r3, lsl #5
	bl	gf_mul_to_acc_inner
	ldrd	r3, r1, [sp, #480]
	adds	r1, #96
	add	r1, r1, r3, lsl #6
	ACC_STORE  r1
	adds	r3, #1
	cmp	r3, #7
	bne	.Lwindow_fill_8_to_xu_affine_inner_loop2

	@ All done, exit.
	add	sp, #492
	add	sp, #40
	pop	{ pc }
	.size	window_fill_8_to_xu_affine_inner, .-window_fill_8_to_xu_affine_inner

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_split_scalar(i128 *k0, i128 *k1, const void *k)
@ =======================================================================

	.align	1
	.global	CN(split_scalar)
	.thumb
	.thumb_func
	.type	CN(split_scalar), %function
CN(split_scalar):
	push	{ r0, r1, r4, r5, r6, r7, r8, r10, r11, lr }
	sub	sp, #72

	@ Decode source scalar; scalar_decode_inner() also reduces
	@ it modulo r.
	mov	r0, sp
	movs	r1, r2
	bl	scalar_decode_inner

	@ Do the split into the stack buffer.
	bl	split_inner

	@ Copy k0 and k1 to output.
	ldrd	r0, r1, [sp, #72]
	ACC_LOAD  sp
	stm	r0!, { r4, r5, r6, r7 }
	stm	r1!, { r8, r10, r11, r12 }

	add	sp, #80
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(split_scalar), .-CN(split_scalar)

@ =======================================================================
@ uint32_t CURVE_split_recode4_scalar(uint8_t *sd0, uint8_t *sd1,
@                                     const void *k)
@ =======================================================================

	.align	1
	.global	CN(split_recode4_scalar)
	.thumb
	.thumb_func
	.type	CN(split_recode4_scalar), %function
CN(split_recode4_scalar):
	push	{ r0, r1, r4, r5, r6, r7, r8, r10, r11, lr }
	sub	sp, #72

	@ Decode source scalar, with reduction modulo r. The value
	@ is now at most 254 bits.
	mov	r0, sp
	movs	r1, r2
	bl	scalar_decode_inner

	@ Do the split into the stack buffer.
	bl	split_inner

	@ Recode k0 and k1.
	@ recode4_small_inner() does not modify r12.
	ldr	r0, [sp, #72]
	mov	r1, sp
	bl	recode4_small_inner
	mov	r12, r0
	ldr	r0, [sp, #76]
	add	r1, sp, #16
	bl	recode4_small_inner
	add	r0, r12, r0, lsl #1

	@ Release stack buffers.
	add	sp, #80
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(split_recode4_scalar), .-CN(split_recode4_scalar)

@ =======================================================================
@ void CURVE_window_fill_8_to_xu_affine(CURVE_point_affine_xu *win,
@                                       const CURVE_point *P)
@ =======================================================================

	.align	1
	.global	CN(window_fill_8_to_xu_affine)
	.thumb
	.thumb_func
	.type	CN(window_fill_8_to_xu_affine), %function
CN(window_fill_8_to_xu_affine):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	window_fill_8_to_xu_affine_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(window_fill_8_to_xu_affine), .-CN(window_fill_8_to_xu_affine)
