@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It implements Lagrange's algorithm for lattice basis reduction in
@ dimenstion two. The following macro must have been defined:
@   CN(name)   makes a symbol name suitable for exporting
@ =======================================================================

@ =======================================================================
@ int bitlength_inner(const uint32_t *x, int len)
@
@ Compute the bit length of a signed integer (len words, with len > 0).
@
@ ABI: registers r0-r3 are consumed; other registers are untouched.
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	bitlength_inner, %function
bitlength_inner:
	@ Make r1 the offset of the last word (in bytes).
	subs	r1, #1
	lsls	r1, #2

	@ Get sign mask.
	ldr	r3, [r0, r1]
	asrs	r3, r3, #31

	@ Find the top word not equal to the sign mask.
.Lbitlength__1:
	ldr	r2, [r0, r1]
	eors	r2, r3
	bne	.Lbitlength__2
	subs	r1, #4
	bpl	.Lbitlength__1

	@ All limbs are equal to the sign mask; bit length is zero.
	movs	r0, #0
	bx	lr

.Lbitlength__2:
	@ Word r2 (sign-adjusted) was found at offset r1 (in bytes).
	lsls	r0, r1, #3
	adds	r0, #1

	lsrs	r1, r2, #16
	beq	.Lbitlength__3
	movs	r2, r1
	adds	r0, #16
.Lbitlength__3:
	lsrs	r1, r2, #8
	beq	.Lbitlength__4
	movs	r2, r1
	adds	r0, #8
.Lbitlength__4:
	lsrs	r1, r2, #4
	beq	.Lbitlength__5
	movs	r2, r1
	adds	r0, #4
.Lbitlength__5:
	lsrs	r1, r2, #2
	beq	.Lbitlength__6
	movs	r2, r1
	adds	r0, #2
.Lbitlength__6:
	lsrs	r1, r2, #1
	beq	.Lbitlength__7
	adds	r0, #1
.Lbitlength__7:
	bx	lr
	.size	bitlength_inner, .-bitlength_inner

@ =======================================================================
@ void add_lshift(uint32_t *a, const uint32_t *b, int len, int s)
@
@ ABI: registers r0-r3 are consumed.
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	add_lshift_inner, %function
add_lshift_inner:
	push	{ r4, r5, r6, r7, lr }

	@ Set actual copy values:
	@   r0   a (adjusted)
	@   r1   b (adjusted)
	@   r2   length of b, in words (adjusted)
	@   r3   shift count (0..31)
	lsrs	r4, r3, #5
	subs	r2, r4
	bls	.Ladd_lshift__exit
	lsls	r4, r4, #2
	adds	r0, r4
	lsls	r3, r3, #27
	lsrs	r3, r3, #27

	@ If there is no shift, use a specialized loop.
	cmp	r3, #0
	beq	.Ladd_lshift_noshift

	eors	r4, r4        @ Storage for carry (0 or 1)
	eors	r5, r5        @ Partial block for shift
	movs	r6, #32
	subs	r6, r3        @ Shift count complement
	subs	r2, #2
	bcc	.Ladd_lshift_cont
	mov	r14, r2       @ Loop counter
.Ladd_lshift_loop1:
	ldm	r1!, { r7 }
	movs	r2, r7
	lsrs	r2, r6
	lsls	r7, r3
	orrs	r7, r5
	ldr	r5, [r0]
	lsrs	r4, r4, #1
	adcs	r5, r7
	adcs	r4, r4
	stm	r0!, { r5 }
	ldm	r1!, { r7 }
	movs	r5, r7
	lsrs	r5, r6
	lsls	r7, r3
	orrs	r7, r2
	ldr	r2, [r0]
	lsrs	r4, r4, #1
	adcs	r2, r7
	adcs	r4, r4
	stm	r0!, { r2 }
	mov	r2, r14
	subs	r2, #2
	mov	r14, r2
	bcs	.Ladd_lshift_loop1
.Ladd_lshift_cont:
	mov	r2, r14
	adds	r2, #1
	bcc	.Ladd_lshift__exit
	ldm	r1!, { r7 }
	lsls	r7, r3
	orrs	r7, r5
	ldr	r5, [r0]
	lsrs	r4, r4, #1
	adcs	r5, r7
	stm	r0!, { r5 }
	b	.Ladd_lshift__exit

.Ladd_lshift_noshift:
	eors	r3, r3        @ Storage for carry (0 or 1)
	subs	r2, #2
	bcc	.Ladd_lshift_noshift_cont
.Ladd_lshift_noshift_loop2:
	ldm	r0!, { r4, r5 }
	ldm	r1!, { r6, r7 }
	lsrs	r3, r3, #1    @ Retrieve carry
	adcs	r4, r6
	adcs	r5, r7
	adcs	r3, r3        @ Store carry
	subs	r0, #8
	stm	r0!, { r4, r5 }
	subs	r2, #2
	bcs	.Ladd_lshift_noshift_loop2
.Ladd_lshift_noshift_cont:
	adds	r2, #1
	bcc	.Ladd_lshift__exit
	ldr	r4, [r0]
	ldm	r1!, { r6 }
	lsrs	r3, r3, #1    @ Retrieve carry
	adcs	r4, r6
	stm	r0!, { r4 }

.Ladd_lshift__exit:
	pop	{ r4, r5, r6, r7, pc }
	.size	add_lshift_inner, .-add_lshift_inner

@ =======================================================================
@ void sub_lshift(uint32_t *a, const uint32_t *b, int len, int s)
@
@ ABI: registers r0-r3 are consumed.
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	sub_lshift_inner, %function
sub_lshift_inner:
	push	{ r4, r5, r6, r7, lr }

	@ Set actual copy values:
	@   r0   a (adjusted)
	@   r1   b (adjusted)
	@   r2   length of b, in words (adjusted)
	@   r3   shift count (0..31)
	lsrs	r4, r3, #5
	subs	r2, r4
	bls	.Lsub_lshift__exit
	lsls	r4, r4, #2
	adds	r0, r4
	lsls	r3, r3, #27
	lsrs	r3, r3, #27

	@ If there is no shift, use a specialized loop.
	cmp	r3, #0
	beq	.Lsub_lshift_noshift

	movs	r4, #1        @ Storage for carry (0 or 1)
	eors	r5, r5        @ Partial block for shift
	movs	r6, #32
	subs	r6, r3        @ Shift count complement
	subs	r2, #2
	bcc	.Lsub_lshift_cont
	mov	r14, r2       @ Loop counter
.Lsub_lshift_loop1:
	ldm	r1!, { r7 }
	movs	r2, r7
	lsrs	r2, r6
	lsls	r7, r3
	orrs	r7, r5
	ldr	r5, [r0]
	lsrs	r4, r4, #1
	sbcs	r5, r7
	adcs	r4, r4
	stm	r0!, { r5 }
	ldm	r1!, { r7 }
	movs	r5, r7
	lsrs	r5, r6
	lsls	r7, r3
	orrs	r7, r2
	ldr	r2, [r0]
	lsrs	r4, r4, #1
	sbcs	r2, r7
	adcs	r4, r4
	stm	r0!, { r2 }
	mov	r2, r14
	subs	r2, #2
	mov	r14, r2
	bcs	.Lsub_lshift_loop1
.Lsub_lshift_cont:
	mov	r2, r14
	adds	r2, #1
	bcc	.Lsub_lshift__exit
	ldm	r1!, { r7 }
	lsls	r7, r3
	orrs	r7, r5
	ldr	r5, [r0]
	lsrs	r4, r4, #1
	sbcs	r5, r7
	stm	r0!, { r5 }
	b	.Lsub_lshift__exit

.Lsub_lshift_noshift:
	movs	r3, #1        @ Storage for carry (0 or 1)
	subs	r2, #2
	bcc	.Lsub_lshift_noshift_cont
.Lsub_lshift_noshift_loop2:
	ldm	r0!, { r4, r5 }
	ldm	r1!, { r6, r7 }
	lsrs	r3, r3, #1    @ Retrieve carry
	sbcs	r4, r6
	sbcs	r5, r7
	adcs	r3, r3        @ Store carry
	subs	r0, #8
	stm	r0!, { r4, r5 }
	subs	r2, #2
	bcs	.Lsub_lshift_noshift_loop2
.Lsub_lshift_noshift_cont:
	adds	r2, #1
	bcc	.Ladd_lshift__exit
	ldr	r4, [r0]
	ldm	r1!, { r6 }
	lsrs	r3, r3, #1    @ Retrieve carry
	sbcs	r4, r6
	stm	r0!, { r4 }

.Lsub_lshift__exit:
	pop	{ r4, r5, r6, r7, pc }
	.size	sub_lshift_inner, .-sub_lshift_inner

@ =======================================================================
@ void CURVE_reduce_basis(uint8_t *c0, uint8_t *c1, const uint8_t *k)
@
@ This function follows the external ABI.
@ =======================================================================

	.align	1
	.global	CN(reduce_basis_vartime)
	.thumb
	.thumb_func
	.type	CN(reduce_basis_vartime), %function
CN(reduce_basis_vartime):
	push	{ r4, r5, r6, r7, lr }
	mov	r4, r8
	mov	r5, r10
	mov	r6, r11
	push	{ r4, r5, r6 }

	@ Data layout on stack, depending on the value of r5:
	@    r5 = 0 or 1:
	@       if r5 == 0 then:
	@            sp+0     u0
	@            sp+20    u1
	@            sp+40    v0
	@            sp+60    v1
	@            sp+80    nu
	@            sp+144   nv
	@       else:
	@            sp+0     v0
	@            sp+20    v1
	@            sp+40    u0
	@            sp+60    u1
	@            sp+80    nv
	@            sp+144   nu
	@ Therefore:
	@    u0   is at address   sp+40*r5
	@    u1   is at address   sp+20+40*r5
	@    v0   is at address   sp+40-40*r5
	@    v1   is at address   sp+60-40*r5
	@    nu   is at address   sp+80+64*r5
	@    nv   is at address   sp+144-64*r5
	@    pp   is at address   sp+208
	@
	@ Extra slots:
	@    272   &c0
	@    276   &c1

	push	{ r0, r1 }
	sub	sp, #272

	@ Decode and reduce input value (k) into stack (0).
	mov	r0, sp
	movs	r1, r2
	bl	scalar_decode_inner

	@ pp <- r*k
	add	r0, sp, #208
	movs	r1, #16
	mov	r2, sp
	movs	r3, #8
	adr	r4, const_reduce_basis_r
	movs	r5, #8
	bl	iXX_mul_inner

	@ nv <- k^2 + 1
	add	r0, sp, #144
	movs	r1, #16
	mov	r2, sp
	movs	r3, #8
	mov	r4, sp
	movs	r5, #8
	bl	iXX_mul_inner
	add	r0, sp, #144
	add	r1, sp, #144
	eors	r7, r7
	ldm	r1!, { r2, r3, r4, r5, r6 }
	adds	r2, #1
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	stm	r0!, { r2, r3, r4, r5, r6 }
	ldm	r1!, { r2, r3, r4, r5, r6 }
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	stm	r0!, { r2, r3, r4, r5, r6 }
	ldm	r1, { r1, r2, r3, r4, r5, r6 }
	adcs	r1, r7
	adcs	r2, r7
	adcs	r3, r7
	adcs	r4, r7
	adcs	r5, r7
	adcs	r6, r7
	stm	r0!, { r1, r2, r3, r4, r5, r6 }

	@ nu <- r^2
	add	r0, sp, #80
	adr	r1, const_reduce_basis_r2
	ldm	r1!, { r2, r3, r4, r5, r6 }
	stm	r0!, { r2, r3, r4, r5, r6 }
	ldm	r1!, { r2, r3, r4, r5, r6 }
	stm	r0!, { r2, r3, r4, r5, r6 }
	ldm	r1, { r1, r2, r3, r4, r5, r6 }
	stm	r0!, { r1, r2, r3, r4, r5, r6 }

	@ v0 <- k (truncated), v1 <- 1
	add	r0, sp, #40
	mov	r1, sp
	ldm	r1!, { r2, r3, r4, r5, r6 }
	stm	r0!, { r2, r3, r4, r5, r6 }
	movs	r2, #1
	eors	r3, r3
	eors	r4, r4
	stm	r0!, { r2, r3, r4 }
	stm	r0!, { r3, r4 }

	@ u0 <- r (truncated), u1 <- 0
	mov	r0, sp
	adr	r1, const_reduce_basis_r
	ldm	r1!, { r2, r3, r4, r5, r6 }
	stm	r0!, { r2, r3, r4, r5, r6 }
	eors	r2, r2
	eors	r3, r3
	eors	r4, r4
	stm	r0!, { r2, r3, r4 }
	stm	r0!, { r3, r4 }

	@ Jump to second part.
	b	.Lreduce_basis_cont2

	.align	2
const_reduce_basis_r:
	.long	0x396152C7, 0xDCF2AC65, 0x912B7F03, 0x2ACF567A
	.long	0x00000000, 0x00000000, 0x00000000, 0x40000000
const_reduce_basis_r2:
	.long	0x739216B1, 0xA31F34E2, 0x835B5211, 0x86A297C9
	.long	0xF04303AD, 0x95DCE66B, 0x2F0F9E3C, 0x8728B04D
	.long	0x9CB0A963, 0xEE795632, 0x4895BF81, 0x1567AB3D
	.long	0x00000000, 0x00000000, 0x00000000, 0x10000000

.Lreduce_basis_cont2:
	@ Initialization is done, here is the actual loop.
	@ Variables:
	@    r5   0 or 1, depending on current swap status
	@    r6   current length of nu, nv and pp
	@ Called functions consume only registers r0-r3.
	eors	r5, r5
	movs	r6, #16

.Lreduce_basis_loop2:
	@ Compare nu and nv; if nu < nv, then swap.
	add	r0, sp, #80
	add	r1, sp, #144
	lsls	r2, r5, #6
	adds	r0, r2
	subs	r1, r2
	lsls	r2, r6, #2
.Lreduce_basis_loop3:
	subs	r2, #4
	bcc	.Lreduce_basis_cont3
	ldr	r3, [r0, r2]
	ldr	r4, [r1, r2]
	cmp	r3, r4
	beq	.Lreduce_basis_loop3
	cmp	r4, r3
	bls	.Lreduce_basis_cont3
	movs	r0, #1
	eors	r5, r0

.Lreduce_basis_cont3:
	@ nu >= nv
	@ Compute current size. We scan for the topmost non-zero word of
	@ nu, and adjust the size accordingly. We must account for a possible
	@ sign bit in pp, i.e. add one to the length if the top word of nu
	@ has its top bit set.
	@ Since nu is always non-zero, we do not need a loop guard here.
	add	r0, sp, #80
	lsls	r2, r5, #6
	adds	r0, r2
	lsls	r2, r6, #2
.Lreduce_basis_loop4:
	subs	r2, #4
	ldr	r3, [r0, r2]
	cmp	r3, #0
	beq	.Lreduce_basis_loop4
	lsrs	r6, r2, #2
	adds	r6, #1
	lsrs	r3, r3, #31
	adds	r6, r3

	@ Get bit length of nv; if 255 or lower, then we exit.
	add	r0, sp, #144
	lsls	r2, r5, #6
	subs	r0, r2
	movs	r1, r6
	bl	bitlength_inner
	cmp	r0, #255
	bhi	.Lreduce_basis_cont4
	@ We copy the low 17 bytes of v0 and v1 into the output buffers.
	ldr	r0, [sp, #272]
	add	r1, sp, #40
	movs	r3, #40
	muls	r3, r5
	subs	r1, r3
	movs	r2, #17
	bl	memcpy
	ldr	r0, [sp, #276]
	add	r1, sp, #60
	movs	r3, #40
	muls	r3, r5
	subs	r1, r3
	movs	r2, #17
	bl	memcpy
	add	sp, #280
	pop	{ r4, r5, r6 }
	mov	r8, r4
	mov	r10, r5
	mov	r11, r6
	pop	{ r4, r5, r6, r7, pc }

.Lreduce_basis_cont4:
	@ Get bit length of pp and compute shift count s (into r4).
	movs	r4, r0
	add	r0, sp, #208
	movs	r1, r6
	bl	bitlength_inner
	subs	r4, r0, r4
	sbcs	r7, r7
	eors	r4, r7

	@ If sp > 0, then:
	@    u <- u - lshift(v, s)
	@    nu <- nu + lshift(nv, 2*s) - lshift(pp, s+1)
	@    pp <- pp - lshift(nv, s)
	@ else:
	@    u <- u + lshift(v, s)
	@    nu <- nu + lshift(nv, 2*s) + lshift(pp, s+1)
	@    pp <- pp + lshift(nv, s)
	add	r0, sp, #204   @ address of pp, minus 4.
	lsls	r1, r6, #2
	ldr	r0, [r0, r1]
	asrs	r0, r0, #31
	bne	.Lreduce_basis_cont5

	@ sp >= 0
	mov	r0, sp
	add	r1, sp, #40
	movs	r3, #40
	muls	r3, r5
	adds	r0, r3
	subs	r1, r3
	movs	r2, #5
	movs	r3, r4
	bl	sub_lshift_inner
	add	r0, sp, #20
	add	r1, sp, #60
	movs	r3, #40
	muls	r3, r5
	adds	r0, r3
	subs	r1, r3
	movs	r2, #5
	movs	r3, r4
	bl	sub_lshift_inner
	add	r0, sp, #80
	add	r1, sp, #144
	lsls	r3, r5, #6
	adds	r0, r3
	subs	r1, r3
	movs	r2, r6
	lsls	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r1, sp, #208
	lsls	r3, r5, #6
	adds	r0, r3
	movs	r2, r6
	adds	r3, r4, #1
	bl	sub_lshift_inner
	add	r0, sp, #208
	add	r1, sp, #144
	lsls	r3, r5, #6
	subs	r1, r3
	movs	r2, r6
	movs	r3, r4
	bl	sub_lshift_inner

	@ Iteration done, we loop.
	b	.Lreduce_basis_loop2

.Lreduce_basis_cont5:
	@ sp < 0
	mov	r0, sp
	add	r1, sp, #40
	movs	r3, #40
	muls	r3, r5
	adds	r0, r3
	subs	r1, r3
	movs	r2, #5
	movs	r3, r4
	bl	add_lshift_inner
	add	r0, sp, #20
	add	r1, sp, #60
	movs	r3, #40
	muls	r3, r5
	adds	r0, r3
	subs	r1, r3
	movs	r2, #5
	movs	r3, r4
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r1, sp, #144
	lsls	r3, r5, #6
	adds	r0, r3
	subs	r1, r3
	movs	r2, r6
	lsls	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r1, sp, #208
	lsls	r3, r5, #6
	adds	r0, r3
	movs	r2, r6
	adds	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #208
	add	r1, sp, #144
	lsls	r3, r5, #6
	subs	r1, r3
	movs	r2, r6
	movs	r3, r4
	bl	add_lshift_inner

	@ Iteration done, we loop.
	b	.Lreduce_basis_loop2

	.size	CN(reduce_basis_vartime), .-CN(reduce_basis_vartime)
