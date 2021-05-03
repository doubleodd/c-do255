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
	@ Find the top word not equal to the sign mask.
	add	r0, r0, r1, lsl #2
	ldr	r3, [r0, #-4]!
	eors	r2, r3, r3, asr #31
	bne	.Lbitlength_cont
	subs	r1, #1
	beq	.Lbitlength_exit_zero
.Lbitlength_loop:
	ldr	r2, [r0, #-4]!
	eors	r2, r2, r3, asr #31
	bne	.Lbitlength_cont
	subs	r1, #1
	bne	.Lbitlength_loop

.Lbitlength_exit_zero:
	@ If we reach that point, then all words are equal to the sign
	@ mask, and the bitlength is zero.
	eors	r0, r0
	bx	lr

.Lbitlength_cont:
	@ r1 contains one more than the index of the top word, which is
	@ in r2 (already sign-adjusted).
	clz	r2, r2
	lsls	r0, r1, #5
	subs	r0, r2
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
	push	{ r4, r5, r6, r7, r8 }

	@ Set actual copy values:
	@   r0   a (adjusted)
	@   r1   b (adjusted)
	@   r2   length of b, in words (adjusted)
	@   r4   shift count (0..31)
	lsrs	r4, r3, #5
	subs	r2, r4
	bls	.Ladd_lshift_exit
	add	r0, r0, r4, lsl #2
	and	r4, r3, #0x1F

	@ If there is only one remaining word, process it now.
	subs	r2, #1
	beq	.Ladd_lshift_oneword

	@ Multiplier (for a left-shift by r4 bits).
	movs	r3, #1
	lsls	r3, r4

	@ We process two words at each iteration.
	@ Carry word for shift+addition is in r8.
	eor	r8, r8
.Ladd_lshift_loop1:
	ldrd	r4, r5, [r0]
	ldm	r1!, { r6, r7 }
	umaal	r4, r8, r3, r6
	umaal	r5, r8, r3, r7
	stm	r0!, { r4, r5 }
	subs	r2, #2
	bhi	.Ladd_lshift_loop1
	@ Process the final word (if any).
	bcc	.Ladd_lshift_exit
	ldr	r4, [r0]
	ldm	r1!, { r6 }
	umaal	r4, r8, r3, r6
	stm	r0!, { r4 }

.Ladd_lshift_exit:
	pop	{ r4, r5, r6, r7, r8 }
	bx	lr

.Ladd_lshift_oneword:
	@ Only one word to process; shift count in r4.
	ldr	r5, [r0]
	ldr	r6, [r1]
	lsls	r6, r4
	adds	r5, r6
	str	r5, [r0]
	pop	{ r4, r5, r6, r7, r8 }
	bx	lr
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
	push	{ r4, r5, r6, r7, r8 }

	@ Set actual copy values:
	@   r0   a (adjusted)
	@   r1   b (adjusted)
	@   r2   length of b, in words (adjusted)
	@   r4   shift count (0..31)
	lsrs	r4, r3, #5
	subs	r2, r4
	bls	.Lsub_lshift_exit
	add	r0, r0, r4, lsl #2
	and	r4, r3, #0x1F

	@ Multiplier (for a left-shift by r4 bits).
	@ We use umaal for the shift, so we need to subtract 1 from the
	@ multiplier.
	movs	r3, #1
	lsls	r3, r4
	subs	r3, #1

	@ Carry word for shift is in r8.
	eor	r8, r8

	@ If there is an odd number of words, process one now.
	@ Since r2 < 2^31, the subs opcode will make a borrow if and
	@ only if r2 is odd.
	subs	r4, r2, r2, lsl #31
	lsr	r2, #1
	bcs	.Lsub_lshift_cont
	ldr	r4, [r0]
	ldm	r1!, { r6 }
	umaal	r6, r8, r3, r6
	subs	r4, r6
	stm	r0!, { r4 }

.Lsub_lshift_cont:
	@ r2 contains the number of iterations (with two words being
	@ processed at each iteration).
	@ Current borrow has been set (if jumping here directly, the
	@ subs opcode set the carry, i.e. cleared the borrow).

	@ We process two words at each iteration.
.Lsub_lshift_loop:
	cbz	r2, .Lsub_lshift_exit
	ldrd	r4, r5, [r0]
	ldm	r1!, { r6, r7 }
	umaal	r6, r8, r3, r6
	umaal	r7, r8, r3, r7
	sbcs	r4, r6
	sbcs	r5, r7
	stm	r0!, { r4, r5 }
	sub	r2, #1
	b	.Lsub_lshift_loop

.Lsub_lshift_exit:
	pop	{ r4, r5, r6, r7, r8 }
	bx	lr
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
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }

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
	mov	r1, sp
	adr	r2, const_reduce_basis_r
	bl	CN(mul256x256)

	@ nv <- k^2 + 1
	add	r0, sp, #144
	mov	r1, sp
	mov	r2, sp
	bl	CN(mul256x256)

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
	add	r0, r0, r5, lsl #6
	sub	r1, r1, r5, lsl #6
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
	eor	r5, r5, #1

.Lreduce_basis_cont3:
	@ nu >= nv
	@ Compute current size. We scan for the topmost non-zero word of
	@ nu, and adjust the size accordingly. We must account for a possible
	@ sign bit in pp, i.e. add one to the length if the top word of nu
	@ has its top bit set.
	@ Since nu is always non-zero, we do not need a loop guard here.
	add	r0, sp, #80
	add	r4, r0, r5, lsl #6
	add	r0, r4, r6, lsl #2
.Lreduce_basis_loop4:
	ldr	r3, [r0, #-4]!
	cmp	r3, #0
	beq	.Lreduce_basis_loop4
	subs	r6, r0, r4
	lsrs	r6, #2
	adds	r6, #1
	add	r6, r6, r3, lsr #31

	@ Get bit length of nv; if 255 or lower, then we exit.
	add	r0, sp, #144
	sub	r0, r0, r5, lsl #6
	movs	r1, r6
	bl	bitlength_inner
	cmp	r0, #255
	bhi	.Lreduce_basis_cont4
	@ We copy the low 17 bytes of v0 and v1 into the output buffers.
	ldr	r0, [sp, #272]
	add	r1, sp, #40
	sub	r1, r1, r5, lsl #5
	sub	r1, r1, r5, lsl #3
	movs	r2, #17
	bl	memcpy
	ldr	r0, [sp, #276]
	add	r1, sp, #60
	sub	r1, r1, r5, lsl #5
	sub	r1, r1, r5, lsl #3
	movs	r2, #17
	bl	memcpy
	add	sp, #280
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }

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
	ldr	r0, [r0, r6, lsl #2]
	asrs	r0, r0, #31
	bne	.Lreduce_basis_cont5

	@ sp >= 0
	add	r7, r5, r5, lsl #2
	add	r0, sp, r7, lsl #3
	add	r1, sp, #40
	sub	r1, r1, r7, lsl #3
	movs	r2, #5
	movs	r3, r4
	bl	sub_lshift_inner
	add	r0, sp, #20
	add	r0, r0, r7, lsl #3
	add	r1, sp, #60
	sub	r1, r1, r7, lsl #3
	movs	r2, #5
	movs	r3, r4
	bl	sub_lshift_inner
	add	r0, sp, #80
	add	r0, r0, r5, lsl #6
	add	r1, sp, #144
	sub	r1, r1, r5, lsl #6
	movs	r2, r6
	lsls	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r0, r0, r5, lsl #6
	add	r1, sp, #208
	movs	r2, r6
	adds	r3, r4, #1
	bl	sub_lshift_inner
	add	r0, sp, #208
	add	r1, sp, #144
	sub	r1, r1, r5, lsl #6
	movs	r2, r6
	movs	r3, r4
	bl	sub_lshift_inner

	@ Iteration done, we loop.
	b	.Lreduce_basis_loop2

.Lreduce_basis_cont5:
	@ sp < 0
	add	r7, r5, r5, lsl #2
	add	r0, sp, r7, lsl #3
	add	r1, sp, #40
	sub	r1, r1, r7, lsl #3
	movs	r2, #5
	movs	r3, r4
	bl	add_lshift_inner
	add	r0, sp, #20
	add	r0, r0, r7, lsl #3
	add	r1, sp, #60
	sub	r1, r1, r7, lsl #3
	movs	r2, #5
	movs	r3, r4
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r0, r0, r5, lsl #6
	add	r1, sp, #144
	sub	r1, r1, r5, lsl #6
	movs	r2, r6
	lsls	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #80
	add	r0, r0, r5, lsl #6
	add	r1, sp, #208
	movs	r2, r6
	adds	r3, r4, #1
	bl	add_lshift_inner
	add	r0, sp, #208
	add	r1, sp, #144
	sub	r1, r1, r5, lsl #6
	movs	r2, r6
	movs	r3, r4
	bl	add_lshift_inner

	@ Iteration done, we loop.
	b	.Lreduce_basis_loop2

	.size	CN(reduce_basis_vartime), .-CN(reduce_basis_vartime)
