@ =======================================================================
@ This file is not meant to be compiled by itself; it is included through
@ the preprocessor.
@ It defines functions for computation in the base field. The following
@ macros must have been defined:
@   CN(name)   makes a symbol name suitable for exporting
@   MQ         literal integer (1..32767 range); field modulus is: 2^255 - MQ
@   INVT510    contents of 1/2^510 in the field (eight 32-bit constants,
@              comma separated, little-endian order)
@ =======================================================================

@ =======================================================================
@ CONVENTIONS: we use the notion of an "accumulator", i.e. a value which
@ is held in registers r4:r5:r6:r7:r8:r10:r11:r12. Core arithmetic
@ functions use the accumulator as one input, and return the value in
@ the accumulator; they are tagged "_acc_inner". Such functions use an
@ internal ABI that does not preserve any other register.
@ =======================================================================

@ debug macros, disabled.

@.macro DEBUG_GEN  skip, num, tt
@	push    { r0, r1, r2, r3, r4, r5, r6, r7 }
@	mov     r0, r8
@	mov     r1, r9
@	mov     r2, r10
@	mov     r3, r11
@	mov     r4, r12
@	mrs     r5, APSR
@	mov     r6, r14
@	movs    r7, #0
@	push    { r0, r1, r2, r3, r4, r5, r6, r7 }
@
@	mov     r0, sp
@	movs    r1, #(\skip)
@	movs    r2, #(\num)
@	movs    r3, #(\tt)
@	bl      print_DEBUG
@
@	pop     { r0, r1, r2, r3, r4, r5, r6, r7 }
@	mov     r8, r0
@	mov     r10, r2
@	mov     r11, r3
@	mov     r12, r4
@	msr     APSR_nzcvq, r5
@	mov     r14, r6
@	pop     { r0, r1, r2, r3, r4, r5, r6, r7 }
@.endm
@
@.macro DEBUG  num
@	DEBUG_GEN  0, \num, 1
@.endm

@ =======================================================================

@ Load the accumulator, using the provided register as source pointer.
@ The source pointer is not modified, unless it is one of the accumulator
@ registers, in which case it is overwritten in the process.
.macro ACC_LOAD  rp
	ldm	\rp, { r4, r5, r6, r7, r8, r10, r11, r12 }
.endm

@ Equivalent to ACC_LOAD, except that the source pointer is incremented
@ by 32 as a side-effect.
.macro ACC_LOAD_UPDATE  rp
	ldm	\rp!, { r4, r5, r6, r7, r8, r10, r11, r12 }
.endm

@ Alternate macro for loading the accumulator, with tolerance for
@ an unaligned input buffer.
.macro ACC_LOAD_UNALIGNED  rp
	ldr	r4, [\rp]
	ldr	r5, [\rp, #4]
	ldr	r6, [\rp, #8]
	ldr	r7, [\rp, #12]
	ldr	r8, [\rp, #16]
	ldr	r10, [\rp, #20]
	ldr	r11, [\rp, #24]
	ldr	r12, [\rp, #28]
.endm

@ Store the accumulator, using the provided register as destination
@ pointer. The destination pointer must be in { r0, r1, r2, r3, sp, r14 };
@ it is not modified.
.macro ACC_STORE  rp
	stm	\rp, { r4, r5, r6, r7, r8, r10, r11, r12 }
.endm

@ Equivalent to ACC_STORE, except that the source pointer is incremented
@ by 32 as a side-effect.
.macro ACC_STORE_UPDATE  rp
	stm	\rp!, { r4, r5, r6, r7, r8, r10, r11, r12 }
.endm

@ Reserve 32 bytes on the stack, pushing the accumulator onto it. The
@ stack pointer is updated (decremented by 32).
.macro ACC_PUSH
	push	{ r4, r5, r6, r7, r8, r10, r11, r12 }
.endm

@ Alternate macro for storing the accumulator, with tolerance for
@ an unaligned output buffer.
.macro ACC_STORE_UNALIGNED  rp
	str	r4, [\rp]
	str	r5, [\rp, #4]
	str	r6, [\rp, #8]
	str	r7, [\rp, #12]
	str	r8, [\rp, #16]
	str	r10, [\rp, #20]
	str	r11, [\rp, #24]
	str	r12, [\rp, #28]
.endm

@ Multiply the accumulator by a small integer. The small integer is
@ provided in register rm and must be lower than 2^16. Two extra
@ scratch registers ru and rv are provided. All three of rm, ru and
@ rv are consumed. rm, ru and rv must be among { r0, r1, r2, r3, r14 }.
@ Flags are not modified.
@ Cost: 13
.macro ACC_MULSMALL  rm, ru, rv
	@ Do multiplication on top word and fold the upper bits.
	@ We compute 2*a7*m into r12:r2; the upper half (rv) contains
	@ the bits that must be folded.
	add	\ru, \rm, \rm
	umull	r12, \rv, \ru, r12
	movw	\ru, #MQ
	mul	\ru, \ru, \rv

	@ Apply multiplication on words a0..a6. Running carry word is
	@ in ru (initialized with the folding result).
	@ We compute x*m as x+ru+(m-1)*x.
	sub	\rm, \rm, #1
	umaal	r4, \ru, \rm, r4
	umaal	r5, \ru, \rm, r5
	umaal	r6, \ru, \rm, r6
	umaal	r7, \ru, \rm, r7
	umaal	r8, \ru, \rm, r8
	umaal	r10, \ru, \rm, r10
	umaal	r11, \ru, \rm, r11

	@ We must add the carry word to the top word which was computed
	@ first. Also, r12 currently contains twice the top word, so it
	@ must be halved prior to the addition. Since the carry word is
	@ small, and r12 is first halved, this addition cannot overflow.
	add	r12, \ru, r12, lsr #1
.endm

@ Multiply the accumulator by a small constant m. The constant must be
@ between 1 and 32767 (inclusive). Two scratch registers \ru and \rv are
@ provided (taken from { r0, r1, r2, r3, r14 }).
@ Flags are not modified.
@ Cost: 13
.macro ACC_MULCONST  m, ru, rv
	@ See ACC_MULSMALL for details.
	movw	\ru, #(2 * (\m))
	umull	r12, \rv, \ru, r12
	movw	\ru, #MQ
	mul	\ru, \ru, \rv
	movw	\rv, #((\m) - 1)
	umaal	r4, \ru, \rv, r4
	umaal	r5, \ru, \rv, r5
	umaal	r6, \ru, \rv, r6
	umaal	r7, \ru, \rv, r7
	umaal	r8, \ru, \rv, r8
	umaal	r10, \ru, \rv, r10
	umaal	r11, \ru, \rv, r11
	add	r12, \ru, r12, lsr #1
.endm

@ =======================================================================
@ void gf_add_acc_inner(gf *b)
@
@ Add field element b to accumulator.
@
@ Cost: 30
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_add_acc_inner, %function
gf_add_acc_inner:
	@ Add top words and truncate to fold extra bits.
	ldr	r1, [r0, #28]
	movs	r3, #2
	mov	r2, r12
	umaal	r12, r2, r3, r1
	movw	r3, #MQ
	muls	r2, r3

	@ Compute the remainder of the addition. Carry word is in r2.
	movs	r3, #1
	ldm	r0!, { r1 }
	umaal	r4, r2, r3, r1
	ldm	r0!, { r1 }
	umaal	r5, r2, r3, r1
	ldm	r0!, { r1 }
	umaal	r6, r2, r3, r1
	ldm	r0!, { r1 }
	umaal	r7, r2, r3, r1
	ldm	r0!, { r1 }
	umaal	r8, r2, r3, r1
	ldm	r0, { r0, r1 }
	umaal	r10, r2, r3, r0
	umaal	r11, r2, r3, r1

	@ Propagate carry into top word (which must also be halved).
	add	r12, r2, r12, lsr #1
	bx	lr
	.size	gf_add_acc_inner, .-gf_add_acc_inner

@ =======================================================================
@ void gf_sub_acc_inner(gf *b)
@
@ Subtract field element b from accumulator.
@
@ Cost: 35
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_sub_acc_inner, %function
gf_sub_acc_inner:
	@ Do subtraction with carry propagation.
	ldm	r0!, { r1, r2, r3 }
	subs	r4, r4, r1
	sbcs	r5, r5, r2
	sbcs	r6, r6, r3
	ldm	r0!, { r1, r2, r3 }
	sbcs	r7, r7, r1
	sbcs	r8, r8, r2
	sbcs	r10, r10, r3
	ldm	r0!, { r1, r2 }
	sbcs	r11, r11, r1
	sbcs	r12, r12, r2

	@ If there is a borrow (br) then we must subtract 2^256. Since
	@ 2^256 = 2*mq mod q, this is equivalent to subtracting 2*mq.
	@ We do the subtraction by adding the 256-bit value
	@ r0:r1:r1:...:r1 = 2^256 - 2*mq  (or 0, if there was no borrow).
	sbcs	r1, r1
	movw	r2, #(2 * MQ)
	muls	r2, r1
	adds	r4, r2
	adcs	r5, r1
	adcs	r6, r1
	adcs	r7, r1
	adcs	r8, r1
	adcs	r10, r1
	adcs	r11, r1
	adcs	r12, r1

	@ Usually, if there was a borrow (r1 == -1), we get a carry here,
	@ which cancels the borrow (and we cannot get a carry if there was
	@ no borrow, because in that case we added zero). However, if the
	@ value after the addition was in the range 2^256..2^256+2*mq-1,
	@ then there might not be a carry, in which case the operations
	@ above yielded a 256-bit value in the 2^256-2*mq..2^256-1. In
	@ that case, we add 2*q again, i.e. 2^256-2*mq. The subtraction of
	@ 2*mq cannot yield a borrow, so only one word needs to be
	@ modified.

	@ If there was initially a borrow, then r2 contains 2^32-2*mq, and
	@ r1 contains -1; otherwise, r2 and r1 both contain 0.
	adcs	r1, #0
	ands	r2, r1
	adds	r4, r2

	bx	lr
	.size	gf_sub_acc_inner, .-gf_sub_acc_inner

@ =======================================================================
@ void gf_rsub_acc_inner(gf *b)
@
@ Subtract accumulator from field element b, result in accumulator
@ (reverse subtraction).
@
@ Cost: 35
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_rsub_acc_inner, %function
gf_rsub_acc_inner:
	@ Do subtraction with carry propagation.
	ldm	r0!, { r1, r2, r3 }
	subs	r4, r1, r4
	sbcs	r5, r2, r5
	sbcs	r6, r3, r6
	ldm	r0!, { r1, r2, r3 }
	sbcs	r7, r1, r7
	sbcs	r8, r2, r8
	sbcs	r10, r3, r10
	ldm	r0!, { r1, r2 }
	sbcs	r11, r1, r11
	sbcs	r12, r2, r12

	@ Rest of the code is identical to gf_sub_acc_inner; see
	@ the comments there.
	sbcs	r1, r1
	movw	r2, #(2 * MQ)
	muls	r2, r1
	adds	r4, r2
	adcs	r5, r1
	adcs	r6, r1
	adcs	r7, r1
	adcs	r8, r1
	adcs	r10, r1
	adcs	r11, r1
	adcs	r12, r1

	adcs	r1, #0
	ands	r2, r1
	adds	r4, r2

	bx	lr
	.size	gf_rsub_acc_inner, .-gf_rsub_acc_inner

@ =======================================================================
@ void gf_neg_acc_inner(void)
@
@ Negate accumulator (in the field).
@
@ Cost: 16
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_neg_acc_inner, %function
gf_neg_acc_inner:
	@ To negate, we first apply a partial reduction, i.e. we get
	@ the high bit and propagate it (with value multiplied by mq);
	@ then the result is subtracted from 2*q = 2^256 - 2*mq. The
	@ two operations are combined:
	@  - we extract the high bit h, and clear that bit in the value;
	@  - we then subtract the result from 2^256 - 2*mq - h*mq.

	@ Get h*mq into r1, and truncate the high word.
	movw	r0, #MQ
	and	r1, r0, r12, asr #31
	bic	r12, r12, #0x80000000

	@ Compute (-2*mq - h*mq) mod 2^32 into r0, and 0xFFFFFFFF into r1.
	add	r0, r1, r0, lsl #1
	rsbs	r0, r0, #0
	orn	r1, r1

	@ Do the subtraction. Since the top word was truncated, this cannot
	@ overflow.
	subs	r4, r0, r4
	sbcs	r5, r1, r5
	sbcs	r6, r1, r6
	sbcs	r7, r1, r7
	sbcs	r8, r1, r8
	sbcs	r10, r1, r10
	sbcs	r11, r1, r11
	sbcs	r12, r1, r12

	bx	lr
	.size	gf_neg_acc_inner, .-gf_neg_acc_inner

@ =======================================================================
@ void gf_condneg_acc_inner(uint32_t ctl)
@
@ Negate accumulator (in the field) if ctl == 1.
@ Do nothing if ctl == 0.
@ ctl MUST be 0 or 1.
@
@ Cost: 26
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_condneg_acc_inner, %function
gf_condneg_acc_inner:
	@ Load GE flags with ctl.
	rsbs	r0, r0, #0
	uadd8	r0, r0, r0

	@ Compute the negation (as in gf_neg_acc_inner()) but use the
	@ sel opcode to keep the result only if ctl == 1.
	movw	r0, #MQ
	and	r1, r0, r12, asr #31
	bic	r3, r12, #0x80000000
	add	r0, r1, r0, lsl #1
	rsbs	r0, r0, #0
	orn	r1, r1
	subs	r2, r0, r4
	sel	r4, r2, r4
	sbcs	r2, r1, r5
	sel	r5, r2, r5
	sbcs	r2, r1, r6
	sel	r6, r2, r6
	sbcs	r2, r1, r7
	sel	r7, r2, r7
	sbcs	r2, r1, r8
	sel	r8, r2, r8
	sbcs	r2, r1, r10
	sel	r10, r2, r10
	sbcs	r2, r1, r11
	sel	r11, r2, r11
	sbcs	r2, r1, r3
	sel	r12, r2, r12

	bx	lr
	.size	gf_condneg_acc_inner, .-gf_condneg_acc_inner

@ =======================================================================
@ void gf_normalize_acc_inner(void)
@
@ Normalize accumulator contents into 0..q-1 range.
@
@ Cost: 25
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_normalize_acc_inner, %function
gf_normalize_acc_inner:
	@ Fold the top bit and add mq at the same time.
	lsrs	r1, r12, #31
	movw	r0, #MQ
	umlal	r0, r1, r0, r1
	bic	r12, r12, #0x80000000
	adds	r4, r0
	adcs	r5, #0
	adcs	r6, #0
	adcs	r7, #0
	adcs	r8, #0
	adcs	r10, #0
	adcs	r11, #0
	adc	r12, #0

	@ We have a value equal to a+mq (modulo q), in the mq..2^255+2*mq
	@ range. If its top bit is set, then we want to subtract 2^255
	@ from it (since a + mq - 2^255 = a mod (2^255-mq)); this will not
	@ overflow since in that case the value is at least equal to 2^255.
	@ Otherwise, the intermediate value is in mq..2^255-1 and we want
	@ to subtract mq from it, to obtain a result in 0..q-1 and equal
	@ to the original value modulo q.
	movw	r0, #MQ
	bic	r0, r0, r12, asr #31
	bic	r12, r12, #0x80000000
	subs	r4, r0
	sbcs	r5, #0
	sbcs	r6, #0
	sbcs	r7, #0
	sbcs	r8, #0
	sbcs	r10, #0
	sbcs	r11, #0
	sbcs	r12, #0

	bx	lr
	.size	gf_normalize_acc_inner, .-gf_normalize_acc_inner

@ =======================================================================
@ int gf_iszero_acc_inner(void)
@
@ Set r0 to 1 if the value in the accumulator is 0, or to 0 otherwise.
@
@ Cost: 31
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_iszero_acc_inner, %function
gf_iszero_acc_inner:
	@ There are three possible representations for 0: 0, q and 2*q.
	@ We accumulate three masks in r1, r2 and r3; at the end:
	@  - if value is 0, then r1 = 0; otherwise, r1 != 0
	@  - if value is q, then r6 = -1; otherwise, r6 != -1
	@  - if value is 2*q, then r7 = -1; otherwise, r7 != -1

	@ Set r1 to 0 is the value is all-zeros.
	orr	r1, r4, r5
	orr	r1, r1, r6
	orr	r1, r1, r7
	orr	r1, r1, r8
	orr	r1, r1, r10
	orr	r1, r1, r11
	orr	r1, r1, r12

	@ Set r2 to 0xFFFFFFFF if a1..a6 are all-ones.
	and	r2, r5, r6
	and	r2, r2, r7
	and	r2, r2, r8
	and	r2, r2, r10
	and	r2, r2, r11

	@ Return value is 1 if any of the following holds:
	@   r1 == 0
	@   r2 == 0xFFFFFFFF and r4 == -mq and r12 == 0x7FFFFFFF
	@   r2 == 0xFFFFFFFF and r4 == -2*mq and r12 == 0xFFFFFFFF

	@ r1 == 0xFFFFFFFF iff a == 0.
	mvn	r1, r1

	@ r3 == 0xFFFFFFFF iff a == 2*q.
	movw	r0, #(MQ - 1)
	add	r3, r4, r0, lsl #1
	adds	r3, #1
	and	r3, r3, r2
	and	r3, r3, r12

	@ r2 == 0xFFFFFFFF iff a == q.
	adds	r0, r4
	and	r2, r2, r0
	eor	r0, r12, #0x80000000
	and	r2, r2, r0

	@ Result is 1 iff one of r1, r2 or r3 is equal to 0xFFFFFFFF.
	eors	r0, r0
	adds	r1, #1
	adcs	r0, #0
	adds	r2, #1
	adcs	r0, #0
	adds	r3, #1
	adcs	r0, #0

	bx	lr
	.size	gf_iszero_acc_inner, .-gf_iszero_acc_inner

@ =======================================================================
@ void gf_half_acc_inner(void)
@
@ Halve accumulator (in the field).
@
@ Cost: 21
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_half_acc_inner, %function
gf_half_acc_inner:
	@ If integer is odd, then add q first.
	sbfx	r0, r4, #0, #1
	movw	r1, #(MQ - 1)
	bic	r1, r0, r1
	adds	r4, r1
	adcs	r5, r0
	adcs	r6, r0
	adcs	r7, r0
	adcs	r8, r0
	adcs	r10, r0
	adcs	r11, r0
	adcs	r12, r12, r0, lsr #1

	@ Perform halving. We may have an extra bit in the carry flag;
	@ it will be imported properly with the first rrxs.
	rrxs	r12, r12
	rrxs	r11, r11
	rrxs	r10, r10
	rrxs	r8, r8
	rrxs	r7, r7
	rrxs	r6, r6
	rrxs	r5, r5
	rrxs	r4, r4

	bx	lr
	.size	gf_half_acc_inner, .-gf_half_acc_inner

@ =======================================================================
@ void gf_mul_to_acc_inner(const gf *a, const gf *b)
@
@ Multiplication (in the field). Output goes to the accumulator.
@ WARNING: pointers to inputs a and b are in r1 and r2, respectively.
@
@ Cost: 174 (altentry: 169)
@ =======================================================================

	.align	2    @ 32-bit alignment
	.thumb
	.thumb_func
	.type	gf_mul_to_acc_inner, %function
	nop
gf_mul_to_acc_inner:
	@ IMPLEMENTATION NOTES:
	@ ---------------------
	@
	@ Function principle is about the same as in:
	@   https://github.com/Emill/X25519-Cortex-M4/blob/master/x25519-cortex-m4-gcc.s
	@ The two main differences are:
	@  - We use 2^255-mq with various values of mq (with mq < 2^15),
	@    instead of always 2^255-19.
	@  - For maximum portability, we do not ever touch register r9.
	@ The different value of mq does not impact performance (it just
	@ changes some constants). The abstinence from using r9 forces us
	@ to use slightly more traffic to/from RAM, yielding an overhead
	@ of about 5 cycles. Hence, we have a total cost of 174 cycles,
	@ instead of 170.
	@
	@ A purely quadratic algorithm is applied, with 64 multiplications
	@ (64x64->128). All extra additions are smuggled into the
	@ multiplications by way of using the umaal opcode. We have
	@ in total 4 umull and 60 umaal. We furthermore have to initialize
	@ registers with 0 in 8 places (because umaal needs _two_ words
	@ to add to the product, and there is no opcode that adds a single
	@ word to a 64-bit product and outputs a 64-bit result). The rest
	@ is data traffic.
	@
	@ Input words from the first operand (a) are read once each, in
	@ two chunks of 4 words. Each of these reads is an ldm, and the
	@ values from a are stored in r4..r7. Input words from the second
	@ operand (b) are read twice each; read sequences are 3, 3 and 2
	@ words, for both the first and second passes on b. Extra exchanges
	@ of partial results with the stack are needed to cope with the
	@ register shortage. The pointer to a is reloaded from a stack
	@ slot when reading the second chunk of 4 words; similarly, the
	@ pointer to b must be reloaded from a stack slot for the second
	@ pass.
	@
	@ Once a 512-bit product has been assembled (with 9 words on the
	@ stack and 7 still in registers), modular reduction is applied.
	@ Some pre-folding computations are applied to allow this reduction
	@ to proceed in a single pass. The reduction itself costs only
	@ 28 cycles.
	@
	@ An alternate entry is provided, for the case when one operand
	@ is already in the accumulator. For the alternate entry, the
	@ four low words of the first input (a) should be already loaded
	@ in r4..r7, and register r1 must point to where the four high
	@ words are stored in RAM.

	@ Load a0..a3. We do this first so that:
	@  - The already updated r1 pointer is stored (saves an addition
	@    later on).
	@  - An alternate entry with a0..a3 already loaded can be used.
	ldm	r1!, { r4, r5, r6, r7 }

gf_mul_to_acc_inner_altentry:
	@ For the alternate entry, words a0..a3 MUST have been already
	@ read in r4..r7, and register r1 MUST point to a4.

	@ We need to save some values (9 words) on the stack in addition
	@ to the input pointers, and return address. In the course of the
	@ function, we will also expand the stack by an extra 4 words, but
	@ they will be popped before reaching the end of the function.
	push	{ r1, r2, lr }
	sub	sp, #36

	@ Load b0..b2.
	ldm	r2!, { r10, r11, r12 }

	@ a0*b0
	@ Write out c0 (complete).
	umull	r1, r0, r4, r10
	str	r1, [sp]

	@ a0*b1 and a1*b0
	umull	r3, r8, r4, r11
	umaal	r0, r3, r5, r10

	@ Write out c1 (complete).
	@ NB: we could have delayed the writing out of c0 and merge it
	@ here into a single strd, but it would be slower by 1 cycle
	@ (strd is 3 cycles, but lone str surrounded by non-memory-access
	@ opcodes are just 1 cycle each thanks to the write buffer).
	str	r0, [sp, #4]

	@ a0*b2, a1*b1 and a2*b0
	umaal	r3, r8, r4, r12
	umull	r1, r0, r5, r11
	umaal	r1, r3, r6, r10

	@ Write out c2 (complete).
	str	r1, [sp, #8]

	@ a1*b2, a2*b1 and a3*b0
	umaal	r0, r3, r5, r12
	umaal	r0, r8, r6, r11

	@ a3*b0 is computed independently, addition to r1 is delayed.
	umull	r10, r14, r7, r10

	@ a2*b2 and a3*b1
	umaal	r3, r8, r6, r12
	umaal	r3, r14, r7, r11

	@ a3*b2
	umaal	r8, r14, r7, r12

	@ We have _almost_ computed a0:a1:a2:a3 * b0:b1:b2; the three low
	@ words have been written out to the stack; the four other words
	@ are in r0:r3:r8:r14, except that we must still add r10 to r0.

	@ Load b3..b5.
	ldm	r2!, { r1, r11, r12 }

	@ a0*b3; also add the delayed addition from a3*b0.
	umaal	r0, r10, r4, r1

	@ Write out c3 (complete).
	str	r0, [sp, #12]

	@ a0*b4 and a1*b3
	movs	r0, #0
	umaal	r0, r10, r4, r11
	umaal	r0, r3, r5, r1

	@ a0*b5, a1*b4 and a2*b3
	@ also a3*b3 (in advance)
	umaal	r3, r10, r6, r1
	umaal	r10, r14, r7, r1
	umaal	r3, r8, r5, r11
	movs	r1, #0
	umaal	r1, r3, r4, r12

	@ a1*b5 and a2*b4
	umaal	r10, r3, r5, r12
	umaal	r10, r8, r6, r11

	@ a2*b5 and a3*b4
	umaal	r3, r8, r6, r12
	umaal	r3, r14, r7, r11

	@ a3*b5
	umaal	r8, r14, r7, r12

	@ We have computed a0..a3 * b0..b5; the four low words have
	@ been written out to the stack; the six other words are in
	@ r0:r1:r10:r3:r8:r14.

	@ Read the last two words of b (b6 and b7) into r2:r11.
	ldm	r2, { r2, r11 }

	@ a0*b6
	mov	r12, #0
	umaal	r10, r12, r4, r2

	@ a0*b7 and a1*b6
	@ consume a0
	umaal	r3, r12, r4, r11
	movs	r4, #0
	umaal	r3, r4, r5, r2

	@ a1*b7 and a2*b6
	umaal	r4, r8, r5, r11
	umaal	r4, r12, r6, r2

	@ a2*b7 and a3*b6
	umaal	r8, r14, r6, r11
	umaal	r8, r12, r7, r2

	@ a3*b7
	umaal	r12, r14, r7, r11

	@ We now have computed (a0..a3)*b (12 words); low 4 words are
	@ in the stack, and top 8 words are in r0:r1:r10:r3:r4:r8:r12:r14.

	@ Write c8..c11 (partial) to the stack.
	@ We _expand_ the stack so that we may write them all in a
	@ single push, and read them back with pop.
	push	{ r4, r8, r12, r14 }

	@ Load a4..a7.
	ldr	r4, [sp, #52]
	ldm	r4, { r4, r5, r6, r7 }

	@ Load b0..b2. We also keep the pointer to b3 in r14.
	ldr	r14, [sp, #56]
	ldm	r14!, { r2, r8, r11 }

	@ a4*b0
	mov	r12, #0
	umaal	r0, r12, r4, r2

	@ Store c4 (completed).
	str	r0, [sp, #32]

	@ a4*b1 and a5*b0
	movs	r0, #0
	umaal	r1, r0, r4, r8
	umaal	r1, r12, r5, r2

	@ Store c5 (complete).
	str	r1, [sp, #36]

	@ a4*b2, a5*b1 and a6*b0
	movs	r1, #0
	umaal	r1, r10, r4, r11
	umaal	r1, r12, r5, r8
	umaal	r1, r0, r6, r2

	@ Store c6 (complete).
	str	r1, [sp, #40]

	@ a5*b2, a6*b1 and a7*b0
	umaal	r10, r3, r5, r11
	umaal	r10, r12, r6, r8
	umaal	r10, r0, r7, r2

	@ a6*b2 and a7*b1
	umaal	r3, r12, r6, r11
	umaal	r3, r0, r7, r8

	@ a7*b2
	umaal	r12, r0, r7, r11

	@ We have now processed b0..b2 against a4..a7.
	@  - c6 is in r1 and is complete.
	@  - c7 is in r10 and still needs a4*b3
	@  - c8 is in r3 and needs the stored partial, and upper products
	@  - c9 is in r12 and needs the stored partial, and upper products
	@  - c10 is in r0 and needs the stored partial, and upper products

	@ Load b3, b4 and b5. We still have the pointer to b3 in r14.
	ldm	r14!, { r1, r2, r11 }

	@ a4*b3
	mov	r8, #0
	umaal	r8, r10, r4, r1

	@ Store c7 (complete).
	str	r8, [sp, #44]

	@ Load partial c8 from stack.
	pop	{ r8 }

	@ a4*b4 and a5*b3
	umaal	r8, r10, r4, r2
	umaal	r8, r3, r5, r1

	@ Write out c8 (complete).
	str	r8, [sp, #44]

	@ Load partial c9 from stack.
	pop	{ r8 }

	@ a4*b5, a5*b4 and a6*b3
	umaal	r8, r10, r4, r11
	umaal	r8, r3, r5, r2
	umaal	r8, r12, r6, r1

	@ a5*b5, a6*b4 and a7*b3
	umaal	r10, r3, r5, r11
	umaal	r10, r12, r6, r2
	umaal	r10, r0, r7, r1

	@ a6*b5 and a7*b4
	umaal	r3, r12, r6, r11
	umaal	r3, r0, r7, r2

	@ a7*b5
	umaal	r12, r0, r7, r11

	@ We have completed the product words c0..c9. c9 is still in r8.
	@ We also have some partial values:
	@  - c10: in r10 (needs a4*b6 and partial from stack)
	@  - c11: in r3 (needs a4*b7, a5*b6 and partial from stack)
	@  - c12: in r12 (needs a5*b7, a6*b6)
	@  - c13: in r0 (needs a6*b7 and a7*b6)

	@ Load b6 and b7. r14 contains the pointer to b6 at that point.
	ldm	r14, { r2, r11 }

	@ Load partial c10 and c11 from stack.
	ldrd	r14, r1, [sp], #8

	@ a4*b6
	umaal	r14, r10, r4, r2

	@ a4*b7 and a5*b6
	umaal	r1, r10, r4, r11
	umaal	r1, r3, r5, r2

	@ a5*b7 and a6*b6
	umaal	r10, r3, r5, r11
	umaal	r10, r12, r6, r2

	@ a6*b7 and a7*b6
	umaal	r3, r12, r6, r11
	umaal	r3, r0, r7, r2

	@ a7*b7
	umaal	r12, r0, r7, r11

	@ --------------------------------------------------------
	@ We now have the complete product c0..c15:
	@  - c0..c8 are in the stack
	@  - top words are in r8:r14:r1:r10:r3:r12:r0
	@ We need to apply the reduction, with an output in the
	@ conventional accumulator.

	@ We split c into the following parts:
	@   u = c0..c6    (224 bits)
	@   v = c7        (32 bits)
	@   w = c8..c14   (224 bits)
	@   x = c15       (32 bits)
	@ The value is equal, modulo q, to:
	@   u + 2*mq*w + 2^224*(v + 2*mq*x)
	@ We then write v + 2*mq*x = y + z*2^31, with y < 2^31. The value
	@ is then equal to:
	@   u + 2*mq*w + 2^224*y + mq*z
	@ Since mq < 2^15, we have:
	@   u < 2^224
	@   2*mq*w < 2^240
	@   2^224*y < 2^255
	@   mq*z < 2^32
	@ Therefore, the sum cannot overflow.
	@
	@ (v + 2*mq*x <= (2^32-1)*(2^16-1), thus z <= 131069)

	@ Read c7:c8 from stack.
	ldrd	r2, r7, [sp, #28]

	@ Compute z (in r11) and 2*y (in r2).
	movw	r6, #MQ
	lsls	r5, r6, #2
	mov	r11, r2
	umaal	r2, r11, r5, r0

	@ Compute z*mq (into r11, always fits in 32 bits).
	mul	r11, r11, r6

	@ We now add 2*mq*(c8..c14) to c0..c6, with r11 as input carry.

	@ Set r0 to 2*MQ.
	lsls	r0, r6, #1

	@ Read c0:c1:c2 from stack, and process them.
	pop	{ r4, r5, r6 }
	umaal	r4, r11, r7, r0
	umaal	r5, r11, r8, r0
	umaal	r6, r11, r14, r0

	@ Read c3:c4 from stack, and process them.
	pop	{ r7, r8 }
	umaal	r7, r11, r1, r0
	umaal	r8, r11, r10, r0

	@ Read c5:c6 from stack and process it.
	ldrd	r10, r1, [sp]
	umaal	r10, r11, r3, r0
	umaal	r11, r1, r12, r0

	@ Compute output word d7. It's the carry word from the previous
	@ umaal (in r1) added with y. We still have 2*y in r2.
	add	r12, r1, r2, lsr #1

	add	sp, #24
	pop	{ pc }
	.size	gf_mul_to_acc_inner, .-gf_mul_to_acc_inner

@ =======================================================================
@ void gf_sqr_acc_inner(void)
@
@ Squaring (in the field); input is the accumulator, output is written
@ back into the accumulator.
@
@ Cost: 117
@ =======================================================================

	.align	2   @ 32-bit alignment
	.thumb
	.thumb_func
	.type	gf_sqr_acc_inner, %function
gf_sqr_acc_inner:
	push	{ lr }
	sub	sp, #28

	@ We need an extra register, so we push away a7 here.
	push	{ r12 }

	@ We keep a value 0 in r14; this allows doing 32+32->64 in a
	@ single cycle with umaal, without impacting the carry flag.
	movw	r14, #0

	@ All non-square products are doubled with a sequence of adcs;
	@ this is interleaved with the whole computation, and all other
	@ operations must leave the carry flag untouched. At any time,
	@ we have:
	@  - carry flag: carry for the running doublings
	@  - pre-doubling carries: registers with values to add to
	@    the non-square products (before doubling)
	@  - post-doubling carry: single register with a value to add
	@    to the squares (after doubling)

	@ a0*a0
	@ write c0
	umull	r0, r1, r4, r4
	str	r0, [sp, #4]

	@ a0*a1
	umull	r2, r3, r4, r5
	adds	r2, r2
	umaal	r1, r2, r14, r14

	@ write c1
	str	r1, [sp, #8]

	@ pre: r3
	@ post: r2

	@ a0*a2
	eors	r0, r0
	umaal	r0, r3, r4, r6
	adcs	r12, r0, r0

	@ a1*a1
	umaal	r12, r2, r5, r5

	@ pre: r3
	@ post: r2

	@ a0*a3 and a1*a2
	@ write c2:c3
	umull	r0, r1, r4, r7
	umaal	r0, r3, r5, r6
	adcs.w	r0, r0
	str	r12, [sp, #12]
	umaal	r0, r2, r14, r14
	str	r0, [sp, #16]

	@ pre: r1, r3
	@ post: r2

	@ a0*a4 and a1*a3
	eors	r0, r0
	umaal	r0, r1, r4, r8
	umaal	r0, r3, r5, r7
	adcs	r12, r0, r0

	@ a2*a2
	umaal	r12, r2, r6, r6

	@ write c4
	str	r12, [sp, #20]

	@ pre: r1, r3
	@ post: r2

	@ a0*a5, a1*a4 and a2*a3
	umull	r0, r12, r4, r10
	umaal	r0, r1, r5, r8
	umaal	r0, r3, r6, r7
	adcs	r0, r0
	umaal	r0, r2, r14, r14

	@ write c5
	str	r0, [sp, #24]

	@ pre: r1, r3, r12
	@ post: r2

	@ a0*a6, a1*a5 and a2*a4
	eors	r0, r0
	umaal	r0, r1, r4, r11
	umaal	r0, r3, r5, r10
	umaal	r0, r12, r6, r8
	adcs	r0, r0

	@ a3*a3
	umaal	r0, r2, r7, r7

	@ write c6
	str	r0, [sp, #28]

	@ pre: r1, r3, r12
	@ post: r2

	@ a0*a7, a1*a6, a2*a5 and a3*a4
	eors	r0, r0
	umaal	r0, r1, r5, r11
	umaal	r0, r3, r6, r10
	umaal	r0, r12, r7, r8
	@ We need to concentrate values a bit, because we are out of
	@ registers. We have three "pre" carries in r1, r3 and r12;
	@ we turn them into one "pre" and one "pre2" (carry for the
	@ round after).
	umaal	r1, r3, r14, r14
	umaal	r1, r12, r14, r14
	add	r3, r12
	@ r12 is free, we can load a7 back.
	pop	{ r12 }
	umaal	r0, r14, r4, r12
	adcs	r0, r0
	@ We used our zero register (r14) but now r4 is free.
	eors	r4, r4
	umaal	r0, r2, r4, r4

	@ r0 contains c7

	@ pre: r1, r14
	@ pre2: r3
	@ post: r2

	@ a1*a7, a2*a6 and a3*a5
	@ consume a1
	umaal	r4, r1, r5, r12
	umaal	r14, r4, r6, r11
	eors	r5, r5
	umaal	r14, r5, r7, r10
	adcs	r14, r14

	@ a4*a4
	umaal	r14, r2, r8, r8

	@ r14 contains c8

	@ pre: r0, r1, r3, r5
	@ post: r2

	@ a2*a7, a3*a6 and a4*a5
	@ consume a2
	umaal	r4, r1, r6, r12
	umaal	r4, r3, r7, r11
	umaal	r4, r5, r8, r10

	@ (processing of r4 delayed)

	@ pre: r1, r3, r5
	@ post: r2

	@ a3*a7 and a4*a6
	@ consume a3
	umaal	r1, r3, r7, r12
	umaal	r1, r5, r8, r11

	@ a4*a7
	umaal	r3, r5, r8, r12

	@ (resume processing of r4)
	adcs	r7, r4, r4
	adcs	r8, r1, r1
	@ Here we make r4 our new register with the value zero, for
	@ additions with umaal.
	eors	r4, r4
	umaal	r7, r2, r4, r4

	@ a5*a5
	umaal	r8, r2, r10, r10

	@ r7 contains c9
	@ r8 contains c10

	@ pre: r3, r5
	@ post: r2

	@ Set both r1 and r6 to zero, using r4.
	umull	r1, r6, r4, r4

	@ a5*a6
	umaal	r1, r3, r10, r11
	adcs	r1, r1
	umaal	r1, r2, r4, r4

	@ r1 contains c11

	@ pre: r3, r5
	@ post: r2

	@ a5*a7
	umaal	r5, r3, r10, r12
	adcs	r10, r5, r5

	@ a6*a6
	umaal	r10, r2, r11, r11

	@ r10 contains c12

	@ pre: r3
	@ post: r2

	@ a7*a6
	umaal	r6, r3, r11, r12
	adcs	r6, r6
	umaal	r2, r6, r4, r4

	@ r2 contains c13

	@ pre: r3
	@ post: r6

	@ Continue doubling.
	adcs	r3, r3

	@ a7*a7
	umaal	r3, r6, r12, r12
	adc	r11, r6, #0

	@ r3 contains c14
	@ r11 contains c15

	@ --------------------------------------------------------
	@ We now have the complete product c0..c15:
	@  - c0..c6 are in the stack
	@  - top words are in r0:r14:r7:r8:r1:r10:r2:r3:r11
	@ We need to apply the reduction, with an output in the
	@ conventional accumulator.

	@ We split c into the following parts:
	@   u = c0..c6    (224 bits)
	@   v = c7        (32 bits)
	@   w = c8..c14   (224 bits)
	@   x = c15       (32 bits)
	@ The value is equal, modulo q, to:
	@   u + 2*mq*w + 2^224*(v + 2*mq*x)
	@ We then write v + 2*mq*x = y + z*2^31, with y < 2^31. The value
	@ is then equal to:
	@   u + 2*mq*w + 2^224*y + mq*z
	@ Since mq < 2^15, we have:
	@   u < 2^224
	@   2*mq*w < 2^240
	@   2^224*y < 2^255
	@   mq*z < 2^32
	@ Therefore, the sum cannot overflow.
	@
	@ (v + 2*mq*x <= (2^32-1)*(2^16-1), thus z <= 131069)

	@ Compute z (in r0) and 2*y (in r12).
	movw	r4, #MQ
	lsls	r5, r4, #2
	mov	r12, r0
	umaal	r12, r0, r5, r11

	@ Compute z*mq (into r11, always fits in 32 bits).
	mul	r11, r0, r4

	@ We now add 2*mq*(c8..c14) to c0..c6, with r11 as input carry.

	@ Set r0 to 2*MQ.
	lsls	r0, r4, #1

	@ Read c0:c1:c2 from stack, and process them.
	pop	{ r4, r5, r6 }
	umaal	r4, r11, r14, r0
	umaal	r5, r11, r7, r0
	umaal	r6, r11, r8, r0

	@ Read c3:c4 from stack, and process them.
	pop	{ r7, r8 }
	umaal	r7, r11, r1, r0
	umaal	r8, r11, r10, r0

	@ Read c5:c6 from stack, and process them.
	pop	{ r10, r14 }
	umaal	r10, r11, r2, r0
	umaal	r11, r14, r3, r0

	@ Compute output word d7. It's the carry word from the previous
	@ umaal (in r14) added with y. We still have 2*y in r12.
	add	r12, r14, r12, lsr #1

	@ We removed all the stack contents, we can return right away.
	pop	{ pc }
	.size	gf_sqr_acc_inner, .-gf_sqr_acc_inner

@ =======================================================================
@ void gf_sqr_x_acc_inner(uint32_t num)
@
@ Multiple squaring (in the field). Input and output use the accumulator.
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_sqr_x_acc_inner, %function
gf_sqr_x_acc_inner:
	push	{ lr }
	cmp	r0, #5
	bhi	.Lgf_sqr_x_acc_inner_cont
.Lgf_sqr_x_acc_inner_tail:
	tbb	[pc, r0]
	.byte	13, 11, 9, 7, 5, 3
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	pop	{ pc }

.Lgf_sqr_x_acc_inner_cont:
	subs	r0, #5
	push	{ r0 }
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	bl	gf_sqr_acc_inner
	pop	{ r0 }
	cmp	r0, #5
	bhi	.Lgf_sqr_x_acc_inner_cont
	b	.Lgf_sqr_x_acc_inner_tail
	.size	gf_sqr_x_acc_inner, .-gf_sqr_x_acc_inner

@ =======================================================================
@ void gf_sel3_inner(gf *d, const gf *a, const gf *b, const gf *c,
@                    uint32_t fa, uint32_t fb)
@
@ This function expects fa and fb in registers r4 and r5, respectively.
@
@ All registers are consumed.
@
@ Cost: 73
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_sel3_inner, %function
gf_sel3_inner:
	push	{ r0, r1, r4, lr }

	@ Load GE flags with fb.
	rsbs	r5, r5, #0
	uadd8	r5, r5, r5

	@ Select between b and c.
	ACC_LOAD  r3
	ldm	r2!, { r0, r1, r3, r14 }
	sel	r4, r0, r4
	sel	r5, r1, r5
	sel	r6, r3, r6
	sel	r7, r14, r7
	ldm	r2!, { r0, r1, r3, r14 }
	sel	r8, r0, r8
	sel	r10, r1, r10
	sel	r11, r3, r11
	sel	r12, r14, r12

	@ Load GE flags with fa.
	ldrd	r2, r3, [sp, #4]
	rsbs	r3, r3, #0
	uadd8	r3, r3, r3

	@ Select between a and the current accumulator.
	ldm	r2!, { r0, r1, r3, r14 }
	sel	r4, r0, r4
	sel	r5, r1, r5
	sel	r6, r3, r6
	sel	r7, r14, r7
	ldm	r2!, { r0, r1, r3, r14 }
	sel	r8, r0, r8
	sel	r10, r1, r10
	sel	r11, r3, r11
	sel	r12, r14, r12

	@ Write the result.
	pop	{ r0 }
	ACC_STORE  r0

	add	sp, #8
	pop	{ pc }
	.size	gf_sel3_inner, .-gf_sel3_inner

@ =======================================================================

@ One step of finding the two top words.
@ If r5 == r10 == 0, then:
@     r5 <- r4
@     r4 <- ra
@     r10 <- r8
@     r8 <- rb
@ else:
@     do nothing
@ Input registers ra and rb are not modified.
@ Register r12 is consumed.
@ Cost: 8
.macro APPROX_FIND_STEP  ra, rb
	orr	r12, r5, r10
	cmp	r12, #1
	sbcs	r12, r12
	uadd8	r12, r12, r12
	sel	r5, r4, r5
	sel	r4, \ra, r4
	sel	r10, r8, r10
	sel	r8, \rb, r8
.endm

@ =======================================================================
@ (uint32_t,uint32_t) approximate_ab(const gf *ab)
@
@ This routine approximates the two provided 255-bit nonnegative integers
@ into two 32-bit values (a and b are consecutive in RAM, in that order,
@ starting at address ab; b must be non-zero). If n = max(len(a), len(b)),
@ then, on output:
@
@  - If n <= 32, then r0 contains a[0] and r1 contains b[0]
@  - Otherwise:
@      r0 contains (a mod 2^15) + 2^15 * floor(a / 2^(n-17))
@      r1 contains (b mod 2^15) + 2^15 * floor(b / 2^(n-17))
@
@ All registers are consumed.
@
@ Cost: 85
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	approximate_ab, %function
approximate_ab:
	@ Find the top non-zero words in c = a OR b.
	adds	r0, #16
	add	r1, r0, #32
	ldm	r0!, { r2, r3, r4, r5 }
	ldm	r1!, { r6, r7, r8, r10 }
	APPROX_FIND_STEP  r3, r7
	APPROX_FIND_STEP  r2, r6
	ldrd	r2, r3, [r0, #-24]
	ldrd	r6, r7, [r1, #-24]
	APPROX_FIND_STEP  r3, r7
	APPROX_FIND_STEP  r2, r6
	ldrd	r2, r3, [r0, #-32]
	ldrd	r6, r7, [r1, #-32]
	APPROX_FIND_STEP  r3, r7
	APPROX_FIND_STEP  r2, r6

	@ The top two words for a are r4:r5.
	@ The top two words for b are r8:r10.
	@ r2 contains a[0].
	@ r6 contains b[0].

	@ Shift values to get the top bits.
	orr	r3, r5, r10
	clz	r3, r3
	rsb	r1, r3, #32
	lsl	r0, r5, r3
	lsr	r5, r4, r1
	orrs	r0, r5
	lsl	r10, r10, r3
	lsr	r5, r8, r1
	orr	r1, r10, r5

	@ If len(a) > 32 or len(b) > 32, then we should take the top
	@ 17 bits of r0 and r1, and combine with the low 15 bits of
	@ r2 and r6.
	@ If len(a) <= 32 and len(b) <= 32, then r2 == r0 and r6 == r1
	@ at this point, and we want to use these values; the 15/17
	@ split+assembly combination still yields the expected result.
	bfi	r0, r2, #0, #15
	bfi	r1, r6, #0, #15
	bx	lr
	.size	approximate_ab, .-approximate_ab

@ =======================================================================

@ One step of lin2.
@ Input:
@   r0    pointer to a and b (unmodified)
@   r1    constant 0x10000 (unmodified)
@   r2    update factor f0 + 0x10000 (unmodified)
@   r3    update factor g0 + 0x10000 (unmodified)
@   r4    update factor f1 + 0x10000 (unmodified)
@   r5    update factor g1 + 0x10000 (unmodified)
@   r6    current carry word for a (ignored if i == 0)
@   r7    current carry word for b (ignored if i == 0)
@ Output:
@   r6    new carry word for a
@   r7    new carry word for b
@ Regusters r8, r10, r11, r12 and r14 are consumed.
@ Cost: 18 if i == 0, 20 otherwise
.macro LIN2_STEP  i
	@ The M4 manages to pipeline these two loads, provided that the
	@ whole sequence has a proper alignment. NB: all opcodes in this
	@ macro are 32-bit, which maintains code alignment when multiple
	@ LIN2_STEP are used successively.
	ldr	r11, [r0, #(4 * (\i))]
	ldr	r12, [r0, #(4 * (\i) + 32)]
	.if ((\i) == 0)
	umull	r6, r8, r2, r11
	umull	r7, r10, r4, r11
	.else
	asr	r8, r6, #31
	umlal	r6, r8, r2, r11
	asr	r10, r7, #31
	umlal	r7, r10, r4, r11
	.endif
	umlal	r6, r8, r3, r12
	umlal	r7, r10, r5, r12
	umull	r12, r14, r1, r12
	umlal	r12, r14, r1, r11
	subs	r11, r6, r12
	str	r11, [r0, #(4 * (\i))]
	sbcs	r6, r8, r14
	subs	r12, r7, r12
	str	r12, [r0, #(4 * (\i) + 32)]
	sbcs	r7, r10, r14
.endm

@ =======================================================================
@ (int32_t, int32_t) lin2(gf *ab, int32_t f0, int32_t g0,
@                         int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, a*f0+b*g0 and
@ a*f1+b*g1. f0, g0, f1 and g1 are signed integers, in the -2^15+1..+2^15
@ range. The results are also signed; the two extra words are returned. Values
@ a and b are assumed to be consecutive in RAM.
@
@ ABI:
@   r0 = pointer to ab
@   r2 = f0
@   r3 = g0
@   r4 = f1
@   r5 = g1
@ The two extra words are returned in r6 and r7, respectively.
@ r0 is preserved. All other registers are consumed.
@
@ Cost: 169
@ =======================================================================

	.align	2   @ Enforce 32-bit alignment here, it somehow helps.
	.thumb
	.thumb_func
	.type	lin2, %function
	nop
lin2:
	push	{ lr }
	@ Add 2^16 to all update factors so that they are positive.
	mov	r1, #0x10000
	adds	r2, r1
	adds	r3, r1
	adds	r4, r1
	adds	r5, r1
	LIN2_STEP  0
	LIN2_STEP  1
	LIN2_STEP  2
	LIN2_STEP  3
	LIN2_STEP  4
	LIN2_STEP  5
	LIN2_STEP  6
	LIN2_STEP  7
	pop	{ pc }
	.size	lin2, .-lin2

@ =======================================================================
@ uint32_t s256_lin_div15_abs(gf *ab, int32_t f0, int32_t g0,
@                             int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, (a*f0+b*g0)/2^15 and
@ (a*f1+b*g1)/2^15. The divisions are assumed to be exact (i.e. the 15
@ low bits are dropped. The two 255-bit integers a and b are consecutive
@ in RAM, starting at address 'ab'.
@
@ f0, g0, f1 and g1 are signed integers, in the -2^15+1..+2^15 range.
@ Moreover, it is known that f0+g0 and f1+g1 are also both in the
@ -2^15+1..+2^15 range.
@
@ If the resulting value of a is negative, then it is negated. Similarly,
@ if b is negative, then it is negated. Returned value (in r0) is:
@   0  both values were positive, no negation applied
@   1  a was negative, negation applied on a
@   2  b was negative, negation applied on b
@   3  both values were negative, two negations applied
@
@ ABI:
@   r0 = pointer to ab
@   r2 = f0
@   r3 = g0
@   r4 = f1
@   r5 = g1
@ No register is preserved.
@
@ Cost: 89 + cost(lin2) = 258
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	s256_lin_div15_abs, %function
s256_lin_div15_abs:
	push	{ lr }

	@ Perform linear combination.
	bl	lin2

	@ Load r14 with 0x1FFFF = 2^17 - 1.
	mvn	r14, #0
	lsr	r14, r14, #15

	@ Right-shift value a by 15 bits (signed) and negate it if it
	@ is negative. Each umaal below computes carry + x + (2^17 - 1)*x,
	@ i.e. combines a left shift by 17 bits (with a 64-bit result)
	@ with an addition of a carry word.
	ldm	r0, { r1, r2, r3, r4, r5, r8, r10, r11 }
	asr	r12, r6, #31
	eor	r1, r1, r12
	eor	r2, r2, r12
	eor	r3, r3, r12
	eor	r4, r4, r12
	eor	r5, r5, r12
	eor	r8, r8, r12
	eor	r10, r10, r12
	eor	r11, r11, r12
	eor	r6, r6, r12
	sub	r1, r1, r12
	lsrs	r1, r1, #15
	umaal	r1, r2, r14, r2
	umaal	r2, r3, r14, r3
	umaal	r3, r4, r14, r4
	umaal	r4, r5, r14, r5
	umaal	r5, r8, r14, r8
	umaal	r8, r10, r14, r10
	umaal	r10, r11, r14, r11
	add	r11, r11, r6, lsl #17
	stm	r0!, { r1, r2, r3, r4, r5, r8, r10, r11 }

	@ Right-shift value b by 15 bits (signed) and negate it if it
	@ is negative.
	ldm	r0, { r1, r2, r3, r4, r5, r8, r10, r11 }
	asrs	r6, r7, #31
	eors	r1, r6
	eors	r2, r6
	eors	r3, r6
	eors	r4, r6
	eors	r5, r6
	eor	r8, r8, r6
	eor	r10, r10, r6
	eor	r11, r11, r6
	eor	r7, r7, r6
	sub	r1, r1, r6
	lsrs	r1, r1, #15
	umaal	r1, r2, r14, r2
	umaal	r2, r3, r14, r3
	umaal	r3, r4, r14, r4
	umaal	r4, r5, r14, r5
	umaal	r5, r8, r14, r8
	umaal	r8, r10, r14, r10
	umaal	r10, r11, r14, r11
	add	r11, r11, r7, lsl #17
	stm	r0!, { r1, r2, r3, r4, r5, r8, r10, r11 }

	@ Compute return value.
	@ r12 and r6 contain the original signs of a and b, respectively.
	rsb	r0, r12, #0
	sub	r0, r0, r6, lsl #1

	pop	{ pc }
	.size	s256_lin_div15_abs, .-s256_lin_div15_abs

@ =======================================================================
@ void gf_lin2(gf *ab, int32_t f0, int32_t g0, int32_t f1, int32_t g1)
@
@ This routine replaces a and b with, respectively, a*f0+b*g0 and
@ a*f1+b*g1. Computations are done modulo p.
@
@ f0, g0, f1 and g1 are signed integers, in the -2^15..+2^15 range.
@ Moreover, it is known that f0+g0 and f1+g1 are also both in the
@ -2^15..+2^15 range.
@
@ ABI:
@   r0 = pointer to ab
@   r2 = f0
@   r3 = g0
@   r4 = f1
@   r5 = g1
@ No register is preserved.
@
@ Cost: 73 + cost(lin2) = 242
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_lin2_inner, %function
gf_lin2_inner:
	push	{ lr }

	@ Perform linear combination.
	bl	lin2

	movw	r12, #MQ

	@ Apply modular reduction on value a. Extra word is in r6; since
	@ source values were less than 2^256, and |f0+g0| <= 2^15, we
	@ know that |r6| <= 2^15-1.
	@
	@ We first augment r6 with the top bit; now, |r6| <= 2^16-1, and
	@ the rest of the value is in 0..2^255-1.
	@
	@ Reduction is adding r6*mq. If r6 < 0, we furthermore add q,
	@ i.e. we add 2^255, and subtract mq.
	@
	@ Since |r6| <= 2^16-1, the values r6*mq and r6*mq - mq both fit
	@ on 32 bits (signed).
	ldm	r0, { r1, r2, r3, r4, r5, r8, r10, r11 }
	adds	r11, r11
	adcs	r6, r6
	rrx	r11, r11
	sub	r6, r6, r6, lsr #31
	mul	r6, r6, r12
	adds	r1, r6
	adcs	r2, r2, r6, asr #31
	adcs	r3, r3, r6, asr #31
	adcs	r4, r4, r6, asr #31
	adcs	r5, r5, r6, asr #31
	adcs	r8, r8, r6, asr #31
	adcs	r10, r10, r6, asr #31
	adcs	r11, r11, r6, asr #31
	stm	r0!, { r1, r2, r3, r4, r5, r8, r10, r11 }

	@ Apply modular reduction on value b. Extra word is in r7. This
	@ is the same treatment as above.
	ldm	r0, { r1, r2, r3, r4, r5, r8, r10, r11 }
	adds	r11, r11
	adcs	r7, r7
	rrx	r11, r11
	sub	r7, r7, r7, lsr #31
	mul	r7, r7, r12
	adds	r1, r7
	adcs	r2, r2, r7, asr #31
	adcs	r3, r3, r7, asr #31
	adcs	r4, r4, r7, asr #31
	adcs	r5, r5, r7, asr #31
	adcs	r8, r8, r7, asr #31
	adcs	r10, r10, r7, asr #31
	adcs	r11, r11, r7, asr #31
	stm	r0!, { r1, r2, r3, r4, r5, r8, r10, r11 }

	pop	{ pc }
	.size	gf_lin2_inner, .-gf_lin2_inner

@ =======================================================================

@ One step of binary GCD.
@ Input:
@   ri = xa
@   r1 = xb
@   rk = fg0
@   r3 = fg1
@ Output:
@   rj = xa
@   r1 = xb
@   rl = fg0
@   r3 = fg1
@ ri, rj, rk and rl must be distinct from r1, r3, r4 and r5.
@ Registers ri, rj, rk, rl, r1, r3, r4 and r5 are updated/consumed.
@ Cost: 12
.macro BINGCD_STEP  ri, rj, rk, rl
	@ Set r4 = -lsb(xa).
	sbfx	r4, \ri, #0, #1

	@ Set GE flags to 1 iff lsb(xa) == 1 and xa < xb.
	cmp	\ri, r1
	adcs	r5, r4, #0
	uadd8	r5, r5, r5

	@ Conditional swap.
	sel	\rj, r1, \ri
	sel	r1, \ri, r1
	sel	\rl, r3, \rk
	sel	r3, \rk, r3

	@ Now xa is in rj, and fg0 is in rl.

	@ Conditional subtraction, if xa is odd.
	umlal	\rj, \ri, r4, r1
	umlal	\rl, \ri, r4, r3

	@ Divide xa by 2; multiply f1 and g1 by 2
	lsrs	\rj, \rj, #1
	lsls	r3, r3, #1
.endm

@ Two steps of binary GCD.
@   r0 = xa
@   r1 = xb
@   r2 = fg0
@   r3 = fg1
@ r0, r1, r2 and r3 are updated. r4, r5, r6 and r7 are scratch.
@ Cost: 24
.macro BINGCD_STEP_X2
	BINGCD_STEP  r0, r6, r2, r7
	BINGCD_STEP  r6, r0, r7, r2
.endm

@ One step of binary Kronecker symbol.
@ This macro works only for the first 13 iterations, because it requires
@ exact knowledge of the three least significant bits of a and b.
@ Input:
@   ri = xa
@   r1 = xb
@   rk = fg0
@   r3 = fg1
@   r8 = current symbol value (in lsb)
@ Output:
@   rj = xa
@   r1 = xb
@   rl = fg0
@   r3 = fg1
@   r8 = new symbol value (in lsb)
@ ri, rj, rk and rl must be distinct from r1, r3, r4 and r5.
@ Registers ri, rj, rk, rl, r1, r3, r4, r5 and r8 are updated/consumed.
@ Cost: 17
.macro BINKRONECKER_STEP  ri, rj, rk, rl
	@ Set r4 = -lsb(xa).
	sbfx	r4, \ri, #0, #1

	@ Set GE flags to 1 iff lsb(xa) == 1 and xa < xb.
	cmp	\ri, r1
	adcs	r5, r4, #0
	uadd8	r5, r5, r5

	@ If swapping, then flip the Kronecker symbol if xa and xb are
	@ both equal to 3 modulo 4. This relies on the following facts:
	@  - If swapping, then both a and b are odd.
	@  - It never happens that a < 0 and b < 0 at the same time.
	and	r5, \ri, r5, lsr #4
	ands	r5, r1
	add	r8, r8, r5, lsr #1

	@ Conditional swap.
	sel	\rj, r1, \ri
	sel	r1, \ri, r1
	sel	\rl, r3, \rk
	sel	r3, \rk, r3

	@ Now xa is in rj, and fg0 is in rl.

	@ Conditional subtraction, if xa is odd.
	@ This cannot change the Kronecker symbol, because:
	@  - if b > 0, then (a|b) = (a+k*b|b) for any k \in Z
	@  - if b < 0, then a > 0 and a-b > 0, thus (a|b) = (a-b|b)
	umlal	\rj, \ri, r4, r1
	umlal	\rl, \ri, r4, r3

	@ Divide xa by 2; multiply f1 and g1 by 2
	lsrs	\rj, \rj, #1
	lsls	r3, r3, #1

	@ Kronecker symbol update: flip if and only if b = 3 or 5 mod 8
	@ at this point.
	adds	r4, r1, #2
	add	r8, r8, r4, lsr #2
.endm

@ Two steps of Kronecker symbol.
@ This macro works only for the first 13 iterations, because it requires
@ exact knowledge of the three least significant bits of a and b.
@   r0 = xa
@   r1 = xb
@   r2 = fg0
@   r3 = fg1
@   r8 = current symbol value (in lsb)
@ r0, r1, r2, r3 and r8 are updated. r4, r5, r6 and r7 are scratch.
@ Cost: 34
.macro BINKRONECKER_STEP_X2
	BINKRONECKER_STEP  r0, r6, r2, r7
	BINKRONECKER_STEP  r6, r0, r7, r2
.endm

@ Final step of Kronecker symbol.
@ The low words of the updated a and b values are also provided so that
@ more low bits may be obtained.
@ Input:
@   ri = xa
@   r1 = xb
@   rk = fg0
@   r3 = fg1
@   r8 = current symbol value (in lsb)
@   r12 = current a0
@   r11 = current b0
@ Output:
@   rj = xa
@   r1 = xb
@   rl = fg0
@   r3 = fg1
@   r8 = new symbol value (in lsb)
@   r10 = current a0 (updated)
@   r11 = current b0 (updated)
@ r0, r1, r2, r3, r8, r10 and r11 are updated.
@ r4, r5, r6, r7 and r12 are scratch.
@ Cost: 21
.macro BINKRONECKER_STEP_LAST  ri, rj, rk, rl
	@ Set r4 = -lsb(xa).
	sbfx	r4, \ri, #0, #1

	@ Set GE flags to 1 iff lsb(xa) == 1 and xa < xb.
	cmp	\ri, r1
	adcs	r5, r4, #0
	uadd8	r5, r5, r5

	@ If swapping, then flip the Kronecker symbol if xa and xb are
	@ both equal to 3 modulo 4. This relies on the following facts:
	@  - If swapping, then both a and b are odd.
	@  - It never happens that a < 0 and b < 0 at the same time.
	and	r5, r10, r5, lsr #4
	and	r5, r5, r11
	add	r8, r8, r5, lsr #1

	@ Conditional swap. Also apply it on a0 and b0.
	sel	\rj, r1, \ri
	sel	r1, \ri, r1
	sel	\rl, r3, \rk
	sel	r3, \rk, r3
	sel	r12, r11, r10
	sel	r11, r10, r11

	@ Now:
	@   xa is in rj
	@   fg0 is in rl
	@   a0 is in r12

	@ Conditional subtraction, if xa is odd.
	@ This cannot change the Kronecker symbol, because:
	@  - if b > 0, then (a|b) = (a+k*b|b) for any k \in Z
	@  - if b < 0, then a > 0 and a-b > 0, thus (a|b) = (a-b|b)
	umlal	\rj, \ri, r4, r1
	umlal	\rl, \ri, r4, r3
	umlal	r12, \ri, r4, r11

	@ Divide xa by 2; multiply f1 and g1 by 2
	@ Also divide a0 by 2.
	lsrs	\rj, \rj, #1
	lsr	r10, r12, #1
	lsls	r3, r3, #1

	@ Kronecker symbol update: flip if and only if b = 3 or 5 mod 8
	@ at this point.
	add	r4, r11, #2
	add	r8, r8, r4, lsr #2
.endm

@ =======================================================================
@ bingcd15()
@
@ Inputs:
@   r0 = xa
@   r1 = xb
@
@ Outputs:
@   r1 = final value of xb
@   r2 = f0
@   r3 = f1
@   r4 = g0
@   r5 = g1
@   r6 = final value of xa
@
@ r7, r11 and r14 are consumed.
@
@ Cost: 193
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	bingcd15, %function
bingcd15:
	movs	r2, #1
	lsls	r3, r2, #16

	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP_X2
	BINGCD_STEP  r0, r6, r2, r7

	@ Unpack the update factors.
	mvn	r11, #0x80008000
	add	r2, r7, r11
	add	r3, r11
	uxth	r0, r11
	rsbs	r0, r0, #0
	uxtah	r5, r0, r3, ror #16
	uxtah	r4, r0, r3
	uxtah	r3, r0, r2, ror #16
	uxtah	r2, r0, r2

	bx	lr
	.size	bingcd15, .-bingcd15

@ =======================================================================
@ binkronecker15()
@
@ Inputs:
@   r0 = xa
@   r1 = xb
@   r10 = a0
@   r11 = b0
@
@ Outputs:
@   r2 = f0
@   r3 = f1
@   r4 = g0
@   r5 = g1
@   r8 = Kronecker symbol update (lsb)
@
@ All registers are updated/consumed.
@
@ Cost: 308
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	binkronecker15, %function
binkronecker15:
	movs	r2, #1
	lsls	r3, r2, #16
	eors	r4, r4
	mov	r8, r4

	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP_X2
	BINKRONECKER_STEP  r0, r6, r2, r7

	@ Current state:
	@   r6   xa
	@   r1   xb
	@   r7   fg0
	@   r3   fg1
	@   r8   current symbol value (lsb)
	@   r10  a0 (original)
	@   r11  b0 (original)

	@ Unpack the update factors to compute updated low words.
	mvn	r0, #0x80008000
	movw	r2, #0x7FFF
	rsbs	r2, r2, #0
	
	add	r4, r7, r0
	uxtah	r5, r2, r4, ror #16
	uxtah	r4, r2, r4
	mul	r12, r4, r10
	umlal	r12, r4, r5, r11

	add	r4, r3, r0
	uxtah	r5, r2, r4, ror #16
	uxtah	r4, r2, r4
	mul	r11, r5, r11
	umlal	r11, r4, r4, r10

	lsr	r10, r12, #13
	lsr	r11, r11, #13

	@ r10 = a0 (updated)
	@ r11 = b0 (updated)

	BINKRONECKER_STEP_LAST  r6, r0, r7, r2
	BINKRONECKER_STEP_LAST  r0, r6, r2, r7

	@ Unpack the update factors.
	mvn	r11, #0x80008000
	add	r2, r7, r11
	add	r3, r11
	movw	r6, #0x7FFF
	rsbs	r6, r6, #0
	uxtah	r5, r6, r3, ror #16
	uxtah	r4, r6, r3
	uxtah	r3, r6, r2, ror #16
	uxtah	r2, r6, r2

	bx	lr
	.size	binkronecker15, .-binkronecker15

@ =======================================================================
@ uint32_t gf_inv(gf *d, const gf *y)
@
@ ABI: all registers are consumed.
@
@ Cost: 1497 + cost(gf_normalize_acc_inner) + 33*cost(approximate_ab)
@       + 34*cost(bingcd15) + 33*cost(s256_lin_div15_abs)
@       + 34*cost(gf_lin2_inner) + cost(gf_mul_to_acc_inner)
@       = 28447
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_inv_inner, %function
gf_inv_inner:
	push	{ r0, lr }
	sub	sp, #136

	@ Initial values:
	@ name  index  value
	@   a      0     y (normalized)
	@   b      8     p (modulus)
	@   u     16     1
	@   v     24     0

	@ b <- 2^255-MQ
	add	r0, sp, #32
	movw	r2, #MQ
	rsbs	r2, r2, #0
	asrs	r3, r2, #31
	movs	r4, r3
	movs	r5, r3
	stm	r0!, { r2, r3, r4, r5 }
	movs	r2, r3
	lsrs	r5, r5, #1
	stm	r0!, { r2, r3, r4, r5 }

	@ u <- 1
	@ v <- 0
	@ These initializations are implicit (values u and v are set
	@ after some initial iterations).

	@ v <- 0
	stm	r0!, { r2, r3, r4, r5 }
	stm	r0!, { r2, r3, r4, r5 }

	@ a <- y (normalized)
	ACC_LOAD  r1
	bl	gf_normalize_acc_inner
	ACC_STORE  sp

	@ Perform 33*15 = 495 iterations. Loop counter is kept on
	@ the stack.

	@ We specialize the first two outer iterations.

	@ Approximations of a and b in r0 and r1. Since b is known, we just
	@ take the top bits. Value a is still in the accumulator.
	lsl	r0, r12, #1
	bfi	r0, r4, #0, #15
	movw	r1, #MQ
	rsbs	r1, r1, #0
	@ Run 15 inner iterations of the binary GCD
	bl	bingcd15
	@ Apply update factors on a and b. We save update factors f0 and
	@ g0 (we don't need f1 and g1, since v = 0 at this point).
	mov	r0, sp
	push	{ r2, r4 }
	bl	s256_lin_div15_abs
	@ Restore the update factors. If negation was applied on
	@ a and/or b, then propagate it to the relevant update
	@ factors.
	@ Since u = 1 and v = 0 at this point, the new u and v are:
	@   u' = u*f0 + v*g0 = f0
	@   v' = u*f1 + v*g1 = f1
	@ For now, we simply store these values into the low words of u.
	pop	{ r2, r4 }
	lsls	r6, r0, #31
	movs	r1, #1
	sub	r1, r1, r6, lsr #30
	muls	r2, r1
	str	r2, [sp, #64]
	and	r1, r0, #2
	rsb	r1, r1, #1
	muls	r4, r1
	str	r4, [sp, #68]

	@ Get approximations of a and b in r0 and r1
	mov	r0, sp
	bl	approximate_ab
	@ Run 15 inner iterations of the binary GCD
	bl	bingcd15
	@ Apply update factors on a and b. Take care to save the
	@ update factors.
	mov	r0, sp
	push	{ r2, r3, r4, r5 }
	bl	s256_lin_div15_abs
	@ Restore the update factors. If negation was applied on
	@ a and/or b, then propagate it to the relevant update
	@ factors.
	pop	{ r2, r3, r4, r5 }
	lsls	r6, r0, #31
	movs	r1, #1
	sub	r1, r1, r6, lsr #30
	muls	r2, r1
	muls	r3, r1
	and	r1, r0, #2
	rsb	r1, r1, #1
	muls	r4, r1
	muls	r5, r1
	@ Apply the update factors on u and v. At that point, u and v
	@ are small signed values, currently saved on the stack. We know
	@ that |u| <= 2^15 and |v| <= 2^15; this also applies to |f0|,
	@ |g0|, |f1| and |g1|. Moreover, |f0| and |g0| cannot be both
	@ equal to 2^15. Thus, |u*f0 + v*g0| < 2^31, and we can compute
	@ the new values of u and v over 32-bit signed integers, and only
	@ then convert them to integers modulo p.
	ldrd	r6, r7, [sp, #64]
	muls	r2, r6
	umlal	r2, r1, r3, r7
	muls	r4, r6
	umlal	r4, r1, r5, r7
	@ We now have the new u and v in r2 and r4 (signed), we reduce them
	@ modulo p and write them in the stack buffers.
	add	r0, sp, #64
	movw	r1, #MQ
	asrs	r3, r2, #31
	umlal	r2, r10, r3, r1
	movs	r5, r3
	movs	r6, r3
	movs	r7, r3
	stm	r0!, { r2, r3, r5, r6, r7 }
	lsrs	r7, r7, #1
	stm	r0!, { r5, r6, r7 }
	asrs	r5, r4, #31
	umlal	r4, r10, r5, r1
	movs	r6, r5
	movs	r7, r5
	mov	r8, r5
	stm	r0!, { r4, r5, r6, r7, r8 }
	lsrs	r7, r7, #1
	stm	r0!, { r5, r6, r7 }

	@ The next 31 outer iterations.
	movs	r0, #31
	str	r0, [sp, #128]
.Lgf_inv_inner:
	@ Get approximations of a and b in r0 and r1
	mov	r0, sp
	bl	approximate_ab

	@ Run 15 inner iterations of the binary GCD
	bl	bingcd15

	@ Apply update factors on a and b. Take care to save the
	@ update factors.
	mov	r0, sp
	push	{ r2, r3, r4, r5 }
	bl	s256_lin_div15_abs

	@ Restore the update factors. If negation was applied on
	@ a and/or b, then propagate it to the relevant update
	@ factors.
	pop	{ r2, r3, r4, r5 }
	lsls	r6, r0, #31
	movs	r1, #1
	sub	r1, r1, r6, lsr #30
	muls	r2, r1
	muls	r3, r1
	and	r1, r0, #2
	rsb	r1, r1, #1
	muls	r4, r1
	muls	r5, r1

	@ Apply update factors to u and v.
	add	r0, sp, #64
	bl	gf_lin2_inner

	@ Loop control.
	ldr	r7, [sp, #128]
	subs	r7, #1
	str	r7, [sp, #128]
	bne	.Lgf_inv_inner

	@ We have finished the main loop. 495 inner iterations have
	@ been performed. Values are now at most 15 bits in length;
	@ we can perform the final 15 iterations without any
	@ approximation.
	ldr	r0, [sp]
	ldr	r1, [sp, #32]
	bl	bingcd15

	@ We work modulo a prime p. If the input was zero, then b was
	@ unchanged and is still equal to p. If the input was not zero,
	@ then b contains the GCD, which is 1. Thus, the current value
	@ of xb (in r1) tells us whether the input was zero or not:
	@   input is zero       ->  xb == -mq
	@   input is non-zero   ->  xb == 1
	@ We save the output status on the stack.
	asrs	r1, r1, #31
	adds	r1, #1
	str	r1, [sp]

	@ Final update to u and v.
	@ TODO: we only need to update v here, not u.
	add	r0, sp, #64
	bl	gf_lin2_inner

	@ We have 510 delayed halvings to perform.
	add	r1, sp, #96
	adr	r2, const_gf_inv_invt510
	bl	gf_mul_to_acc_inner
	ldr	r0, [sp, #136]
	ACC_STORE  r0

	@ If the input was zero, then the output has been set to zero,
	@ which is what the API asks for. The returned value has been
	@ computed previously and saved on the stack, we retrieve it
	@ here.
	ldr	r0, [sp]

	add	sp, #140
	pop	{ pc }
	.align	2
const_gf_inv_invt510:
	.long	INVT510
	.size	gf_inv_inner, .-gf_inv_inner

@ =======================================================================
@ int32_t gf_legendre(const gf *y)
@
@ ABI: all registers are consumed.
@
@ Cost: 1269 + cost(gf_normalize_acc_inner)
@       + 33*cost(approximate_ab) + 34*cost(binkronecker15)
@       + 33*cost(s256_lin_div15_abs)
@       = 23217
@ =======================================================================

	.align	1
	.thumb
	.thumb_func
	.type	gf_legendre_inner, %function
gf_legendre_inner:
	push	{ lr }
	sub	sp, #76

	@ Initial values:
	@ name  index  value
	@   a      0     y (normalized)
	@   b      8     p (modulus)

	@ a <- y (normalized)
	ACC_LOAD  r0
	bl	gf_normalize_acc_inner
	ACC_STORE  sp

	@ Check if the input value is zero (needed at the end).
	@ Since it is normalized, this is a simple check that all words
	@ are zero.
	orrs	r4, r5
	orrs	r4, r6
	orrs	r4, r7
	orr	r4, r8
	orr	r4, r10
	orr	r4, r11
	orr	r4, r12
	cmp	r4, #1
	sbcs	r4, r4
	str	r4, [sp, #72]

	@ b <- 2^255-MQ
	add	r0, sp, #32
	movw	r2, #MQ
	rsbs	r2, r2, #0
	asrs	r3, r2, #31
	movs	r4, r3
	movs	r5, r3
	stm	r0!, { r2, r3, r4, r5 }
	movs	r2, r3
	lsrs	r5, r5, #1
	stm	r0!, { r2, r3, r4, r5 }

	@ Initialize Kronecker symbol accumulator.
	eors	r0, r0
	str	r0, [sp, #68]

	@ Perform 33*15 = 495 iterations. Loop counter is kept on
	@ the stack.
	movs	r0, #33
	str	r0, [sp, #64]
.Lgf_legendre_inner:
	@ Get approximations of a and b in r0 and r1
	mov	r0, sp
	bl	approximate_ab

	@ Run 15 inner iterations of the binary Kronecker symbol
	ldr	r10, [sp]
	ldr	r11, [sp, #32]
	bl	binkronecker15
	ldr	r7, [sp, #68]
	add	r7, r8
	str	r7, [sp, #68]

	@ Apply update factors on a and b.
	mov	r0, sp
	bl	s256_lin_div15_abs

	@ We may had to negate a or b (never both, though). Negating
	@ b does not change the result (because, in that case, a > 0).
	@ However, negating a flips the symbol if and only if b = 3 mod 4.
	ldr	r1, [sp, #32]
	and	r0, r0, r1, lsr #1
	ldr	r2, [sp, #68]
	eors	r2, r0
	str	r2, [sp, #68]

	@ Loop control.
	ldr	r7, [sp, #64]
	subs	r7, #1
	str	r7, [sp, #64]
	bne	.Lgf_legendre_inner

	@ We have finished the main loop. 495 inner iterations have
	@ been performed. Values are now at most 15 bits in length;
	@ we can perform the final 15 iterations without any
	@ approximation.
	ldr	r0, [sp]
	ldr	r1, [sp, #32]
	mov	r10, r0
	mov	r11, r1
	bl	binkronecker15

	@ Since the values were exact, there will be no negation on
	@ a and b.

	@ Convert the Legendre symbol to +1/-1 convention.
	ldr	r1, [sp, #68]
	add	r1, r8
	sbfx	r0, r1, #0, #1
	lsls	r0, #1
	adds	r0, #1

	@ Value is correct, unless the source operand was 0, in which
	@ case we need to return 0 and not the value computed above
	@ (as per the API).
	ldr	r1, [sp, #72]
	bics	r0, r1

	add	sp, #76
	pop	{ pc }
	.size	gf_legendre_inner, .-gf_legendre_inner

@ =======================================================================
@ =======================================================================
@ Below are public wrappers for the functions defined above. The wrappers
@ make them callable from C code, by saving all required registers as per
@ the ABI. In practice, most of these wrappers are used only for tests.
@ =======================================================================
@ =======================================================================

@ =======================================================================
@ void CURVE_gf_add(gf *d, const gf *a, const gf *b)
@ =======================================================================

	.align	1
	.global	CN(gf_add)
	.thumb
	.thumb_func
	.type	CN(gf_add), %function
CN(gf_add):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r2
	movs	r0, r1
	bl	gf_add_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_add), .-CN(gf_add)

@ =======================================================================
@ void CURVE_gf_sub(gf *d, const gf *a, const gf *b)
@ =======================================================================

	.align	1
	.global	CN(gf_sub)
	.thumb
	.thumb_func
	.type	CN(gf_sub), %function
CN(gf_sub):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	movs	r0, r2
	bl	gf_sub_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_sub), .-CN(gf_sub)

@ =======================================================================
@ void CURVE_gf_neg(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_neg)
	.thumb
	.thumb_func
	.type	CN(gf_neg), %function
CN(gf_neg):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_neg_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_neg), .-CN(gf_neg)

@ =======================================================================
@ void CURVE_gf_condneg(gf *d, const gf *a, uint32_t ctl)
@ =======================================================================

	.align	1
	.global	CN(gf_condneg)
	.thumb
	.thumb_func
	.type	CN(gf_condneg), %function
CN(gf_condneg):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	movs	r0, r2
	bl	gf_condneg_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_condneg), .-CN(gf_condneg)

@ =======================================================================
@ void CURVE_gf_sub2(gf *d, const gf *a, const gf *b, const gf *c)
@ =======================================================================

	.align	1
	.global	CN(gf_sub2)
	.thumb
	.thumb_func
	.type	CN(gf_sub2), %function
CN(gf_sub2):
	push	{ r0, r3, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	movs	r0, r2
	bl	gf_sub_acc_inner
	ldr	r0, [sp, #4]
	bl	gf_sub_acc_inner
	pop	{ r0, r3 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_sub2), .-CN(gf_sub2)

@ =======================================================================
@ void CURVE_gf_half(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_half)
	.thumb
	.thumb_func
	.type	CN(gf_half), %function
CN(gf_half):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_half_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_half), .-CN(gf_half)

@ =======================================================================
@ void CURVE_gf_sel3(gf *d, const gf *a, const gf *b, const gf *c,
@                    uint32_t fa, uint32_t fb)
@ =======================================================================

	.align	1
	.global	CN(gf_sel3)
	.thumb
	.thumb_func
	.type	CN(gf_sel3), %function
CN(gf_sel3):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ldr	r4, [sp, #32]
	ldr	r5, [sp, #36]
	bl	gf_sel3_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_sel3), .-CN(gf_sel3)

@ =======================================================================
@ void CURVE_gf_mul2(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_mul2)
	.thumb
	.thumb_func
	.type	CN(gf_mul2), %function
CN(gf_mul2):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	ACC_MULCONST  2, r2, r3
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul2), .-CN(gf_mul2)

@ =======================================================================
@ void CURVE_gf_mul4(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_mul4)
	.thumb
	.thumb_func
	.type	CN(gf_mul4), %function
CN(gf_mul4):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	ACC_MULCONST  4, r2, r3
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul4), .-CN(gf_mul4)

@ =======================================================================
@ void CURVE_gf_mul8(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_mul8)
	.thumb
	.thumb_func
	.type	CN(gf_mul8), %function
CN(gf_mul8):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	ACC_MULCONST  8, r2, r3
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul8), .-CN(gf_mul8)

@ =======================================================================
@ void CURVE_gf_mul32(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_mul32)
	.thumb
	.thumb_func
	.type	CN(gf_mul32), %function
CN(gf_mul32):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	ACC_MULCONST  32, r2, r3
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul32), .-CN(gf_mul32)

@ =======================================================================
@ void CURVE_gf_mul_small(gf *d, const gf *a, uint32_t b)
@ =======================================================================

	.align	1
	.global	CN(gf_mul_small)
	.thumb
	.thumb_func
	.type	CN(gf_mul_small), %function
CN(gf_mul_small):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	ACC_MULSMALL  r2, r1, r3
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul_small), .-CN(gf_mul_small)

@ =======================================================================
@ void CURVE_gf_mul(gf *d, const gf *a, const gf *b)
@ =======================================================================

	.align	1
	.global	CN(gf_mul)
	.thumb
	.thumb_func
	.type	CN(gf_mul), %function
CN(gf_mul):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	bl	gf_mul_to_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_mul), .-CN(gf_mul)

@ =======================================================================
@ void CURVE_gf_sqr(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_sqr)
	.thumb
	.thumb_func
	.type	CN(gf_sqr), %function
CN(gf_sqr):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_sqr_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_sqr), .-CN(gf_sqr)

@ =======================================================================
@ void CURVE_gf_sqr_x(gf *d, const gf *a, uint32_t num)
@ =======================================================================

	.align	1
	.global	CN(gf_sqr_x)
	.thumb
	.thumb_func
	.type	CN(gf_sqr_x), %function
CN(gf_sqr_x):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	movs	r0, r2
	bl	gf_sqr_x_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_sqr_x), .-CN(gf_sqr_x)

@ =======================================================================
@ void CURVE_gf_normalize(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_normalize)
	.thumb
	.thumb_func
	.type	CN(gf_normalize), %function
CN(gf_normalize):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_normalize_acc_inner
	pop	{ r0 }
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_normalize), .-CN(gf_normalize)

@ =======================================================================
@ uint32_t CURVE_gf_iszero(const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_iszero)
	.thumb
	.thumb_func
	.type	CN(gf_iszero), %function
CN(gf_iszero):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r0
	bl	gf_iszero_acc_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_iszero), .-CN(gf_iszero)

@ =======================================================================
@ uint32_t CURVE_gf_eq(const gf *a, const gf *b)
@ =======================================================================

	.align	1
	.global	CN(gf_eq)
	.thumb
	.thumb_func
	.type	CN(gf_eq), %function
CN(gf_eq):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_sub_acc_inner
	bl	gf_iszero_acc_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_eq), .-CN(gf_eq)

@ =======================================================================
@ uint32_t CURVE_gf_decode(gf *d, const void *src)
@ =======================================================================

	.align	1
	.global	CN(gf_decode)
	.thumb
	.thumb_func
	.type	CN(gf_decode), %function
CN(gf_decode):
	push	{ r4, r5, r6, r7, r8, r10, r11 }

	@ Load the value (source may be unaligned).
	ACC_LOAD_UNALIGNED  r1

	@ Check that the value is in 0..p-1: subtracting p must yield
	@ a borrow.
	movw	r2, #MQ
	rsbs	r2, r2, #0
	asrs	r3, r2, #31
	subs	r2, r4, r2
	sbcs	r2, r5, r3
	sbcs	r2, r6, r3
	sbcs	r2, r7, r3
	sbcs	r2, r8, r3
	sbcs	r2, r10, r3
	sbcs	r2, r11, r3
	sbcs	r2, r12, r3, lsr #1
	sbcs	r3, r3

	@ Store the decoded value.
	ACC_STORE  r0

	lsrs	r0, r3, #31
	pop	{ r4, r5, r6, r7, r8, r10, r11 }
	bx	lr
	.size	CN(gf_decode), .-CN(gf_decode)

@ =======================================================================
@ uint32_t CURVE_gf_decode_reduce(gf *d, const void *src, size_t len)
@ =======================================================================

	.align	1
	.global	CN(gf_decode_reduce)
	.thumb
	.thumb_func
	.type	CN(gf_decode_reduce), %function
CN(gf_decode_reduce):
	@ If source has length 0, then value is zero.
	cmp	r2, #0
	bne	.Lgf_decode_reduce_cont1
	eors	r1, r1
	umull	r2, r3, r1, r1
	stm	r0!, { r1, r2, r3 }
	stm	r0!, { r1, r2, r3 }
	stm	r0!, { r1, r2 }
	bx	lr

.Lgf_decode_reduce_cont1:
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }

	@ Count number of top bytes (1..32) so that will remain only an
	@ integral number of 32-byte chunks.
	subs	r3, r2, #1
	and	r3, r3, #31
	adds	r3, #1
	subs	r2, r3
	adds	r1, r2

	@ r1 = pointer to the top bytes
	@ r2 = number of non-top bytes (multiple of 32)
	@ r3 = number of top bytes (in 1..32)

	@ Optimize common case of decoding a multiple of 32 bytes.
	cmp	r3, #32
	beq	.Lgf_decode_reduce_cont4

	@ Read r3 bytes into output buffer, and read back value from there.
	eors	r4, r4
	eors	r6, r6
.Lgf_decode_reduce_loop1:
	cmp	r3, r4
	beq	.Lgf_decode_reduce_loop2
	ldrb	r5, [r1, r4]
	strb	r5, [r0, r4]
	adds	r4, #1
	b	.Lgf_decode_reduce_loop1
.Lgf_decode_reduce_loop2:
	cmp	r4, #31
	bhi	.Lgf_decode_reduce_cont3
	strb	r6, [r0, r4]
	adds	r4, #1
	b	.Lgf_decode_reduce_loop2
.Lgf_decode_reduce_cont3:
	ACC_LOAD  r0
	b	.Lgf_decode_reduce_loop5

.Lgf_decode_reduce_cont4:
	ACC_LOAD_UNALIGNED  r1

.Lgf_decode_reduce_loop5:
	@ We now load subsequent chunks. For each input chunk:
	@ 1. multiply the current value by 2^256 = 2*mq mod q
	@ 2. add the new chunk

	@ See if we still have bytes to load.
	cmp	r2, #0
	beq	.Lgf_decode_reduce_cont6
	subs	r2, #32
	subs	r1, #32

	@ Save some values.
	push	{ r0, r1, r2 }

	@ Multiply current value by 2*mq.
	movw	r14, #(2 * MQ)
	ACC_MULSMALL  r14, r2, r3

	@ Write current value in output buffer (we use it as a temporary,
	@ because it's properly aligned).
	ACC_STORE  r0

	@ Load new chunk (possibly unaligned).
	ACC_LOAD_UNALIGNED  r1

	@ Compute addition with the current value (in output buffer,
	@ pointer is already in r0).
	bl	gf_add_acc_inner

	@ Restore saved values and loop.
	pop	{ r0, r1, r2 }
	b	.Lgf_decode_reduce_loop5

.Lgf_decode_reduce_cont6:
	ACC_STORE  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_decode_reduce), .-CN(gf_decode_reduce)

@ =======================================================================
@ void CURVE_gf_encode(void *dst, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_encode)
	.thumb
	.thumb_func
	.type	CN(gf_encode), %function
CN(gf_encode):
	push	{ r0, r4, r5, r6, r7, r8, r10, r11, lr }
	ACC_LOAD  r1
	bl	gf_normalize_acc_inner
	pop	{ r0 }
	ACC_STORE_UNALIGNED  r0
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_encode), .-CN(gf_encode)

@ =======================================================================
@ uint32_t CURVE_gf_inv(gf *d, const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_inv)
	.thumb
	.thumb_func
	.type	CN(gf_inv), %function
CN(gf_inv):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	gf_inv_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_inv), .-CN(gf_inv)

@ =======================================================================
@ int32_t CURVE_gf_legendre(const gf *a)
@ =======================================================================

	.align	1
	.global	CN(gf_legendre)
	.thumb
	.thumb_func
	.type	CN(gf_legendre), %function
CN(gf_legendre):
	push	{ r4, r5, r6, r7, r8, r10, r11, lr }
	bl	gf_legendre_inner
	pop	{ r4, r5, r6, r7, r8, r10, r11, pc }
	.size	CN(gf_legendre), .-CN(gf_legendre)
