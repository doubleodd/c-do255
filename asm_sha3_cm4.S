@ =======================================================================
@ Assembly implementation of the SHA3 / SHAKE core for ARM Cortex M4.
@ =======================================================================

	.syntax	unified
	.cpu	cortex-m4
	.file	"asm_sha3_cm4.S"
	.text

@ =======================================================================

.macro INVERT_WORDS
	@ Invert A[1] and A[2].
	adds	r1, r0, #8
	ldm	r1, { r2, r3, r4, r5 }
	mvns	r2, r2
	mvns	r3, r3
	mvns	r4, r4
	mvns	r5, r5
	stm	r1!, { r2, r3, r4, r5 }
	@ Invert A[8]
	ldrd	r2, r3, [r0, #64]
	mvns	r2, r2
	str	r2, [r0, #64]
	mvns	r3, r3
	str	r3, [r0, #68]
	@ Invert A[12]
	ldrd	r2, r3, [r0, #96]
	mvns	r2, r2
	str	r2, [r0, #96]
	mvns	r3, r3
	str	r3, [r0, #100]
	@ Invert A[17]
	ldrd	r2, r3, [r0, #136]
	mvns	r2, r2
	str	r2, [r0, #136]
	mvns	r3, r3
	str	r3, [r0, #140]
	@ Invert A[20]
	ldrd	r2, r3, [r0, #160]
	mvns	r2, r2
	str	r2, [r0, #160]
	mvns	r3, r3
	str	r3, [r0, #164]
.endm

.macro KHI_LOAD  s0, s1, s2, s3, s4
	ldrd	r12, r14, [sp, #(32 + 8 * (\s0))]
	ldrd	r3, r4, [sp, #(32 + 8 * (\s1))]
	ldrd	r5, r6, [sp, #(32 + 8 * (\s2))]
	ldrd	r7, r8, [sp, #(32 + 8 * (\s3))]
	ldrd	r10, r11, [sp, #(32 + 8 * (\s4))]
.endm

.macro KHI_STEP  op, x0, x1, x2, x3, x4, x5, d
	\op	r1, \x0, \x2
	eors	r1, \x4
	str	r1, [r0, #(8 * (\d))]
	\op	r2, \x1, \x3
	eors	r2, \x5
	str	r2, [r0, #(8 * (\d) + 4)]
	@strd	r1, r2, [r0, #(8 * (\d))]
.endm

.macro XOR_ROL_1  s1, s2
	ldm	r0!, { r5, r6, r7, r8 }
	eors	r3, r5, r10
	eors	r4, r6, r11
	.if (\s1) > 32
	mov	r14, #(1 << ((\s1) - 32))
	umull	r5, r6, r4, r14
	umlal	r6, r5, r3, r14
	.else
	mov	r14, #(1 << (\s1))
	umull	r5, r6, r3, r14
	umlal	r6, r5, r4, r14
	.endif
	eors	r3, r7, r1
	eors	r4, r8, r2
	.if (\s2) > 32
	mov	r14, #(1 << ((\s2) - 32))
	umull	r7, r8, r4, r14
	umlal	r8, r7, r3, r14
	.else
	mov	r14, #(1 << (\s2))
	umull	r7, r8, r3, r14
	umlal	r8, r7, r4, r14
	.endif
	stm	r12!, { r5, r6, r7, r8 }
.endm

.macro XOR_ROL_2  s1, s2
	ldm	r0!, { r5, r6, r7, r8 }
	eors	r10, r5, r1
	eors	r11, r6, r2
	.if (\s1) > 32
	mov	r14, #(1 << ((\s1) - 32))
	umull	r5, r6, r11, r14
	umlal	r6, r5, r10, r14
	.else
	mov	r14, #(1 << (\s1))
	umull	r5, r6, r10, r14
	umlal	r6, r5, r11, r14
	.endif
	eors	r10, r7, r3
	eors	r11, r8, r4
	.if (\s2) > 32
	mov	r14, #(1 << ((\s2) - 32))
	umull	r7, r8, r11, r14
	umlal	r8, r7, r10, r14
	.else
	mov	r14, #(1 << (\s2))
	umull	r7, r8, r10, r14
	umlal	r8, r7, r11, r14
	.endif
	stm	r12!, { r5, r6, r7, r8 }
.endm

.macro XOR_ROL_3  s
	eors	r3, r5, r1
	eors	r4, r6, r2
	.if (\s) > 32
	mov	r14, #(1 << ((\s) - 32))
	umull	r6, r5, r3, r14
	umlal	r5, r6, r4, r14
	.else
	mov	r14, #(1 << (\s))
	umull	r5, r6, r3, r14
	umlal	r6, r5, r4, r14
	.endif
	stm	r12!, { r5, r6 }
.endm

@ =======================================================================
@ void do255_sha3_process_block(shake_context *sc)
@ =======================================================================

	.align	2
	.global	do255_sha3_process_block
	.thumb
	.thumb_func
	.type	do255_sha3_process_block, %function
do255_sha3_process_block:
	push	{ r1, r2, r3, r4, r5, r6, r7, r8, r10, r11, r12, lr }
	sub	sp, #232

	INVERT_WORDS

	@ Do 24 rounds. Each loop iteration performs one rounds. We
	@ keep eight times the current round counter in [sp+28] (i.e.
	@ a multiple of 8, from 0 to 184).
	eors	r1, r1
	str	r1, [sp, #28]
.Ldo255_sha3_process_block_loop:

	@ xor(A[5*i+0]) -> r1:r2
	@ xor(A[5*i+1]) -> r3:r4
	@ xor(A[5*i+2]) -> r5:r6
	@ xor(A[5*i+3]) -> r7:r8
	@ xor(A[5*i+4]) -> r10:r11
	ldm	r0!, { r1, r2, r3, r4, r5, r6, r7, r8 }
	adds	r0, #8
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r1, r10
	eors	r2, r11
	eors	r3, r12
	eors	r4, r14
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r5, r10
	eors	r6, r11
	eors	r7, r12
	eors	r8, r14
	adds	r0, #8
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r1, r10
	eors	r2, r11
	eors	r3, r12
	eors	r4, r14
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r5, r10
	eors	r6, r11
	eors	r7, r12
	eors	r8, r14
	adds	r0, #8
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r1, r10
	eors	r2, r11
	eors	r3, r12
	eors	r4, r14
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r5, r10
	eors	r6, r11
	eors	r7, r12
	eors	r8, r14
	adds	r0, #8
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r1, r10
	eors	r2, r11
	eors	r3, r12
	eors	r4, r14
	ldm	r0!, { r10, r11, r12, r14 }
	eors	r5, r10
	eors	r6, r11
	eors	r7, r12
	eors	r8, r14
	ldm	r0!, { r10, r11 }
	subs	r0, #200
	ldrd	r12, r14, [r0, #32]
	eors	r10, r12
	eors	r11, r14
	ldrd	r12, r14, [r0, #72]
	eors	r10, r12
	eors	r11, r14
	ldrd	r12, r14, [r0, #112]
	eors	r10, r12
	eors	r11, r14
	ldrd	r12, r14, [r0, #152]
	eors	r10, r12
	eors	r11, r14

	@ t0 = xor(A[5*i+4]) ^ rotl1(xor(A[5*i+1])) -> r10:r11
	@ t1 = xor(A[5*i+0]) ^ rotl1(xor(A[5*i+2])) -> r1:r2
	@ t2 = xor(A[5*i+1]) ^ rotl1(xor(A[5*i+3])) -> r3:r4
	@ t3 = xor(A[5*i+2]) ^ rotl1(xor(A[5*i+4])) -> r5:r6
	@ t4 = xor(A[5*i+3]) ^ rotl1(xor(A[5*i+0])) -> r7:r8
	mov	r12, r10
	mov	r14, r11
	eors	r10, r10, r3, lsl #1
	eors	r10, r10, r4, lsr #31
	eors	r11, r11, r4, lsl #1
	eors	r11, r11, r3, lsr #31
	eors	r3, r3, r7, lsl #1
	eors	r3, r3, r8, lsr #31
	eors	r4, r4, r8, lsl #1
	eors	r4, r4, r7, lsr #31
	eors	r7, r7, r1, lsl #1
	eors	r7, r7, r2, lsr #31
	eors	r8, r8, r2, lsl #1
	eors	r8, r8, r1, lsr #31
	eors	r1, r1, r5, lsl #1
	eors	r1, r1, r6, lsr #31
	eors	r2, r2, r6, lsl #1
	eors	r2, r2, r5, lsr #31
	eors	r5, r5, r12, lsl #1
	eors	r6, r6, r12, lsr #31
	eors	r5, r5, r14, lsr #31
	eors	r6, r6, r14, lsl #1

	@ Save t2, t3 and t4 on the stack.
	stm	sp, { r3, r4, r5, r6, r7, r8 }

	@ We XOR one of the t0..t4 values into each A[] word, and
	@ rotate the result by some amount (each word has its own
	@ amount). The results are written back into a stack buffer
	@ that starts at sp+32
	add	r12, sp, #32

	@ XOR t0 into A[5*i+0] and t1 into A[5*i+1]; each A[i] is also
	@ rotated left by some amount.

	@ A[0] and A[1]
	ldm	r0!, { r5, r6, r7, r8 }
	eors	r5, r10
	eors	r6, r11
	eors	r3, r7, r1
	eors	r4, r8, r2
	mov	r14, #(1 << 1)
	umull	r7, r8, r3, r14
	umlal	r8, r7, r4, r14
	stm	r12!, { r5, r6, r7, r8 }

	@ A[5] and A[6]
	adds	r0, #24
	XOR_ROL_1  36, 44

	@ A[10] and A[11]
	adds	r0, #24
	XOR_ROL_1  3, 10

	@ A[15] and A[16]
	adds	r0, #24
	XOR_ROL_1  41, 45

	@ A[20] and A[21]
	adds	r0, #24
	XOR_ROL_1  18, 2

	@ XOR t2 into A[5*i+2] and t3 into A[5*i+3]; each A[i] is also
	@ rotated left by some amount. We reload t2 into r1:r2 and t3
	@ into r3:r4.
	ldm	sp, { r1, r2, r3, r4 }

	@ A[2] and A[3]
	subs	r0, #160
	XOR_ROL_2  62, 28

	@ A[7] and A[8]
	adds	r0, #24
	XOR_ROL_2  6, 55

	@ A[12] and A[13]
	adds	r0, #24
	XOR_ROL_2  43, 25

	@ A[17] and A[18]
	adds	r0, #24
	XOR_ROL_2  15, 21

	@ A[22] and A[23]
	adds	r0, #24
	XOR_ROL_2  61, 56

	@ XOR t4 into A[5*i+4]; each A[i] is also rotated left by some
	@ amount. We reload t4 into r1:r2.
	ldrd	r1, r2, [sp, #16]

	@ r0 currently points to A[24]; reset it to point to A[0].
	subs	r0, #192

	@ A[4]
	ldrd	r5, r6, [r0, #32]
	XOR_ROL_3  27

	@ A[9]
	ldrd	r5, r6, [r0, #72]
	XOR_ROL_3  20

	@ A[14]
	ldrd	r5, r6, [r0, #112]
	XOR_ROL_3  39

	@ A[19]
	ldrd	r5, r6, [r0, #152]
	XOR_ROL_3  8

	@ A[24]
	ldrd	r5, r6, [r0, #192]
	XOR_ROL_3  14

	@ At that point, the stack buffer at sp+32 contains the words
	@ at the following indexes (0 to 24) and offsets (from sp)
	@   A[ 0]    0      32
	@   A[ 1]    1      40
	@   A[ 2]   10     112
	@   A[ 3]   11     120
	@   A[ 4]   20     192
	@   A[ 5]    2      48
	@   A[ 6]    3      56
	@   A[ 7]   12     128
	@   A[ 8]   13     136
	@   A[ 9]   21     200
	@   A[10]    4      64
	@   A[11]    5      72
	@   A[12]   14     144
	@   A[13]   15     152
	@   A[14]   22     208
	@   A[15]    6      80
	@   A[16]    7      88
	@   A[17]   16     160
	@   A[18]   17     168
	@   A[19]   23     216
	@   A[20]    8      96
	@   A[21]    9     104
	@   A[22]   18     176
	@   A[23]   19     184
	@   A[24]   24     224

	@ A[0], A[6], A[12], A[18] and A[24]
	KHI_LOAD  0, 3, 14, 17, 24
	KHI_STEP  orrs, r3, r4, r5, r6, r12, r14, 0
	KHI_STEP  orns, r7, r8, r5, r6, r3, r4, 1
	KHI_STEP  ands, r7, r8, r10, r11, r5, r6, 2
	KHI_STEP  orrs, r12, r14, r10, r11, r7, r8, 3
	KHI_STEP  ands, r12, r14, r3, r4, r10, r11, 4

	@ A[3], A[9], A[10], A[16] and A[22]
	KHI_LOAD  11, 21, 4, 7, 18
	KHI_STEP  orrs, r3, r4, r5, r6, r12, r14, 5
	KHI_STEP  ands, r7, r8, r5, r6, r3, r4, 6
	KHI_STEP  orns, r7, r8, r10, r11, r5, r6, 7
	KHI_STEP  orrs, r12, r14, r10, r11, r7, r8, 8
	KHI_STEP  ands, r12, r14, r3, r4, r10, r11, 9

	@ A[1], A[7], A[13], A[19] and A[20]
	KHI_LOAD  1, 12, 15, 23, 8
	KHI_STEP  orrs, r3, r4, r5, r6, r12, r14, 10
	KHI_STEP  ands, r7, r8, r5, r6, r3, r4, 11
	KHI_STEP  bics, r10, r11, r7, r8, r5, r6, 12
	mvns	r7, r7
	mvns	r8, r8
	KHI_STEP  orrs, r12, r14, r10, r11, r7, r8, 13
	KHI_STEP  ands, r12, r14, r3, r4, r10, r11, 14

	@ A[4], A[5], A[11], A[17] and A[23]
	KHI_LOAD  20, 2, 5, 16, 19
	KHI_STEP  ands, r3, r4, r5, r6, r12, r14, 15
	KHI_STEP  orrs, r7, r8, r5, r6, r3, r4, 16
	KHI_STEP  orns, r10, r11, r7, r8, r5, r6, 17
	mvns	r7, r7
	mvns	r8, r8
	KHI_STEP  ands, r12, r14, r10, r11, r7, r8, 18
	KHI_STEP  orrs, r12, r14, r3, r4, r10, r11, 19

	@ A[2], A[8], A[14], A[15] and A[21]
	KHI_LOAD  10, 13, 22, 6, 9
	KHI_STEP  bics, r5, r6, r3, r4, r12, r14, 20
	KHI_STEP  ands, r12, r14, r3, r4, r10, r11, 24
	mvns	r3, r3
	mvns	r4, r4
	KHI_STEP  orrs, r7, r8, r5, r6, r3, r4, 21
	KHI_STEP  ands, r7, r8, r10, r11, r5, r6, 22
	KHI_STEP  orrs, r12, r14, r10, r11, r7, r8, 23

	@ Get round counter XOR round constant into A[0]
	ldr	r1, [sp, #28]
	adr	r2, .const_do255_sha3_process_block_RC
	adds	r2, r1
	ldrd	r3, r4, [r2]
	ldrd	r5, r6, [r0]
	eors	r5, r3
	str	r5, [r0]
	eors	r6, r4
	str	r6, [r0, #4]

	@ Increment round counter, loop until all 24 rounds are done.
	adds	r1, #8
	str	r1, [sp, #28]
	cmp	r1, #192
	blo	.Ldo255_sha3_process_block_loop

	INVERT_WORDS

	add	sp, #232
	pop	{ r1, r2, r3, r4, r5, r6, r7, r8, r10, r11, r12, pc }

.const_do255_sha3_process_block_RC:
	.word	0x00000001
	.word	0x00000000
	.word	0x00008082
	.word	0x00000000
	.word	0x0000808A
	.word	0x80000000
	.word	0x80008000
	.word	0x80000000
	.word	0x0000808B
	.word	0x00000000
	.word	0x80000001
	.word	0x00000000
	.word	0x80008081
	.word	0x80000000
	.word	0x00008009
	.word	0x80000000
	.word	0x0000008A
	.word	0x00000000
	.word	0x00000088
	.word	0x00000000
	.word	0x80008009
	.word	0x00000000
	.word	0x8000000A
	.word	0x00000000
	.word	0x8000808B
	.word	0x00000000
	.word	0x0000008B
	.word	0x80000000
	.word	0x00008089
	.word	0x80000000
	.word	0x00008003
	.word	0x80000000
	.word	0x00008002
	.word	0x80000000
	.word	0x00000080
	.word	0x80000000
	.word	0x0000800A
	.word	0x00000000
	.word	0x8000000A
	.word	0x80000000
	.word	0x80008081
	.word	0x80000000
	.word	0x00008080
	.word	0x80000000
	.word	0x80000001
	.word	0x00000000
	.word	0x80008008
	.word	0x80000000

	.size	do255_sha3_process_block, .-do255_sha3_process_block
