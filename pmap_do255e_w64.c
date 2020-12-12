/*
 * This file is meant to be included, not compiled by itself.
 *
 * This file implements CURVE_map_to_curve(); it works with any finite
 * field implementation with 64-bit limbs.
 */

/* see do255.h */
void
do255e_map_to_curve(do255e_point *P, const void *src, size_t len)
{
	/*
	 * The input is first mapped into a field element through
	 * reduction.
	 *
	 * Map definition: we map onto a point on the dual curve
	 * y^2 = x^3 + bb*x with bb = -4*b.
	 *
	 * Given a non-zero field element e, we define:
	 *   x1 = e + (1-bb)/4*e
	 *   x2 = d*(e - (1-bb)/4*e)
	 * where d is a square root of -1. Then:
	 *   (x1^3 + bb*x1)*(x2^3 + bb*x2) = (x1*x2)^3 + bb*(x1*x2)
	 * Therefore, at least one of the following is a square:
	 *   x1^3 + bb*x1
	 *   x2^3 + bb*x2
	 *   (x1*x2)^3 + bb*(x1*x2)
	 * For proper interoperability, we need square roots to be
	 * normalized, and all implementations to compute square roots
	 * on the same values. We therefore write:
	 *
	 *   x1^3 + bb*x1 = yy1num / (y12den^2)
	 *   x2^3 + bb*x2 = yy2num / (y12den^2)
	 *
	 * with:
	 *
	 *   yy1num = 64*e^7 + 16*(3+bb)*e^5
	 *            + 4*(3-2*bb-bb^2)*e^3 + (1-bb)^3*e
         *   yy2num = -d*(64*e^7 - 16*(3+bb)*e^5
	 *                + 4*(3-2*bb-bb^2)*e^3 - (1-bb)^3*e)
	 *   y12den = 8*e^2
	 *
	 * Then:
	 *   - If yy1num is a square, then we set:
	 *        x = x1
	 *        y = sqrt(yy1num) / y12den
	 *   - Otherwise, if yy2num is a square, then we set:
	 *        x = x2
	 *        y = sqrt(yy2num) / y12den
	 *   - Otherwise, we set:
	 *        x = x1*x2
	 *        y = sqrt(yy1num*yy2num) / (y12den^2)
	 *
	 * Once x and y have been chosen, we convert the point to (x,w)
	 * and apply the theta_{1/2} isogeny to get a point on the proper
	 * group on the right curve:
	 *   x' = 4*b/w^2
	 *   w' = (x-b/x)/(2*w)
	 *
	 * If e == 0, then we ignore all of the above, and set the result
	 * to the neutral N.
	 */
	gf e, e2, t1, t2, t3, t4;
	uint64_t qr1, qr2;

	/* Decode the source into a field element. */
	gf_decode_reduce(&e, src, len);

	/* We for now assume that e != 0. */

	gf_sqr(&e2, &e);

	/*
	 * yy1num = 64*e^7 + 176*e^5 - 308*e^3 - 343*e       (into t3)
	 * yy2num = -d*(64*e^7 - 176*e^5 - 308*e^3 + 343*e)  (into t4)
	 */
	gf_mul_small(&t4, &e, 343);
	gf_neg(&t3, &t4);
	gf_mul(&t1, &e, &e2);
	gf_mul_small(&t2, &t1, 308);
	gf_sub(&t3, &t3, &t2);
	gf_sub(&t4, &t4, &t2);
	gf_mul(&t1, &t1, &e2);
	gf_mul_small(&t2, &t1, 176);
	gf_add(&t3, &t3, &t2);
	gf_sub(&t4, &t4, &t2);
	gf_mul(&t1, &t1, &e2);
	gf_mul_small(&t2, &t1, 64);
	gf_add(&t3, &t3, &t2);
	gf_add(&t4, &t4, &t2);
	gf_mul(&t4, &t4, &MINUS_ETA);

	/*
	 * Test QR status of yy1num and yy2num, and set ynum in t1.
	 */
	qr1 = gf_issquare(&t3);
	qr2 = gf_issquare(&t4);
	gf_mul(&t2, &t3, &t4);
	gf_sel3(&t1, &t3, &t4, &t2, qr1, qr2);
	gf_sqrt(&t1, &t1);

	/*
	 * x1num = 4*e^2 - 7
	 * x2num = d*(4*e^2 + 7)
	 * x3num = x1num*x2num
	 * We set xnum into t2.
	 * Also, e2 is set 4*e^2.
	 */
	gf_mul4(&e2, &e2);
	gf_sub(&t3, &e2, &GF_SEVEN);
	gf_add(&t4, &e2, &GF_SEVEN);
	gf_mul(&t4, &t4, &ETA);
	gf_mul(&t2, &t3, &t4);
	gf_sel3(&t2, &t3, &t4, &t2, qr1, qr2);

	/*
	 * Set xden (in t3); it is 4*e if x1 or x2 is chosen, of 16*e^2 if
	 * x1*x2 is chosen.
	 */
	gf_sel2(&t3, &e, &e2, qr1 | qr2);
	gf_mul4(&t3, &t3);

	/*
	 * If x1 or x2 is chosen, then yden = 8*e^2 = (xden^2)/2.
	 * If x1*x2 is chosen, then yden = 64*e^4 = (xden^2)/4.
	 *
	 * We want w = y/x = (ynum*xden) / (xnum*yden). Thus:
	 *
	 *   - If x1 or x2 is chosen, then:
	 *       wnum = 2*ynum
	 *       wden = xnum*xden
	 *   - If x1*x2 is chosen, then:
	 *       wnum = 4*ynum
	 *       wden = xnum*xden
	 *
	 * Thus, we simply double ynum conditionally. Afterwards, we
	 * get the point in Jacobian coordinates:
	 *
	 *    X = xnum^3*xden        (in t2)
	 *    W = 2*ynum or 4*ynum   (in t1)
	 *    Z = xnum*xden          (in t3)
	 */
	gf_mul4(&t4, &t1);
	gf_mul2(&t1, &t1);
	gf_sel2(&t1, &t1, &t4, qr1 | qr2);
	gf_mul(&t3, &t2, &t3);
	gf_sqr(&t2, &t2);
	gf_mul(&t2, &t2, &t3);

	/*
	 * Apply the isogeny theta_{1/2} to get a point in the proper
	 * group on the right curve:
	 *   x' = 4*b/w^2 = -8 / w^2
	 *   w' = (x-bb/x)/(2*w)
	 * In Jacobian coordinates:
	 *   X' = -32*Z^4
	 *   W' = 2*X - W^2
	 *   Z' = 2*W*Z
	 */
	gf_add(&t4, &t1, &t3);
	gf_sqr(&t1, &t1);
	gf_sqr(&t3, &t3);
	gf_sqr(&t4, &t4);
	gf_sub2(&P->Z.w64, &t4, &t1, &t3);
	gf_mul2(&t2, &t2);
	gf_sub(&P->W.w64, &t2, &t1);
	gf_sqr(&t3, &t3);
	gf_mul32(&t3, &t3);
	gf_neg(&P->X.w64, &t3);

	/*
	 * If the starting value was e = 0, then all coordinates are zero
	 * at this point. This is almost correct: we want the neutral
	 * element N = (0:W:0) for a non-zero W. We thus fix W in that
	 * case to get a valid representation of the point N.
	 */
	P->W.w64.v0 |= gf_iszero(&e);
}
