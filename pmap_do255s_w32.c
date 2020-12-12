/*
 * This file is meant to be included, not compiled by itself.
 *
 * This file implements CURVE_map_to_curve(); it works with any finite
 * field implementation with 32-bit limbs.
 */

/* see do255.h */
void
do255s_map_to_curve(do255s_point *P, const void *src, size_t len)
{
	/*
	 * The input is first mapped into a field element through
	 * reduction.
	 *
	 * Map definition: we map onto a point on the dual curve
	 * y^2 = x^3 + aa*x^2 + bb*x with aa = -2*a and bb = a^2-4*b.
	 * The map is Elligator2, from:
	 *   https://elligator.cr.yp.to/elligator-20130828.pdf
	 * (section 5.2)
	 *
	 * Given a non-zero field element e, we define:
	 *   x1 = -aa / (1 + d*e^2)
	 *   x2 = -x1 - aa
	 *      = -aa*d*e^2 / (1 + d*e^2)
	 * where d is a fixed non-quadratic residue (in do255s, we
	 * use d = -1).
	 *
	 * The candidates for y^2 are then yy1num/(yden^2) and
	 * yy2num/(yden^2), with:
	 *   yy1num = -aa*bb*d^3*e^6 + (aa^3-3*aa*bb)*d^2*e^4
	 *            + (aa^3-3*aa*bb)*d*e^2 - aa*bb
	 *   yy2num = -aa*bb*d^4*e^8 + (aa^3-3*aa*bb)*d^3*e^6
	 *            + (aa^3-3*aa*bb)*d^2*e^4 - aa*bb*d*e^2
	 *          = yy1num*d*e^2
	 *   yden = (1 + d*e^2)^2
	 * If yy1num is a square, then we use y = sqrt(yy1num) / yden.
	 * Otherwise, we use y = -sqrt(yy2num) / yden.
	 *
	 * Once x and y have been chosen, we convert the point to (x,w)
	 * and apply the theta_{1/2} isogeny to get a point on the proper
	 * group on the right curve:
	 *   x' = 4*b/w^2
	 *   w' = (x-b/x)/(2*w)
	 *
	 * If e == 1 or -1, then we ignore all of the above, and set the
	 * result to the neutral N. If e == 0, the formulas above also
	 * yield N.
	 */
	gf e, e2, t1, t2, t3, t4;
	uint32_t qr;

	/* Decode the source into a field element. */
	gf_decode_reduce(&e, src, len);

	/* We for now assume that e != 1 or -1. */

	gf_sqr(&e2, &e);

	/*
	 * yy1num = -2*e^6 + 14*e^4 - 14*e^2 + 2      (into t1)
	 * yy2num =  2*e^8 - 14*e^6 + 14*e^4 - 2*e^2  (into t2)
	 * Note that: yy2num = -yy1num*e^2
	 */
	gf_sub(&t1, &GF_SEVEN, &e2);
	gf_mul(&t1, &t1, &e2);
	gf_sub(&t1, &t1, &GF_SEVEN);
	gf_mul(&t1, &t1, &e2);
	gf_add(&t1, &t1, &GF_ONE);
	gf_mul2(&t1, &t1);
	gf_neg(&t2, &t1);
	gf_mul(&t2, &t2, &e2);

	/*
	 * Test QR status of yy1num, and set ynum in t1.
	 * Take care of setting the sign of ynum properly (in Elligator2,
	 * if yy2num is chosen, then we use the opposite of the square root).
	 */
	qr = gf_issquare(&t1);
	gf_sel2(&t1, &t1, &t2, qr);
	gf_sqrt(&t1, &t1);
	gf_condneg(&t1, &t1, 1 - qr);

	/*
	 * x1num = -2
	 * x2num = 2*e^2
	 * Set xnum in t2.
	 */
	gf_neg(&t2, &GF_TWO);
	gf_mul2(&t3, &e2);
	gf_sel2(&t2, &t2, &t3, qr);

	/*
	 * xden = 1 - e^2  (in t3)
	 */
	gf_sub(&t3, &GF_ONE, &e2);

	/*
	 * yden = (1 - e^2)^2 = xden^2
	 *
	 * We want w = y/x = (ynum*xden) / (xnum*yden). Thus:
	 *
	 *   wnum = ynum
	 *   wden = xnum*xden
	 *
	 * We want the point in Jacobian coordinates:
	 *
	 *    X = xnum^3*xden   (in t2)
	 *    W = ynum          (in t1)
	 *    Z = xnum*xden     (in t3)
	 */
	gf_mul(&t3, &t2, &t3);
	gf_sqr(&t2, &t2);
	gf_mul(&t2, &t2, &t3);

	/*
	 * Apply the isogeny theta_{1/2} to get a point in the proper
	 * group on the right curve:
	 *   x' = 4*b/w^2
	 *   w' = (x-bb/x)/(2*w) = (2*x + aa - w^2) / (2*w)
	 * In Jacobian coordinates:
	 *   X' = 8*Z^4
	 *   W' = 2*X + 2*Z^2 - W^2
	 *   Z' = 2*W*Z
	 */
	gf_add(&t4, &t1, &t3);
	gf_sqr(&t1, &t1);
	gf_sqr(&t3, &t3);
	gf_sqr(&t4, &t4);
	gf_sub2(&P->Z.w32, &t4, &t1, &t3);
	gf_add(&t2, &t2, &t3);
	gf_mul2(&t2, &t2);
	gf_sub(&P->W.w32, &t2, &t1);
	gf_sqr(&t3, &t3);
	gf_mul8(&P->X.w32, &t3);

	/*
	 * If the starting value was e = 1 or e = -1, then what we
	 * computed above is (0:0:0); this is the special case of
	 * Elligator2, and we want to map to the neutral N in that case;
	 * we must thus set W to a non-zero value.
	 *
	 * Similarly, if e = 0, then Elligator2, for our target curve,
	 * selects N as output, but the equations above yielded (0:0:0),
	 * which is not a valid representation of N. In that case too,
	 * we must set W to a non-zero value.
	 */
	P->W.w32.v[0] |= gf_iszero(&P->Z.w32);
}
