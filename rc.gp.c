/******************************************************************************/
/**													Real and complex functions											**/
/******************************************************************************/


// Checked up to 5000 with the 64-bit version
GEN
Bell(long n)
{
	static const long BellMod32[] = {
		1,1,2,5,15,20,11,13,12,27,7,10,29,21,26,17,3,12,15,1,12,15,27,26,9,25,18,13,
		7,4,3,5,12,19,31,10,5,13,10,25,27,28,7,25,12,7,19,26,17,17,2,21,31,20,27,29,
		12,11,23,10,13,5,26,1,19,12,31,17,12,31,11,26,25,9,18,29,23,4,19,21,12,3,15,
		10,21,29,10,9,11,28,23,9,12,23,3,26
	};
	if (n < 2) {
		if (n < 0)
			pari_err(talker, "negative argument to Bell");
		return gen_1;
	}
	pari_sp ltop = avma;
	GEN B, Br, f, t;
	
	double logn = log(n);
	double w = W_small(n);
	long sz = n * (logn - log(w) - 1 + 1/w) + logn;	// - log(w)/2 - 1;
	sz /= BITS_IN_LONG * log(2);
	sz += 5;	// Guard digits

	// This while re-does the calculation if the precision was too low.
	while (1) {
		long k = 1;
		f = real_1(sz);
		B = real_0(sz);
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		
		// This is the hot loop where all calculations are done.
		while (1) {
			f = divrs(f, k);
			t = mulir(powis(stoi(k), n), f);
			if (expo(t) < -25)
				break;
			B = addrr(B, t);
			++k;
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 2, &f, &B);
		}
		B = mulrr(B, gexp(gen_m1, sz));
		if (nbits2prec(expo(B)+1) > lg(B))
			goto FIX_PRECISION;
		
		Br = ceilr(B); // calculate answer
		
		// Check if answer is within a millionth of an integer
		if (cmprr(absr(mpsub(B, Br)), real2n(-20, DEFAULTPREC)) < 0) {
			Br = gerepilecopy(ltop, Br);
			
			// Sanity check on the floating-point arithmetic
			if (mod32(Br) != BellMod32[n % 96])
				pari_err(talker, "Bell(%d) failed verification mod 32, please report", n);
			if (DEBUGLEVEL > 4) {
				pari_printf("Precision for Bell(%d): %d (%d digits), an excess of %d (%d digits)\n", n, sz - 2, prec2ndec(sz), sz - expi(Br) / BITS_IN_LONG - 3, prec2ndec(sz - expi(Br) / BITS_IN_LONG - 1));
			}
			return Br;
		}
pari_printf("Not accurate enough: error %f\n", rtodbl(absr(mpsub(B, Br))));
		
FIX_PRECISION:
		sz += logf(sz);	// increase precion slightly
		pari_warn(warnprec, "Bell", sz);
		avma = ltop;
	}
}


GEN
deBruijnXi(GEN x)
{
	double xx = rtodbl(x), left, right;
	if (xx < 1)
		pari_err(talker, "deBruijnXi: Can't find a xi given x < 1.");
	if (xx > 1)
		left = log(xx);
	else
		left = DBL_EPSILON;
	right = 1.35 * log(xx) + 1;	// Heuristic

	// Bisection
	while (right - left > left * DBL_EPSILON) {
		double m = (left + right) / 2;
		if (expm1(m) > xx * m)
			right = m;
		else
			left = m;
	}
	return dbltor((left + right) / 2);
}

GEN
rhoest(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN xi, ret;
	x = gtor(x, "rhoest", prec);
	xi = deBruijnXi(x);
	ret = gexp(gsub(gneg(veceint1(negr(xi), NULL, prec)), gmul(x, xi)), prec);
	ret = gdiv(gdiv(ret, sqrtr(mulrr(mulsr(2, mppi(prec)), x))), xi);
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
DickmanRho(GEN x, long prec)
{
	static const double rhoTable[] = {
		NEVER_USED, 1, 3.068528194e-1, 4.860838829e-2, 4.910925648e-3,
		3.547247005e-4, 1.964969635e-5, 8.745669953e-7, 3.232069304e-8,
		1.016248283e-9,	2.770171838e-11, 6.644809070e-13, 1.419713165e-14,
		2.729189030e-16, 4.760639989e-18, 7.589908004e-20
	};
	static const int rhoTableLen = 15;
	static const double rhoScale = 1.130709295873035782;
	// Last table entry, divided by rhoest at that point

	pari_sp ltop = avma;
	GEN ret, left, right, scale;
	x = gtor(x, "DickmanRho", prec);
	if (cmprs(x, 2) <= 0) {
		ret = gsubsg(1, glog(gmaxgs(x, 1), prec));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	if (gcmpgs(x, 3) <= 0) {
		ret = gadd(gadd(gsubsg(1, mulrr(subsr(1, mplog(subrs(x, 1))), mplog(x))), greal(dilog(subsr(1, x), prec))), divrs(sqrr(mppi(prec)), 12));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
  
	double xx = rtodbl(x);
  
	// Asymptotic estimate (scaled for continuity)
	if (xx > rhoTableLen) {
		double scale = rhoScale;
		scale = (scale - 1) * sqrt(sqrt(rhoTableLen / xx)) + 1;
		/* Let the scale factor dwindle away, since the estimate is (presumably) */
		/* better in the long run than any scaled version of it.  The exponent */
		/* of 0.25 has been chosen to give the best results for 10 < x < 100 */
		/* with a table size of 10. */

		ret = precision0(mulrr(rhoest(x, prec), dbltor(scale)), 9);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
  
	// Scaling factors: the factor by which the true value of rho differs from
	// the estimates at the endpoints.
	left = divrr(dbltor(rhoTable[(int)floor(xx)]), rhoest(floorr(x), prec));
	right = divrr(dbltor(rhoTable[(int)ceil(xx)]), rhoest(ceilr(x), prec));
	
	// Linear interpolation on the scale factors.
	scale = gadd(left, gmul(gsub(right, left), mpsub(x, floorr(x))));
	
	// Return a result based on the scale factor and the asymptotic formula.
	ret = precision0(gmul(rhoest(x, prec), scale), 9);
	ret = gerepileupto(ltop, ret);
	return ret;
}


// Convenience function: binary logarithm of x
GEN
log_2(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN ret = NEVER_USED;	// to silence compiler, which doesn't know that pari_err never returns
	switch(typ(x)) {
		case t_INT:
			ret = mplog(itor(x, prec));
			break;
		case t_REAL:
			ret = mplog(x);
			break;
		case t_FRAC:
		case t_COMPLEX:
			ret = mplog(cxcompotor(x, prec));
			break;
		default:
			pari_err(typeer, "lg");
	}
	ret = divrr(ret, mplog2(prec));
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
contfracback(GEN v, GEN terms)
{
	pari_sp ltop = avma;
	GEN x = gen_0;
	long tterms = NEVER_USED;
	if (!terms)
		tterms = glength(v) - 1;
	else if (typ(terms) == t_INT)
		tterms = itos(terms);
	else
		pari_err(typeer, "contfracback");
	x = gel(v, tterms + 1);
	if (tterms == 1)
		x = gcopy(x);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i = 0;
	for (i = tterms; i >= 1; i--)
	{
		x = gadd(gel(v, i), ginv(x));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &x);
	}
	x = gerepileupto(ltop, x);
	return x;
}


double
W_small(double x)
{
	if (x <= 0)
	{
		if (!x)
			return 0.0;
		if (x == -exp(-1))
			return -1.0;	// otherwise, sometimes sqrt becomes complex?
		if (x < -exp(-1))
			pari_err(talker, "out of range");
	}
	double e, w, t = 1.0;
	
	// Initial approximation for iteration
	if (x < 1)
	{
		double tmp = sqrt(2 * exp(1) * x + 2);
		w = (1 - (11/72*tmp + 1/3)*tmp)*tmp - 1;
	} else {
		w = log(x);
	}
	if (x > 3)
		w -= log(w);

	double ep = DBL_EPSILON * (1 + fabs(w));
	while (fabs(t) > ep)
	{
		// Halley loop
		e = exp(w);
		t = w*e - x;
		t /= e*(w+1) - 0.5*(w+2)*t/(w+1);
		w -= t;
	}
	return w;
}


GEN
W(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN e, t, w, ep, tmp;
	
	if (typ(x) == t_INT)
		x = itor(x, prec);
	else if (typ(x) != t_REAL)
		pari_err(typeer, "W");
	prec = precision(x);
	
	if (signe(x) <= 0)
	{
		if (!signe(x))
		{
			avma = ltop;
			return real_0(prec);
		}
		pari_sp btop = avma;
		GEN oneOverE = mpexp(stor(-1, prec));
		long c = absr_cmp(x, oneOverE);
		avma = btop;
		if (!c)
			return real_m1(prec);	// otherwise, sometimes sqrt becomes complex
		if (c > 0)	// x < -1/e
			pari_err(talker, "out of range");
		// TODO: Parabolic iteration rather than Halley in this case?
	}
	t = real_1(prec);
	
	// Initial approximation for iteration
	if (cmprs(x, 1) < 0)
	{
		// t = 1 already, might as well use it for this calculation
		tmp = sqrtr(addrs(mulrr(shiftr(mpexp(t), 1), x), 2));
		w = subrs(mulrr(tmp, subir(gen_1, mulrr(tmp, addrr(invr(stor(3, prec)), divrs(mulrs(tmp, 11), 72))))), 1);
		// tmp = sqrt(2e * x + 2)
		// w = tmp * (1 - tmp * (1/3 + 11/72 * tmp)) - 1
		if (precision(w) < prec) {
			w = rtor(w, prec);
		}
	} else {
		w = mplog(x);
	}
	if (cmprs(x, 3) > 0)
		w = subrr(w, mplog(w));

	ep = mulrr(eps(prec), addsr(1, absr(w)));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);

	while (cmprr(absr(t), ep) > 0)
	{
		// Halley loop
		e = mpexp(w);
		t = subrr(mulrr(w, e), x);
//pari_printf("	%Ps (off by %Ps)\n", w, t);
		if (cmprr(absr(t), ep) <= 0)
			break;	// Stops calculation when answer is already very close, to avoid division by 0
		tmp = addrs(w, 1);
		t = divrr(t, subrr(mulrr(e, tmp), divrr(mulrr(shiftr(addrs(w, 2), -1), t), tmp)));
		w = subrr(w, t);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &w);
	}

	w = gerepileupto(ltop, w);
	if (precision(w) < prec)
		pari_warn(warner, "precision loss");
	return w;
}

/******************************************************************************/
/**															 Statistics																 **/
/******************************************************************************/

long
infinite(GEN x)
{
	if (typ(x) != t_VEC || glength(x) != 1)
		return 0;
	GEN e = gel(x, 1);	// Nothing is created, so no garbage... right?
	
	// If e is gen_0, then is_pm1 is unpredicatable, but that doesn't matter
	// because then both branches return 0.
	return (typ(e) == t_INT && is_pm1(e)) ? signe(e) : 0;
}


long
isExtendedReal(GEN x)
{
	long t = typ(x);
	if (t == t_INT || t == t_FRAC || t == t_REAL)
		return 1;
	return infinite(x);
}


// FIXME: Infinities are broken?
GEN
normd(GEN a, GEN b, long prec)
{
	if (!isExtendedReal(a) || !isExtendedReal(b))
		pari_err(talker, "incorrect endpoint in normd");
	pari_sp ltop = avma;
	long tmp;
	GEN ret = NEVER_USED;
	
	/* Infinities */
	if ((tmp = infinite(a)))	// Assignment and test-if-0
	{
pari_warn(warner, "Doesn't work properly with infinities");
		if (tmp < 0)	// (-oo, b)
		{
			tmp = infinite(b);
			if (tmp > 0)
				ret = gen_1;
			else if (tmp < 0)
				ret = gen_0;
			else
				ret = gdivgs(mpneg(gerfc(mpdiv(b, gsqrt(gen_2, prec)), prec)), 2);
		} else {		// (oo, b)
			if (infinite(b) == 1)
				ret = gen_1;
			else
				pari_err(talker, "incorrect endpoint in normd");
		}
	} else if ((tmp = infinite(b))) {	// Assignment and test-if-0
pari_warn(warner, "Doesn't work properly with infinities");
		if (tmp < 0)	// (a, -oo)
			pari_err(talker, "incorrect endpoint in normd");
		ret = gdivgs(gerfc(mpdiv(a, gsqrt(gen_2, prec)), prec), 2);
	} else {
		GEN root2 = gsqrt(gen_2, prec);
		ret = gdivgs(gsub(gerfc(gdiv(a, root2), prec), gerfc(gdiv(b, root2), prec)), 2);
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


// Use the Box-Muller transform to generate random normal variables. Caches
// values, so multiple calls at the same precision are fast.
GEN
rnormal(long prec)
{
	if (rnormal_cached) {
		if (precision(rnormal_cached) != prec) {
			gunclone(rnormal_cached);
		} else {
			GEN ret = gcopy(rnormal_cached);
			gunclone(rnormal_cached);
			rnormal_cached = 0;
			return ret;
		}
	}
	pari_sp ltop = avma;
	GEN u1, u2, ret, outside, inside, cos_inside;
	u1 = randomr(prec);
	u2 = randomr(prec);
	outside = sqrtr_abs(shiftr(mplog(u1), 1));
	inside = mulrr(shiftr(mppi(prec), 1), u2);
	cos_inside = mpcos(inside);
	
	ret = mulrr(outside, cos_inside);
	rnormal_cached = gclone(ret);	// Cache for later use
	ret = mulrr(outside, cos_inside);
	ret = gerepileupto(ltop, ret);
		return ret;
}
