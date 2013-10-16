/******************************************************************************/
/**												 Other number theory															**/
/******************************************************************************/

long
Collatz(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("Collatz",n);
	if (signe(n) < 1)
		pari_err_DOMAIN("Collatz", "n", "<", gen_1, n);
	pari_sp ltop = avma;

	long v = vali(n);
	if (v) n = shifti(n, -v);	// Possible stack garbage
	ulong nn = itou_or_0(n);
#ifdef LONG_IS_64BIT
	if (nn && nn < 23035537407UL)
#else
	if (nn && nn < 159487UL)
#endif
		return Collatz_tiny(nn);
	
	long iterations = 0;
#ifdef LONG_IS_64BIT
	while (cmpiu(n, 23035537407UL) >= 0) {
#else
	while (cmpiu(n, 159487UL) >= 0) {
#endif
		iterations++;
		n = addii(n, addis(shifti(n, -1), 1));
		v = vali(n);
		if (v) n = shifti(n, -v);
		// Not sure if it's worthwhile to collect garbage in this loop.
	}
	
	avma = ltop;
	return iterations + Collatz_tiny(itou_or_0(n));
}

long
Collatz_tiny(ulong n)
{
	long iterations = 0;
	//n >>= vals(n);	// It is the caller's responsibility to pass an odd number
	while (n > 1) {
		iterations++;
		n += (n>>1)+1;
		n >>= vals(n);
	}
	return iterations;
}


long
isfactorial(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("isfactorial",n);
	if (signe(n) < 1)
		return 0;
	if(mod2(n))
		return isint1(n);
	
	// Remove small factorials with a combined binary-linear search
	ulong nn = itou_or_0(n);
	if (nn) {
		return
#ifdef LONG_IS_64BIT
		nn < 6227020800 ? (
#endif
		nn < 40320 ?
			nn==2 || nn==6 || nn==24 || nn==120 || nn==720 || nn==5040
		:
			nn==40320 || nn==362880 || nn==3628800 || nn==39916800 || nn==479001600
#ifdef LONG_IS_64BIT
		) : nn < 355687428096000 ?
			nn==6227020800 || nn==87178291200 || nn==1307674368000 || nn==20922789888000
		:
			nn==355687428096000 || nn==6402373705728000 || nn==121645100408832000 || nn==2432902008176640000
#endif
		;
	}
	
	// 
	pari_sp ltop = avma;
	long v2 = vali(n);	// 2-adic valuation of n: 2^k || n.
#ifdef LONG_IS_64BIT
	if (v2 < 18)
#else
	if (v2 < 10)
#endif
		return 0;

	long mn = v2 + 1, mx = mn + __builtin_clzl(v2), t, c;
	
	// Find the 
	while (mx - mn > 1) {
		t = mn + (mx - mn)/2;
		c = factorial_lval(t, 2);
		if (c < v2)
			mn = t+1;
		else if (c > v2)
			mx = t-1;
		else {
			mx = t|1;
			mn = mx==mn ? mx : mx-1;
		}
	}
	if (mn < mx) {
		long p = itos(lpf(stoi(mx)));	// Doesn't really need to be the smallest, just some prime factor
		avma = ltop;
		t = Z_lval(n, p);
		c = factorial_lval(mx, p);
		if (t > c)
			return 0;
		
		if (t == c)
			mn = mx;
		else if (t == c - u_lval(mx, p))
			mx = mn;
		else
			return 0;
	}
	//c = equalii(mpfact(mn), n);
	c = equalii(mulu_interval(2UL, (ulong)mn), n);
	avma = ltop;
	return c;
}


GEN
Faulhaber(long e, GEN a)
{
	pari_sp ltop = avma;
	GEN ret = gen_0;
	GEN x = pol_x(0);
	if (!a)
		a = x;
	else if (!gcmpX(a))
		pari_err_TYPE("Faulhaber", a);
	pari_sp btop = avma;
	long i = 0;
	for (; i <= e; ++i)
	{
		// TODO: Form the polynomial in a sensible way
		ret = gadd(ret, gmul(gmul(binomialuu(e + 1, i), bernfrac(i)), gpowgs(x, e + 1 - i)));
		ret = gerepileupto(btop, ret);
	}
	ret = gsubstpol(gadd(gdivgs(ret, e + 1), gpowgs(x, e)), x, a);
	ret = gerepileupto(ltop, ret);
	return ret;
}


ulong
cuberoot(ulong n)
{
	int i;
	ulong m = 0, b;

	// Another possibility, if needed to eke out more speed, would be to use
	// __builtin_clzl and a static array to find the first bit.  This
	// would help a little for 32-bit calculations and a lot for 64-bit.
#ifdef LONG_IS_64BIT
	for (i = 63; i; i -= 3) {
#else
	for (i = 30; i; i -= 3) {
#endif
		b = 3*m*(m + 1) | 1;
		if ((n >> i) >= b) {
			n -= b << i;
			m |= 1;
		}
		m += m;
	}
	b = 3*m*(m + 1) | 1;
	return n < b ? m : m | 1;
}


GEN
cuberootint(GEN x)
{
	pari_sp ltop = avma;
	long t = typ(x);
	GEN ret = NEVER_USED;
	if (t == t_INT) {
		ulong n = itou_or_0(x);
		if (n)
			return utoipos(cuberoot(n));
		ret = gfloor(powrfrac(itor(x, lgefint(x)), 1, 3));
	} else if (t == t_REAL) {
		ret = gfloor(powrfrac(x, 1, 3));
		x = gfloor(x);
	} else {
		pari_err_TYPE("cuberootint", x);
	}

	if (cmpii(powis(ret, 3), x) == 1) {
		ret = subis(ret, 1);
		while (cmpii(powis(ret, 3), x) == 1)
			ret = subis(ret, 1);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	while (cmpii(powis(addis(ret, 1), 3), x) < 1)
		ret = addis(ret, 1);
	ret = gerepileupto(ltop, ret);
	return ret;
}


long
issquarefree_small(ulong n)
{
#define CUTOFF 1627UL
	long tmp = n&3;
	if (!tmp)
		return 0;
	if (tmp == 2)
		n >>= 1;
	
	long p;
    forprime_t primepointer;
    u_forprime_init(&primepointer, 3, CUTOFF);

	// First loop: remove tiny primes, don't calculate cube roots
	// 99.7% of non-squarefree numbers are detected by this loop (or the above)
    while ((p = u_forprime_next(&primepointer)) < 98)
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
		}
	}
	
	long last = (long)cuberoot(n);
	long last1 = (long)minuu(last, CUTOFF);

	// Beyond this point, 99.89% of numbers are squarefree.
    while ((p = u_forprime_next(&primepointer)))
	{
		if (p > last1)
			break;
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = (long)cuberoot(n);
			last1 = (long)minuu(last, CUTOFF);
		}
	}

#ifdef LONG_IS_64BIT
	if (n < CUTOFF * CUTOFF * CUTOFF)
#elif CUTOFF <= 1621
	pari_warn(warner, "issquarefree: cutoff marginal, performace suffers");
	if (n < CUTOFF * CUTOFF * CUTOFF)
#endif
	return n == 1 || !uissquare(n);

    u_forprime_init(&primepointer, p+2, last);

	// n is at least CUTOFF^3 and is not divisible by any prime under CUTOFF
	if (last < 65536) {	// maxprime() > 65536
		while ((p = u_forprime_next(&primepointer)))
		{
			if (n%p == 0) {
				n /= p;
				if (n%p == 0)
					return 0;
				last = (long)cuberoot(n);
			}
			if (p > last)
				return !uissquare(n);
		}
	}
	
	// n is at least 49 bits
	pari_sp ltop = avma;
	GEN f = pollardbrent(utoipos(n));
	if (f == NULL) {
		// Do nothing
	} else if (typ(f) == t_INT) {
		ulong ff = itos(f);
		n /= ff;
		if (!issquarefree_small(ff) || ugcd(ff, n % ff) > 1) {
			avma = ltop;
			return 0;
		}
		last = cuberoot(n);
	} else {	// typ(f) == t_VEC
		long l = lg(f);
		int i = 1;
		while (i < l) {
			/*
			 // Exponent from pollardbrent is guaranteed to equal 1
			if (gel(f[i], i+1) != gen_1) {
				avma = ltop;
				return 0;
			}*/
			ulong ff = itos(gel(f[i], i));
			n /= ff;
			if (!issquarefree_small(ff) || ugcd(ff, n % ff) > 1) {
				avma = ltop;
				return 0;
			}
			i += 3;	// [factor, exponent, class, factor, exponent, class, ...]
		}
		last = cuberoot(n);
	}
	avma = ltop;
	
	// cP(2^80)
	// 750000	72.81
	// 700000	72.26
	// 600000	72.76
	// 500000	72.50 (was: 71.86)
	// 400000	72.41
	// 300000	72.06
	// 250000	72.33
	// 100000	72.51
	if (last > 300000 || (ulong)last > maxprime())	// TODO: Find good breakover point here
	{
		long ret = Z_issquarefree(stoi(n));
		avma = ltop;
		return ret;
	}
	while ((p = u_forprime_next(&primepointer)))
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = (long)cuberoot(n);
		}
		
		if (p > last)
			break;
	}
	
	return !uissquare(n);
#undef CUTOFF
}


ulong
ucountPowerfulu(ulong n)
{
	if (n < 4)	/* avoid log(0) */
		return n > 0;
	int i, crossover = (int)pow(n, .14);
	ulong res = 0, k, breakpoint = cuberoot(n / (crossover * crossover));
/*
pari_printf("Counting squarefrees until %d, values ranging from %lu to %lu;\n",
crossover, breakpoint, cuberoot(n));
pari_printf("Counting directly to %lu.\n", breakpoint);
*/
	for (i = 1; i <= crossover; i++)
		res += ucountSquarefree(cuberoot(n / (i * i)));
	res -= crossover * ucountSquarefree(breakpoint);
	for (k = 1; k <= breakpoint; k++)
		if (issquarefree_small(k))
			res += usqrt(n / k) / k;
	return res;
}


// Assumption: n > 1 and cuberoot(n) fits into a word.
// Output: the *least significant word* of the number of powerful numbers up to
// n.  If there is overflow, let the calling function handle it if the actual
// number (not mod 2^32 or 2^64) is desired.
// TODO: Use a sieve to determine if k is squarefree.
ulong
ucountPowerfuli(GEN n)
{
#define uCBRTis(n,k) itou(cuberootint(divii((n), mulss(k,k))))
	pari_sp ltop = avma;
	int i, crossover = (int)exp(dbllog2r(itor(n, MEDDEFAULTPREC)) * 0.09);
	ulong k, res = 0;
	ulong breakpoint = uCBRTis(n, crossover);
pari_printf("Counting squarefrees until %d, values ranging from %lu to %lu, ",crossover, breakpoint, itou(cuberootint(n)));
pari_printf("then counting directly to %lu.\n", breakpoint);
	avma = ltop;
pari_timer T;
timer_start(&T);
	for (i = 1; i <= crossover; i++) {
		res += ucountSquarefree(uCBRTis(n, i));
		avma = ltop;
	}
	res -= crossover * ucountSquarefree(breakpoint);
GEN est = divrr(czeta(gdiv(stoi(3), gen_2), MEDDEFAULTPREC),czeta(stoi(3), MEDDEFAULTPREC));
est = gmul(est, sqrtr(itor(n, MEDDEFAULTPREC)));
GEN correction = divrr(czeta(gdiv(gen_2, stoi(3)), MEDDEFAULTPREC),czeta(gen_2, MEDDEFAULTPREC));
correction = gmul(correction, powrfrac(itor(n, MEDDEFAULTPREC), 1, 3));
est = gadd(est, correction);
est = gerepileupto(ltop, est);
err_printf("Estimate: %Ps\n", gfloor(est));
err_printf("Squarefree time:   %6ld ms; partial sum = %lu (%.2f%% of estimated total)\n",
timer_delay(&T), res, 100.0 * res / rtodbl(est));
// zeta(3/2)/zeta(3)*sqrt(n)+zeta(2/3)/zeta(2)*n^(1/3)
timer_start(&T);
	if (breakpoint > LONG_MAX)
		pari_err_OVERFLOW("Breakpoint much too large, something funny is going on.");
	for (k = 1; k <= breakpoint; k++) {
		if (issquarefree_small(k)) {
			res += itos(divis(sqrti(divis(n, k)), k));
			avma = ltop;
		}
	}
err_printf("Direct count time: %6ld ms\n", timer_delay(&T));
	return res;
#undef uCBRTis
}


// Varies as zeta(3/2)/zeta(3) n^1/2 + zeta(2/3)/zeta(2) n^1/3 + o(n^1/6)
// estPowerful(n)=zeta(3/2)/zeta(3)*sqrt(n) + zeta(2/3)/zeta(2)*n^(1/3)
GEN
countPowerful(GEN n)
{
	if (signe(n) < 1)
		return gen_0;
	if (typ(n) != t_INT && typ(n) != t_REAL)
		pari_err_TYPE("countPowerful", n);
	
	GEN ret;
	pari_sp ltop = avma;
	long sz = lgefint(n = gfloor(n));
	if (sz <= 4) {
		if (sz <= 3) {
			ret = utoipos(ucountPowerfulu(itou(n)));
			ret = gerepileupto(ltop, ret);
			return ret;
		}
		// TODO: use bounds on countSquarefree to use this function even
		// when the return value exceeds wordsize.  Count using
		// ucountPowerfuli, allowing it to overflow, and calculate the
		// high-order bits through the asymptotic formula.  Good bounds are
		// needed to make this work, but these can be found.
#ifdef LONG_IS_64BIT
		// FIXME: Find appropriate breakpoint here.
		// around 72047453657149422928610771848621443808
		// Slightly-conservative estimate below.
		else if (1)	//(cmpii(n, mkintn(4, 909366712, 745454542, 1225794890, 2766209793)) < 0)
#else
		// This is the last number such that the result fits in a 32-bit ulong.
		else if (cmpii(n, uu32toi(910359010, 4006567939)) <= 0)
#endif
		{
			ret = utoipos(ucountPowerfuli(n));
			ret = gerepileupto(ltop, ret);
			return ret;
		}
	}
	
	GEN breakpoint = cuberootint(shifti(n, -2));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k = gen_0;
	ret = subii(countSquarefree(cuberootint(n)), countSquarefree(breakpoint));
	for (k = gen_1; cmpii(k, breakpoint) <= 0; k = addis(k, 1))
	{
		if (Z_issquarefree(k))
			ret = addii(ret, sqrti(gdivent(n, powis(k, 3))));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &ret, &k);
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


// TODO: Doesn't really save much time vs. the original.  To improve, a sieve
// would be needed (to calculate the values of moebius faster).
ulong
ucountSquarefree(ulong lim)
{
	ulong b = usqrt(lim >> 1);
	ulong k;
	ulong ret = 0;
	for (k = 1; k <= b; k++)
		ret += moebiusu(k) * (lim / (k * k));
	ulong p3 = usqrt(lim);
	for (k = b + 1; k <= p3; k++)
		ret += moebiusu(k);
	return ret;
}


GEN
countSquarefree(GEN lim)
{
	pari_sp ltop = avma;

	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err_TYPE("countSquarefree", lim);

	GEN b, ret = gen_0, p3 = gen_0, p4 = gen_0;
	b = sqrti(shifti(lim, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k = gen_0;
	for (k = gen_1; cmpii(k, b) <= 0; k = addis(k, 1))
	{
		// Clean up operations
		ret = addii(ret, mulis(gdivent(lim, gsqr(k)), moebius(k)));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &ret, &k);
	}
	p3 = sqrti(lim);
	btop = avma;
	st_lim = stack_lim(btop, 1);
	p4 = gen_0;
	for (k = addis(b, 1); cmpii(k, p3) <= 0; k = addis(k, 1))
	{
		p4 = addis(p4, moebius(k));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &p4, &k);
	}
	ret = addii(ret, p4);
	ret = gerepileupto(ltop, ret);
	return ret;
}

/******************************************************************************/
/**															 Factoring																	**/
/******************************************************************************/

GEN
Mfactor(GEN p, GEN lim, GEN start)
{
	pari_sp ltop = avma;

	// Check types
	if (typ(p) != t_INT)
		pari_err_TYPE("Mfactor", p);
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err_TYPE("Mfactor", lim);
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err_TYPE("Mfactor", start);

	GEN v, k, p1, p2;
	v = cgetg(1, t_VEC);
	if (signe(p) < 1 || !(isprime(p)))
		pari_err_PRIME("Mfactor", p);
	if (mod4(p) != 3)
		pari_warn(warner, "p must be a Mersenne exponent equal to 3 mod 4... I think");

	/* Really, only k in [0, 5, 8, 9] mod 12 need to be checked (at last for Mersenne prime exponents?). */
	/* So this should use a quick check, then a full loop with step [10p, 6p, 2p, 6p]. */
	/* Check for a factor of 3 also before starting loop. */
	p2 = shifti(p, 1);
	k = gceil(gdiv(subis(start, 1), p2));
	p1 = addis(mulii(p2, k), 1);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN q = gen_0;
	for (q = p1; cmpii(q, lim) <= 0; q = addii(q, p2))
	{
		if (mod8(q) != 1 && mod8(q) != 7)
			continue;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &q, &v);
		if (!gequal1(powgi(gmodulsg(2, q), modii(p, subis(q, 1)))))
			continue;
		v = concat(v, q);
		long i = 2;
		while (gequal1(powgi(gmodulsg(2, powis(q, i)), modii(p, mulii(subis(q, 1), powis(q, i - 1))))))
		{
			v = concat(v, q);
			i++;
		}
	}
	v = gerepileupto(ltop, v);
	return v;
}


// FIXME: Gives spurious factors of 2 sometimes
GEN
bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start)
{
	pari_sp ltop = avma;

	GEN v = cgetg(1, t_VEC);
	GEN p1 = gen_0, p2 = gen_0;
	if (typ(a) != t_INT)
		pari_err_TYPE("bigfactor", a);
	if (typ(b) != t_INT)
		pari_err_TYPE("bigfactor", b);
	if (typ(c) != t_INT)
		pari_err_TYPE("bigfactor", c);
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err_TYPE("bigfactor", start);
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err_TYPE("bigfactor", lim);
	ulong lm = itou_or_0(lim);
	if (lm > maxprime())
		pari_err_MAXPRIME(lm);
	if (!lm)
		pari_err(e_MISC, "bigfactor: Not enough primes on your architechture for that!");

	if (signe(b) < 0)
	{
		// TODO: These two should have their formats changed to match that given below.
		if (equali1(a))
		{
			p1 = Z_factor(subsi(1, c));
			p1 = gerepileupto(ltop, p1);
			return p1;
		}
		if (equalim1(a))
		{
			p2 = Z_factor(subii(mod2(b) ? gen_m1 : gen_1, c));
			p2 = gerepileupto(ltop, p2);
			return p2;
		}
		/* a^b not in Z */
		pari_err_OP("bigfactor", a, b);
	}
	ulong p3 = minuu(itou(a), lm);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN p5 = gen_0;		/* int */

	long p;
    forprime_t primepointer;
    u_forprime_init(&primepointer, 3, lm);
    while ((ulong)(p = u_forprime_next(&primepointer)) <= p3)
	{
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &v);
		if (cmpis(gcdii(a, stoi(p)), 1) > 0)
			p5 = b;	// What's the right way to handle this case?	It seems that something more efficient could be done.
		else
			p5 = stoi(smodis(b, p - 1));
		if (!gequal(powgi(gmodulo(a, stoi(p)), p5), c))
			continue;
		v = concat(v, stoi(p));
		long i = 2;
		pari_sp ctop = avma, c_lim = stack_lim(btop, 1);
		GEN p6 = gen_0;
		for(;;)
		{
			if (cmpis(gcdii(a, stoi(p)), 1) > 0)
				p6 = b;
			else
				p6 = modii(b, mulis(powuu(p, i - 1), p - 1));
			if (!gequal(powgi(gmodulo(a, powuu(p, i)), p6), c))
				break;
			v = concat(v, stoi(p));
			i++;
			if (low_stack(c_lim, stack_lim(ctop, 1)))
				v = gerepileupto(ctop, v);
		}
	}
	v = gerepileupto(btop, v);
	btop = avma;
	st_lim = stack_lim(btop, 1);

	// Second loop -- most of the work is done here
    while ((p = u_forprime_next(&primepointer)))
	{
		if (cmpsi(p, lim) > 0)
			break;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &v);
		if (!gequal(gpowgs(gmodulo(a, stoi(p)), smodis(b, p - 1)), c))
			continue;
		v = concat(v, stoi(p));
		long i = 2;
		while (gequal(powgi(gmodulo(a, powis(stoi(p), i)), modii(b, mulis(powis(stoi(p), i - 1), p - 1))), c))
		{
			v = concat(v, stoi(p));
			i++;
		}
	}
	v = gerepileupto(ltop, v);
	return v;
}


// Does d divide a^b - c?
long
bigdiv(GEN a, GEN b, GEN c, GEN d)
{
	pari_sp ltop = avma;
	long ret;
	if (typ(a) != t_INT)
		pari_err_TYPE("bigdiv", a);
	if (typ(b) != t_INT)
		pari_err_TYPE("bigdiv", b);
	if (typ(c) != t_INT)
		pari_err_TYPE("bigdiv", c);
	if (typ(d) != t_INT)
		pari_err_TYPE("bigdiv", d);
	
	if (signe(b) < 0)
	{
		if (equali1(a))
		{
			// Does d divide 1 - c?
			ret = !signe(modii(subsi(1, c), d));
			avma = ltop;
			return ret;
		}
		if (equalim1(a))
		{
			// Does d divide (-1)^b - c?
			ret = !signe(modii(subsi(mpodd(b) ? -1 : 1, c), d));
			avma = ltop;
			return ret;
		}
		/* a^b not in Z */
		pari_err_OP("bigdiv", a, b);
	} else if (!signe(b)) {
		// Does d divide a^0 - c?
		ret = !signe(modii(signe(a) ? subii(a, c) : subsi(1, c), d));
		avma = ltop;
		return ret;
	}
	
	if (signe(d) < 0)
		d = negi(d);
	if (cmpis(d, 2) <= 0)
	{
		if (cmpis(d, 0) == 0)
			pari_err(e_INV);
		if (equali1(d))
		{
			// Does 1 divide a^b - c?
			avma = ltop;
			return 1;
		}
		
		// Does 2 divide a^b - c?
		ret = !mpodd(subii(a, c));
		avma = ltop;
		return ret;
	}
	
	if (cmpis(gcdii(a, d), 1) > 0)
		ret = gequal(powgi(gmodulo(a, d), b), c);	// Not as slow as it looks
	else
		ret = gequal(powgi(gmodulo(a, d), modii(b, eulerphi(d))), c);
	avma = ltop;
	return ret;
}


GEN
solvePell(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("solvePell", n);
	if (signe(n) < 1)
		pari_err_DOMAIN("solvePell", "n", "<=", gen_0, n);
	if (Z_issquare(n))
		pari_err_IMPL("square Pell solutions");
	pari_sp ltop = avma;
	GEN C, t, x, y;
	long k = 1;
	long myprec = 125;
	pari_sp btop = avma;
	while (1) {
		C = contfrac0(gsqrt(n, myprec), NULL, 0);
		pari_sp ctop = avma;
		GEN p1, p2;	  /* vec */
		while (k++ <= glength(C)) {
			long i;
			p1 = cgetg(k, t_VEC);
			for (i = 1; i <= k; ++i)
				gel(p1, i) = gel(C, i);
			t = contfracback(p1, NULL);
			x = numer(t);
			y = denom(t);
			if (equali1(subii(sqri(x), mulii(n, sqri(y))))) {
				p2 = cgetg(3, t_VEC);
				gel(p2, 1) = gcopy(x);
				gel(p2, 2) = gcopy(y);
				p2 = gerepileupto(ltop, p2);
				return p2;
			}
			avma = ctop;
		}
		myprec <<= 1;
		avma = btop;
	}
	avma = ltop;
	return gen_0;
}


GEN
tetrMod(GEN a, GEN b, GEN M)
{
	pari_err_IMPL("tetrMod");
	// FIXME: Handle the case where gcd(a, M) > 1.
	if (typ(a) != t_INT)
		pari_err_TYPE("tetrMod", a);
	if (typ(b) != t_INT)
		pari_err_TYPE("tetrMod", b);
	if (typ(M) != t_INT)
		pari_err_TYPE("tetrMod", M);
	switch (signe(b)) {
		case -1:
			pari_err_OP("tetrMod", a, b);
		case 0:
			return gen_1;
	}
	pari_sp ltop = avma;
	GEN e = icopy(a), v;
	long vlen = itos_or_0(b)-1;
	if (vlen < 0)
		vlen = expi(subis(M, 1)) + 1;	// Upper bound on A003434
	else
		vlen = minuu(vlen, expi(subis(M, 1)) + 1);
	v = cgetg(vlen+1, t_VEC);
	gel(v, 1) = M;
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 2; i <= vlen; ++i)
	{
		gel(v, i) = geulerphi(gel(v, i - 1));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			v = gerepilecopy(btop, v);
	}
	for (i = glength(v); i >= 1; i--)
	{
		e = lift(powgi(gmodulo(a, gel(v, i)), e));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			e = gerepilecopy(btop, e);
	}
	e = gerepileuptoint(ltop, e);
	return e;
}


ulong
tetrMod_tiny(ulong a, ulong b, ulong M)
{
	(void)a;(void)b;(void)M;
	/*
	if (typ(a) != t_INT || typ(b) != t_INT || typ(M) != t_INT)
		pari_err(typeer, "tetrMod");
	pari_sp ltop = avma;
	GEN e = icopy(a), v, p1;
	long l3, l4;
	p1 = subis(b, 1);
	v = cgetg(itos(p1)+1, t_VEC);
	gel(v, 1) = M;
	l3 = glength(v);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 2; i <= l3; ++i)
	{
		gel(v, i) = geulerphi(gel(v, i - 1));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			v = gerepilecopy(btop, v);
	}
	l4 = glength(v);
	btop = avma, st_lim = stack_lim(btop, 1);
	GEN ii = gen_0;
	long l6 = -1 > 0;
	for (ii = stoi(l4); l6?gcmpgs(ii, 1) <= 0:gcmpgs(ii, 1) >= 0; ii = gaddgs(ii, -1))
	{
		e = lift(gpow(gmodulo(a, gel(v, gtos(ii))), e, prec));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &i, &e);
	}
	e = gerepileuptoint(ltop, e);
	return e;
	*/
	return NEVER_USED;
}


long
checkmult(GEN v, long verbose)
{
	pari_sp ltop = avma;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("checkmult", v);
	RgV_check_ZV(v, "checkmult");
	long n, l1 = lg(v), l2;

	// I suppose [] is a multiplicative sequence... if nothing else guard
	// against checking v[1] below.
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	GEN f, target, pr, ex;

	for (n = 6; n < l1; ++n) {
		if (uisprimepower(n, NULL))
			continue;
		f = factoru(n);
		pr = gel(f, 1);
		ex = gel(f, 2);
		l2 = lg(pr);
		long i;
		
		// Set target = prod v[p^e] for each p^e || n
		target = gen_1;
		for (i = 1; i < l2; ++i)
			target = mulii(target, gel(v, itos(powuu(pr[i], ex[i]))));
		
		// If v[n] is not equal to the target, the sequence is not multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not multiplicative at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	return 1;
}


long
checkcmult(GEN v, long verbose)
{
	pari_sp ltop = avma;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("checkcmult", v);
	RgV_check_ZV(v, "checkcmult");
	long n, l1 = lg(v), l2;

	// I suppose [] is a (completely) multiplicative sequence... if nothing else
	// guard against checking v[1] below.
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	GEN f, target, pr, ex;

	for (n = 4; n < l1; ++n) {
		if (uisprime(n))
			continue;
		f = factoru(n);
		pr = gel(f, 1);
		ex = gel(f, 2);
		l2 = lg(pr);
		long i;
		
		// Set target = prod v[p]^e for each p^e || n
		target = gen_1;
		for (i = 1; i < l2; ++i)
			target = mulii(target, powiu(gel(v, pr[i]), ex[i]));
		
		// If v[n] is not equal to the target, the sequence is not
		// completely multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not completely multiplicative at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	return 1;
}


long
checkadd(GEN v, long verbose)
{
	pari_sp ltop = avma;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("checkadd", v);
	RgV_check_ZV(v, "checkadd");
	long n, l1 = lg(v), l2;

	// I suppose [] is an additive sequence... if nothing else guard
	// against checking v[1] below.
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 0; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 0))
		return 0;
		
	GEN f, target, pr, ex;

	for (n = 6; n < l1; ++n) {
		if (uisprimepower(n, NULL))
			continue;
		f = factoru(n);
		pr = gel(f, 1);
		ex = gel(f, 2);
		l2 = lg(pr);
		long i;
		
		// Set target = sum v[p^e] for each p^e || n
		target = gen_0;
		for (i = 1; i < l2; ++i)
			target = addii(target, gel(v, itos(powuu(pr[i], ex[i]))));
		
		// If v[n] is not equal to the target, the sequence is not multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not additive at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	return 1;
}


long
checkcadd(GEN v, long verbose)
{
	pari_sp ltop = avma;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("checkcadd", v);
	RgV_check_ZV(v, "checkcadd");
	long n, l1 = lg(v), l2;

	// I suppose [] is a (completely) additive sequence... if nothing else
	// guard against checking v[1] below.
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 0; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 0))
		return 0;
		
	GEN f, target, pr, ex;

	for (n = 4; n < l1; ++n) {
		if (uisprime(n))
			continue;
		f = factoru(n);
		pr = gel(f, 1);
		ex = gel(f, 2);
		l2 = lg(pr);
		long i;
		
		// Set target = sum v[p]*e for each p^e || n
		target = gen_0;
		for (i = 1; i < l2; ++i)
			target = addii(target, muliu(gel(v, pr[i]), ex[i]));
		
		// If v[n] is not equal to the target, the sequence is not
		// completely multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not completely additive at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	return 1;
}


long
checkdiv(GEN v, long verbose/*=1*/)
{
	pari_sp ltop = avma;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("checkdiv", v);
	RgV_check_ZV(v, "checkdiv");
	long n, l1 = lg(v);
	for (n = 1; n < l1; ++n) {
		GEN vn = gel(v, n);
		long i;
		for (i = n + n; i < l1; i += n) {
			if (!dvdii(gel(v, i), vn)) {
				if (verbose)
					pariprintf("Not a divisibility sequence: a(%ld) = %Ps does not divide a(%ld) = %Ps.\n", n, gel(v, n), i, gel(v, i));
				avma = ltop;
				return 0;
			}
		}
		avma = ltop;
	}
	avma = ltop;
	return 1;
}


GEN
HurwitzClassNumber_small(ulong n)
{
	pari_sp ltop = avma;
	ulong c = coreu(n);
	ulong sq = usqrt(n/c);	// n = c * sq^2 with c squarefree
	GEN div = divisorsu(sq);
	pari_sp btop = avma;
	long i, last = lg(div), sum = 0;
	if (c == 3)
		last--;
	else if (c == 1 && (n&3) == 0)
		last -= 2;
	
	for (i = 1; i < last; ++i)
	{
		ulong D = n/(div[i]*div[i]);
		if ((1<<(D&3)) & 0x9) {
			sum += itos(classno(utoineg(D)));
			avma = btop;
		}
	}
	
	GEN ret = stoi(sum);
	
	// This loop handles the special case where the result may be a fraction.
	for (i = last; i < lg(div); ++i)
	{
		ulong D = n/(div[i]*div[i]);
		if ((1<<(D&3)) & 0x9)
			ret = gadd(ret, gdivgs(classno(utoineg(D)), D > 4 ? 1 : 6-D));
	}
	
	return gerepileupto(ltop, ret);
}


// Class number of n. Not stack clean.
INLINE GEN
classno_fast(GEN n)
{
	if (lgefint(n) == 3)
		return classno(n);
	return member_no(quadclassunit0(n, 0, NULL, DEFAULTPREC));
}


// The Hurwitz class number of n. An integer unless n is of the form 3k^2 or 4k^2.
GEN
HurwitzClassNumber(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("HurwitzClassNumber",n);
	if (signe(n) < 1) {
		if (isintzero(n))
			return gdiv(gen_m1, stoi(12));
		pari_err(e_IMPL,"HurwitzClassNumber(n) for n < 0.");
	}
	ulong nn = itou_or_0(n);
#ifdef LONG_IS_64BIT
	if (nn && nn < 500000000000000000UL)
#else
	if (nn)
#endif
		return HurwitzClassNumber_small(nn);
	
	pari_sp ltop = avma;
	GEN div = divisors(gel(core0(n, 1), 2)), ret = gen_0;
	pari_sp btop = avma;
	long i, mx = lg(div);
	for (i = 1; i < mx; ++i) {
		GEN D = diviiexact(negi(n), sqri(gel(div, i)));
		if ((1<<mod4(D)) & 0x9)	{
			// mod4(D) = |D| mod 4, so this selects |D| in {0, 3} mod 4 thus D in {0, 1} mod 4
			ret = gadd(ret, gdiv(classno_fast(D), gmaxsg(1, addis(D, 6))));
			ret = gerepileupto(btop, ret);	// perhaps do this only when the stack is low?
		}
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


// Computes tau(p) assuming p is prime.
GEN
taup_small(ulong p)
{
	ulong k;
	pari_sp ltop = avma, st_lim = stack_lim(ltop, 1);
	GEN ret = gen_0;
	for (k = 1; k < p; ++k) {
		ret = addii(ret, mulii(sumdivk(utoipos(k), 5), sumdivk(utoipos(p - k), 5)));
		if (low_stack(st_lim, stack_lim(ltop, 1)))
			ret = gerepileuptoint(ltop, ret);
	}
	ret = diviuexact(subii(addii(mulsi(65, sumdivk(utoipos(p), 11)), mulsi(691, sumdivk(utoipos(p), 5))), mulsi(174132, ret)), 756);
	return gerepileuptoint(ltop, ret);
}


// Computes tau(p) assuming p > 3 is prime.
// Algorithm based on the Selberg trace formula, as analyzed by
// Denis Xavier Charles, Computing the Ramanujan tau function, The Ramanujan
// Journal 11:2 (April 2006), pp. 221-224.
GEN
taup_big(GEN p)
{
	pari_sp ltop = avma;
	GEN P = cgetg(7, t_VEC);
	gel(P, 1) = negi(powiu(p, 5));
	gel(P, 2) = mulsi(15, powiu(p, 4));
	gel(P, 3) = mulsi(-35, powiu(p, 3));
	gel(P, 4) = mulsi(28, sqri(p));
	gel(P, 5) = mulsi(-9, p);
	gel(P, 6) = gen_1;
	//P = gtopolyrev(P, -1);
	P = gerepileupto(ltop, P);
	GEN p4 = mulsi(4, p), lim = sqrtint(p4);
	
	pari_sp btop = avma, st_lim = stack_lim(btop, 2);
	GEN t, ret = gen_0;
	for (t = gen_1; cmpii(t, lim) <= 0; t = addis(t, 1)) {
		GEN t2 = sqri(t);
		GEN tmp = gmul(poleval_denseint(P, t2), HurwitzClassNumber(subii(p4, t2)));
		ret = gadd(ret, tmp);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &t, &ret);
	}
	ret = gsub(gsubgs(gdivgs(mulii(powiu(p, 5), HurwitzClassNumber(p4)), 2), 1), ret);
	return gerepileupto(ltop, ret);
}


GEN
taup(GEN p, long e)
{
	if (typ(p) != t_INT)
		pari_err_TYPE("taup",p);
	
	ulong pp = itou_or_0(p);
	pari_sp ltop = avma;
	GEN t = pp && pp < 220000 ? taup_small(pp) : taup_big(p);
	
	// Special case for common low exponents.
	if (e == 1)
		return t;
	if (e < 4) {
		GEN ret;
		if (e == 2)
			ret = subii(sqri(t), powiu(p, 11));
		else
			ret = mulii(subii(sqri(t), shifti(powiu(p, 11), 1)), t);
		return gerepileupto(ltop, ret);
	}
	
	// Exponent calculation in full generality
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN ret = gen_0;
	long j, l2 = e/2;
	for (j = 0; j <= l2; ++j) {
		// (-1)^j*binomial(e-j, e-2*j)*p^(11*j)*t^(e-2*j)
		GEN tmp = mulii(mulii(binomialuu(e - j, e - 2*j), powiu(t, e-2*j)), powiu(p, 11*j));
		if (j&1)
			ret = subii(ret, tmp);
		else
			ret = addii(ret, tmp);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			ret = gerepileupto(btop, ret);
	}
	return gerepileupto(ltop, ret);
}


GEN
tau(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("tau",n);
	pari_sp ltop = avma;
	GEN ret = gen_1, f = Z_factor(n);
	long i, l1 = glength(gel(f, 1));
	for (i = 1; i <= l1; ++i)
		ret = mulii(ret, taup(gcoeff(f, i, 1), itos(gcoeff(f, i, 2))));
	return gerepileupto(ltop, ret);
}
