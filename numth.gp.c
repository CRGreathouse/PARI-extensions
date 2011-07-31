/******************************************************************************/
/**												 Other number theory															**/
/******************************************************************************/

long
istotient(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "istotient");
	if (signe(n) < 1)
		return 0;
	if(mod2(n))
		return isint1(n);

	pari_sp ltop = avma;
	GEN k, p, d, p2;
	k = n;
	while (1) {
		if (totientHelper(k, gen_2)) {
			avma = ltop;
			return 1;
		}
		if (mod2(k))
			break;
		k = shifti(k, -1);
		k = gerepileuptoint(ltop, k);
	}
	p2 = divisors(shifti(n, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 1; i < lg(p2); ++i) {
		d = shifti(gel(p2, i), 1);
		if (!(isprime(p = addis(d, 1))))
			continue;
		k = diviiexact(n, d);
		while (1) {
			if (totientHelper(k, p)) {
				avma = ltop;
				return 1;
			}
			if (!dvdii(k, p))
				break;
			k = diviiexact(k, p);
		}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &k);
	}
	avma = ltop;
	return 0;
}


/* This function should be called only internally.  n and m are positive t_INT
 * values.  Test with
#
sum(n=1,35214/2,istotient(2*n)==0)
forstep(n=2,1e6,2,istotient(2*n))
 */
long
totientHelper(GEN n, GEN m)
{
	if (mod2(n))
		return equali1(n);
	pari_sp ltop = avma;
	GEN k, p, d, p1;
	p1 = divisors(shifti(n, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long l2;
	for (l2 = 1; l2 < lg(p1); ++l2) {
		d = shifti(gel(p1, l2), 1);
		if ((cmpii(d, m) < 0) || !(isprime(p = addis(d, 1))))
			continue;
		k = diviiexact(n, d);
		while (1) {
			if (totientHelper(k, p)) {
				avma = ltop;
				return 1;
			}
			if (!dvdii(k, p))
				break;
			k = diviiexact(k, p);
			/* // Slower:
			GEN quotient, remainder;
			quotient = dvmdii(k, p, &remainder);
			if (remainder != gen_0)	// documentation guarantees that gen_0, not a copy, is returned if p|k.
				break;
			k = quotient;
			*/
		}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &k);
	}
	avma = ltop;
	return 0;
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
		pari_err(typeer, "Faulhaber; must be a variable");
	pari_sp btop = avma;
	long i = 0;
	for (; i <= e; ++i)
	{
		ret = gadd(ret, gmul(gmul(binomialuu(e + 1, i), bernfrac(i)), gpowgs(x, e + 1 - i)));
		ret = gerepileupto(btop, ret);
	}
	ret = gsubstpol(gadd(gdivgs(ret, e + 1), gpowgs(x, e)), x, a);
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
rp(long b)
{
	pari_sp ltop = avma;
	GEN ret;
	GEN B = int2n(b - 1);
	ret = gnextprime(addii(B, randomi(B)));
	ret = gerepileupto(ltop, ret);
	return ret;
}


// FIXME: Take pity and do some numerical analysis here!
ulong
cuberoot(ulong n)
{
	ulong ret = pow(n + 0.5, 1.0/3);
#ifdef LONG_IS_64BIT
	if (n < 100000000000001UL)
#endif
		return ret;
	if (n >= 18446724184312856125UL)
		return 2642245UL;
	if (ret * ret * ret > n) {
		ret--;
		while (ret * ret * ret > n)
			ret--;
		return ret;
	}
	while ((ret + 1) * (ret + 1) * (ret + 1) <= n)
		ret++;
	return ret;
}


GEN
cuberootint(GEN x)
{
	pari_sp ltop = avma;
	long t = typ(x);
	GEN ret = NEVER_USED;
	if (t == t_INT) {
		ret = gfloor(powrfrac(itor(x, lgefint(x)), 1, 3));
	} else if (t == t_REAL) {
		ret = gfloor(powrfrac(x, 1, 3));
		x = gfloor(x);
	} else {
		pari_err(typeer, "cuberootint");
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
	
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2

	// First loop: remove tiny primes, don't calculate cube roots
	// 99.7% of non-squarefree numbers are detected by this loop (or the above)
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > 97)
			break;
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
		}
	}
	
	// Beyond this point, 99.89% of numbers are squarefree.
	ulong last = cuberoot(n);
	ulong last1 = minuu(last, CUTOFF);
	for (;;)
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = cuberoot(n);
			last1 = minuu(last, CUTOFF);
		}
		
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > last1)
			break;
	}

#ifdef LONG_IS_64BIT
	if (n < CUTOFF * CUTOFF * CUTOFF)
#elif CUTOFF <= 1621
	pari_warn(warner, "issquarefree: cutoff marginal, performace suffers");
	if (n < CUTOFF * CUTOFF * CUTOFF)
#endif
	return n == 1 || !uissquare(n);

	// n is at least CUTOFF^3 and is not divisible by any prime under CUTOFF
	if (last < 65536) {	// maxprime() > 65536
		for (;;)
		{
			if (n%p == 0) {
				n /= p;
				if (n%p == 0)
					return 0;
				last = cuberoot(n);
			}
			
			NEXT_PRIME_VIADIFF(p, primepointer);
			if (p > last)
				break;
		}
		return !uissquare(n);
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
	if (last > 300000 || last > maxprime())	// TODO: Find good breakover point here
	{
		long ret = Z_issquarefree(stoi(n));
		avma = ltop;
		return ret;
	}
	for (;;)
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = cuberoot(n);
		}
		
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > last)
			break;
	}
	
	return !uissquare(n);
#undef CUTOFF
}


// TODO: All of the countPowerful functions could be improved greatly
// by using countSquarefree to check the # of squarefree numbers at the point
// that dividing will give exactly 1.	Cost: two invocations of countSquarefree
// at size ~ cuberoot(n).	Savings: ~0.63 cuberoot(n) invocations of
// issquarefree and ~0.38 cuberoot(n) powerings, square roots, and divisions.
// Big win!
ulong
ucountPowerfulu(ulong n)
{
#if 1
	// About 33% faster
	ulong k, breakpoint = cuberoot(n >> 2);
	ulong res = ucountSquarefree(cuberoot(n)) - ucountSquarefree(breakpoint);
#else
	ulong k, breakpoint = cuberoot(n), res = 0;
#endif
	for (k = 1; k <= breakpoint; k++)
		if (issquarefree_small(k))
			res += usqrtsafe(n / k) / k;
	return res;
}


ulong
ucountPowerfuli(GEN n)
{
	pari_sp ltop = avma;
	ulong cube_root = itou(cuberootint(n));
	ulong res = 0, k;
	for (k = 1; k <= cube_root; k++)
		if (issquarefree_small(k)) {
			res += itos(divis(sqrti(divis(n, k)), k));
			//res += itos(sqrti(divii(n, powuu(k, 3))));    // About 35% slower
			avma = ltop;
		}
	return res;
}


/*
// FIXME: Write this, I apparently needed it sometime before but never
// finished coding it.
GEN
listPowerful(GEN lim)
{
	// Should check for memer before attempting further...
	pari_sp ltop = avma;
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "listPowerful");

#ifdef LONG_IS_64BIT
	if (lim > mkintn(4, 227191940, 3022159750, 2605788960, 0))
		pari_err(overflower, "listPowerful");	// around 2^64 or more entries;
												// surely not enough memory
#else
	// 977971126747547528
	if (lim > uu32toi(227701646, 3932178312))
		pari_err(overflower, "listPowerful");	// 2^32 or more entries
#endif
	
}
*/


// FIXME: seems broken, at least for small values, unlike ucountPowerfuli;
// c.f. isPowerful
// Varies as zeta(3/2)/zeta(3) n^1/2 + zeta(2/3)/zeta(2) n^1/3 + o(n^1/6)
// estPowerful(n)=zeta(3/2)/zeta(3)*sqrt(n) + zeta(2/3)/zeta(2)*n^(1/3)
GEN
countPowerful(GEN n)
{
	pari_sp ltop = avma;
	GEN p1, ret;
pari_warn(warner, "countPowerful is possibly broken, FIXME");
	
	if (lgefint(p1 = gfloor(n)) <= 3) {
		ret = utoi(ucountPowerfulu(itou(p1)));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	// ui, 32-bit: 3909962179575504899
	// ui, 64-bit: ~72047453657149422936422171552392422209
	if (typ(n) == t_INT)
		n = itor(n, lgefint(n));	// Is this the right amount of precision?
	else if (typ(n) != t_REAL)
		pari_err(typeer, "countPowerful");

	//p1 = gfloor(gpow(gadd(n, ghalf), ginv(stoi(3)), FAKE_PREC));
	p1 = gfloor(powrfrac(gadd(n, ghalf), 1, 3));
	n = gfloor(n);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k = gen_0;
	ret = gen_0;
	for (k = gen_1; cmpii(k, p1) <= 0; k = addis(k, 1))
	{
		if (Z_issquarefree(k))
			ret = addii(ret, sqrti(gdivent(n, powis(k, 3))));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &ret, &k);
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


INLINE long
moebiusu(ulong n)
{
	pari_sp ltop = avma;
	long ret = moebius(utoi(n));
	avma = ltop;
	return ret;
}


// TODO: Doesn't really save much time vs. the original.	To improve, a sieve
// would be needed (to calculate the values of moebius faster).
ulong
ucountSquarefree(ulong lim)
{
	ulong b = usqrtsafe(lim >> 1);
	ulong k;
	ulong ret = 0;
	for (k = 1; k <= b; k++)
		ret += moebiusu(k) * (lim / (k * k));
	ulong p3 = usqrtsafe(lim);
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
		pari_err(typeer, "countSquarefree");

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
		pari_err(arither1, "Mfactor");
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "Mfactor");
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err(arither1, "Mfactor");

	GEN v, k, p1, p2;
	v = cgetg(1, t_VEC);
	if (signe(p) < 1)
		pari_err(talker, "p must be positive");
	if (!(isprime(p)))
		pari_err(talker, "p must be prime");
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
		if (!gequalgs(gpow(gmodulsg(2, q), modii(p, subis(q, 1)), FAKE_PREC), 1))
			continue;
		v = concat(v, q);
		long i = 2;
		while (gequal1(gpow(gmodulsg(2, powis(q, i)), modii(p, mulii(subis(q, 1), powis(q, i - 1))), FAKE_PREC)))
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
	if (typ(a) != t_INT || typ(b) != t_INT || typ(c) != t_INT)
		pari_err(arither1, "bigfactor");
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err(arither1, "bigfactor");
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "bigfactor");
	long lm = itos(lim);
	if (lm > maxprime())
		pari_err(primer1, lim);

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
			p2 = Z_factor(subii(stoi(1 - 2 * mpodd(b)), c));
			p2 = gerepileupto(ltop, p2);
			return p2;
		}
		/* a^b not in Z */
		pari_err(talker, "not an integer power in bigfactor");
	}
	long p3 = minss(itos(a), lm);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long p = 0;
	byteptr primepointer = diffptr;
	GEN p5 = gen_0;		/* int */
	if (p3 > maxprime())
		pari_err(primer1, stoi(p3));

	// First loop -- deal with small primes
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > p3)
			break;
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
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN p6 = gen_0;
		for(;;)
		{
			if (cmpis(gcdii(a, stoi(p)), 1) > 0)
				p6 = b;
			else
				p6 = modii(b, mulis(powuu(p, i - 1), p - 1));
			if (!gequal(gpow(gmodulo(a, powuu(p, i)), p6, FAKE_PREC), c))
				break;
			v = concat(v, stoi(p));
			i++;
			if (low_stack(st_lim, stack_lim(btop, 1)))
				v = gerepileupto(btop, v);
		}
	}
	//long p4 = p;//itos(gprecprime(gmaxgs(start, p3))) + 1;
	v = gerepileupto(btop, v);
	btop = avma;	// Should this be reset (as currently) or removed?
	st_lim = stack_lim(btop, 1);	// Ditto
	p = 0;
	primepointer = diffptr;		/* bptr */

	//while (p < p4)
	//	NEXT_PRIME_VIADIFF(p, primepointer);

	// Second loop -- most of the work is done here
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		//if (p < p4)
		//	continue;
		if (cmpsi(p, lim) > 0)
			break;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &v);
		if (!gequal(gpowgs(gmodulo(a, stoi(p)), smodis(b, p - 1)), c))
			continue;
		v = concat(v, stoi(p));
		long i = 2;
		while (gequal(gpow(gmodulo(a, powis(stoi(p), i)), modii(b, mulis(powis(stoi(p), i - 1), p - 1)), FAKE_PREC), c))
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
	if (typ(a) != t_INT || typ(b) != t_INT || typ(c) != t_INT || typ(d) != t_INT)
		pari_err(arither1, "bigdiv");
	
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
		pari_err(arither1, "bigdiv");
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
			pari_err(gdiver, "bigdiv");
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


// Assumes that f is a nonzero t_POL with t_INT coefficients.  Returns a
// polynomial with roots that are precisely the squares of the roots of f.
GEN
graeffe(GEN f)
{
	pari_sp ltop = avma;
	GEN g, h, ret, x = mkpoln(2, gen_1, gen_0);
	if (typ(f) != t_POL)
		pari_err(typeer, "graeffe");
		
	long d = degpol(f);
	long gsize = (d >> 1) + 1, hsize = (d + 1) >> 1;
	long i;
	
	g = cgetg(gsize+2, t_POL);
	g[1] = evalvarn(0);
	for (i = 1; i <= gsize; ++i)
		gel(g, i + 1) = polcoeff0(f, (2*i) - 2, -1);
	g = normalizepol_lg(g, gsize + 2);
	
	h = cgetg(hsize+2, t_POL);
	h[1] = evalvarn(0);
	for (i = 1; i <= hsize; ++i)
		gel(h, i + 1) = polcoeff0(f, (2*i) - 1, -1);
	h = normalizepol_lg(h, hsize + 2);
	
	ret = gsub(gsqr(g), gmul(x, gsqr(h)));
	ret = gerepileupto(ltop, ret);
	setvarn(ret, varn(f));
	return ret;
}


long
poliscyclo(GEN f)
{
	if (typ(f) != t_POL)
		pari_err(typeer, "poliscyclo");
	if (!isint1(leading_term(f)))
		return 0;
	long degree = degpol(f);
	if (degree < 2)
		return degree == 1 && is_pm1(constant_term(f));
	
	return RgX_is_ZX(f) && BradfordDavenport(f) && gisirreducible(f) == gen_1;
}


long
poliscycloproduct(GEN f, long flag)
{
	if (typ(f) != t_POL)
		pari_err(typeer, "poliscyclo");
	if (!isint1(leading_term(f)))
		return 0;
	long degree = degpol(f);
	if (degree < 2)
		return degree == 1 && is_pm1(constant_term(f));
	
	if (!RgX_is_ZX(f) || !BradfordDavenportProduct(f))
		return 0;
	if (flag == 0)
		return 1;
	
	// Determine degree
	pari_err(impl, "finding the degree");
	return NEVER_USED;
}


// Checks if f is a cyclotomic polynomial.  Assumes f is an irreducible t_POL
// of t_INTS with degree > 1.
long
BradfordDavenport(GEN f) {
	pari_sp ltop = avma;
	GEN f1, f2, fn, mx;
	
	f1 = graeffe(f);
	if (polequal(f, f1)) {
		avma = ltop;
		return 1;
	}

	// Set up variables
	long var = varn(f);	// Variable number in polynomial
	mx = mkpoln(2, gen_m1, gen_0);	// -x
	setvarn(mx, var);
	
	fn = gsubst(f, var, mx);
	if (ZX_equal(f1, fn) && BradfordDavenport(fn)) {
		avma = ltop;
		return 1;
	}
	long ret = polissquareall(f1, &f2) && BradfordDavenport(f2);
	avma = ltop;
	return ret;
}


// Checks if f is a product of distinct cyclotomic polynomials.  Assumes f is
// a t_POL of t_INTS with degree > 1.
long
BradfordDavenportProduct(GEN f) {
	pari_sp ltop = avma;
	GEN f1, f2, fn, mx;
	
	f1 = graeffe(f);
	if (polequal(f, f1)) {
		avma = ltop;
		return 1;
	}

	// Set up variables
	long var = varn(f);	// Variable number in polynomial
	mx = mkpoln(2, gen_m1, gen_0);	// -x
	setvarn(mx, var);
	
	fn = gsubst(f, var, mx);
	if (ZX_equal(f1, fn) && BradfordDavenport(fn)) {
		avma = ltop;
		return 1;
	}
	long ret = polissquareall(f1, &f2) && BradfordDavenport(f2);
	avma = ltop;
	return ret;
}

// TODO: Better algorithm than continued fractions?  Failing that, at least
// control number of digits?
GEN
solvePell(GEN n, long prec)
{
	pari_sp ltop = avma;
	GEN myprec = gen_0, C = gen_0, k = gen_1, t = gen_0, x = gen_0, y = gen_0;
	if (typ(n) != t_INT)
		pari_err(typeer, "solvePell");
	myprec = stoi(getrealprecision());
	setrealprecision(125, &prec);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while (1) {
		C = contfrac0(gsqrt(n, prec), NULL, 0);
{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN p1 = gen_0, p2 = gen_0;	  /* vec */
		while (gcmpgs(k, glength(C)) <= 0) {
			long i;
			p1 = cgetg(gtos(k)+1, t_VEC);
			for (i = 1; gcmpsg(i, k) <= 0; ++i)
				gel(p1, i) = gcopy(gel(C, i));
			t = contfracback(p1, NULL);
			x = numer(t);
			y = denom(t);
			if (gequal1(gsub(gsqr(x), gmul(n, gsqr(y))))) {
				setrealprecision(gtos(myprec), &prec);
				p2 = cgetg(3, t_VEC);
				gel(p2, 1) = gcopy(x);
				gel(p2, 2) = gcopy(y);
				p2 = gerepileupto(ltop, p2);
				return p2;
			}
			k = gaddgs(k, 1);
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 6, &p1, &t, &x, &y, &p2, &k);
		}
}
		setrealprecision(2*getrealprecision(), &prec);
		if (low_stack(st_lim, stack_lim(btop, 1)))
		gerepileall(btop, 5, &C, &t, &x, &y, &k);
	}
	avma = ltop;
	return gen_0;
}


GEN
tetrMod(GEN a, GEN b, GEN M)
{
	pari_err(talker, "busted");
	// FIXME: Handle the case where gcd(a, M) > 1.
	if (typ(a) != t_INT || typ(b) != t_INT || typ(M) != t_INT)
		pari_err(typeer, "tetrMod");
	switch (signe(b)) {
		case -1:
			pari_err(talker, "negative argument");
		case 0:
			return gen_1;
	}
	pari_sp ltop = avma;
	GEN e = icopy(a), v;
	// __builtin_ffsl or expi
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
		e = lift(gpow(gmodulo(a, gel(v, i)), e, FAKE_PREC));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			e = gerepilecopy(btop, e);
	}
	e = gerepileuptoint(ltop, e);
	return e;
}


GEN
tetrMod_tiny(ulong a, ulong b, ulong M)
{
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
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkmult");

	// I suppose [] is a multiplicative sequence... if nothing else guard
	// against checking v[1] below.
	long n, l1 = lg(v), l2;
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	pari_sp ltop = avma;
	GEN f, target;

	// Require all arguments to be integers
	for (n = 2; n < 6 && n < l1; ++n)
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
			
	for (n = 6; n < l1; ++n) {
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
		if (uisprimepower(n))
			continue;
		f = Z_factor(stoi(n));
		l2 = glength(gel(f, 1));
		long i;
		
		// Set target = prod v[p^e] for each p^e || n
		target = gen_1;
		for (i = 1; i <= l2; ++i)
			target = mulii(target, gel(v, itos(powii(gcoeff(f, i, 1), gcoeff(f, i, 2)))));
		
		// If v[n] is not equal to the target, the sequence is not multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not multiplicative at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	avma = ltop;
	return 1;
}


long
checkcmult(GEN v, long verbose)
{
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkmult");

	// I suppose [] is a (completely) multiplicative sequence... if nothing else
	// guard against checking v[1] below.
	long n, l1 = lg(v), l2;
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	pari_sp ltop = avma;
	GEN f, target;

	// Require all arguments to be integers
	for (n = 2; n < 4 && n < l1; ++n)
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
			
	for (n = 4; n < l1; ++n) {
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
		if (uisprime(n))
			continue;
		f = Z_factor(stoi(n));
		l2 = glength(gel(f, 1));
		long i;
		
		// Set target = prod v[p^e] for each p^e || n
		target = gen_1;
		for (i = 1; i <= l2; ++i)
			target = mulii(target, gel(v, itos(powii(gcoeff(f, i, 1), gcoeff(f, i, 2)))));
		
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
	avma = ltop;
	return 1;
}


long
checkdiv(GEN v, long verbose/*=1*/)
{
	pari_sp ltop = avma;
	long l1;	  /* lg */
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkdiv");

	l1 = lg(v);
	pari_sp btop = avma;
	long n, l2;
	long l3;	  /* lg */
	for (n = 1; n < l1; ++n) {
		l2 = n + n;
		l3 = lg(v);
{
		pari_sp btop = avma;
		GEN i = gen_0;
		long l4 = n > 0;	  /* bool */
		for (i = stoi(l2); l4?gcmpgs(i, l3-1) <= 0:gcmpgs(i, l3-1) >= 0; i = gaddgs(i, n)) {
			if (!gequal0(gmod(gel(v, gtos(i)), gel(v, n)))) {
				if (verbose)
					pariprintf("Not a divisibility sequence: a(%ld) = %Ps does not divide a(%Ps) = %Ps.\n", n, gel(v, n), i, gel(v, gtos(i)));
				avma = ltop;
				return 0;
			}
			i = gerepileupto(btop, i);
		}
}
		avma = btop;
	}
	avma = ltop;
	return 1;
}

