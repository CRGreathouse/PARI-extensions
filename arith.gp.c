/******************************************************************************/
/**											 Other arithmetic functions													*/
/******************************************************************************/

GEN
oddres(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("oddres", n);
	long v = vali(n);
	return v ? shifti(n, -v) : icopy(n);
}


// Faster than
//	 !cmpii(n, int2n(vali(n)))
//	 !cmpis(shifti(n, -vali(n)), 1)
//	 expi(n) == vali(n)
//	 hamming(n) == 1
// even for single-word values, and much faster for multiword values.
// If all you have is a word, you can just use n & !(n & (n-1)).
long
ispow2(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("ispow2", n);
	if (signe(n) < 1)
		return 0;
	
	GEN xp = int_LSW(n);
	long lx = lgefint(n);
	ulong u = *xp;
	long i = 3;
	for (; i < lx; ++i)
	{
		if (u != 0)
			return 0;
		xp = int_nextW(xp);
		u = *xp;
	}
	//return hamming_word(u) == 1;	// 14% slower
	return !(u & (u-1));
}


long
ispow3(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("ispow3", n);

	// Remove negatives and 25% of random numbers
	if (signe(n) < 1 || !(10 & (1 << mod8(n))))
		return 0;

	ulong nn = itou_or_0(n);
	if (nn)
		return ispow3_tiny(nn);

	// Remove almost all random numbers
#ifdef LONG_IS_64BIT
	if (!ispow3_tiny(smodis(n, 18236498188585393201UL)))
#else
	if (!ispow3_tiny(smodis(n, 1743392200UL)))
#endif
		return 0;
	
	pari_sp ltop = avma;
	double sz3 = dbllog2r(itor(n, DEFAULTPREC)) / log2(3);
	int e = (int)(sz3 + 0.5);
	long ret = equalii(n, powis(stoi(3), e));
	avma = ltop;
	return ret;
}


long
issm3(long n)
{
	static const long sm3table[] = {
#ifdef LONG_IS_64BIT
		0, 0, 4052555153018976267, 1350851717672992089, 0, 450283905890997363, 150094635296999121, 0, 50031545098999707, 0, 16677181699666569, 5559060566555523, 0, 1853020188851841, 617673396283947, 0, 205891132094649, 0, 68630377364883, 22876792454961, 0, 7625597484987, 2541865828329, 0, 847288609443, 282429536481, 0, 94143178827, 0, 31381059609, 10460353203, 0,
#endif
		3486784401, 1162261467, 0, 387420489, 0, 129140163, 43046721, 0, 14348907, 4782969, 0, 1594323, 531441, 0, 177147, 0, 59049, 19683, 0, 6561, 2187, 0, 729, 0, 243, 81, 0, 27, 9, 0, 3, 1
	};
	
	n >>= __builtin_ctzl(n);
	return n > 0 && n == sm3table[__builtin_clzl(n)];
}


long ispow3_tiny(ulong n)
{
	static const ulong pow3table[] = {
#ifdef LONG_IS_64BIT
		12157665459056928801UL, 0, 4052555153018976267UL, 1350851717672992089UL, 0, 450283905890997363UL, 150094635296999121UL, 0, 50031545098999707UL, 0, 16677181699666569UL, 5559060566555523UL, 0, 1853020188851841UL, 617673396283947UL, 0, 205891132094649UL, 0, 68630377364883UL, 22876792454961UL, 0, 7625597484987UL, 2541865828329UL, 0, 847288609443UL, 282429536481UL, 0, 94143178827UL, 0, 31381059609UL, 10460353203UL, 0,
#endif
		3486784401, 1162261467, 0, 387420489, 0, 129140163, 43046721, 0, 14348907, 4782969, 0, 1594323, 531441, 0, 177147, 0, 59049, 19683, 0, 6561, 2187, 0, 729, 0, 243, 81, 0, 27, 9, 0, 3, 1
	};
	
	return n == pow3table[__builtin_clzl(n)];
}


// Approximate speed on a Phenom II at various input sizes:
// size	 cycles		time
// 1e00			110	 40 ns
// 1e10			150	 50 ns
// 1e20			268	100 ns
// 1e30			282	100 ns
// 1e50			297	110 ns
// 1e75			316	110 ns
// 1e100		 322	120 ns
// 1e200		 385	140 ns
// 1e500		 569	200 ns
// 1e1000		941	300 ns
// 1e2000	 1710	600 ns
// 1e5000	 3910	1.4 µs
// 1e10000	7540		3 µs
// 1e20000 15200		5 µs
// 1e50000 37000	 13 µs
long
isFibonacci(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("isFibonacci", n);
	if (!is_bigint(n))
		return isSmallFib(itos(n));
	pari_sp ltop = avma;

	// Good residue classes: A189761 (discovery predates sequence!)
#ifdef LONG_IS_64BIT
	// Does not catch enough non-Fibonacci numbers to save time.
	// if ((1 << smodis(n, 55)) & 0x2F7FFBFFDFDED0) return 0;
#endif
	long rem = smodis(n, 17711);
	if (rem & 1) {
		if (rem > 5 && rem != 13 && rem != 21 && rem != 55 && rem != 89 && rem != 233 && rem != 377 && rem != 987 && rem != 1597 && rem != 4181 && rem != 6765 && rem != 15127 && rem != 17567 && rem != 17703)
			return 0;
	} else {
		if (rem > 2 && rem < 17708 && rem != 8 && rem != 34 && rem != 144 && rem != 610 && rem != 2584 && rem != 10946 && rem != 16724 && rem != 17334 && rem != 17656 && rem != 17690)
			return 0;
	}
	GEN k = sqri(n);
	k = addis(mulis(k, 5), 4);	// Multiplication via shifts is slower here (by ~4 cycles)
	long res = Z_issquare(k) || (signe(n) > 0 && Z_issquare(subis(k, 8)));
	avma = ltop;
	return res;
}


static const long smallfib[] = {
#ifdef LONG_IS_64BIT
	-7540113804746346429,-2880067194370816120,-1100087778366101931,
	-420196140727489673,-160500643816367088,-61305790721611591,
	-23416728348467685,-8944394323791464,-3416454622906707,
	-1304969544928657,-498454011879264,-190392490709135,-72723460248141,
	-27777890035288,-10610209857723,-4052739537881,-1548008755920,
	-591286729879,-225851433717,-86267571272,-32951280099,-12586269025,
	-4807526976,-1836311903,-701408733,-267914296,-102334155,-39088169,
	-14930352,-5702887,-2178309,-832040,-317811,-121393,-46368,-17711,
	-6765,-2584,-987,-377,-144,-55,-21,-8,-3,-1,0,1,1,2,3,5,8,13,21,34,
	55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657,
	46368,75025,121393,196418,317811,514229,832040,1346269,2178309,
	3524578,5702887,9227465,14930352,24157817,39088169,63245986,
	102334155,165580141,267914296,433494437,701408733,1134903170,
	1836311903,2971215073,4807526976,7778742049,12586269025,20365011074,
	32951280099,53316291173,86267571272,139583862445,225851433717,
	365435296162,591286729879,956722026041,1548008755920,2504730781961,
	4052739537881,6557470319842,10610209857723,17167680177565,
	27777890035288,44945570212853,72723460248141,117669030460994,
	190392490709135,308061521170129,498454011879264,806515533049393,
	1304969544928657,2111485077978050,3416454622906707,5527939700884757,
	8944394323791464,14472334024676221,23416728348467685,
	37889062373143906,61305790721611591,99194853094755497,
	160500643816367088,259695496911122585,420196140727489673,
	679891637638612258,1100087778366101931,1779979416004714189,
	2880067194370816120,4660046610375530309,7540113804746346429
#else
	-1836311903,-701408733,-267914296,-102334155,-39088169,-14930352,
	-5702887,-2178309,-832040,-317811,-121393,-46368,-17711,-6765,-2584,
	-987,-377,-144,-55,-21,-8,-3,-1,0,1,1,2,3,5,8,13,21,34,55,89,144,
	233,377,610,987,1597,2584,4181,6765,10946,17711,28657,46368,75025,
	121393,196418,317811,514229,832040,1346269,2178309,3524578,5702887,
	9227465,14930352,24157817,39088169,63245986,102334155,165580141,
	267914296,433494437,701408733,1134903170,1836311903
#endif
};


long
isSmallFib(long n)
{
	int left = 0;
#ifdef LONG_IS_64BIT
	int right = 139;
#else
	int right = 70;
#endif
	while (left + 1 < right) {
		int mid = (left + right) >> 1;
		if (n < smallfib[mid])
			right = mid;
		else
			left = mid;
	}
	return n == smallfib[left] || n == smallfib[right];
}


static void
lucasmod(ulong n, GEN m, GEN *a, GEN *b)
{
	GEN z, t, zt;
	if (!n) { *a = gen_2; *b = gen_1; return; }
	lucasmod(n >> 1, m, &z, &t);
	zt = mulii(z, t);
	switch(n & 3) {
		case  0: *a = addsi(-2,sqri(z)); *b = addsi(-1,zt); break;
		case  1: *a = addsi(-1,zt);      *b = addsi(2,sqri(t)); break;
		case  2: *a = addsi(2,sqri(z));  *b = addsi(1,zt); break;
		case  3: *a = addsi(1,zt);       *b = addsi(-2,sqri(t));
	}
	*a = modii(*a, m);
	*b = modii(*b, m);
}


GEN
fibomod(long n, GEN m)
{
	if (typ(m) != t_INT)
		pari_err_TYPE("fibomod", m);
	pari_sp av = avma;
	GEN a, b;
	if (!n) return gen_0;
	lucasmod((ulong)(labs(n)-1), mulis(m, 5), &a, &b);
	a = modii(diviuexact(addii(shifti(a,1),b), 5), m);
	if (n < 0 && !odd(n)) setsigne(a, -1);
	return gerepileuptoint(av, a);
}


static void
lucasmod_tiny(ulong n, ulong m, ulong *a, ulong *b)
{
	ulong z, t;
	if (!n) { *a = 2; *b = 1; return; }
	lucasmod_tiny(n >> 1, m, &z, &t);
	switch(n & 3) {
		case  0: *a = z*z - 2; *b = z*t - 1; break;
		case  1: *a = z*t - 1; *b = t*t + 2; break;
		case  2: *a = z*z + 2; *b = z*t + 1; break;
		case  3: *a = z*t + 1; *b = t*t - 2;
	}
	*a = *a % m;
	*b = *b % m;
}


ulong
fibomod_tiny(long n, ulong m)
{
	ulong a, b;
	if (!n) return 0;
	lucasmod_tiny((ulong)(labs(n)-1), (ulong)(5*m), &a, &b);
	a = (((a << 1) + b) / 5) % m;
	
	return n < 0 && !odd(n) ? m-a : a;
}


GEN
fibmod(GEN n, GEN m)
{
	static const ulong fmod3[] = {0,1,1,2,0,2,2,1};
 	static const ulong fmod4[] = {0,1,1,2,3,1};
	static const ulong fmod5[] = {0,1,1,2,3,0,3,3,1,4,0,4,4,3,2,0,2,2,4,1};
	if (typ(n) != t_INT)
		pari_err_TYPE("fibmod", n);
	if (typ(m) != t_INT)
		pari_err_TYPE("fibmod", m);
	long nn = itos_or_0(n);
	if (nn) {
		ulong mm = itou_or_0(m);
#ifdef LONG_IS_64BIT
		if (mm && mm <= 858993459UL)
#else
		if (mm && mm <= 13107UL)
#endif
			return utoi(fibomod_tiny(nn, mm));
		return fibomod(nn, m);
	}
	if (!signe(n))
		return gen_0;

	GEN f, t, res;
	long l;
	pari_sp ltop = avma;
	if (cmpis(m, 6) < 0)
	{
		if (signe(m) < 0)
			m = negi(m);
		if (!signe(m))
			pari_err(e_INV);
			
		switch (itos_or_0(m)) {
		case 1:	
			avma = ltop;
			return gen_0;
		case 2:
			l = smodis(n, 3) > 0;
			avma = ltop;
			return utoipos(l);
		case 3:
			l = fmod3[mod8(n)];
			avma = ltop;
			return utoipos(l);
		case 4:
			l = fmod4[smodis(n, 6)];
			avma = ltop;
			return utoipos(l);
		case 5:
			l = fmod5[smodis(n, 20)];
			avma = ltop;
			return utoipos(l);
		}
	}
	// Rough estimate of the break-even point.	Yes, it's quite large.
	/*
	if (cmpis(n, 1000) < 0)
	{
		res = modii(fibo(itos(n)), m);
		res = gerepileupto(ltop, res);
		return res;
	}*/
	
	f = Z_factor(m);
	long e;
	t = gmodulo(gen_0, gen_1);
	l = glength(gel(f, 1));
	pari_sp btop = avma;
	long i;
	for (i = 1; i <= l; ++i)
	{
		e = itos(gcoeff(f, i, 2));
		t = chinese(t, gmodulo(fibo(smodis(n, Pisano(itos(gcoeff(f, i, 1)), e))), powis(gcoeff(f, i, 1), e)));
		t = gerepileupto(btop, t);
	}
	res = lift(t);
	res = gerepileupto(ltop, res);
	return res;
}


// Pisano periods of odd primes (with 0s elsewhere)
static const long PisanoArr[] = {
	0,8,20,16,0,10,28,0,36,18,0,48,0,0,14,30,0,0,76,0,40,88,0,32,0,0,108,0,0,58,
	60,0,0,136,0,70,148,0,0,78,0,168,0,0,44,0,0,0,196,0,50,208,0,72,108,0,76,0,
	0,0,0,0,0,256,0,130,0,0,276,46,0,0,0,0,148,50,0,0,316,0,0,328,0,336,0,0,348,
	0,0,178,90,0,0,0,0,190,388,0,396,22,0,0,0,0,0,42,0,0,0,0,0,448,0,456,114,0,
	52,0,0,238,240,0,0,0,0,250,0,0,516,0,0,176,0,0,268,270,0,0,556,0,56,568,0,0,
	0,0,588,0,0,0,0,0,0,88,0,310,628,0,636,0,0,0,0,0,0,110,0,0,676,0,0,0,0,232,
	174,0,236,0,0,358,0,0,0,736,0,0,748,0,0,378,0,768,0,0,388,0,0,0,796,0,200,0,
	0,0,408,0,0,0,0,418,84,0,0,0,0,430,868,0,0,438,0,888,0,0,448,0,0,0,916,0,46,
	928,0,936,0,0,0,0,0,478,0,0,0,976,0,490,0,0,0,498,0,1008,0,0,254,0,0,0,0,0,
	26,1048,0,0,0,0,0,0,0,0,90,0,0,1096,0,0,0,0,124,0,0,376,0,0,568,570,0,0,
	1156,0,0,0,0,1176,0,0,1188,0,0,598,600,0,0,1216,0,0,1228,0,1236,206,0,0,0,0,
	0,630,0,0,0,0,640,1288,0,1296,0,0,1308,0,0,658,220,0,0,0,0,0,1348,0,452,0,0,
	1368,0,0,0,138,0,0,0,0,700,0,0,0,118,0,0,0,0,718,0,0,0,1456,0,0,1468,0,0,
	738,0,496,0,0,0,750,0,0,1516,0,380,0,0,0,192,0,1548,0,0,0,0,0,0,1576,0,0,0,
	0,228,0,0,0,0,0,202,270,0,0,0,0,820,1648,0,1656,276,0,0,0,0,838,0,0,0,0,0,0,
	1708,0,1716,78,0,1728,0,0,0,0,0,0,1756,0,176,1768,0,1776,0,0,0,0,0,0,0,0,0,
	1816,0,70,0,0,0,102,0,0,0,0,928,0,0,0,1876,0,470,0,0,1896,0,0,212,0,0,0,0,0,
	0,176,0,970,0,0,652,0,0,1968,0,0,0,198,0,0,1996
};


/* Helper function for fibmod; returns a multiple of the period of Fibonacci numbers mod p^e. */
/* Assumes p is prime and e > 0. */
long
Pisano(long p, long e)
{
	long p1;
	if (p < 999) {
		if (p == 5)
			return (long)upowuu(5, (ulong)e) << 2;
		e--;
		if (p == 2)
			return 3 << e;
		p1 = PisanoArr[p>>1];
		return e ? p1 * (long)upowuu(p, e) : p1;
	}
	long t = p%5;
	if (t == 2 || t == 3)
		p1 = (p + 1) << 1;
	else
		p1 = p - 1;
	return e == 1 ? p1 : p1 * (long)upowuu(p, e - 1);
}


INLINE long
hamming_word(ulong w)
{
	return __builtin_popcountl(w);
#if 0
static const long byte_weight[] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
	long sum = 0;
	while (w) {
		//sum += w & 1;
		//w >>= 1;

		//sum++;
		//w &= w - 1;

		sum += byte_weight[w & 255];
		w >>= 8;
	}
	return sum;
#endif
}


long
hamming(GEN n)
{
	pari_sp ltop = avma;
	if (typ(n) != t_INT)
		pari_err_TYPE("hamming", n);
	if (!signe(n))
		return 0;

	GEN xp = int_MSW(n);
	long lx = lgefint(n);
	ulong u = *xp;
	long sum = 0;
	long i = 3;
	for (; i < lx; ++i)
	{
		sum += hamming_word(u);
		xp=int_precW(xp);
		u = *xp;
	}
	avma = ltop;
	return sum + hamming_word(u);
}


long
istwo(GEN n)
{
	// TODO: Serious improvements available with incremental factoring
	pari_sp ltop = avma;
	GEN f = gen_0;
	long ret, l;
	if (typ(n) != t_INT)
		pari_err_TYPE("istwo", n);
	
	// Remove factors of 2
	long tmp = vali(n);
	if (tmp)
		n = shifti(n, -tmp);
	
	// Remove negatives
	if (cmpis(n, 3) < 0)
	{
		ret = signe(n) >= 0;
		avma = ltop;
		return ret;
	}

	// End early if possible: numbers that are 3 mod 4 always have some prime factor 3 mod 4 raised to an odd power.
	if (mod4(n) == 3)
		return 0;
	
	f = Z_factor(n);
	l = glength(gel(f, 1));
	long i;
	for (i = 1; i <= l; ++i)
	{
		if (mod4(gcoeff(f, i, 1)) == 3 && mod2(gcoeff(f, i, 2)))
		{
			avma = ltop;
			return 0;
		}
	}
	avma = ltop;
	return 1;
}


GEN
ways2(GEN n)
{
	// TODO: Serious improvements available with incremental factoring
	if (typ(n) != t_INT)
		pari_err_TYPE("ways2", n);

	pari_sp ltop = avma;
	GEN res = gen_1, f = factor(oddres(n));
	long l1 = glength(gel(f, 1));
	pari_sp btop = avma;
	long i;
	for (i = 1; i <= l1; ++i)
	{
		if (mod4(gcoeff(f, i, 1)) == 3)
		{
			if (mod2(gcoeff(f, i, 2)))
			{
				avma = ltop;
				return gen_0;
			}
		} else {
			res = mulii(res, addis(gcoeff(f, i, 2), 1));
		}
		res = gerepileupto(btop, res);
	}
	res = shifti(addis(res, 1), -1);
	res = gerepileuptoleaf(ltop, res);
	return res;
}


long
isthree(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("isthree", n);
	if (signe(n) < 0)
		return 0;

	long tmp = vali(n);
	if (tmp & 1)
		return 1;
	pari_sp ltop = avma;
	long l1 = mod8(shifti(n, -tmp)) != 7;
	avma = ltop;
	return l1;
}


long
sways3s(ulong n)
{
	pari_sp ltop = avma;
	long p1 = (long)usqrt(n);
	long res = 0;
	long k, t, m, p2, tmod4;
	for (k = 0; k <= p1; k++)
	{
		tmod4 = (n - (k&1))&3;
		if (tmod4 == 3)
			continue;
		t = n - k * k;
		if (!istwo(stoi(t)))
			continue;
		avma = ltop;
		p2 = (long)usqrt(t>>1);
		switch (tmod4) {
			case 0:
				for (m = ((k + 1)|1) - 1; m <= p2; m += 2)
					if (uissquare(t - m * m))
						res++;
				break;
			case 1:
				for (m = k; m <= p2; m++)
					if (uissquare(t - m * m))
						res++;
				break;
			case 2:
				for (m = k|1; m <= p2; m += 2)
					if (uissquare(t - m * m))
						res++;
				break;
		}
	}
	return res;
}


GEN
ways3(GEN n)
{
	pari_sp ltop = avma;
	GEN p1, t, res;
	if (typ(n) != t_INT)
		pari_err_TYPE("ways3", n);
	if (signe(n) < 1)
		return signe(n) ? gen_0 : gen_1;
	
	// ways3(4n) = ways3(n); ways3(4^k(8n + 7)) = 0
	if (!mod4(n)) {
		n = shifti(n, -(vali(n) >> 1) << 1);
	}
	if (mod8(n) == 7) {
		avma = ltop;
		return gen_0;
	}
	ulong small = itou_or_0(n);
	if (small) {
		avma = ltop;
		return stoi(sways3s(small));
	}
	
	p1 = sqrtint(truedivis(n, 3));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k, p3, p4;
	res = gen_0;
	for (k = gen_0; cmpii(k, p1) <= 0; k = addis(k, 1))
	{
		t = subii(n, sqri(k));
		p3 = sqrtint(shifti(t, -1));
		pari_sp ctop = avma, c_lim = stack_lim(ctop, 1);
		GEN m;
		p4 = gen_0;
		for (m = k; cmpii(m, p3) <= 0; m = addis(m, 1))
		{
			if (Z_issquare(subii(t, gsqr(m))))
				p4 = addis(p4, 1);
			if (low_stack(c_lim, stack_lim(ctop, 1)))
				gerepileall(ctop, 2, &p4, &m);
		}
		res = addii(res, p4);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 4, &res, &k, &p3, &p4);
			//gerepileall(btop, 2, &res, &k);	// This is significantly (15%+) slower -- I don't know why.
	}
	res = gerepileupto(ltop, res);
	return res;
}


// Looking for A000523?  Try expi.  A029837 and A070939 are similar.
// See also __builtin_clzl (also __builtin_ffsl for the lsb)
GEN
msb(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("msb", n);
	if (signe(n) < 1) {
		if (signe(n))
			pari_err_DOMAIN("msb", "n", "<", gen_0, n);
		return gen_0;	// Convention from A053644
	}

	return int2n(expi(n));
}


GEN
fusc(GEN n)
{
	if (typ(n) != t_INT)
		pari_err_TYPE("fusc", n);
	if (signe(n) < 1)
		return gen_0;
	pari_sp ltop = avma;
	GEN ret;
	
	// Note: Different constants depending on word size!
#ifdef LONG_IS_64BIT
	long l = lgefint(n);
	if (l < 4 || (l == 4 && (ulong)(*int_MSW(n)) <= 13153337344UL))
#else
	if (cmpii(n, u2toi(9386, 2863311530UL)) <= 0)
#endif
		ret = utoi(fusc_small(n));
	else
		ret = fusc_large(n);
	ret = gerepileupto(ltop, ret);
	return ret;
}


static const ulong fuscAA[] ={
	1,8,7,13,6,17,11,16,5,19,14,23,9,22,13,17,4,19,15,26,11,29,18,25,7,24,17,27,10,23,13,16,3,17,14,25,11,30,19,27,8,29,21,34,13,31,18,23,5,22,17,29,12,31,19,26,7,23,16,25,9,20,11,13,2,13,11,20,9,25,16,23,7,26,19,31,12,29,17,22,5,23,18,31,13,34,21,29,8,27,19,30,11,25,14,17,3,16,13,23,10,27,17,24,7,25,18,29,11,26,15,19,4,17,13,22,9,23,14,19,5,16,11,17,6,13,7,8,1,7,6,11,5,14,9,13,4,15,11,18,7,17,10,13,3,14,11,19,8,21,13,18,5,17,12,19,7,16,9,11,2,11,9,16,7,19,12,17,5,18,13,21,8,19,11,14,3,13,10,17,7,18,11,15,4,13,9,14,5,11,6,7,1,6,5,9,4,11,7,10,3,11,8,13,5,12,7,9,2,9,7,12,5,13,8,11,3,10,7,11,4,9,5,6,1,5,4,7,3,8,5,7,2,7,5,8,3,7,4,5,1,4,3,5,2,5,3,4,1,3,2,3,1,2,1,1
};
static const ulong fuscAB[] ={
	8,7,13,6,17,11,16,5,19,14,23,9,22,13,17,4,19,15,26,11,29,18,25,7,24,17,27,10,23,13,16,3,17,14,25,11,30,19,27,8,29,21,34,13,31,18,23,5,22,17,29,12,31,19,26,7,23,16,25,9,20,11,13,2,13,11,20,9,25,16,23,7,26,19,31,12,29,17,22,5,23,18,31,13,34,21,29,8,27,19,30,11,25,14,17,3,16,13,23,10,27,17,24,7,25,18,29,11,26,15,19,4,17,13,22,9,23,14,19,5,16,11,17,6,13,7,8,1,7,6,11,5,14,9,13,4,15,11,18,7,17,10,13,3,14,11,19,8,21,13,18,5,17,12,19,7,16,9,11,2,11,9,16,7,19,12,17,5,18,13,21,8,19,11,14,3,13,10,17,7,18,11,15,4,13,9,14,5,11,6,7,1,6,5,9,4,11,7,10,3,11,8,13,5,12,7,9,2,9,7,12,5,13,8,11,3,10,7,11,4,9,5,6,1,5,4,7,3,8,5,7,2,7,5,8,3,7,4,5,1,4,3,5,2,5,3,4,1,3,2,3,1,2,1,1,0
};
static const ulong fuscBA[] ={
	0,1,1,2,1,3,2,3,1,4,3,5,2,5,3,4,1,5,4,7,3,8,5,7,2,7,5,8,3,7,4,5,1,6,5,9,4,11,7,10,3,11,8,13,5,12,7,9,2,9,7,12,5,13,8,11,3,10,7,11,4,9,5,6,1,7,6,11,5,14,9,13,4,15,11,18,7,17,10,13,3,14,11,19,8,21,13,18,5,17,12,19,7,16,9,11,2,11,9,16,7,19,12,17,5,18,13,21,8,19,11,14,3,13,10,17,7,18,11,15,4,13,9,14,5,11,6,7,1,8,7,13,6,17,11,16,5,19,14,23,9,22,13,17,4,19,15,26,11,29,18,25,7,24,17,27,10,23,13,16,3,17,14,25,11,30,19,27,8,29,21,34,13,31,18,23,5,22,17,29,12,31,19,26,7,23,16,25,9,20,11,13,2,13,11,20,9,25,16,23,7,26,19,31,12,29,17,22,5,23,18,31,13,34,21,29,8,27,19,30,11,25,14,17,3,16,13,23,10,27,17,24,7,25,18,29,11,26,15,19,4,17,13,22,9,23,14,19,5,16,11,17,6,13,7,8
};
static const ulong fuscBB[] ={
	1,1,2,1,3,2,3,1,4,3,5,2,5,3,4,1,5,4,7,3,8,5,7,2,7,5,8,3,7,4,5,1,6,5,9,4,11,7,10,3,11,8,13,5,12,7,9,2,9,7,12,5,13,8,11,3,10,7,11,4,9,5,6,1,7,6,11,5,14,9,13,4,15,11,18,7,17,10,13,3,14,11,19,8,21,13,18,5,17,12,19,7,16,9,11,2,11,9,16,7,19,12,17,5,18,13,21,8,19,11,14,3,13,10,17,7,18,11,15,4,13,9,14,5,11,6,7,1,8,7,13,6,17,11,16,5,19,14,23,9,22,13,17,4,19,15,26,11,29,18,25,7,24,17,27,10,23,13,16,3,17,14,25,11,30,19,27,8,29,21,34,13,31,18,23,5,22,17,29,12,31,19,26,7,23,16,25,9,20,11,13,2,13,11,20,9,25,16,23,7,26,19,31,12,29,17,22,5,23,18,31,13,34,21,29,8,27,19,30,11,25,14,17,3,16,13,23,10,27,17,24,7,25,18,29,11,26,15,19,4,17,13,22,9,23,14,19,5,16,11,17,6,13,7,8,1
};


#define fusc8bits(a, b, idx) {int i = (idx) & 0xFF; int newA = a*fuscAA[i] + b*fuscAB[i]; b = a*fuscBA[i] + b*fuscBB[i]; a = newA;}
static void
fusc_word(ulong u, ulong *a, ulong *b) {
	*a = fuscAA[u & 0xFF];
	*b = 0;
	fusc8bits(*a, *b, u >> 8)
	fusc8bits(*a, *b, u >> 16)
	fusc8bits(*a, *b, u >> 24)
#ifdef LONG_IS_64BIT
	fusc8bits(*a, *b, u >> 32)
	fusc8bits(*a, *b, u >> 40)
	fusc8bits(*a, *b, u >> 48)
	fusc8bits(*a, *b, u >> 56)
#endif
}


ulong
fusc_small(GEN n)
{
	ulong a = 1, b = 0;
	pari_sp ltop = avma;
	
	GEN xp = int_LSW(n);
	
	// If n has two words, handle the least-significant one (otherwise it has
	// only one word).
	if (lgefint(n) > 3)
	{
		fusc_word(*xp, &a, &b);
		xp=int_nextW(xp);
	}
	
	ulong tmp = *xp;
	while (tmp)
	{
		fusc8bits(a, b, tmp)
		tmp >>= 8;
	}
	
	avma = ltop;
	return b;
	#undef fusc8bits
}


GEN
fusc_large(GEN n)
{
	pari_sp ltop = avma;
	GEN a = gen_1, b = gen_0;
	#define fusc8bits(a, b, idx) {ulong i = (idx) & 0xFF; GEN newA = addii(muliu(a, fuscAA[i]), muliu(b, fuscAB[i])); b = addii(muliu(a, fuscBA[i]), muliu(b, fuscBB[i])); a = newA;}
	GEN xp = int_LSW(n);
	int len = lgefint(n);
	
	while(1)
	{
		ulong u = *xp;
		// TODO: Rather than work with large numbers directly, compute the
		// (wordsize) coefficients for 1 word's worth of bits then do the
		// GEN multiplications all at once.
		fusc8bits(a, b, u)
		fusc8bits(a, b, u >> 8)
		fusc8bits(a, b, u >> 16)
		fusc8bits(a, b, u >> 24)
#ifdef LONG_IS_64BIT
		fusc8bits(a, b, u >> 32)
		fusc8bits(a, b, u >> 40)
		fusc8bits(a, b, u >> 48)
		fusc8bits(a, b, u >> 56)
#endif
		if (--len < 3)
			break;
		xp=int_nextW(xp);
	}
	b = gerepileupto(ltop, b);
	return b;
	#undef fusc8bits
}


// Evaluate the integer polynomial x at integer y.
// Also accepts a t_VEC or t_COL as Polrev(x).
// Optimized for dense polynomials using Horner's rule.
GEN
poleval_denseint(GEN x, GEN y)
{
	long i = lg(x)-1, imin = (typ(x) == t_POL) ? 2 : 1;
	if (i<=imin)
		return (i==imin) ? icopy(gel(x,imin)) : gen_0;

	pari_sp ltop = avma;
	GEN ret = gel(x,i);
	for (i-- ; i>=imin; i--)
		ret = addii(mulii(ret,y),gel(x,i));
	return gerepileupto(ltop,ret);
}
