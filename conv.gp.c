/******************************************************************************/
/**														 Convenience																	**/
/******************************************************************************/

INLINE GEN
gtor(GEN x, const char* funcName, long prec)
{
	switch (typ(x)) {
		case t_REAL:
			return x;	// x, not a copy of x
		case t_INT:
		case t_FRAC:
			return cxcompotor(x, prec);
		default:
			pari_err(typeer, funcName);
	}
	return NEVER_USED;
}


GEN
vecsum(GEN v)
{
	pari_sp ltop = avma;
	GEN p2 = gen_0;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "vecsum");
	long l1 = lg(v);
	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		p2 = gadd(p2, gel(v, i));
		p2 = gerepileupto(btop, p2);
	}
	p2 = gerepileupto(ltop, p2);
	return p2;
}


// TODO: Binary splitting for smaller subproducts.
GEN
vecprod(GEN v)
{
	pari_sp ltop = avma;
	GEN p2 = gen_1;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "vecprod");
	long l1 = lg(v);

	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		p2 = gmul(p2, gel(v, i));
		p2 = gerepileupto(btop, p2);
	}
	p2 = gerepileupto(ltop, p2);
	return p2;
}


GEN
vecgcd(GEN v)
{
	pari_sp ltop = avma;
	GEN l = gen_0;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "veclcm");
	long l1 = lg(v);
	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		l = ggcd(l, gel(v, i));
		l = gerepileupto(btop, l);
	}
	l = gerepileupto(ltop, l);
	return l;
}


GEN
veclcm(GEN v)
{
	pari_sp ltop = avma;
	GEN l = gen_1;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "veclcm");
	long l1 = lg(v);
	long i;
	for (i = 1; i < l1; ++i)
	{
		l = glcm(l, gel(v, i));
		l = gerepileupto(ltop, l);
	}
	return l;
}


// TODO (possibly): give 32- and 64-bit specific code if it would differ, since
// sometimes code needs to be split anyway (or only one is relevant, etc.).
// TODO: Handle negatives, and possibly non-integer types?
void
toC(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(typeer, "toC");
	if (cmpis(n, 3) < 0)
	{
		if (cmpis(n, 2) == 0)
			pari_printf("gen_2\n");
		else if (equali1(n))
			pari_printf("gen_1\n");
		else if (signe(n) == 0)
			pari_printf("gen_0\n");
		else if (equalim1(n))
			pari_printf("gen_m1\n");
		else
			pari_err(alarmer, "can't handle negatives yet"); // white lie
		return;
	}
	if (ispow2(n)) {
		pari_printf("int2n(%ld)\n", expi(n));
		return;
	}
	
	pari_sp ltop = avma;
	GEN t;
	// words: number of 32-bit words in n.
#ifdef LONG_IS_64BIT
	long words = (lgefint(n) - 2) << 1;
	if ((ulong)*int_MSW(n) <= 0xFFFFFFFF)
		words--;
#else
	long words = lgefint(n) - 2;
#endif

	if (words == 1)
	{
		pari_printf("utoipos(%Ps)\n", n);
		avma = ltop;
		return;
	}
	if (words == 2)
	{
		pari_printf("uu32toi(%Ps, %Ps)\n", shifti(n, -32), remi2n(n, 32));
		avma = ltop;
		return;
	}

	// Large numbers
	// If efficiency mattered, walking through the binary representation
	// would be far more efficient.
	pari_printf("mkintn(%Ps", stoi(words));
	long i = words - 1;
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	for (; i >= 0; i--)
	{
		t = shifti(n, -(i * 32));
		pari_printf(", %Ps", t);
		n = subii(n, shifti(t, i * 32));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &n);
	}
	pari_printf(")\n");
	avma = ltop;
	return;
}


long
digits(GEN x)
{
	pari_sp av = avma;
	long s = sizedigit(x) - 1;
	if (gcmp(x, powis(stoi(10), s)) >= 0)
		s++;
	avma = av;
	return s;
}


GEN
eps(long prec)
{
	GEN ret = real_1(DEFAULTPREC);
	setexpo(ret, 1 - bit_accuracy(prec));
	return ret;
}
