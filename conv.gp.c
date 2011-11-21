/******************************************************************************/
/**														 Convenience																	**/
/******************************************************************************/

GEN
facmod(long n, long m)
{
  pari_sp ltop = avma;
  long t;
  GEN s = gen_0, s1 = gen_0;
  long e, sq, l1, l2;
  GEN p3 = gen_0;
  s = gmodulss(1, m);
  sq = itos(sqrtint(stoi(n)));
  l1 = n/sq;
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    long p = 0;
    byteptr primepointer = diffptr;	  /* bptr */
    if (l1 > maxprime())
      pari_err_MAXPRIME(l1);
    for (;;)
    {
      NEXT_PRIME_VIADIFF(p, primepointer);
      if (p > l1)
        break;
      t = n;
      e = 0;
      {
        pari_sp btop = avma;
        while (t)
        {
          e += (t /= p);
          avma = btop;
        }
      }
      s = gmul(s, gpowgs(gmodulss(p, m), e));
      if (low_stack(st_lim, stack_lim(btop, 1)))
        s = gerepilecopy(btop, s);
    }
  }
  l2 = sq - 1;
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    long j;
    long l4 = -1 > 0;	  /* bool */
    for (j = l2; l4?j <= 1:j >= 1; j += -1)
    {
      {
        long l5, l6;
        s1 = gmodulss(1, m);
        l5 = (n/(j + 1)) + 1;
        l6 = n/j;
        {
          pari_sp btop = avma, st_lim = stack_lim(btop, 1);
          long p = 0;
          byteptr primepointer = diffptr;	  /* bptr */
          if (l6 > maxprime())
            pari_err_MAXPRIME(l6);
          for (;;)
          {
            NEXT_PRIME_VIADIFF(p, primepointer);
            if (p > l6)
              break;
            if (p < l5)
              continue;
            s1 = gmulgs(s1, p);
            if (low_stack(st_lim, stack_lim(btop, 1)))
              s1 = gerepilecopy(btop, s1);
          }
        }
        s = gmul(s, gpowgs(s1, j));
      }
      if (low_stack(st_lim, stack_lim(btop, 1)))
        gerepileall(btop, 2, &s, &s1);
    }
  }
  p3 = lift(s);
  p3 = gerepilecopy(ltop, p3);
  return p3;
}


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
			pari_err_TYPE(funcName, x);
	}
	return NEVER_USED;
}


GEN
vecsum(GEN v)
{
	pari_sp ltop = avma;
	GEN p2 = gen_0;
	if (!is_matvec_t(typ(v)))
		pari_err_TYPE("vecsum", v);
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
		pari_err_TYPE("vecprod", v);
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
		pari_err_TYPE("veclcm", v);
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
		pari_err_TYPE("veclcm", v);
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
		pari_err_TYPE("toC", n);
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
			pari_err_IMPL("negatives in toC"); // white lie
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
