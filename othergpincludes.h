////////////////////////////////////////////////////////////////////////////////
#ifndef S_SPLINT_S
// A private function from language/sumiter.c
static byteptr
prime_loop_init(GEN ga, GEN gb, ulong *a, ulong *b, ulong *p)
{
  byteptr d = diffptr;

  ga = gceil(ga); gb = gfloor(gb);
  if (typ(ga) != t_INT || typ(gb) != t_INT)
    pari_err(typeer,"prime_loop_init");
  if (signe(gb) < 0) return NULL;
  if (signe(ga) < 0) ga = gen_1;
  if (lgefint(ga)>3 || lgefint(gb)>3)
  {
    if (cmpii(ga, gb) > 0) return NULL;
    pari_err(primer1, 0);
  }
  *a = itou(ga);
  *b = itou(gb); if (*a > *b) return NULL;
  maxprime_check(*b);
  *p = init_primepointer(*a, 0, &d); return d;
}


// A private function from basemath/prime.c
static ulong
u_LucasMod(ulong n, ulong P, ulong N)
{
  long j = 1 + bfffo(n);
  ulong v = P, v1 = P*P - 2, mP = N - P, m2 = N - 2, m = n << j;

  j = BITS_IN_LONG - j;
  for (; j; m<<=1,j--)
  { /* v = v_k, v1 = v_{k+1} */
    if (((long)m) < 0)
    { /* set v = v_{2k+1}, v1 = v_{2k+2} */
      v = Fl_add(Fl_mul(v,v1,N), mP, N);
      v1= Fl_add(Fl_mul(v1,v1,N),m2, N);
    }
    else
    {/* set v = v_{2k}, v1 = v_{2k+1} */
      v1= Fl_add(Fl_mul(v,v1,N),mP, N);
      v = Fl_add(Fl_mul(v,v,N), m2, N);
    }
  }
  return v;
}


// A private function from basemath/prime.c
static int
u_IsLucasPsP(ulong n)
{
  long i, v;
  ulong b, z, m2, m = n + 1;

  for (b=3, i=0;; b+=2, i++)
  {
    ulong c = b*b - 4; /* = 1 mod 4 */
    if (krouu(n % c, c) < 0) break;
    if (i == 64 && uissquareall(n, &c)) return 0; /* oo loop if N = m^2 */
  }
  if (!m) return 0; /* neither 2^32-1 nor 2^64-1 are Lucas-pp */
  v = vals(m); m >>= v;
  z = u_LucasMod(m, b, n);
  if (z == 2) return 1;
  m2 = n - 2;
  if (z == m2) return 1;
  for (i=1; i<v; i++)
  {
    if (!z) return 1;
    z = Fl_add(Fl_mul(z,z, n), m2, n);
    if (z == 2) return 0;
  }
  return 0;
}


// A private function from gen2.c
/* x,y t_POL */
static int
polequal(GEN x, GEN y)
{
  long lx, ly;
  if (x[1] != y[1]) return 0;
  lx = lg(x); ly = lg(y);
  while (lx > ly) if (!gequal0(gel(x,--lx))) return 0;
  while (ly > lx) if (!gequal0(gel(y,--ly))) return 0;
  for (lx--; lx >= 2; lx--) if (!gequal(gel(x,lx), gel(y,lx))) return 0;
  return 1;
}


// A private function from arith1.c
static int
is_char_2(GEN a)
{
  long j;
  GEN b;
  switch(typ(a))
  {
  case t_INTMOD:
    b = gel(a,1);
    if (!mod2(b))
    {
      if (!equaliu(b, 2)) pari_err(impl, "issquare for this input");
      return 1;
    }
    return 0;
  case t_FFELT:
    if (equaliu(FF_p_i(a), 2)) return 1;
    return 0;
  case t_POLMOD:
    if (is_char_2(gel(a,1)) || is_char_2(gel(a,2))) return 1;
    return 0;
  case t_POL:
    for (j = 2; j < lg(a); j++)
      if (is_char_2(gel(a,j))) return 1;
    return 0;
  }
  return 0;
}


// A private function from arith1.c
static long
polissquareall(GEN x, GEN *pt)
{
  pari_sp av;
  long v, l = degpol(x);
  GEN y, a, b;

  if (!signe(x))
  {
    if (pt) *pt = gcopy(x);
    return 1;
  }
  if (pt) *pt = gen_0;
  if (l&1) return 0; /* odd degree */
  av = avma;
  v = RgX_valrem(x, &x);
  if (v) {
    l = degpol(x);
    if (l&1) return 0;
  }
  a = gel(x,2);
  switch (typ(a))
  {
    case t_INT: y =  Z_issquareall(a,&b)? gen_1: gen_0; break;
    case t_POL: y = polissquareall(a,&b)? gen_1: gen_0; break;
    default: y = gissquare(a); b = NULL; break;
  }
  if (y == gen_0) { avma = av; return 0; }
  if (!l) {
    if (!pt) { avma = av; return 1; }
    if (!b) b = gsqrt(a,DEFAULTPREC);
    y = scalarpol(b, varn(x)); goto END;
  }
  if (is_char_2(x))
  {
    long i, lx;
    x = gmul(x, mkintmod(gen_1, gen_2));
    lx = lg(x);
    if ((lx-3) & 1) { avma = av; return 0; }
    for (i = 3; i < lx; i+=2)
      if (!gequal0(gel(x,i))) { avma = av; return 0; }
    if (pt) {
      y = cgetg((lx+3) / 2, t_POL);
      for (i = 2; i < lx; i+=2)
        if (!gissquareall(gel(x,i), &gel(y, (i+2)>>1))) { avma = av; return 0; }
      y[1] = evalsigne(1) | evalvarn(varn(x));
      goto END;
    } else {
      for (i = 2; i < lx; i+=2)
        if (!gissquare(gel(x,i))) { avma = av; return 0; }
      avma = av; return 1;
    }
  }
  else
  {
    x = RgX_Rg_div(x,a);
    y = gtrunc(gsqrt(RgX_to_ser(x,2+l),0));
    if (!RgX_equal(gsqr(y), x)) { avma = av; return 0; }
    if (!pt) { avma = av; return 1; }
    if (!gequal1(a))
    {
      if (!b) b = gsqrt(a,DEFAULTPREC);
      y = gmul(b, y);
    }
  }
END:
  *pt = v? gerepilecopy(av, RgX_shift_shallow(y, v >> 1)): gerepileupto(av, y);
  return 1;
}
#endif
////////////////////////////////////////////////////////////////////////////////
