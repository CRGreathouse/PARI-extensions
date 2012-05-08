////////////////////////////////////////////////////////////////////////////////
#ifndef S_SPLINT_S
// A private function from language/sumiter.c
static byteptr
prime_loop_init(GEN ga, GEN gb, ulong *a, ulong *b, ulong *p)
{
  byteptr d = diffptr;

  ga = gceil(ga); gb = gfloor(gb);
  if (typ(ga) != t_INT) pari_err_TYPE("prime_loop_init",ga);
  if (typ(gb) != t_INT) pari_err_TYPE("prime_loop_init",gb);
  if (signe(gb) < 0) return NULL;
  if (signe(ga) < 0) ga = gen_1;
  if (lgefint(ga)>3 || lgefint(gb)>3)
  {
    if (cmpii(ga, gb) > 0) return NULL;
    pari_err_MAXPRIME(0);
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

#endif
////////////////////////////////////////////////////////////////////////////////
