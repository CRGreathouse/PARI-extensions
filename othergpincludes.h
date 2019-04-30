////////////////////////////////////////////////////////////////////////////////
#define VALUE(x) gel(x,0)
#define EXPON(x) gel(x,1)
#define CLASS(x) gel(x,2)

INLINE long hamming_word(ulong w) __attribute__ ((pure));
static void ensure_nb(GEN L, long l);
static GEN pollardbrent(GEN n);
static ulong utridiv_bound(ulong n);
INLINE ulong u_forprime_next_fast(forprime_t *T);
static void one_iter(GEN *x, GEN *P, GEN x1, GEN n, long delta);
static void rho_dbg(pari_timer *T, long c, long msg_mask);
INLINE void INIT(GEN x, GEN v, GEN e, GEN c);


INLINE void
INIT(GEN x, GEN v, GEN e, GEN c) {
  VALUE(x) = v;
  EXPON(x) = e;
  CLASS(x) = c;
}


static void
rho_dbg(pari_timer *T, long c, long msg_mask)
{
  if (c & msg_mask) return;
  err_printf("Rho: time = %6ld ms,\t%3ld round%s\n",
             timer_delay(T), c, (c==1?"":"s"));
  err_flush();
}


static void
one_iter(GEN *x, GEN *P, GEN x1, GEN n, long delta)
{
  *x = addis(remii(sqri(*x), n), delta);
  *P = modii(mulii(*P, subii(x1, *x)), n);
}


static GEN
pollardbrent(GEN n)
{
  const long tune_pb_min = 14; /* even 15 seems too much. */
  long tf = lgefint(n), size = 0, delta, retries = 0, msg_mask;
  long c0, c, k, k1, l;
  pari_sp av;
  GEN x, x1, y, P, g, g1, res;
  pari_timer T;

  if (DEBUGLEVEL >= 4) timer_start(&T);

  if (tf >= 4)
    size = expi(n) + 1;
  else if (tf == 3) /* keep purify happy */
    size = 1 + expu(uel(n,2));

  if (size <= 28)
    c0 = 32;/* amounts very nearly to 'insist'. Now that we have squfof(), we
             * don't insist any more when input is 2^29 ... 2^32 */
  else if (size <= 42)
    c0 = tune_pb_min;
  else if (size <= 59) /* match squfof() cutoff point */
    c0 = tune_pb_min + ((size - 42)<<1);
  else if (size <= 72)
    c0 = tune_pb_min + size - 24;
  else if (size <= 301)
    /* nonlinear increase in effort, kicking in around 80 bits */
    /* 301 gives 48121 + tune_pb_min */
    c0 = tune_pb_min + size - 60 +
      ((size-73)>>1)*((size-70)>>3)*((size-56)>>4);
  else
    c0 = 49152; /* ECM is faster when it'd take longer */

  c = c0 << 5; /* 2^5 iterations per round */
  msg_mask = (size >= 448? 0x1fff:
                           (size >= 192? (256L<<((size-128)>>6))-1: 0xff));
  y = cgeti(tf);
  x1= cgeti(tf);
  av = avma;

PB_RETRY:
 /* trick to make a 'random' choice determined by n.  Don't use x^2+0 or
  * x^2-2, ever.  Don't use x^2-3 or x^2-7 with a starting value of 2.
  * x^2+4, x^2+9 are affine conjugate to x^2+1, so don't use them either.
  *
  * (the point being that when we get called again on a composite cofactor
  * of something we've already seen, we had better avoid the same delta) */
  switch ((size + retries) & 7)
  {
    case 0:  delta=  1; break;
    case 1:  delta= -1; break;
    case 2:  delta=  3; break;
    case 3:  delta=  5; break;
    case 4:  delta= -5; break;
    case 5:  delta=  7; break;
    case 6:  delta= 11; break;
    /* case 7: */
    default: delta=-11; break;
  }
  if (DEBUGLEVEL >= 4)
  {
    if (!retries)
      err_printf("Rho: searching small factor of %ld-bit integer\n", size);
    else
      err_printf("Rho: restarting for remaining rounds...\n");
    err_printf("Rho: using X^2%+1ld for up to %ld rounds of 32 iterations\n",
               delta, c >> 5);
    err_flush();
  }
  x = gen_2; P = gen_1; g1 = NULL; k = 1; l = 1;
  affui(2, y);
  affui(2, x1);
  for (;;) /* terminated under the control of c */
  { /* use the polynomial  x^2 + delta */
    one_iter(&x, &P, x1, n, delta);

    if ((--c & 0x1f)==0)
    { /* one round complete */
      g = gcdii(n, P); if (!is_pm1(g)) goto fin;
      if (c <= 0)
      { /* getting bored */
        if (DEBUGLEVEL >= 4)
        {
          err_printf("Rho: time = %6ld ms,\tPollard-Brent giving up.\n",
                     timer_delay(&T));
          err_flush();
        }
        return NULL;
      }
      P = gen_1;
      if (DEBUGLEVEL >= 4) rho_dbg(&T, c0-(c>>5), msg_mask);
      affii(x,y); x = y; avma = av;
    }

    if (--k) continue; /* normal end of loop body */

    if (c & 0x1f) /* otherwise, we already checked */
    {
      g = gcdii(n, P); if (!is_pm1(g)) goto fin;
      P = gen_1;
    }

   /* Fast forward phase, doing l inner iterations without computing gcds.
    * Check first whether it would take us beyond the alloted time.
    * Fast forward rounds count only half (although they're taking
    * more like 2/3 the time of normal rounds).  This to counteract the
    * nuisance that all c0 between 4096 and 6144 would act exactly as
    * 4096;  with the halving trick only the range 4096..5120 collapses
    * (similarly for all other powers of two) */
    if ((c -= (l>>1)) <= 0)
    { /* got bored */
      if (DEBUGLEVEL >= 4)
      {
        err_printf("Rho: time = %6ld ms,\tPollard-Brent giving up.\n",
                   timer_delay(&T));
        err_flush();
      }
      return NULL;
    }
    c &= ~0x1f; /* keep it on multiples of 32 */

    /* Fast forward loop */
    affii(x, x1); avma = av; x = x1;
    k = l; l <<= 1;
    /* don't show this for the first several (short) fast forward phases. */
    if (DEBUGLEVEL >= 4 && (l>>7) > msg_mask)
    {
      err_printf("Rho: fast forward phase (%ld rounds of 64)...\n", l>>7);
      err_flush();
    }
    for (k1=k; k1; k1--)
    {
      one_iter(&x, &P, x1, n, delta);
      if ((k1 & 0x1f) == 0) gerepileall(av, 2, &x, &P);
    }
    if (DEBUGLEVEL >= 4 && (l>>7) > msg_mask)
    {
      err_printf("Rho: time = %6ld ms,\t%3ld rounds, back to normal mode\n",
                 timer_delay(&T), c0-(c>>5));
      err_flush();
    }
    affii(x,y); avma = av; x = y;
  } /* forever */

fin:
  /* An accumulated gcd was > 1 */
  if  (!equalii(g,n))
  { /* if it isn't n, and looks prime, return it */
    if (MR_Jaeschke(g))
    {
      if (DEBUGLEVEL >= 4)
      {
        rho_dbg(&T, c0-(c>>5), 0);
        err_printf("\tfound factor = %Ps\n",g);
        err_flush();
      }
      return g;
    }
    avma = av; g1 = icopy(g);  /* known composite, keep it safe */
    av = avma;
  }
  else g1 = n; /* and work modulo g1 for backtracking */

  /* Here g1 is known composite */
  if (DEBUGLEVEL >= 4 && size > 192)
  {
    err_printf("Rho: hang on a second, we got something here...\n");
    err_flush();
  }
  x = y;
  for(;;)
  { /* backtrack until period recovered. Must terminate */
    x = addis(remii(sqri(x), g1), delta);
    g = gcdii(subii(x1, x), g1); if (!is_pm1(g)) break;

    if (DEBUGLEVEL >= 4 && (--c & 0x1f) == 0) rho_dbg(&T, c0-(c>>5), msg_mask);
  }

  if (g1 == n || equalii(g,g1))
  {
    if (g1 == n && equalii(g,g1))
    { /* out of luck */
      if (DEBUGLEVEL >= 4)
      {
        rho_dbg(&T, c0-(c>>5), 0);
        err_printf("\tPollard-Brent failed.\n"); err_flush();
      }
      if (++retries >= 4) return NULL;
      goto PB_RETRY;
    }
    /* half lucky: we've split n, but g1 equals either g or n */
    if (DEBUGLEVEL >= 4)
    {
      rho_dbg(&T, c0-(c>>5), 0);
      err_printf("\tfound %sfactor = %Ps\n", (g1!=n ? "composite " : ""), g);
      err_flush();
    }
    res = cgetg(7, t_VEC);
    /* g^1: known composite when g1!=n */
    INIT(res+1, g, gen_1, (g1!=n? gen_0: NULL));
    /* cofactor^1: status unknown */
    INIT(res+4, diviiexact(n,g), gen_1, NULL);
    return res;
  }
  /* g < g1 < n : our lucky day -- we've split g1, too */
  res = cgetg(10, t_VEC);
  /* unknown status for all three factors */
  INIT(res+1, g,                gen_1, NULL);
  INIT(res+4, diviiexact(g1,g), gen_1, NULL);
  INIT(res+7, diviiexact(n,g1), gen_1, NULL);
  if (DEBUGLEVEL >= 4)
  {
    rho_dbg(&T, c0-(c>>5), 0);
    err_printf("\tfound factors = %Ps, %Ps,\n\tand %Ps\n", res[1], res[4], res[7]);
    err_flush();
  }
  return res;
}

INLINE long
hamming_word(ulong w)
{
#if 1
  return __builtin_popcountl(w);
#else
  static long byte_weight[] = {
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
  while (w) { sum += byte_weight[w & 255]; w >>= 8; }
  return sum;
#endif
}

#ifndef S_SPLINT_S

static void
ensure_nb(GEN L, long l)
{
  long nmax = list_nmax(L);
  GEN v;
  if (l <= nmax) return;
  if (nmax)
  {
    nmax <<= 1;
    if (l > nmax) nmax = l;
    v = (GEN)pari_realloc(list_data(L), (nmax+1) * sizeof(long));
  }
  else /* unallocated */
  {
    nmax = 32;
    if (list_data(L))
      pari_err(e_MISC, "store list in variable before appending elements");
    v = (GEN)pari_malloc((nmax+1) * sizeof(long));
    v[0] = evaltyp(t_VEC) | _evallg(1);
  }
  list_data(L) = v;
  L[1] = evaltyp(list_typ(L))|evallg(nmax);
}


/* return a value <= (48 << 10) = 49152 < primelinit */
static ulong
utridiv_bound(ulong n)
{
#ifdef LONG_IS_64BIT
  if (n & HIGHMASK)
    return ((ulong)expu(n) + 1 - 16) << 10;
#else
  (void)n;
#endif
  return 1UL<<14;
}

INLINE ulong
u_forprime_next_fast(forprime_t *T)
{
  if (*(T->d))
  {
    NEXT_PRIME_VIADIFF(T->p, T->d);
    return T->p > T->b ? 0: T->p;
  }
  return u_forprime_next(T);
}

#endif
////////////////////////////////////////////////////////////////////////////////
