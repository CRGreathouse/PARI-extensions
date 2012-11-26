////////////////////////////////////////////////////////////////////////////////
#ifndef S_SPLINT_S

/* Where to stop trial dividing in factorization. Guaranteed >= 2^14 */
static ulong
tridiv_bound(GEN n)
{
  ulong l = (ulong)expi(n) + 1;
  if (l <= 32)  return 1UL<<14;
  if (l <= 512) return (l-16) << 10;
  return 1UL<<19; /* Rho is generally faster above this */
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
