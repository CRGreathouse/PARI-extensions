////////////////////////////////////////////////////////////////////////////////
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
  list_nmax(L) = nmax;
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
