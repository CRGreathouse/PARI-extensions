#include "ext.h"
#include "extprv.h"

static GEN mulNii(void* a, GEN x, GEN y);
static GEN prodtree(GEN A, long start, long stop);
static GEN prodtree_small(GEN A, long start, long stop);
static ulong ucomposite(long n);

/******************************************************************************/
/**                     Prime-related arithmetic functions                   **/
/******************************************************************************/

/*
GP;install("issemiprime","lG","issemi","./auto.gp.so");
GP;addhelp(issemi, "issemi(n): Is n a semiprime? Sloane's A064911; characteristic function of Sloane's A001358.");
*/
/**
 * @brief Test if the number is semiprime.
 * 
 * @param n Number to test
 * @return long 1 if the number is semiprime and 0 otherwise
 */
long
issemiprime(GEN n)
{
  // TODO: Optimize multi-word inputs
  if (typ(n) != t_INT) pari_err_TYPE("issemiprime", n);
  if (signe(n) <= 0) return 0;

  ulong nn = itou_or_0(n);
  if (nn) return uissemiprime(nn);

  pari_sp ltop = avma;
  if (!mpodd(n))
  {
    return gc_bool(ltop, mod4(n) && isprime(shifti(n, -1)));
  }


  long p;
  forprime_t primepointer;
  u_forprime_init(&primepointer, 3, 997); // What is a good breakpoint here?
  while ((p = u_forprime_next(&primepointer)))
  {
    if (dvdis(n, p))
    {
      return gc_bool(ltop, isprime(diviuexact(n, p)));
    }
  }

  if (isprime(n)) return 0;

  if (DEBUGLEVEL > 3)
    pari_printf(
      "issemi: Number is a composite with no small prime factors; using general factoring mechanisms.");

  GEN fac = Z_factor_until(
    n, shifti(
         n, -1)); // Find a nontrivial factor -- returns just the factored part
  GEN expo = gel(fac, 2);
  GEN pr = gel(fac, 1);
  long len = glength(expo);
  if (len > 2)
  {
    return gc_bool(ltop, 0);
  }
  if (len == 2)
  {
    if (cmpis(gel(expo, 1), 1) > 0 || cmpis(gel(expo, 2), 1) > 0)
    {
      return gc_bool(ltop, 0);
    }
    GEN P = gel(pr, 1);
    GEN Q = gel(pr, 2);
    return gc_bool(ltop, isprime(P) && isprime(Q) && equalii(mulii(P, Q), n));
  }
  if (len == 1)
  {
    long e = itos(gel(expo, 1));
    if (e == 2)
    {
      GEN P = gel(pr, 1);
      return gc_bool(ltop, isprime(P) && equalii(sqri(P), n));
    }
    else if (e > 2)
    {
      return gc_bool(ltop, 0);
    }
    GEN P = gel(pr, 1);
    return gc_bool(ltop, isprime(P) && isprime(diviiexact(n, P)));
  }

  pari_err_BUG(pari_sprintf(
    "Z_factor_until returned an unexpected value %Ps at n = %Ps, exiting...",
    fac, n));
  __builtin_unreachable();
  set_avma(ltop);
  return NEVER_USED;
}


// FIXME: Handle the case of small primelimit
/*
P=primorial(661)/2;
v=vector(10^4);i=0;forstep(n=2^64-1,1,-2,if(1==gcd(n,P),v[i++]=n;if(i==#v,return(#v))))
#
sum(i=1,#v,issemi(v[i]))
v=vector(10^5);i=0;forstep(n=2^64-1,1,-2,if(1==gcd(n,P),v[i++]=n;if(i==#v,return(#v))))
sum(i=1,#v,issemi(v[i]))

sum(i=1,#v,print(v[i]);issemi(v[i]))

No rho, crossover  661: 27.68s for 10k, 4:38.21s for 100k
   Rho, crossover  661: 28.06s for 10k, 4:45.50s for 100k
   Rho, crossover 1000: 28.05s for 10k, 4:43.82s for 100k
   Rho, crossover 3000: 27.98s for 10k, 4:44.21s for 100k

With prime test fronting:
   Rho, crossover  661: 13.60s for 10k, 2:12.36s for 100k
   Rho, crossover 1000: 13.72s for 10k, 2:12.53s for 100k
   Rho, crossover 1500: 13.36s for 10k, 2:12.98s for 100k
   Rho, crossover 2000: 13.65s for 10k, 2:12.80s for 100k
   Rho, crossover 2500: 13.57s for 10k, 2:12.98s for 100k
   Rho, crossover 3000: 13.50s for 10k, 2:13.30s for 100k
   Rho, crossover 4000: 13.48s for 10k, 2:12.98s for 100k
   Rho, crossover 9000: 13.54s for 10k, 2:14.17s for 100k

With prime test fronting, before even rho:
   Rho, crossover  661: 12.43s for 10k, 2:01.02s for 100k
   Rho, crossover  750: 12.50s for 10k, 2:01.03s for 100k
   Rho, crossover 1000: 12.38s for 10k, 2:00.59s for 100k
   Rho, crossover 1250: 12.40s for 10k, 2:01.12s for 100k
   Rho, crossover 1500: 12.38s for 10k, 2:00.77s for 100k
   Rho, crossover 2000: 12.58s for 10k, 2:01.75s for 100k
   Rho, crossover 3000: 12.42s for 10k, 2:01.08s for 100k
*/
long
uissemiprime(ulong n)
{
#define CUTOFF 1000
#ifdef LONG_IS_64BIT
#if CUTOFF <= 2642245
#define CUTOFF_CUBE CUTOFF* CUTOFF* CUTOFF
#else
#define CUTOFF_CUBE 18446744073709551615UL
#endif
#else
#if CUTOFF <= 1625UL
#define CUTOFF_CUBE CUTOFF* CUTOFF* CUTOFF
#else
#define CUTOFF_CUBE 4294967295UL
#endif
#endif
#if CUTOFF < 661
#error uissemiprime misconfigured, needs more primes to use uisprime_661.
#endif

  // Remove even numbers. Half of random inputs are caught here.
  if (!(n & 1)) return uisprime(n >> 1);

  // If n is small, simply test up to the cube root; no need for fancy stuff
  long lim;
  if (n <= CUTOFF_CUBE)
  {
    if (n < 9) return 0;
    lim = (long)cuberoot(n);

    long p = 0;
    byteptr primepointer = diffptr;
    NEXT_PRIME_VIADIFF(p, primepointer); // Skip 2
    for (;;)
    {
      NEXT_PRIME_VIADIFF(p, primepointer);
      if (p > lim) break;
      if (n % p == 0) return uisprime(n / p);
    }

    return !uisprime(n);
  }

  // First trial division loop, catches 'easy' numbers.
  // 83% of 'random' odd numbers trapped by this loop for CUTOFF = 661,
  // or more for larger values.
  long p = 0;
  byteptr primepointer = diffptr;
  NEXT_PRIME_VIADIFF(p, primepointer); // Skip 2
  for (;;)
  {
    NEXT_PRIME_VIADIFF(p, primepointer);
    if (p > CUTOFF) break;
    if (n % p == 0) return uisprime(n / p);
  }

  // Test for primality. About 27% of 661-rough numbers are caught here.
  if (uisprime_661(n)) return 0;

  // Check for a small prime factor with rho. Catches about 70% of remaining
  // composites, based on testing with CUTOFF = 1000.
  pari_sp ltop = avma;
  GEN fac = pollardbrent(utoipos(n));
  if (fac == NULL)
  {
    set_avma(ltop);
  }
  else if (typ(fac) == t_INT)
  {
    ulong f = itou(fac);
    return gc_bool(ltop, uisprime_661(f) && uisprime_661(n / f));
  }
  else if (typ(fac) == t_VEC)
  {
    // TODO: Slight speedup possible by paying attention to format instead
    // of just taking first factor:
    //   "a vector of t_INTs, each triple of successive entries containing
    //   a factor, an exponent (equal to one),  and a factor class (NULL
    //   for unknown or zero for known composite)"
    ulong f = itou(gel(fac, 1));
    return gc_bool(ltop, uisprime_661(f) && uisprime_661(n / f));
  }

  // Second part of trial division loop: avoids the cube root calculation
  // for numbers with a tiny prime divisor, and allows the use of
  // uisprime_nosmalldiv instead of uisprime.  Neither really matter for
  // hard numbers, but for 'average' numbers the first, at least, is
  // worthwhile.
  lim = cuberoot(n);
  for (;;)
  {
    if (p > lim) break;
    if (n % p == 0) return uisprime_661(n / p);
    NEXT_PRIME_VIADIFF(p, primepointer);
  }

  return 1;
#undef CUTOFF_CUBE
#undef CUTOFF
}


static GEN
mulNii(void* a, GEN x, GEN y)
{
  (void)a;
  return mulii(x, y);
} // for rad


/*
GP;install("rad","G","rad","./auto.gp.so");
GP;addhelp(rad, "rad(n): Radical of n, the largest squarefree number dividing n. Sloane's A007947.");
*/
/**
 * @brief Radical: the largest squarefree number dividing the given number.
 * 
 * @param n The number of which to take the radical
 * @return GEN The radical of n
 */
GEN
rad(GEN n)
{
  pari_sp ltop = avma;
  GEN f, ret;
  f = check_arith_all(n, "rad");
  if (!f) f = Z_factor(n);
  ret = gen_product(gel(f, 1), NULL, &mulNii);
  if (signe(ret) < 0) setabssign(ret);
  ret = gerepileupto(ltop, ret);
  return ret;
}


// 2-adic valuation of n
long
valu(ulong n)
{
#if 1
  return n ? __builtin_ctzl(n) : -1;
#else
  if (n == 0) return -1;
  long count = 0;
  while (!(n & 1))
  {
    n >>= 1;
    count++;
  }
  return count;
#endif
}


/*
GP;install("prp","lGDG","prp","./auto.gp.so");
GP;addhelp(prp, "prp(n,b=2): Is n a b-probable prime?");
*/
/**
 * @brief Tests if a number is a probable prime.
 * 
 * @param n Number to test
 * @param b Base of test (defaults to 2)
 * @return long 1 if the number is a probable prime or 0 otherwise.
 */
long
prp(GEN n, GEN b)
{
  pari_sp ltop = avma;
  if (typ(n) != t_INT) pari_err_TYPE("prp", n);
  if (!b)
    b = gen_2;
  else if (typ(b) != t_INT)
    pari_err_TYPE("prp", b);
  return gc_bool(ltop, gequal1(powgi(gmodulo(b, n), subis(n, 1))));
}


/*
GP;install("sprp","lGDG","sprp","./auto.gp.so");
GP;addhelp(sprp, "sprp(n,b=2): Is n a b-strong probable prime?");
*/
/**
 * @brief Tests if a number is a strong probable prime.
 * 
 * @param n Number to test
 * @param b Base of test (defaults to 2)
 * @return long 1 if the number is a strong probable prime or 0 otherwise.
 */
long
sprp(GEN n, GEN b)
{
  pari_sp ltop = avma;
  GEN d;
  if (typ(n) != t_INT)
    pari_err_TYPE("sprp", n);
  else if (cmpis(n, 3) < 0)
    return cmpis(n, 1) > 0; // Doesn't like even primes
  if (!b)
    b = gen_2;
  else if (typ(b) != t_INT)
    pari_err_TYPE("sprp", b);

  d = subis(n, 1);
  long s = vali(d);
  d = shifti(d, -s);
  d = powgi(gmodulo(b, n), d);
  if (gequal1(d))
  {
    return gc_bool(ltop, 1);
  }

  pari_sp btop = avma, st_lim = stack_lim(btop, 1);
  long i;
  for (i = 1; i <= s; i++)
  {
    if (gequalm1(d))
    {
      return gc_bool(ltop, 1);
    }
    d = gsqr(d);
    if (low_stack(st_lim, stack_lim(btop, 1))) gerepileall(btop, 1, &d);
  }
  return gc_bool(ltop, gequalm1(d));
}


/*
GP;install("sopf","G","sopf","./auto.gp.so");
GP;addhelp(sopf, "sopf(n): Sum of distinct prime factors of n. Sloane's A008472.");
*/
/**
 * @brief Sum of the distinct prime factors of the number.
 * 
 * @param n The input number
 * @return GEN The sum of the distinct prime factors of n
*/
GEN
sopf(GEN n)
{
  pari_sp ltop = avma;
  GEN f, ret = gen_0;
  f = check_arith_all(n, "sopf");
  if (!f) f = Z_factor(n);
  long l1 = glength(gel(f, 1));
  long i;
  for (i = 1; i <= l1; ++i)
  {
    ret = addii(ret, gcoeff(f, i, 1));
  }
  ret = gerepileupto(ltop, ret);
  return ret;
}


/*
GP;install("sopfr","G","sopfr","./auto.gp.so");
GP;addhelp(sopfr, "sopfr(n): Sum of prime factors of n (with multiplicity). Sloane's A001414.");
*/
/**
 * @brief Sum of the prime factors of the number, with repetition.
 * 
 * @param n The input number
 * @return GEN The sum of the prime factors of n with repetition
*/
GEN
sopfr(GEN n)
{
  pari_sp ltop = avma;
  GEN f, ret = gen_0;
  f = check_arith_all(n, "sopfr");
  if (!f) f = Z_factor(n);
  long l1 = glength(gel(f, 1));
  long i;
  for (i = 1; i <= l1; ++i)
  {
    ret = addii(ret, mulii(gcoeff(f, i, 1), gcoeff(f, i, 2)));
  }
  ret = gerepileupto(ltop, ret);
  return ret;
}


/*
GP;install("gpf","G","gpf","./auto.gp.so");
GP;addhelp(gpf, "gpf(n): The greatest prime factor of a number. Sloane's A006530.");
*/
/**
 * @brief Greatest prime factor of a given number.
 * 
 * @param n Number of which to find the greatest prime factor
 * @return GEN Greatest prime factor
 */
GEN
gpf(GEN n)
{
  GEN f, ret;
  pari_sp ltop = avma;
  f = check_arith_all(n, "gpf");
  if (!f) f = Z_factor(n);

  f = gel(f, 1);
  long len = glength(f); // Number of distinct prime factors

  if (len == 0) return gen_1;
  if (equalii(gel(f, 1), gen_0)) return gen_0;
  // My choice of convention: gpf(0) = 0

  ret = gel(f, len);
  if (equalii(ret, gen_m1)) ret = gen_1;
  ret = gerepileupto(ltop, ret);
  return ret;
}


static GEN
prodtree(GEN A, long start, long stop)
{
  pari_sp ltop = avma, st_lim = stack_lim(ltop, 1);
  long diff = stop - start;
  if (diff >= 8)
  {
    diff >>= 1;
    GEN leftprod = prodtree(A, start, start + diff);
    if (low_stack(st_lim, stack_lim(ltop, 1)))
      leftprod = gerepileupto(ltop, leftprod);
    pari_sp btop = avma;
    GEN rightprod = prodtree(A, start + diff + 1, stop);
    // if (low_stack(st_lim, stack_lim(ltop, 1)))
    rightprod = gerepileupto(btop, rightprod);
    GEN ret = mulii(leftprod, rightprod);
    ret = gerepileupto(ltop, ret);
    return ret;
  }

  GEN ret, a, b, c, d;
  switch (diff)
  {
    case 7:
      a = mulss(A[start], A[start + 7]);
      b = mulss(A[start + 1], A[start + 6]);
      c = mulss(A[start + 2], A[start + 5]);
      d = mulss(A[start + 3], A[start + 4]);
      ret = mulii(mulii(a, b), mulii(c, d));
      break;
    case 6:
      a = mulss(A[start], A[start + 3]);
      b = mulss(A[start + 1], A[start + 2]);
      c = mulss(A[start + 4], A[start + 6]);
      ret = mulii(mulii(a, b), mulis(c, A[start + 5]));
      break;
    case 5:
      a = mulss(A[start], A[start + 5]);
      b = mulss(A[start + 1], A[start + 4]);
      ret = mulii(mulis(a, A[start + 2]), mulis(b, A[start + 3]));
      break;
    case 4:
      a = mulss(A[start], A[start + 2]);
      b = mulss(A[start + 3], A[start + 4]);
      ret = mulii(mulis(a, A[start + 1]), b);
      break;
    case 3:
      a = mulss(A[start], A[start + 3]);
      b = mulss(A[start + 1], A[start + 2]);
      ret = mulii(a, b);
      break;
    default:
      pari_err_BUG("prodtree passed small argument");
      __builtin_unreachable();
  }
  ret = gerepileupto(ltop, ret);
  return ret;
}


static GEN
prodtree_small(GEN A, long start, long stop)
{
  pari_sp ltop = avma, st_lim = stack_lim(ltop, 1);
  long diff = stop - start;
  if (diff >= 8)
  {
    diff >>= 1;
    GEN leftprod = prodtree(A, start, start + diff);
    if (low_stack(st_lim, stack_lim(ltop, 1)))
      leftprod = gerepileupto(ltop, leftprod);
    pari_sp btop = avma;
    GEN rightprod = prodtree(A, start + diff + 1, stop);
    // if (low_stack(st_lim, stack_lim(ltop, 1)))
    rightprod = gerepileupto(btop, rightprod);
    GEN ret = mulii(leftprod, rightprod);
    ret = gerepileupto(ltop, ret);
    return ret;
  }

  GEN ret;
  long a, b, c, d;
  switch (diff)
  {
    case 7:
      a = A[start] * A[start + 7];
      b = A[start + 1] * A[start + 6];
      c = A[start + 2] * A[start + 5];
      d = A[start + 3] * A[start + 4];
      ret = mulii(mulss(a, b), mulss(c, d));
      break;
    case 6:
      a = A[start] * A[start + 3];
      b = A[start + 1] * A[start + 2];
      c = A[start + 4] * A[start + 6];
      ret = mulii(mulss(a, b), mulss(c, A[start + 5]));
      break;
    case 5:
      a = A[start] * A[start + 5];
      b = A[start + 1] * A[start + 4];
      ret = mulii(mulss(a, A[start + 2]), mulss(b, A[start + 3]));
      break;
    case 4:
      a = A[start] * A[start + 2];
      b = A[start + 3] * A[start + 4];
      ret = mulis(mulss(a, A[start + 1]), b);
      break;
    case 3:
      a = A[start] * A[start + 3];
      b = A[start + 1] * A[start + 2];
      ret = mulss(a, b);
      break;
    default:
      pari_err_BUG("prodtree_small passed small argument");
      __builtin_unreachable();
  }
  ret = gerepileupto(ltop, ret);
  return ret;
}


/*
GP;install("primorial","G","primorial","./auto.gp.so");
GP;addhelp(primorial, "primorial(n): Returns the product of each prime less than or equal to n. Sloane's A034386.");
*/
/**
 * @brief Primorial, that is, product of the primes up to the given number.
 * 
 * @param n Number up to which primes are multiplied
 * @return GEN Product of primes
*/
GEN
primorial(GEN n)
{
  static const ulong smallpr[] = {
    1UL, 1UL, 2UL, 6UL, 6UL, 30UL, 30UL, 210UL, 210UL, 210UL, 210UL, 2310UL,
    2310UL, 30030UL, 30030UL, 30030UL, 30030UL, 510510UL, 510510UL, 9699690UL,
    9699690UL, 9699690UL, 9699690UL, 223092870UL, 223092870UL, 223092870UL,
    223092870UL, 223092870UL, 223092870UL
#ifdef LONG_IS_64BIT
    , 6469693230UL, 6469693230UL, 200560490130UL, 200560490130UL,
    200560490130UL, 200560490130UL, 200560490130UL, 200560490130UL,
    7420738134810UL, 7420738134810UL, 7420738134810UL, 7420738134810UL,
    304250263527210UL, 304250263527210UL, 13082761331670030UL,
    13082761331670030UL, 13082761331670030UL, 13082761331670030UL,
    614889782588491410UL, 614889782588491410UL, 614889782588491410UL,
    614889782588491410UL, 614889782588491410UL, 614889782588491410UL
#endif
  };

  pari_sp ltop = avma;
  ulong nn;
  GEN ret;
  if (typ(n) == t_REAL)
  {
    if (signe(n) < 1) return gen_1;
    nn = itou_or_0(floorr(n));
    set_avma(ltop);
  }
  else if (typ(n) == t_INT)
  {
    if (signe(n) < 1) return gen_1;
    nn = itou_or_0(n);
  }
  else
  {
    pari_err_TYPE("primorial", n);
    __builtin_unreachable();
  }

  if (nn > maxprime() || nn == 0) // nn == 0 if n didn't fit into a word
    pari_err_MAXPRIME(nn);
#ifdef LONG_IS_64BIT
  if (nn < 53)
#else
  if (nn < 29)
#endif
    return utoipos(smallpr[nn]);

  ulong primeCount = uprimepi(nn);
  GEN pr = primes_zv(primeCount);
#ifdef LONG_IS_64BIT
  if (primeCount < 146144319)
#else
  if (primeCount < 4793)
#endif
  {
    ret = prodtree_small(pr, 1, primeCount);
  }
  else
  {
    /*
    pari_sp btop;
    GEN right = prodtree(pr, primeCount>>1 + 1, primeCount);
    // how to mark latter half of pr as unallocated?
    right = gerepileupto(btop, right);
    ret = mulii(prodtree(pr, 1, primeCount>>1), right);
    */
    ret = prodtree(pr, 1,
                   primeCount); // Possible TODO: Free memory from latter half
                                // of array once its product is calculated?
  }
  ret = gerepileupto(ltop, ret);
  return ret;
}


ulong
lpfu(ulong n)
{
  if (n < 2) return n; // Convention: lpf(1) = lpf(-1) = 1, lpf(0) = 0
  if (!(n & 1)) return 2;

  ulong p;
  forprime_t S;
  u_forprime_init(&S, 3, utridiv_bound(n)); // At least 2^14
  while ((p = u_forprime_next_fast(&S)) <= 661)
  {
    if (n % p == 0) return p;
  }
  if (uisprime_661(n)) return n;
  while (p)
  {
    if (n % p == 0) return p;
    p = u_forprime_next_fast(&S);
  }

  // No small prime factors, factor using general mechanisms.
  // This is a performance disaster, since factoru repeats trial division.
  // ifac_factoru would be better.
  pari_sp ltop = avma;
  return gc_ulong(ltop, vecsmall_min(gel(factoru(n), 1)));
}


/*
GP;install("lpf","G","lpf","./auto.gp.so");
GP;addhelp(lpf, "lpf(n): The least prime factor of a number. Sloane's A020639.");
*/
/**
 * @brief Least prime factor of a given number.
 * 
 * @param n Number of which to find the least prime factor
 * @return GEN Least prime factor
 */
GEN
lpf(GEN n)
{
  pari_sp ltop = avma;
  GEN res, f;
  if ((f = check_arith_all(n, "lpf")))
  {
    f = gel(f, 1);
    long len = glength(f); // Number of distinct prime factors

    // Convention: lpf(1) = lpf(-1) = 1, lpf(0) = 0
    if (len == 0 || (len == 1 && equalii(gel(f, 1), gen_m1)))
    {
      set_avma(ltop);
      return gen_1;
    }
    if (len == 1 && equalii(gel(f, 1), gen_0))
    {
      set_avma(ltop);
      return gen_0;
    }

    // Skip the -1 if present
    if (equalii(gel(f, 1), gen_m1))
    {
      res = gel(f, 2);
    }
    else
    {
      res = gel(f, 1);
    }
    res = gerepileupto(ltop, res);
    return res;
  }

  ulong nn = itou_or_0(n);
  if (nn) return utoipos(lpfu(nn));

  if (!signe(n)) return gen_0; // My choice of convention: lpf(0) = 0

  /* Trial division to 661 */
  if (!mod2(n)) return gen_2;
#ifdef LONG_IS_64BIT
  static const ulong G = 16294579238595022365UL;
#else
  static const ulong G = 3234846615UL;
#endif
  ulong g = ugcd(umodiu(n, G), G);
  if (g > 1) return utoipos(lpfu(g));
  forprime_t S;
#ifdef LONG_IS_64BIT
  u_forprime_init(&S, 59, 661);
#else
  u_forprime_init(&S, 31, 661);
#endif
  long p = 0;
  while ((p = u_forprime_next_fast(&S)))
  {
    if (dvdis(n, p))
    {
      return utoipos(gc_ulong(ltop, p));
    }
  }

  /* Check if we have a prime and can stop. */
  if (BPSW_psp_nosmalldiv(n) && BPSW_isprime(n))
  {
    set_avma(ltop);
    return icopy(n);
  }

  /* Check table of special primes */
  long l = lg(primetab), i;
  for (i = 1; i < l; i++)
  {
    GEN q = gel(primetab, i);
    if (dvdii(n, q))
    {
      return gerepileuptoint(ltop, q);
    }
  }


  /* How far should we search? If the number is small, just factor it. */
  ulong lim = tridiv_bound(n);
#ifdef LONG_IS_64BIT
  if (expi(n) > 132)
  {
    lim = 10000000000UL;
#else
  if (expi(n) > 95)
  {
    lim = 0xFFFFFFFFUL;
#endif
  }
  else if (expi(n) > 62)
  {
    /* Big enough that factoring the whole number might be impractical */
    lim = itou(sqrtnint(n, 4));
  }

  u_forprime_init(&S, 673, lim);
  p = 0;
  while ((p = u_forprime_next_fast(&S)))
  {
    if (dvdis(n, p))
    {
      return utoipos(gc_ulong(ltop, p));
    }
  }
  res = gcoeff(Z_factor(n), 1,
               1); // TODO: Partial factorization?  Tricky to do right...
  res = gerepileupto(ltop, res);
  return res;
}


// TODO: Would be much more efficient to just walk through the primelist...
static ulong
ucomposite(long n)
{
  if (n < 2)
  {
    if (n < 1) pari_err_DOMAIN("ucomposite", "n", "<", gen_1, utoi(n));
    return 4;
  }
  double l = 1.0 / log(n);

  // Series apparently due to Bojarincev 1967; more terms are known if desired
  ulong c = (ulong)(n * (1 + l * (1 + 2 * l * (1 + 2 * l))));
  int i = 0;
  int limit = (int)(n - c + uprimepi(c));
  for (; i <= limit; i++)
    if (uisprime(++c)) i--;
  return c;
}


/*
GP;install("composite","L","composite","./auto.gp.so");
GP;addhelp(composite, "composite(n): Returns the n-th composite. Sloane's A002808.");
*/
/**
 * @brief Returns the n-th composite number.
 * 
 * @param n Index
 * @return GEN Composite number
 */
GEN
composite(long n)
{
  return utoipos(ucomposite(n));
}
