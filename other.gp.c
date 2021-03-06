#include "ext.h"
#include "extprv.h"

static GEN sumset_self(GEN a);
static GEN Eng_small(long n);
static GEN Eng_tiny(long n);
static GEN Edigit(long n);

/******************************************************************************/
/* General stuff */
/******************************************************************************/

void
init_auto(void)
{
    rnormal_cached = 0;
}


void
assume(int expr)
{
#ifdef DEBUG
  if (!expr) pari_err_BUG("assume");
#else
  if (!expr) __builtin_unreachable();
#endif
}


/******************************************************************************/
/* Set stuff */
/******************************************************************************/

GEN
vecsum_zv(GEN x)
{
  if (typ(x) != t_VECSMALL) pari_err_TYPE("vecsum_zv", x);
  /*
   * TODO: handle likely special case where some top bits are unused and so
   * many additions can be handled in long at a time.
  long i, minval = 0, maxval = 0, lx = lg(x);
  for (i=1; i<lx; i++) {
      if (x[i] < minval) minval = x[i];
      if (x[i] > maxval) maxval = x[i];
  }
  */
  pari_sp av = avma;
  GEN s = gen_0;
  long i, lx = lg(x);
  for (i = 1; i < lx; i++)
  {
    s = addis(s, x[i]);
    if (gc_needed(av, 6)) s = gerepileuptoint(av, s); // memory is 97% used
    // if (gc_needed(av, 3)) s = gerepileuptoint(av, s);   // memory is 80% used
  }
  return gerepileuptoint(av, s);
}


static GEN
sumset_self(GEN a)
{
  pari_sp av = avma;
  long i, j, k = 1, len = lg(a);
  GEN z = cgetg(len * (len - 1) / 2 + 1, t_VEC);
  for (i = 1; i < len; i++)
    for (j = 1; j <= i; j++)
      gel(z, k++) = gadd(gel(a, i), gel(a, j));
  return gerepileupto(av, gtoset(z));
}


/*
GP;install("sumset","GDG","sumset","./auto.gp.so");
GP;addhelp(sumset, "sumset(A, B): Set of all numbers of the form a+b, a in A, b in B. If B is omitted, assume A = B.");
*/
/**
 * @brief Given sets A and B, returns the set of all numbers of the form a+b, a in A, b in B.
 * 
 * @param a First input set
 * @param b Second input set
 * @return GEN Sumset of a and b
 */
GEN
sumset(GEN a, GEN b)
{
  if (typ(a) != t_VEC) pari_err_TYPE("sumset", a);
  if (b == NULL || gequal(a, b)) return sumset_self(a);
  if (typ(b) != t_VEC) pari_err_TYPE("sumset", b);
  // if (RgV_is_ZV(a) && RgV_is_ZV(b)) return sumset_ZV(a, b);
  pari_sp av = avma;
  long i, j, k = 1, lenA = lg(a), lenB = lg(b);
  GEN z = cgetg((lenA - 1) * (lenB - 1) + 1, t_VEC);
  for (i = 1; i < lenA; i++)
    for (j = 1; j < lenB; j++)
      gel(z, k++) = gadd(gel(a, i), gel(b, j));
  return gerepileupto(av, gtoset(z));
}


/*
GP;install("sumset_lim","GGDG","sumset_lim","./auto.gp.so");
GP;addhelp(sumset_lim, "sumset_lim(A, B, lim): Set of all numbers of the form a+b <= lim, a in A, b in B.");
*/
/**
 * @brief Given sets A and B, returns the set of all numbers of the form a+b <= lim, a in A, b in B.
 * 
 * @param a First input set
 * @param b Second input set
 * @param lim Highest sum to consider
 * @return GEN Sumset of a and b up to a limit of lim
 */
GEN
sumset_lim(GEN a, GEN b, GEN lim)
{
  if (typ(a) != t_VEC) pari_err_TYPE("sumset_lim", a);
  if (typ(b) != t_VEC) pari_err_TYPE("sumset_lim", b);
  if (!is_intreal_t(typ(lim))) pari_err_TYPE("sumset_lim", lim);
  pari_sp ltop = avma;

  long stack_bits = 8 * pari_mainstack->size;
  if (gcmpgs(lim, stack_bits) >= 0)
    pari_err(e_MEM); // morally: we'd overflow as soon as we'd try. Revise once
                     // we try sumset() first for appropriate inputs.
  if (signe(lim) != 1) return cgetg(1, t_VEC); // return empty vector

  long i, k = 0,
          slim = itos(gfloor(lim)); // fits in long because gcmpgs above is 0

  if ((gcmpgs(gel(a, 1), 1) < 0) || (gcmpgs(gel(b, 1), 1) < 0))
    pari_err_IMPL("not implemented: nonpositive in sumset_lim");
  RgV_check_ZV(a, "sumset_lim");
  RgV_check_ZV(b, "sumset_lim");

  if (!setisset(a)) a = gtoset(a);
  if (!setisset(b)) b = gtoset(b);
  long lenA = lg(a);
  long lenB = lg(b);

  pari_sp vtop = avma;
  GEN u, v = cgetg(slim + 1, t_VECSMALL);
  for (i = 1; i <= slim; ++i)
    v[i] = 0; // TODO: There are better ways to do this...

  for (i = 1; i < lenA; ++i)
  {
    long j;
    for (j = 1; j < lenB; ++j)
    {
      long t = itos(addii(gel(a, i), gel(b, j)));
      if (t > slim) break;
      v[t] = 1;
    }
    if (gc_needed(vtop, 2)) v = gerepileupto(vtop, v);
  }

  u = cgetg(itos(vecsum_zv(v)) + 1, t_VEC);
  pari_sp btop = avma;
  for (i = 1; i <= slim; ++i)
  {
    if (v[i])
    {
      gel(u, ++k) = stoi(i);
      if (gc_needed(btop, 1)) u = gerepilecopy(btop, u);
    }
  }
  return gerepilecopy(ltop, u);
}


/*
GP;install("diffset","D0,G,D0,G,","diffset","./auto.gp.so");
GP;addhelp(diffset, "diffset(A, B): Set of all numbers of the form a-b, a in A, b in B.");
*/
/**
 * @brief Given sets A and B, returns the set of all numbers of the form a-b, a in A, b in B.
 * 
 * @param a First set
 * @param b Second set
 * @return GEN Difference set of a and b
 */
GEN
diffset(GEN a, GEN b)
{
  pari_sp ltop = avma;
  long lenA = glength(a), lenB = glength(b);
  GEN c = cgetg(lenA * lenB + 1, t_VEC);
  pari_sp btop = avma, st_lim = stack_lim(btop, 1);
  long i, k = 0;
  for (i = 1; i <= lenA; ++i)
  {
    long j;
    for (j = 1; j <= lenB; ++j)
    {
      gel(c, ++k) = gsub(gel(a, i), gel(b, j));
    }
    if (low_stack(st_lim, stack_lim(btop, 1))) c = gerepileupto(btop, c);
  }
  return gerepileupto(ltop, gtoset(c));
}


/******************************************************************************/
/* Verbose monstrosities */
/******************************************************************************/

void
pBounds(GEN n, GEN verbose, long prec)
{
  pari_sp ltop = avma;
  GEN lower, upper, appx, l, ll;
  if (!verbose) verbose = gen_0;
  if (gcmpgs(n, 6548) < 0)
  {
    if (gcmpgs(n, 1) < 0)
      pari_printf("There are no negative primes.\n");
    else
    {
      n = gfloor(n);
      pari_printf("p_%Ps = %Ps (exactly)\n", n, prime(gtos(n)));
    }
    set_avma(ltop);
    return;
  }
  n = gfloor(n);
  l = glog(n, prec);
  ll = glog(l, prec);
  lower = gmul(n, gsubgs(gadd(l, ll), 1));
  /* Dusart, n >= 2 */
  if (gcmpgs(n, 13196) > 0)
    lower = gmul(n, gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)),
                         gdiv(strtor("2.25", prec), l)));
  /* Dusart, n >= 2 */
  appx = gmul(n, gadd(gadd(gsub(gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)),
                                     gdivsg(2, l)),
                                gdiv(gdivgs(gsqr(ll), 2), gsqr(l))),
                           gdiv(gmulsg(3, ll), gsqr(l))),
                      gdiv(gdivgs(utoipos(11), 2), gsqr(l))));
  /* + O(ll^3/l^3) */
  upper = gmul(n, gadd(l, ll));
  /* ?, n >= 6 */
  if (gcmpgs(n, 27076) >= 0)
    upper = gmul(n, gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)),
                         gdiv(strtor("1.8", prec), l)));
  /* Dusart, n >= 27076 */
  if (gcmpgs(n, 39017) >= 0)
    upper = gmin(upper, gmul(n, gsub(gadd(l, ll), strtor(".9484", prec))));
  /* Dusart, n >= 39017 */
  lower = gceil(lower);
  upper = gfloor(upper);
  pari_printf("%Ps (lower bound)\n", lower);
  if ((gcmp(lower, appx) < 0) && (gcmp(appx, upper) < 0))
    pari_printf("%Ps (approximate)\n", appx);
  pari_printf("%Ps (upper bound)\n", upper);
  if (!gequal0(verbose))
  {
    pari_printf("\nPierre Dusart, 'Autour de la fonction qui compte le nombre "
                "de nombres\n");
    pari_printf(
      "premiers', doctoral thesis for l'Université de Limoges (1998).\n");
    if ((gcmp(lower, appx) < 0) && (gcmp(appx, upper) < 0))
    {
      pari_printf("Ernest Cesàro (1894). \"Sur une formule empirique de M. "
                  "Pervouchine\". Comptes\n");
      pari_printf("rendus hebdomadaires des séances de l'Académie des sciences "
                  "119, pp. 848-849.\n");
    }
  }
  set_avma(ltop);
  return;
}


/******************************************************************************/
/* Works-in-progress / limited utility */
/******************************************************************************/

/*
GP;install("checkVDW","lD0,G,DG","checkVDW","./auto.gp.so");
GP;addhelp(checkVDW, "checkVDW(vv): Given a partition vv = [p1, p2, ...] with union(p1, p2, ...) = [1, 2, ..., n], finds a lower-bound proof for van der Waerden numbers based on the partition. Returns 0 if vv is not a partition of any initial segment, and k if vv proves that W(#vv, k) > n.");
*/
/**
 * @brief Given a partition vv = [p1, p2, ...] with union(p1, p2, ...) = [1, 2, ..., n], finds a lower-bound proof for van der Waerden numbers based on the partition. Returns 0 if vv is not a partition of any initial segment, and k if vv proves that W(#vv, k) > n.
 * 
 * @param vv A vector of disjoint vectors vv = [p1, p2, ...] with union(p1, p2, ...) = [1, 2, ..., n]
 * @param verbose 1 to print additional details or 0 to omit
 * @return GEN Returns 0 if vv is not a partition of any initial segment, and k if vv proves that W(#vv, k) > n.
 */
long
checkVDW(GEN vv, GEN verbose)
{
  pari_sp ltop = avma;
  GEN r, s = gen_0, c;
  long k = 0;
  if (!verbose) verbose = gen_1;
  r = stoi(glength(vv));
  c = cgetg(1, t_VEC);
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    GEN i = gen_0;
    for (i = gen_1; gcmp(i, r) <= 0; i = gaddgs(i, 1))
    {
      if (!gequal(gel(vv, gtos(i)), vecsort0(gel(vv, gtos(i)), NULL, 8)))
      {
        if (!gequal0(verbose))
          pari_printf("Not a partition: numbers repeated in color %Ps.\n", i);
        set_avma(ltop);
        return 0;
      }
      s = gaddgs(s, glength(gel(vv, gtos(i))));
      c = concat(c, gel(vv, gtos(i)));
      if (low_stack(st_lim, stack_lim(btop, 1)))
        gerepileall(btop, 3, &i, &s, &c);
    }
  }
  c = vecsort0(c, NULL, 8);
  if (gcmpgs(gel(c, 1), 1) < 0)
  {
    if (!gequal0(verbose))
      pari_printf(
        "Not a natural number partition: negative numbers in a member array.\n");
    return gc_long(ltop, 0);
  }
  if ((gcmpgs(gel(c, 1), 1) > 0) ||
      (gcmpgs(gel(c, glength(c)), glength(c)) > 0))
  {
    if (!gequal0(verbose))
      pari_printf(
        "Not a partition of an initial segment: not all numbers {1, 2, ..., %Ps} appear.\n",
        gel(c, glength(c)));
    return gc_long(ltop, 0);
  }
  k = longestProgression(gel(vv, 1)) + 1;
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    GEN i = gen_0;
    for (i = gen_2; gcmp(i, r) <= 0; i = gaddgs(i, 1))
    {
      k = maxss(k, longestProgression(gel(vv, gtos(i))) + 1);
      if (low_stack(st_lim, stack_lim(btop, 1))) i = gerepileuptoint(btop, i);
    }
  }
  if (!gequal0(verbose)) pari_printf("W(%Ps, %d) > %ld\n", r, k, glength(c));
  return k;
}


/*
GP;install("longestProgression","lD0,G,","longestProgression","./auto.gp.so");
GP;addhelp(longestProgression, "longestProgression(s): Finds the length of the longest arithmetic progression in s. Assumes that s is a set of integers, that is, is sorted from smallest to largest. Uses a space-efficient naive algorithm.");
*/
/**
 * @brief Finds the length of the longest arithmetic progression in the given set s of integers.
 * 
 * @param v A vector
 * @return GEN Length of the longest arithmetic progression
 */
long
longestProgression(GEN v)
{
  if(!is_vec_t(typ(v))) pari_err_TYPE("longestProgression", v);
  pari_sp ltop = avma;
  if (!setisset(v)) v = gtoset(v);
  long len = glength(v), i, r = 0;
  if (len < 3)
  {
    return gc_long(ltop, len);
  }

  for (i = 1; i < len; ++i)
  {
    long j;
    pari_sp btop = avma;
    for (j = i+1; j <= len; ++j)
    {
      long t = 2;
      GEN d = subii(gel(v, j), gel(v, i));
      pari_sp ctop = avma;
      while (setsearch(v, addii(gel(v, i), mulis(d, t)), 0))
      {
        ++t;
      }
      r = maxss(r, t);
      set_avma(ctop);
    }
    set_avma(btop);
  }
  return gc_long(ltop, r);
}

/*
GP;install("longestProgression1","lD0,G,","longestProgression1","./auto.gp.so");
GP;addhelp(longestProgression1, "longestProgression1(v): Uses a quadratic algorithm of Jeff Erickson, which is worst-case optimal; better algorithms are available when there are long progressions (> lg #v lg lg #v).");
*/
/**
 * @brief Finds the length of the longest arithmetic progression in the given vector. Uses a quadratic algorithm of Jeff Erickson, which is worst-case optimal; better algorithms are available when there are long progressions (> lg #v lg lg #v).
 * 
 * @param v A vector
 * @return GEN Length of the longest arithmetic progression
 */
long
longestProgression1(GEN v)
{
  if(!is_vec_t(typ(v))) pari_err_TYPE("longestProgression1", v);
  pari_sp ltop = avma;
  GEN L;
  long Lstar = 2, len = glength(v), len1 = len + 1, l6, l7;
  L = cgetg(len1, t_MAT);
  for (l7 = 1; l7 <= len; ++l7)
  {
    gel(L, l7) = cgetg(len1, t_COL);
    for (l6 = 1; l6 <= len; ++l6)
      gcoeff(L, l6, l7) = gen_0;
  }
  if (len < 3)
  {
    return gc_long(ltop, len);
  }
  pari_sp btop = avma, st_lim = stack_lim(btop, 1);
  long j;
  for (j = len - 1; j >= 1; --j)
  {
    long i = j - 1, k = j + 1;
    pari_sp ctop = avma, c_lim = stack_lim(ctop, 1);
    while (i > 0 && k <= len)
    {
      GEN tmp = subii(addii(gel(v, i), gel(v, k)), mulis(gel(v, j), 2));
      if (signe(tmp) < 0)
      {
        ++k;
      }
      else if (signe(tmp) > 0)
      {
        gcoeff(L, i, j) = gen_2;
        --i;
      }
      else
      {
        gcoeff(L, i, j) = addis(gcoeff(L, j, k), 1);
        Lstar = maxss(Lstar, itos(gcoeff(L, i, j)));
        --i;
        ++k;
      }
      if (low_stack(c_lim, stack_lim(ctop, 1))) L = gerepileupto(ctop, L);
    }
    while (i-- > 0)
    {
      gcoeff(L, i, j) = gen_2;
    }
    if (low_stack(st_lim, stack_lim(btop, 1))) L = gerepileupto(btop, L);
  }
  set_avma(ltop);
  return Lstar;
}


/*
Continued fraction:
        av = avma; lx = lg(x);
        e = bit_accuracy(lx)-1-expo(x);
        if (e < 0) pari_err(talker,"integral part not significant in gboundcf");
        c = trunc2nr_lg(x,lx,0);
        y = int2n(e);
        a = Qsfcont(c,y, NULL, 0);
        b = addsi(signe(x), c);
        return gerepilecopy(av, Qsfcont(b,y, a, 0));*/

// FIXME: Needs numerical analysis to determine stopping point.  Also needs to
// handle rational numbers and intgers.
/*
GP;install("Engel","G","Engel","./auto.gp.so");
GP;addhelp(Engel, "Engel(x): Engel expansion of x.");
*/
/**
 * @brief Engel expansion of the given number.
 * 
 * @param x Number to find Engel expansion for.
 * @return GEN Engel expansion
 */
GEN
Engel(GEN x)
{
  GEN v, t, ret;
  switch (typ(x))
  {
    case t_INT:
      if (signe(x) < 0) pari_err_DOMAIN("Engel", "x", "<", gen_0, x);
      long n = itos(x);
      long i = 1;
      v = cgetg(n + 1, t_VEC);
      for (; i <= n; ++i)
        gel(v, i) = gen_1;
      return v;

    case t_FRAC:
      pari_err_IMPL("fractions in Engel");

    case t_REAL:
      break;

    default:
      pari_err_TYPE("Engel", x);
  }

  pari_sp ltop = avma;
  v = listcreate();
  pari_sp btop = avma, st_lim = stack_lim(btop, 1);
  while (1)
  {
    if (!signe(x))
    {
      ret = gtovec(v);
      ret = gerepileupto(ltop, ret);
      return ret;
    }
    else
    {
      t = ceilr(invr(x));
    }
    listput_shallow(v, t);
    x = subrs(mulri(x, t), 1);
    if (low_stack(st_lim, stack_lim(btop, 1))) x = gerepileupto(btop, x);
  }
  __builtin_unreachable();
  return NEVER_USED;
}


/*
GP;install("Eng","G","Eng","./auto.gp.so");
GP;addhelp(Eng, "Eng(n): English name of the number n.");
*/
/**
 * @brief English name of the number.
 * 
 * @param n Number to name
 * @return GEN English name
 */
GEN
Eng(GEN n)
{
  if (typ(n) != t_INT) pari_err_TYPE("Eng", n);

  if (cmpis(n, 1000000000) >= 0)
  {
    pari_sp av = avma;
    GEN tmp = truedivis(n, 1000000000);
    GEN s = Str(mkvec2(Eng(tmp), strtoGENstr(" billion")));
    n = subii(n, mulis(tmp, 1000000000));
    if (signe(n)) s = Str(mkvec3(s, strtoGENstr(" "), Eng_small(itos(n))));
    s = gerepileupto(av, s);
    return s;
  }
  else
  {
    return Eng_small(itos(n));
  }
}


static GEN
Eng_small(long n)
{
  pari_sp av = avma;
  GEN s = strtoGENstr("");
  GEN space = strtoGENstr(" ");
  if (n >= 1000000)
  {
    long tmp = n / 1000000;
    n -= 1000000 * tmp;
    if (n)
      s = Str(mkvec4(s, Eng_tiny(tmp), strtoGENstr(" million"), space));
    else
      s = Str(mkvec3(s, Eng_tiny(tmp), strtoGENstr(" million")));
  }
  if (n >= 1000)
  {
    long tmp = n / 1000;
    s = Str(mkvec3(s, Eng_tiny(tmp), strtoGENstr(" thousand")));
    n -= 1000 * tmp;
    if (n) s = Str(mkvec3(s, space, Eng_tiny(n)));
    s = gerepileupto(av, s);
    return s;
  }
  s = Eng_tiny(n);
  s = gerepileupto(av, s);
  return s;
}


static GEN
Eng_tiny(long n)
{
  GEN s;
  pari_sp av = avma;
  if (n >= 100)
  {
    long tmp = n / 100;
    s = Str(mkvec2(Edigit(tmp), strtoGENstr(" hundred")));
    n -= 100 * tmp;
    if (!n)
    {
      s = gerepileupto(av, s);
      return s;
    }
    s = Str(mkvec2(s, strtoGENstr(" ")));
  }
  else
  {
    s = strtoGENstr("");
  }

  if (n < 20)
  {
    static const char* lookup[] = {
      "",        "one",     "two",       "three",    "four",
      "five",    "six",     "seven",     "eight",    "nine",
      "ten",     "eleven",  "twelve",    "thirteen", "fourteen",
      "fifteen", "sixteen", "seventeen", "eighteen", "ninteen"
    };
    GEN ret = Str(mkvec2(s, strtoGENstr(lookup[n])));
    ret = gerepileupto(av, ret);
    return ret;
  }

  long tmp = n / 10;
  static const char* lookup[] = { "",       "",      "twenty", "thirty",
                                  "forty",  "fifty", "sixty",  "seventy",
                                  "eighty", "ninety" };
  s = Str(mkvec2(s, strtoGENstr(lookup[tmp])));
  n -= 10 * tmp;
  if (n) s = Str(mkvec3(s, strtoGENstr("-"), Edigit(n)));
  s = gerepileupto(av, s);
  return s;
}


static GEN
Edigit(long n)
{
  static const char* lookup[] = { "",     "one", "two",   "three", "four",
                                  "five", "six", "seven", "eight", "nine" };
  return strtoGENstr(lookup[n]);
}

