/* Error checking */

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
vecsum_zv(GEN x);
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

GEN
sumset_self(GEN a);
GEN
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

GEN
diffset(GEN a, GEN b) /* vecsmall */
{
  pari_sp ltop = avma;
  GEN c;
  long l1;
  GEN p2; /* vec */
  long l3;
  GEN p4; /* vecsmall */
  l1 = glength(a) * glength(b);
  {
    long l5;
    p2 = cgetg(l1 + 1, t_VEC);
    for (l5 = 1; l5 <= l1; ++l5)
      gel(p2, l5) = gen_0;
  }
  c = p2;
  l3 = glength(a);
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    long i, l6;
    for (i = 1; i <= l3; ++i)
    {
      l6 = glength(b);
      pari_sp ctop = avma, c_lim = stack_lim(ctop, 1);
      long j;
      for (j = 1; j <= l6; ++j)
      {
        gel(c, ((i - 1) * glength(b)) + j) = gsub(gel(a, i), gel(b, j));
        if (low_stack(c_lim, stack_lim(ctop, 1))) c = gerepilecopy(ctop, c);
      }
      if (low_stack(st_lim, stack_lim(btop, 1))) c = gerepilecopy(btop, c);
    }
  }
  p4 = vecsort0(geval(gtoset(c)), NULL, 0);
  p4 = gerepileuptoleaf(ltop, p4);
  return p4;
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

GEN
checkVDW(GEN vv, GEN verbose)
{
  pari_sp ltop = avma;
  GEN r, k = gen_0, s = gen_0, c;
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
        return gen_0;
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
    set_avma(ltop);
    return gen_0;
  }
  if ((gcmpgs(gel(c, 1), 1) > 0) ||
      (gcmpgs(gel(c, glength(c)), glength(c)) > 0))
  {
    if (!gequal0(verbose))
      pari_printf(
        "Not a partition of an initial segment: not all numbers {1, 2, ..., %Ps} appear.\n",
        gel(c, glength(c)));
    set_avma(ltop);
    return gen_0;
  }
  k = gaddgs(longestProgression(gel(vv, 1)), 1);
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    GEN i = gen_0;
    for (i = gen_2; gcmp(i, r) <= 0; i = gaddgs(i, 1))
    {
      k = gmax(k, gaddgs(longestProgression(gel(vv, gtos(i))), 1));
      if (low_stack(st_lim, stack_lim(btop, 1))) gerepileall(btop, 2, &i, &k);
    }
  }
  if (!gequal0(verbose)) pari_printf("W(%Ps, %Ps) > %ld\n", r, k, glength(c));
  k = gerepileupto(ltop, k);
  return k;
}

GEN
longestProgression(GEN v)
{
  pari_sp ltop = avma;
  GEN r = gen_0, s, t = gen_0, d = gen_0;
  long l2;
  if (glength(v) < 3)
  {
    return stoi(gc_long(ltop, glength(v)));
  }
  s = gtoset(v);
  l2 = glength(v) - 1;
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    GEN i = gen_0, p3 = gen_0;
    long l4;
    for (i = gen_1; gcmpgs(i, l2) <= 0; i = gaddgs(i, 1))
    {
      p3 = gaddgs(i, 1);
      l4 = glength(v);
      pari_sp ctop = avma, c_lim = stack_lim(ctop, 1);
      GEN j = gen_0;
      for (j = p3; gcmpgs(j, l4) <= 0; j = gaddgs(j, 1))
      {
        t = gen_2;
        d = gsub(gel(v, gtos(j)), gel(v, gtos(i)));
        pari_sp dtop = avma;
        while (setsearch(s, gadd(gel(v, gtos(i)), gmul(d, t)), 0))
        {
          t = gaddgs(t, 1);
          t = gerepileupto(dtop, t);
        }
        r = gmax(r, t);
        if (low_stack(c_lim, stack_lim(ctop, 1)))
          gerepileall(ctop, 4, &j, &t, &d, &r);
      }
      if (low_stack(st_lim, stack_lim(btop, 1)))
        gerepileall(btop, 5, &i, &p3, &t, &d, &r);
    }
  }
  r = gerepileupto(ltop, r);
  return r;
}

GEN
longestProgression1(GEN v)
{
  pari_sp ltop = avma;
  GEN Lstar = gen_2, L = gen_0, i = gen_0, k = gen_0, tmp = gen_0;
  long l1, l2;
  GEN p3; /* vec */
  long l4, l5;
  l1 = glength(v);
  l2 = glength(v);
  {
    long l6, l7;
    p3 = cgetg(l1 + 1, t_MAT);
    for (l7 = 1; l7 <= l1; ++l7)
    {
      gel(p3, l7) = cgetg(l2 + 1, t_COL);
      for (l6 = 1; l6 <= l2; ++l6)
        gcoeff(p3, l6, l7) = gen_0;
    }
  }
  L = p3;
  if (glength(v) < 3)
  {
    return stoi(gc_long(ltop, glength(v)));
    avma = ltop;
    return stoi(l4);
  }
  l5 = glength(v) - 1;
  {
    pari_sp btop = avma, st_lim = stack_lim(btop, 1);
    GEN j = gen_0;
    long l8 = -1 > 0; /* bool */
    for (j = stoi(l5); l8 ? gcmpgs(j, 1) <= 0 : gcmpgs(j, 1) >= 0;
         j = gaddgs(j, -1))
    {
      i = gsubgs(j, 1);
      k = gaddgs(j, 1);
      pari_sp ctop = avma, c_lim = stack_lim(ctop, 1);
      while ((gcmpgs(i, 0) > 0) && (gcmpgs(k, glength(v)) <= 0))
      {
        tmp = gsub(gadd(gel(v, gtos(i)), gel(v, gtos(k))),
                   gmulsg(2, gel(v, gtos(j))));
        if (gcmpgs(tmp, 0) < 0)
          k = gaddgs(k, 1);
        else
        {
          if (gcmpgs(tmp, 0) > 0)
          {
            gcoeff(L, gtos(i), gtos(j)) = gen_2;
            i = gsubgs(i, 1);
          }
          else
          {
            gcoeff(L, gtos(i), gtos(j)) =
              gaddgs(gcoeff(L, gtos(j), gtos(k)), 1);
            Lstar = gmax(Lstar, gcoeff(L, gtos(i), gtos(j)));
            i = gsubgs(i, 1);
            k = gaddgs(k, 1);
          }
        }
        if (low_stack(c_lim, stack_lim(ctop, 1)))
          gerepileall(ctop, 5, &tmp, &k, &L, &i, &Lstar);
      }
      ctop = avma;
      c_lim = stack_lim(ctop, 1);
      while (gcmpgs(i, 0) > 0)
      {
        gcoeff(L, gtos(i), gtos(j)) = gen_2;
        i = gsubgs(i, 1);
        if (low_stack(c_lim, stack_lim(ctop, 1))) gerepileall(ctop, 2, &L, &i);
      }
      if (low_stack(st_lim, stack_lim(btop, 1)))
        gerepileall(btop, 6, &j, &i, &k, &tmp, &L, &Lstar);
    }
  }
  Lstar = gerepileupto(ltop, Lstar);
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

GEN
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

GEN
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

GEN
Edigit(long n)
{
  static const char* lookup[] = { "",     "one", "two",   "three", "four",
                                  "five", "six", "seven", "eight", "nine" };
  return strtoGENstr(lookup[n]);
}
