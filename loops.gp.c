/******************************************************************************/
/* Looping constructs */
/******************************************************************************/

// Digit reversal
GEN
rev(GEN n, long B)
{
  pari_sp av = avma;
  if (typ(n) != t_INT) pari_err_TYPE("rev", n);
  GEN m = modis(n, B);
  n = divis(n, B);

  pari_sp btop = avma, st_lim = stack_lim(btop, 1);
  while (signe(n))
  {
    m = addis(mulis(m, B), smodis(n, B));
    n = divis(n, B);
    if (low_stack(st_lim, stack_lim(btop, 1))) gerepileall(btop, 2, &m, &n);
  }
  m = gerepilecopy(av, m);
  return m;
}


// Return value: Did the user break out of the loop?
// Not stack clean.
int
palhelper(long digits, GEN a, GEN b, GEN code)
{
  GEN p10 = powuu(10, (digits + 1) >> 1);
  GEN aLeft = divii(a, p10);
  GEN bLeft = divii(b, p10);
  GEN cur;

  // TODO: Handle case of digits odd (middle digit)

  pari_sp btop = avma, lim = stack_lim(btop, 1);
  cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
  if (cmpii(cur, a) < 0)
  {
    aLeft = addis(aLeft, 1);
    cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
  }

  push_lex(cur, code);
  while (cmpii(aLeft, bLeft) < 0)
  {
    closure_evalvoid(code);
    if (loop_break())
    {
      pop_lex(1);
      return (1);
    }
    // cur = get_lex(-1);	// Allow the user to modify the variable
    aLeft = addis(aLeft, 1);
    cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
    set_lex(-1, cur); // Set the variable atop the stack to the current value

    if (low_stack(lim, stack_lim(btop, 1)))
    {
      if (DEBUGMEM > 1) pari_warn(warnmem, "forpal");
      cur = gerepileupto(btop, cur);
    }
  }
  // TODO: Handle final few numbers
  pop_lex(1);
  return 0;
}


void
forpal(GEN a, GEN b, GEN code)
{
  pari_sp av = avma;

  if (typ(a) == t_REAL)
    a = gceil(a);
  else if (typ(a) != t_INT)
    pari_err_TYPE("forpal", a);
  if (typ(b) == t_REAL)
    b = gfloor(b);
  else if (typ(b) != t_INT)
    pari_err_TYPE("forpal", b);

  if (cmpii(a, b) > 0) return;

  long lower_digits = countdigits(a);
  long upper_digits = countdigits(b);
  if (lower_digits == upper_digits)
  {
    palhelper(lower_digits, a, b, code);
    set_avma(av);
    return;
  }

  pari_sp btop = avma;
  if (palhelper(lower_digits, a, powuu(10, lower_digits), code))
  {
    set_avma(av);
    return;
  }
  set_avma(btop);
  long d = lower_digits + 1;
  for (; d < upper_digits; d++)
  {
    if (palhelper(d, powuu(10, d - 1), powuu(10, d), code))
    {
      set_avma(av);
      return;
    }
    set_avma(btop);
  }
  palhelper(upper_digits, powuu(10, upper_digits - 1), b, code);
  set_avma(av);
}
