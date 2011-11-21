/******************************************************************************/
/**												Looping constructs																**/
/******************************************************************************/

void
forodd(GEN a, GEN b, GEN code)
{
	pari_sp av, av0 = avma, lim;

	if (typ(a) == t_REAL)
		a = gceil(a);
	else if (typ(a) != t_INT)
		pari_err_TYPE("forodd", a);
	if (typ(b) == t_REAL)
		b = gfloor(b);
	else if (typ(b) != t_INT)
		pari_err_TYPE("forodd", b);
	
	if (!mpodd(a))
		a = gbitor(a, gen_1);
	av=avma;
	lim = stack_lim(av,1);
	push_lex(a, code);
	while (cmpii(a,b) <= 0)
	{
		closure_evalvoid(code); if (loop_break()) break;
		//a = get_lex(-1);	// Allow the user to modify the variable
		a = addis(a, 2);
		set_lex(-1, a);	// Set the variable atop the stack to the value of a

		if (low_stack(lim, stack_lim(av,1)))
		{
			if (DEBUGMEM>1) pari_warn(warnmem,"forodd");
			a = gerepileupto(av,a);
		}
	}
	pop_lex(1);
	avma = av0;
}


// Digit reversal
GEN
rev(GEN n, long B)
{
	pari_sp av = avma;
	if (typ(n) != t_INT)
		pari_err_TYPE("rev", n);
	GEN m = modis(n, B);
	n = divis(n, B);
	
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while (signe(n)) {
		m = addis(mulis(m, B), smodis(n, B));
		n = divis(n, B);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &m, &n);
	}
	m = gerepilecopy(av, m);
	return m;
}


// Return value: Did the user break out of the loop?
// Not stack clean.
int palhelper(long digits, GEN a, GEN b, GEN code)
{
	GEN p10 = powuu(10, (digits+1)>>1);
	GEN aLeft = divii(a, p10);
	GEN bLeft = divii(b, p10);
	GEN cur;

	// TODO: Handle case of digits odd (middle digit)
	
	pari_sp btop = avma, lim = stack_lim(btop,1);
	cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
	if (cmpii(cur, a) < 0) {
		aLeft = addis(aLeft, 1);
		cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
	}
	
	push_lex(cur, code);
	while (cmpii(aLeft, bLeft) < 0)
	{
		closure_evalvoid(code);
		if (loop_break()) {
			pop_lex(1);
			return(1);
		}
		//cur = get_lex(-1);	// Allow the user to modify the variable
		aLeft = addis(aLeft, 1);
		cur = addii(mulii(aLeft, p10), rev(aLeft, 10));
		set_lex(-1, cur);	// Set the variable atop the stack to the current value

		if (low_stack(lim, stack_lim(btop,1)))
		{
			if (DEBUGMEM>1) pari_warn(warnmem,"forpal");
			cur = gerepileupto(btop,cur);
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
	
	if (cmpii(a, b) > 0)
		return;
	
	long lower_digits = digits(a);
	long upper_digits = digits(b);
	if (lower_digits == upper_digits) {
		palhelper(lower_digits, a, b, code);
		avma = av;
		return;
	}
	
	pari_sp btop = avma;
	if (palhelper(lower_digits, a, powuu(10, lower_digits), code)) {
		avma = av;
		return;
	}
	avma = btop;
	long d = lower_digits + 1;
	for (; d < upper_digits; d++) {
		if (palhelper(d, powuu(10, d-1), powuu(10, d), code)) {
			avma = av;
			return;
		}
		avma = btop;
	}
	palhelper(upper_digits, powuu(10, upper_digits-1), b, code);
	avma = av;
}


// TODO: Cleanup
// FIXME: Doesn't work if the user gives the 'wrong' variable (other than x).
// http://pari.math.u-bordeaux.fr/archives/pari-dev-1002/msg00025.html
GEN
sumformal(GEN start, GEN end, GEN expr)
{
	pari_printf("%Ps has degree %d with leading coefficient %Ps.\n", expr, degree(expr), truecoeff(expr, degree(expr)));
	pari_sp ltop = avma;
	GEN c, F, res = gen_0, ret;
	GEN x = pol_x(fetch_var());
push_lex(x, expr);
	pari_printf("Main variable: %Ps\n", x);
	long t = typ(expr);
	if (t == t_INT || t == t_REAL || t == t_COMPLEX)
	{
		ret = gmul(gaddgs(gsub(end, start), 1), expr);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	if (t != t_POL)
		pari_err_TYPE("sumformal", expr);

	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long d = degree(expr) + 1;
	while (--d)
	{
		c = truecoeff(expr, d);
		F = Faulhaber(d, x);
		res = gadd(res, gmul(c, gsub(gsubstpol(F, x, end), gsubstpol(F, x, gsubgs(start, 1)))));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &res);
	}
pop_lex(1);
delete_var();
	ret = gadd(res, gmul(gaddgs(gsub(end, start), 1), truecoeff(expr, 0)));
	ret = gerepileupto(ltop, ret);
	return ret;
}

/*
GEN
sumformal(GEN start, GEN end, GEN expr)
{
	pari_sp ltop = avma;
	GEN c, F, res = gen_0, ret;
	long d;
	GEN x = pol_x(fetch_user_var("x"));
	push_lex(x, expr);
	if ((typ(expr) == t_INT) || (typ(expr) == t_REAL))
	{
		ret = gmul(gaddgs(gsub(end, start), 1), expr);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	if (typ(expr) != t_POL)
		pari_err(notpoler, "sumformal; can only handle polynomials, not arbitrary functions");

	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while ((d = poldegree(expr, -1)) > 0)
	{
		c = polcoeff0(expr, d, -1);
		F = Faulhaber(stoi(d), x);
		res = gadd(res, gmul(c, gsub(gsubstpol(F, x, end), gsubstpol(F, x, gsubgs(start, 1)))));
		// TODO: This is treating expr like a sparse poly, but it's not...
		// should just pull out the coefficients one by one.
		expr = gsub(expr, gmul(c, gpow(x, stoi(d), FAKE_PREC)));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 4, &d, &c, &res, &expr);
	}
	pop_lex(1);
	ret = gadd(res, gmul(gaddgs(gsub(end, start), 1), expr));
	ret = gerepileupto(ltop, ret);
	return ret;
}
*/


void
fortwin(GEN ga, GEN gb, GEN code)
{
	long p[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3), 0};	// magic
	ulong *prime = (ulong*)p;
	ulong a, b;
	pari_sp av = avma;
	byteptr d;

	d = prime_loop_init(ga,gb, &a,&b, (ulong*)&prime[2]);
	if (!d) { avma = av; return; }

	avma = av; push_lex((GEN)prime,code);
	while (prime[2] < b)
	{
		if (*d != 2){
			NEXT_PRIME_VIADIFF(prime[2], d);
			continue;
		}
		closure_evalvoid(code); if (loop_break()) break;
		NEXT_PRIME_VIADIFF(prime[2], d);
		avma = av;
	}
	if (prime[2] == b) {
		closure_evalvoid(code);
		(void)loop_break();
		avma = av;
	}
	pop_lex(1);
}


// TODO: Make a function that allows GENs, especially in the specialized case of 64-bit primes in a 32-bit environment
void
forbigprime(GEN ga, GEN gb, GEN code)
{
	pari_sp av = avma;
	long t = typ(ga);
	if (t == t_REAL || t == t_FRAC)
		ga = gceil(ga);
	else if (t != t_INT)
		pari_err_TYPE("forbigprime", ga);

	t = typ(gb);
	if (t == t_REAL || t == t_FRAC)
		gb = gfloor(gb);
	else if (t != t_INT)
		pari_err_TYPE("forbigprime", gb);

	if (signe(gb) < 1)
	{
		avma = av;
		return;
	}
	if (signe(ga) < 1)
		ga = gen_2;
	ulong a = itou_or_0(ga);
	ulong b = itou_or_0(gb);
	
	if (!b) {
		if (cmpii(ga, gb) <= 0)
			pari_err(e_MISC, "Only works for single-word integers.");
		avma = av;
		return;
	}
	if (!a || a > b) {	// ga > gb
		avma = av;
		return;
	}
	
	if (b < maxprime())
	{
		if (DEBUGLEVEL>3) fprintferr("Using precomputed primes up to %Ps...\n", gb);
		forprime(ga, gb, code);
		avma = av;
		return;
	}
	avma = av;
#ifdef LONG_IS_64BIT
	if (maxprime() < 4294967296UL && b > maxprime() * maxprime())
		pari_err_MAXPRIME(usqrtsafe(b));
#else
	// PARI currently guarantees that maxprime() >= 65557, so no check needed
#endif

	if (a < maxprime())
	{
		if (DEBUGLEVEL>3) fprintferr("Using precomputed primes up to %lu...\n", maxprime());
		forprime(ga, utoi(maxprime()), code);
		avma = av;
		if (loop_break())
			return;
		a = maxprime() + 1;
	}
	
	forbigprime_sieve(a|1, b, code);
}


void
sieve_block(ulong a, ulong b, char* sieve)
{
	if (b < a) {
		pari_warn(warner, "sieve_block called needlessly!");
		return;
	}
    if (DEBUGLEVEL>4) fprintferr("Sieving from %lu to %lu", a, b);
	ulong lim = usqrtsafe(b);
	ulong sz = (b - a + 2) >> 1;
	if (DEBUGLEVEL>4) fprintferr("; size = %lu\n", sz);
	long p = 0;
	
	memset(sieve, 0, sz);
	
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > lim)
			break;
		
		// sieve[0] is a, a+2, a+4, a+6, a+8, a+10, a+12, a+14
		// sieve[1] is a+16, ..., a+30
		// pos n <--> a+2n
		// pos n is sieve[n>>3], bit n&7
		int pos = p - a%p;
		if (pos&1) {
			if (pos == p)
				pos = 0;
			else
				pos = (pos + p) >> 1;
		} else
			pos >>= 1;
		
		while (pos <= sz) {
			sieve[pos>>3] |= 1 << (pos&7);
			pos += p;
		}
	}
}

/*
#
default(primelimit,4*10^9)
default(debug,3)
forbigprime(p=1e19,1e19+1e8,if(p%199==1,print1(".")))
*/
void
forbigprime_sieve(ulong a, ulong b, GEN code)
{
	// TODO: Optimize size (surely < 512k to stay in L1 cache, but not so large
	// as to force recalculating too often).
	// Guesstimate: greater of sqrt(n) * lg(n) or 1M
	ulong chunk = maxuu(0x100000, usqrtsafe(b) * __builtin_ffsl(b));
	ulong tmp = (b - a) / chunk + 1;
	pari_sp ltop = avma;

	if (tmp == 1)
		chunk = b - a + 16;
	else
		chunk = (b - a) / tmp + 15;
	chunk = minuu(chunk, avma - stack_lim(avma, 2));	// Don't take up more than 2/3 of the stack
	
	// chunk + 2 should be divisible by 16
	chunk = (((chunk + 2)>>4)<<4) - 2;
	ulong maxpos = (chunk + 2) >> 4;	// Shift by 1 since only odds are considered; shift by 3 to convert from bits to bytes
	
		if (DEBUGLEVEL>2) {
		tmp = (b - a) / chunk + 1;
		fprintferr("Chunk size %lu (%.2f MB of %.2f MB free), splitting the work into ~%lu parts\n", chunk, (float)(maxpos * 9.53674316e-7), (float)((avma - bot) * 9.53674316e-7), tmp);
	}
	
	char* sieve = pari_malloc(chunk);
	
	// A GEN representing a prime to be passed to the code; its value is in p[2].
	long p[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3), 0};
	push_lex((GEN)p,code);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);

	while (a <= b) {
		ulong end = a + chunk;
		if (end > b)
			break;
		sieve_block(a, end, sieve);	// Sieve the interval
		int pos = 0;
		for (; pos < maxpos; pos++) {
			if (sieve[pos] == 0xFF)
				continue;
			if (!(sieve[pos]&1))
			{
				p[2] = a + (pos << 4);
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&2))
			{
				p[2] = a + (pos << 4) + 2;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&4))
			{
				p[2] = a + (pos << 4) + 4;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&8))
			{
				p[2] = a + (pos << 4) + 6;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&16))
			{
				p[2] = a + (pos << 4) + 8;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&32))
			{
				p[2] = a + (pos << 4) + 10;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&64))
			{
				p[2] = a + (pos << 4) + 12;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (!(sieve[pos]&128))
			{
				p[2] = a + (pos << 4) + 14;
				closure_evalvoid(code);
				if (loop_break()) goto CLEANUP;
			}
			if (low_stack(st_lim, stack_lim(btop, 1)))
				avma = btop;
		}
		a = end + 2;
	}
	
	// Handle the last chunk.	This tests the endpoint at every step.
	if (b < a)
		goto CLEANUP;
	if (DEBUGLEVEL>3) fprintferr("Last chunk: ");
	sieve_block(a, b, sieve);	// Sieve the interval

	int pos = 0;
	chunk = b - a + 2;
	for (; pos <= chunk; pos++) {
		if (sieve[pos] == 0xFF)
			continue;
		if (!(sieve[pos]&1))
		{
			p[2] = a + (pos << 4);
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&2))
		{
			p[2] = a + (pos << 4) + 2;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&4))
		{
			p[2] = a + (pos << 4) + 4;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&8))
		{
			p[2] = a + (pos << 4) + 6;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&16))
		{
			p[2] = a + (pos << 4) + 8;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&32))
		{
			p[2] = a + (pos << 4) + 10;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&64))
		{
			p[2] = a + (pos << 4) + 12;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (!(sieve[pos]&128))
		{
			p[2] = a + (pos << 4) + 14;
			if (p[2] > b) goto CLEANUP;
			closure_evalvoid(code);
			if (loop_break()) goto CLEANUP;
		}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			avma = btop;
	}
	
CLEANUP:
	pop_lex(1);
	pari_free(sieve);
	avma = ltop;
}


// Like forprime, but for b-a small. Assumes a is large -- should be above maxprime.
void
forthinprime(GEN ga, GEN gb, GEN code)
{
	pari_sp av = avma;
	long t = typ(ga);
	if (t == t_REAL || t == t_FRAC)
		ga = gceil(ga);
	else if (t != t_INT)
		pari_err_TYPE("forthinprime", ga);

	t = typ(gb);
	if (t == t_REAL || t == t_FRAC)
		gb = gfloor(gb);
	else if (t != t_INT)
		pari_err_TYPE("forthinprime", gb);

	if (cmpii(ga, gb) > 0) {
		avma = av;
		return;
	}
	
	ulong a = itou_or_0(ga);
	ulong b = itou_or_0(gb);
	
	if (b) {
		avma = av;
		forthinprimeuu(a, b, code);
		return;
	}
	
	/////////////////
}


// Like forprime, but for b-a small. Assumes a is large -- should be above maxprime.
void
forthinprimeuu(ulong a, ulong b, GEN code)
{
	pari_sp ltop = avma;
	long p[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3), 0};	// A GEN representing a prime to be passed to the code; its value is in p[2].
	push_lex((GEN)p,code);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);

	a |= 1;	// Make a odd (round up)
	b = (b-1)|1;	// Make b odd (round down)
pari_printf("%d to %d\n", a, b);
	for (; a <= b; a+= 2) {
		if(!uisprime(a))
			continue;
		p[2] = a;
		closure_evalvoid(code);
		if (loop_break())
			break;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			avma = btop;
	}
	pop_lex(1);
	avma = ltop;
}
