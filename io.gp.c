/******************************************************************************/
/**																	 I/O																		**/
/******************************************************************************/

GEN listtovec_shallow(GEN v)
{
	GEN x = list_data(v);
	long i = 1, lx = x ? lg(x): 1;
	GEN y = cgetg(lx, t_VEC);
	for (; i < lx; i++)
		gel(y,i) = gel(x,i);
	return y;
}


// Parse a zero-terminated line to find the second number:
// 20 10485876
// would yield "1048576". Note that the input string is modified in the process.
// TODO: Should check values more carefully -- does end-of-line happen when expected?
// TODO: Should check that the numbers are sequential
char*
getBValue(char* line)
{
	int start = 0;
	while (line[start] == ' ' || line[start] == '\t')
		start++;
	if (line[start] == '#')
		return NULL;	// Comment
	if (line[start] == '-')
		start++;
	while (line[start] >= '0' && line[start] <= '9')
		start++;
	while (line[start] == ' ' || line[start] == '\t')
		start++;
	int end = start;
	if (line[end] == '-')
		end++;
	while (line[end] >= '0' && line[end] <= '9')
		end++;
	if (start == end)
		return NULL;	// Blank line, or no numbers found
	line[end] = '\0';
	
	return line + start;
}

#define MAX_VECLEN 180000
#define MAX_LINELEN 5000
// Returns at most 18,000 values (warning the user if there are more left)
// Throws an error if it encounters an extremely long line -- b-files shouldn't
// have lines this long.
GEN
bfilein(char* name)
{
	FILE *f = fopen(name, "r");
	if (!f) {
#if defined(_WIN32) || defined(__CYGWIN32__)
		pari_err_FILE("input file", name);
#else
		if (strlen(name) == 11 && name[0] == 'b' && name[1] >= '0' && name[1] <= '9' && name[2] >= '0' && name[2] <= '9' && name[3] >= '0' && name[3] <= '9' && name[4] >= '0' && name[4] <= '9' && name[5] >= '0' && name[5] <= '9' && name[6] >= '0' && name[6] <= '9' && name[7] == '.' && name[8] == 't' && name[9] == 'x' && name[10] == 't') {
			char command[65];
			sprintf(command, "wget http://oeis.org/A%c%c%c%c%c%c/%s", name[1],name[2],name[3],name[4],name[5],name[6],name);
			int result = system(command);
			if (result == -1)
				pari_warn(warner, "Download failed.");
		}
		f = fopen(name, "r");
		if (!f)
			pari_err_FILE("input file", name);
#endif
	}
	
	GEN v = vectrunc_init(MAX_VECLEN + 1);
	char line[MAX_LINELEN];
	int i = 0;
	while(fgets(line, MAX_LINELEN, f) != NULL) {
		if (strlen(line) > MAX_LINELEN - 5)
			pari_err(e_MISC, "Maximum line length exceeded; b-file probably not valid");
		char* kept = getBValue(line);
		if (kept == NULL)
			continue;
		GEN value = strtoi(kept);
		if (++i > MAX_VECLEN) {
			pari_warn(warner, "only %d terms used; b-file has unread terms", MAX_VECLEN);
			break;
		}
		vectrunc_append(v, value);
	}
	fclose(f);
	return v;
}


GEN
bfile(GEN name, GEN v, GEN offset)
{
	pari_sp ltop = avma;
	GEN Anum = NEVER_USED;
	long cur = NEVER_USED;

	if (v) {
		switch(typ(v)) {
			case t_VEC:
			case t_COL:
				break;
			case t_LIST:
				v = listtovec_shallow(v);
				break;
			default:
				pari_err_TYPE("bfile", v);
		}
	}
	
	if (!offset)
		cur = 0;
	else if (typ(offset) == t_INT)
		cur = itos(offset) - 1;
	else
		pari_err_TYPE("bfile", offset);
	
	if (typ(name) == t_INT)
	{
		name = gtovec(GENtoGENstr(name));
		GEN p2;
		while (glength(name) < 6)
		{
			p2 = cgetg(2, t_VEC);
			gel(p2, 1) = strtoGENstr("0");
			name = concat(p2, name);
		}
		Anum = concat(name, NULL);	// "0","0","0","0","4","0" -> "000040"
		name = concat(concat(strtoGENstr("b"), Anum), strtoGENstr(".txt"));
	} else if (typ(name) == t_STR) {
		// TODO: Try to extract a reasonable A-number, or just set to blank?
		Anum = strtoGENstr("000000");
		//Anum = concat(extract0(gtovec(name), stoi(126), NULL), NULL);
	} else {
		pari_err_TYPE("bfile", name);
	}
	
	char* filename = GSTR(name);
	if (!v)
		return bfilein(filename);
	bfileout(filename, name, v, Anum, cur+1);
	avma = ltop;
	return gnil;
}


void
bfileout(char* filename, GEN name, GEN v, GEN Anum, long offset)
{
	FILE *f = fopen(filename, "r");
	if (f) {
		pari_warn(warner, "File `%Ps' already exists. Moving to bfile.old...", name);
		fclose(f);
		rename(filename, "bfile.old");
	}
	
	f = fopen(filename, "w+");
	long l1 = lg(v), cur = offset-1, i;
	for (i = 1; i < l1; ++i)
	{
		GEN e = gel(v, i);
		if (typ(e) != t_INT)
			pari_err_TYPE("bfile", e);
		if (digits(e) > 1000)
		{
			pari_warn(warner, "Next term has %ld digits; exiting.\n", digits(e));
			break;
		}

		char *num = GENtostr(e);
		fprintf(f, "%ld %s\n", ++cur, num);
		pari_free(num);
	}
	fclose(f);
	pari_printf("A%Ps: Terms %ld..%ld written to %s\n", Anum, offset, cur, filename);
}


GEN
fnice(GEN n)
{
	pari_sp ltop = avma;
	GEN f, s, s1, p1 = gen_0;
	long l2;
	if (typ(n) != t_INT)
		pari_err_TYPE("fnice", n);
	if (signe(n) < 0)
	{
		s = strtoGENstr("-");
		n = negi(n);
	} else {
		s = strtoGENstr("");
	}
	if (cmpis(n, 4) < 0)
	{
		p1 = concat(s, n);
		p1 = gerepileupto(ltop, p1);
		return p1;
	}
	f = Z_factor(n);
	s = Str(mkvec2(s, gcoeff(f, 1, 1)));
	if (itos(gcoeff(f, 1, 2)) > 1)
		s = Str(mkvec3(s, strtoGENstr("^"), gcoeff(f, 1, 2)));
	l2 = glength(gel(f, 1));

	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 2; i <= l2; ++i)
	{
		s1 = Str(mkvec2(strtoGENstr(" * "), gcoeff(f, i, 1)));
		if (!gequalgs(gcoeff(f, i, 2), 1))
			s1 = Str(mkvec3(s1, strtoGENstr("^"), gcoeff(f, i, 2)));
		s = Str(mkvec2(s, s1));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &s1, &s);
	}
	s = gerepileupto(ltop, s);
	return s;
}


/* Can handle only single-variable polynomials. */
GEN
tonice(GEN o, long prec)
{
	pari_sp ltop = avma;
	GEN s = gen_0, t = gen_0, t1 = gen_0, v = gen_0, p1 = gen_0, p2 = gen_0, p3 = gen_0;
	GEN p4 = gen_0;		/* genstr */
	if (typ(o) == t_POL)
	{
		t = content(o);
		o = gdiv(o, t);
		if (gequal1(denom(t)))
		{
			v = gpolvar(o);
			t = stoi(degree(o));
			t1 = polcoeff0(o, gtos(t), gvar(v));
			if (gcmpgs(t1, 0) < 0)
				p1 = Str(mkvec2(strtoGENstr("-"), monomialnice(gneg(t1), t, v)));
			else
				p1 = monomialnice(t1, t, v);
			s = p1;
			o = gsub(o, gmul(t1, gpow(v, t, prec)));
			t = stoi(poldegree(o, gvar(v)));
			{
				pari_sp btop = avma, st_lim = stack_lim(btop, 1);
				GEN p5 = gen_0;
				while (gcmpgs(t, 0) > 0)
				{
					t1 = polcoeff0(o, gtos(t), gvar(v));
					if (gcmpgs(t1, 0) > 0)
						p5 = strtoGENstr(" + ");
					else
						p5 = strtoGENstr(" - ");
					s = Str(mkvec3(s, p5, monomialnice(gabs(t1, prec), t, v)));
					o = gsub(o, gmul(t1, gpow(v, t, prec)));
					t = stoi(poldegree(o, gvar(v)));
					if (low_stack(st_lim, stack_lim(btop, 1)))
						gerepileall(btop, 5, &p5, &t1, &s, &o, &t);
				}
			}
			o = simplify(o);
			if (typ(o) != t_INT)
				pari_err_IMPL("multivariate polynomials in nice");
			o = polcoeff0(o, 0, -1);
			if (gcmpgs(o, 0) > 0)
				s = Str(mkvec3(s, strtoGENstr(" + "), o));
			if (gcmpgs(o, 0) < 0)
				s = Str(mkvec3(s, strtoGENstr(" - "), gneg(o)));
			s = gerepileupto(ltop, s);
			return s;
		}
		o = tonice(gmul(denom(t), o), prec);
		if (gequal1(numer(t)))
		{
			p2 = Str(mkvec4(strtoGENstr("("), o, strtoGENstr(")/"), ginv(t)));
			p2 = gerepileupto(ltop, p2);
			return p2;
		}
		p3 = Str(mkvecn(5, numer(t), strtoGENstr("("), o, strtoGENstr(")/"), denom(t)));
		p3 = gerepileupto(ltop, p3);
		return p3;
	}
	p4 = gcopy(GENtoGENstr(o));
	p4 = gerepileupto(ltop, p4);
	return p4;
}

GEN
initial(GEN n, char *s)
{
	pari_sp ltop = avma;
	char *l1;		/* str */
	GEN p2 = gen_0, p3 = gen_0;
	if (typ(n) != t_INT)
		pari_err_TYPE("initial", n);
	if (cmpis(n, 0) == 0)
	{
		l1 = "";
		avma = ltop;
		return strtoGENstr(l1);
	}
	if (equalim1(n))
	{
		p2 = Str(mkvec2(strtoGENstr("-"), strtoGENstr(s)));
		p2 = gerepileupto(ltop, p2);
		return p2;
	}
	if (equali1(n))
	{
		avma = ltop;
		return strtoGENstr(s);
	}
	p3 = Str(mkvec2(n, strtoGENstr(s)));
	p3 = gerepileupto(ltop, p3);
	return p3;
}


GEN
medial(GEN n, char *s)
{
	pari_sp ltop = avma;
	char *l1;		/* str */
	GEN p2 = gen_0, p3 = gen_0, p4 = gen_0, p5 = gen_0;
	if (typ(n) != t_INT)
		pari_err_TYPE("medial", n);
	if (cmpis(n, 0) == 0)
	{
		l1 = "";
		avma = ltop;
		return strtoGENstr(l1);
	}
	if (equalim1(n))
	{
		p2 = Str(mkvec2(strtoGENstr(" - "), strtoGENstr(s)));
		p2 = gerepileupto(ltop, p2);
		return p2;
	}
	if (equali1(n))
	{
		p3 = Str(mkvec2(strtoGENstr(" + "), strtoGENstr(s)));
		p3 = gerepileupto(ltop, p3);
		return p3;
	}
	if (cmpis(n, 0) < 0)
		p4 = strtoGENstr(" - ");
	else
		p4 = strtoGENstr(" + ");
	p5 = Str(mkvec3(p4, mpabs(n), strtoGENstr(s)));
	p5 = gerepileupto(ltop, p5);
	return p5;
}


/* Degree assumed to be positive */
GEN
monomialnice(GEN coeff, GEN degree, GEN v)
{
	pari_sp ltop = avma;
	GEN p1 = gen_0, p2 = gen_0;
	if (typ(coeff) != t_INT)
		pari_err_TYPE("monomialnice", coeff);
	if (typ(degree) != t_INT)
		pari_err_TYPE("monomialnice", degree);
	if (!v)
		v = pol_x(fetch_user_var("x"));
	if (equali1(coeff))
	{
		if (equali1(degree))
			p1 = gcopy(GENtoGENstr(v));
		else
			p1 = Str(mkvec3(v, strtoGENstr("^"), degree));
		p1 = gerepileupto(ltop, p1);
		return p1;
	}
	if (equali1(degree))
		p2 = Str(mkvec2(coeff, v));
	else
		p2 = Str(mkvec4(coeff, v, strtoGENstr("^"), degree));
	p2 = gerepileupto(ltop, p2);
	return p2;
}
