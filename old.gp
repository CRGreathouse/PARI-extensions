timervalue = default(timer);
default(timer, 0);


\\ ***************************************************************************************************
\\ *                                Reusable base-dependant functions                                *
\\ ***************************************************************************************************

ispal(n:int,B:int=10)={
	my(tmp:int, d:int);
	if (B == 10,
		d = digits(n) - 2
	,
		d = log(n) \ log(B) - 1
	);
	while(d > 0,
		tmp = n%B;
		n = (n - tmp) / B;
		if (n\B^d != tmp, return(0));
		n = n % (B^d);
		d = d - 2;
	);
	d < 0 || n%(B+1)==0
};
addhelp(ispal, "ispal(n): Is n a palindrome? Sloane's A136522; characteristic function of Sloane's A002113.");


term(n)={
	if(type(n) == "t_INT", return(1));
	if(type(n) != "t_FRAC", error("term: Not a fraction!"));
	n = oddres(denominator(n));
	if (n < 25, return(n == 1 || n == 5));
	ispower(n,, &n);
	n == 5
};
addhelp(term, "term(n): Is n a terminating decimal?");


digitsToNum(v:vec,B:int=10)={
	my(s=v[1]);
	for(i=2,#v,
		s = B*s + v[i]
	);
	s
};
addhelp(digitsToNum, "digitsToNum(v, {B=10}): Transforms a vector of base-B digits into the number it represents.");


numToDigits(n:int, B:int=10)={
	if(B == 10,
		eval(Vec(Str(n)))
	,
		my(v=[]);
		while(n>B,
			v = concat(n%B,v);
			n \= B
		);
		concat(n,v)
	)
};
addhelp(numToDigits, "numToDigits(n): Transforms a number to a vector of the digits.");


glue(a:int,b:int)={
	a*10^digits(b)+b
};
addhelp(glue, "glue(a, b): Returns the (decimal) concatenation of a and b.");


rev(n:int,B=10)={
	my(m=n%B);
	n\=B;
	while(n>0,
		m=m*B+n%B;
		n\=B
	);
	m
};
addhelp(rev, "rev(n, {B=10}): Base-B digit reversal of n. Sloane's A004086.");


hexTable=["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"];
hex(n:int)={
	my(res="");
	if (n <= 0,
		if (n == 0, return("0"));
		n = -n;
		res = "-";
	);
	while (n > 0,
		res = concat(hexTable[bitand(n, 15) + 1], res);
		n >>= 4
	);
	res
};
addhelp(hex, "hex(n): Hexadecimal representation of n.");


ssdTable=[1,4,9,16,25,36,49,64,81,1,2,5,10,17,26,37,50,65,82,4,5,8,13,20,29,40,53,68,85,9,10,13,18,25,34,45,58,73,90,16,17,20,25,32,41,52,65,80,97,25,26,29,34,41,50,61,74,89,106,36,37,40,45,52,61,72,85,100,117,49,50,53,58,65,74,85,98,113,130,64,65,68,73,80,89,100,113,128,145,81,82,85,90,97,106,117,130,145,162];
ssd(n:int)={
	my(s=0,tmp);
	while(n>9,
		tmp = n%100;
		s += if (tmp < 2, tmp, ssdTable[tmp]);
		n \= 100
	);
	s+n^2
};
addhelp(ssd, "ssd(n): Sum of squared digits of n. Sloane's A003132.");


dprod(n:int)={
	my(t=1,k);
	while(n,
		k=n%10;
		if(k,t*=k,return(0));
		n=(n-k)/10
	);
	t
};
addhelp(dprod, "dprod(n): Product of digits of n. Sloane's A007954.");


rcd(n:int)={
	my(d=divisors(n));
	n=1;
	for(i=2,#d,n=glue(d[i],n));
	n
};
addhelp(rcd, "rcd(n): Reverse concatenation of the divisors of n. Sloane's A176558.");


fcd(n:int)={
	my(d=divisors(n));
	n=0;
	for(i=1,#d,n=glue(n,d[i]));
	n
};
addhelp(rcd, "fcd(n): (Forward) concatenation of the divisors of n. Sloane's A037278.");


\\ ***************************************************************************************************
\\ *                                   Obscure base-related functions                                *
\\ ***************************************************************************************************

res(f)={
	if(type(f)!="t_POL",error("need poly"));if(poldegree(f)!=4,error("bad deg"));my(u=polroots(f),r);r=[u[1]*u[2]+u[3]*u[4],u[1]*u[3]+u[2]*u[4],u[1]*u[4]+u[2]*u[3]];print(r[1]"\n"r[2]"\n"r[3]);(x-r[1])*(x-r[2])*(x-r[3])
};
addhelp(res, "Resolvant cubic");


has666(n)={
	my(tmp);
	while(n>665,
		tmp = n%1000;
		if (tmp > 699 || tmp < 600,
			n = (n-tmp)/1000
		,
			if (tmp == 666, return(1));
			n\=10
		)
	);
	0
};
addhelp(has666, "has666(n): Does n contain 666 as a substring? Characteristic function of Sloane's A051003.");


noPalSub(n)={
	my(d);
	local(digit);
	digit=numToDigits(n);
	d = #digit;

	for(len=2,d,
		for(i=1,d-len+1,
			if(isPalSub(i,len), return(0))
		)
	);
	1
};
addhelp(noPalSub, "noPalSub(n): Does n lack a palindromic substring?");


\\ Uses local variable digit, a vector of digits.
isPalSub(start,len)={
	my(b=start-1,e=start+len);
	for(j=1,len>>1,
		if(digit[b+j] != digit[e-j], return(0))
	);
	1
};
addhelp(isPalSub, "isPalSub(start, len): Helper function for noPalSub.");


diffdigits(n:int)={
	my(v=[],sz=#Str(n));
	while(n>9,
		v=concat(v,n%10);
		n\=10
	);
	#vecsort(concat(v,n),,8) == sz
};
addhelp(diffdigits, "diffdigits(n): Are the digits of n distinct? Characteristic function of Sloane's A010784.");


PHLenDigits=3;
PHLen=10^PHLenDigits;
panHelper(n)={
	my(t=0);
	for(i=1,PHLenDigits,
		t=bitor(t,1 << (n%10));
		n\=10
	);
	t
};
if(PH=='PH || #PH < #PHLen, PH=vector(PHLen,i,panHelper(i-1)));
pan(n)={
	my(t=0);
	while(n>=PHLen,
		t=bitor(t,PH[n%PHLen + 1]);
		n \= PHLen
	);
	while(n > 9,
		t=bitor(t,1 << (n%10));
		n \= 10
	);
	bitor(t, 1 << n) == 1023
};
addhelp(pan, "pan(n): Is n pandigital (repeated digits allowed)? Characteristic function of Sloane's A050278.");


pan10(n:int)={
	my(v);
	if(n<1023456789 || n > 9876543210, return(0));
	v=vector(10);
	for(i=1,10,
		v[i]=n%10;
		n\=10
	);
	#vecsort(v,,8) == 10
};
addhelp(pan10, "pan10(n): Is the number a 0-to-9 pandigital (exactly 10 distinct digits)? Characteristic function of a finite subset of Sloane's A050278.");


pan9(n:int)={
	my(v);
	if(n<123456789 || n > 987654321, return(0));
	v=vector(9);
	for(i=1,9,
		v[i]=n%10;
		n\=10
	);
	v=vecsort(v,,8);
	#v == 9 && v[1]==1
};
addhelp(pan9, "pan9(n): Is the number a 1-to-9 pandigital (exactly 9 distinct digits, no 0s)? A finite subset of zeroless pandigital numbers, Sloane's A050289.");


\\ ***************************************************************************************************
\\ *                          Obscure (but otherwise reusable) functions                             *
\\ ***************************************************************************************************

\\ These functions relate to A218459.
docheck(d,p)=for(y=1,sqrtint(p\d),if(issquare(p-d*y^2),return(1)));0;
dochk(p)=my(d);while(!docheck(d++,p),);d;
/*
\\ v=apply(do,primes(1000));
\\ bnfisintnorm(bnfinit('x^2+93),2903)
verify(d,p)=for(y=1,sqrtint(p\d),if(issquare(p-d*y^2),return(1)));0;
dv(d,p)={
	if(#DV<9,DV=vector(1000));
	my(B);
	B=if(d<=#DV,
		if(DV[d],DV[d],DV[d]=bnfinit('x^2+d))
	,
		bnfinit('x^2+d)
	);
	if(#bnfisintnorm(B,p),
		if(verify(d,p),
			1
		,
			\\print("Funny stuff at p = "p": bnfisintnorm sez "d" but it ain't so.");
			0
		)
	,
		0
	)
};
dv_k(d,p)={
	kronecker(-d,p)>0 && dv(d,p)
};
dv_k_s(d,p)={
	kronecker(-d,p)>=0 && issquarefree(d) && dv(d,p)
};
\\ Values known to occur: 
do(p)={
	if(p%24<23,return(if(p%4<3,1,if(p%8==3,2,3))));
	if(kronecker(p,7)>0,return(7));
	if(dv_k(11,p), return(11));
	if(dv_k(19,p), return(19));
	if(kronecker(p,11)>0,return(22));
	if(dv_k(23,p), return(23));
	if(dv_k(26,p), return(26));	\\ probably never happens
	if(dv_k(29,p), return(29)); \\ probably never happens
	if(dv_k(31,p), return(31));
	for(d=34,p,
		if(dv_k_s(d,p), return(d))
	)
};
*/


testCI(n,alpha,b:small)={
	my(t,s);
	s=sum(i=1,n,
		t=rand(b);
		sigfac(t)^2/sigfac(sqfac(t)) > 2
	);
	print(s" found; "precision(100*(1-alpha),9)"% confidence interval:");
	binomialInterval(n\1,s,alpha)
};


\\ Returns a factorization matrix from a vector of prime factors.
list2fac(f:vec)={
	if(!#f,return(factor(1)));
	f=vecsort(f);
	my(p=vecsort(f,,8),o=#p,M=matrix(o,2,i,j,if(j==1,p[i],0)),n=1);
	M[1,2]=1;
	for(i=2,#f,
		if(f[i]!=f[i-1],n++);
		M[n,2]++
	);
	M
};


\\ Returns sigma(factorback(f)).
sigfac(f:vec)={
	my(t);
	prod(i=1,#f[,1],
		if(f[i,2]==1,
			f[i,1]+1
		,
			t=1/f[i,1];
			sum(j=0,f[i,2],t*=f[i,1])
		)
	)
};


\\ Returns factor(factorback(f)^2).
sqfac(f:vec)={
	for(i=1,#f[,1],
		f[i,2]*=2
	);
	f
};


numinvphi(n,m=1)={
	my(s=0,k,p);
	if(n==1,return(1+(m<2)));
	fordiv(n,d,
		if(d<m,next);
		if(!isprime(p=d+1),next);
		k=n\d;
		while(1,
			s+=numinvphi(k,p);
			if(k%p,break);
			k\=p
		)
	);
	s
};
addhelp(numinvphi, "numinvphi(n): Number of m with eulerphi(m) = n. Sloane's A014197.");


\\ Lazily evaluate the number of labeled acyclic digraphs.
dag(t:int)={
	my(n:int);
	if (t<1, return(t==0));
	if (DAG == 'DAG, DAG = [1]);
	while(t > #DAG,
		n = #DAG + 1;
		DAG=concat(DAG,-sum(k=1,n-1,
			(-1)^k * binomial(n, k) * DAG[n-k] << (k * (n-k))
		)-(-1)^n)
	);
	DAG[t]
};
addhelp(dag, "dag(n): Number of labeled acyclic digraphs with n labeled nodes. Sloane's A003024.");


\\ Estimates dag(t) based on R. P. Stanley (2006)
estdag(t)={
	1.74106112529322984034218802800555555204450538663852601463*2^binomial(t,2)*t!*
	1.4880785455997102946562460315823576618935161526029908077497^(-t)
};
addhelp(estdag, "estdag(n): Estimates the number of labeled acyclic digraphs with n labeled nodes. See dag(n).");


sumOfNext(n:int)={
	substpol(Faulhaber(n,x),x,x*(x+1)/2)-substpol(Faulhaber(n,x),x,x*(x-1)/2)
};
addhelp(sumOfNext, "sumOfNext(n): 'Sum of next x nth powers.' Related to Sloane's A006003, A072474, A075664, A075665, A075666, A075667, A075668, A075669, A075670, and A075671.");


sindex(n)={
	2*log(core(n,1)[2])/log(n)
};
addhelp(sindex,"Square index of n: 1 if n is a square, 0 if n is squarefree, otherwise strictly between those. NBRTHRY 2008-06-28.");


\\ Sloane's A002386, primes starting a record prime gap
lgap=[2,3,7,23,89,113,523,887,1129,1327,9551,15683,19609,31397,155921,360653,370261,492113,1349533,1357201,2010733,4652353,17051707,20831323,47326693,122164747,189695659,191912783,387096133,436273009,1294268491,1453168141,2300942549,3842610773,4302407359,10726904659,20678048297,22367084959,25056082087,42652618343,127976334671,182226896239,241160624143,297501075799,303371455241,304599508537,416608695821,461690510011,614487453523,738832927927,1346294310749,1408695493609,1968188556461,2614941710599,7177162611713,13829048559701,19581334192423,42842283925351,90874329411493,171231342420521,218209405436543,1189459969825483,1686994940955803,1693182318746371,43841547845541059,55350776431903243,80873624627234849,203986478517455989,218034721194214273,305405826521087869,352521223451364323,401429925999153707,418032645936712127,804212830686677669];


trinv(n:int)={
	(1+sqrtint(1+n<<3))>>1
};
addhelp(trinv, "trinv(n): Triangle inverse function, used for mapping tables to lists. Sloane's A002024 (with a different offset).");


Bspsp(n:int)={
	my(f=factor(n)[,1], nu=valuation(f[1]-1, 2), nn = star(n));
	for(i=2,#f,
		nu = min(nu, valuation(f[i] - 1, 2));
	);
	(1 + (2^(#f * nu) - 1) / (2^#f - 1)) * prod(i=1, #f, gcd(nn, star(f[i]))) - 1
};
addhelp(Bspsp, "Bspsp(n): To how many bases 1 < a < n is n an a-strong pseudoprime?");


star(n:int)={
	oddres(n--)
};
addhelp(star, "star(n): Helper function for Bspsp.");


Bpsp(n:int)={
	my(f=factor(n)[,1]);
	n--;
	prod(i=1,#f,gcd(n,f[i]-1))-1
};
addhelp(Bpsp, "Bpsp(n): To how many bases 1 < a < n is n an a-pseudoprime?");


normalNum(n:int, b:int=10)={
	my(expected=n/b,err=sqrt(2*n*log(log(n))*(1-1/b)/b));
	[floor(max(0, expected-err)), ceil(expected+err)]
};
addhelp(normalNum, "normalNum(n, b=10): Given n digits of a b-normal number, in what range can each digit be expected to occur? This result is due to Khintchine.");


\\ These functions are for constructing a list of numbers not the sum of a prime and two positive powers of two.
\\ Related to Sloane's A156695.
testA156695(n)={
	my(i=0);
	while(p22[i++]<n,
		if(isprime(n-p22[i]), return(0))
	);
	1
};
list(a,b)={
	if(p22=='p22,
		print1("Initializing vectors...");
		p2=vector(80, n, 2^n);
		p22=sumset(p2, p2);
		p22=vector(2000, i, p22[i]);
		print("done.")
	);
	forstep(n=bitor(floor(a),1),b,2,
		if(testA156695(n), print1(", "n))
	)
};
list1(a,b)={
	list(2,2);	\\ initialize, if needed
	print("WARNING: Using the conjecture that members of this sequence are 255 mod 510.");
	forstep(n=(a\510)*510+255,b,510,
		if(testA156695(n),print1(", "n))
	)
};


RhoBounds(upper, n)={
	local(rhoupper, rholower, s1, s2, s3);
	\\ TODO: Test default(parisize) and default(precision).
	n = floor(n);
	rholower=vector(2*n*upper);
	rhoupper=vector(n*upper);
	for(i=n,3*n,
		rhoupper[i] = DickmanRho(i/n);
		rholower[2*i] = rhoupper[i];
		rholower[2*i+1] = DickmanRho((i + .5)/n);
	);
	
	\\ Running tallies introduce substantial error, but the cost of an extra
	\\ limb is tiny compared to the cost of redoing each sum at every step.
	\\ At (8, 15000) the time savings is > 200-fold.
	s1 = sum(k=2*n+1, 3*n-1, rhoupper[k]);
	s2 = sum(k=2*n+1, 3*n, rholower[2*k - 2]);
	s3 = sum(k=2*n+1, 3*n, rholower[2*k - 1]);
	for(i=3*n+1,n*upper,
		s1 = s1 + rhoupper[i-1] - rhoupper[i-n];
		rhoupper[i] = (rhoupper[i-n] + 2 * s1) / (2 * i - 1);
		
		s2 = s2 + rholower[2*i - 2] - rholower[2*(i - n) - 2];
		rholower[2*i - 1] = s2 / (i - .5);
		
		s3 = s3 + rholower[2*i - 1] - rholower[2*(i - n) - 1];
		rholower[2*i] = s3 / i;
	);
	print(rholower[#rholower]);
	print(rholower[#rholower] * 2/3 + rhoupper[#rhoupper] / 3);
	print(rhoupper[#rhoupper]);
	\\[rholower, rhoupper]
};
addhelp(RhoBounds, "Bounds for Dickman's rho function at the given point and with the given number of intervals. Displays lower and upper bounds, plus a 'most likely' value between the two.  Returns an evenly-spaced array up to the requested value.  WARNING: if precision is too low (below 30?), the answers may lose significance (even all significance, if precision is low enough). Always test results!");


residues(n,e)={
	my(v=vector(n));
	for(i=1,n,v[i]=lift(Mod(i,n)^e));
	vecsort(v,,8)
};
addhelp(residues,"residues(n,e): Lists the e-th powers mod n.");


isA129912(n)={
	local(l,old);
	if (n < 3, return (n > 0));
	old = 0;
	while(bitand(n,1) == 0,
		n >>= 1;
		l = 2;
		forprime(p=3,9999999, \\ Good up to 1e4340851
			if (n%p > 0, break());
			n = n / p;
			l = p;
		);
		if (old == l, return(0));
		old = l;
	);
	n == 1;
};
addhelp(isA129912, "Is the number a product of distinct primorials? Characteristic function of Sloane's A129912.");


cumPrimeFactorSum(n)={
	my(s=primepi(n)-primepi(n\2),t);
	forprime(p=2,sqrtint(n),
		t = n\p;
		while(t>0,
			s += t;
			t \= p
		)
	);
	forprime(p=sqrtint(n)+1,n\2,
		s += n\p;
	);
	s
};
addhelp(cumPrimeFactorSum, "cumPrimeFactorSum(n): The sum of the number of prime factors of 1..n.  A faster version of sum(k=1,n,bigomega(k)). Sloane's A022559; related to A071811.");


\\ ***************************************************************************************************
\\ *                                Functions that aren't really needed                              *
\\ ***************************************************************************************************

secant_root(ff,a,b)={
	e = eps() * 2;
	aval=ff(a);
	bval=ff(b);
	while (abs(bval) > e,
		oldb = b;
		b = b - (b - a)/(bval - aval) * bval;
		aval = bval;
		bval = ff(b);
		a = oldb
	);
	b
};
addhelp(secant_root, "secant_root(ff,a,b): Finds a root of ff between a and b using the secant method.");


bisect(ff,a,b)={
	local(e,aval,bval,m,mval);
	e = eps() * 2;
	aval=ff(a);
	bval=ff(b);
	if (aval == 0, return(a));
	if (bval == 0, return(b));
	if ((aval > 0) == (bval > 0), error("bisect: No sign change!  Choose left and right values so there is a sign change."));
	
	m = (a + b)/2.;
	mval = ff(m);
	if(aval<0,
		while (abs(mval) > e,
			if (mval < 0, a=m;aval=mval, b=m;bval=mval);
			m = (a + b)/2.;
			mval = ff(m);
		)
		,	
		while (abs(mval) > e,
			if (mval > 0, a=m;aval=mval, b=m;bval=mval);
			m = (a + b)/2.;
			mval = ff(m);
		)
	);
	m
};
addhelp(bisect, "bisect(ff,a,b): Finds a root of ff between a and b using bisection.");


\\Euler + log(x) + suminf(n=1,x^n/(n*n!))
\\Euler + x + suminf(n=1,x^n/(n*n!))
Ei(x)={
	-eint1(-x)
\\	local(t,s,o);
\\	s = 1;
\\	t = Euler + x;
\\	for(i=1,9999999,
\\		o = t;
\\		s = s * x / i;
\\		t = t + s / i;
\\		if (o == t, break());
\\	);
\\	t
};
addhelp(Ei, "Exponential integral.");


Ei0(x)={
	local(n,t);
	n = round(x);
	t=1.;
	s=1.;
	for(m=1,n-1,
		s=s*m/x;
		t=t+s
	);
	t*exp(x)/x
};
addhelp(Ei0, "The 'usual' asymptotic approximation to the exponential integral.");


Ei2(x)={
	Ei0(x) - sqrt(2*Pi/x)*(1/3 + round(x) - x)
};
addhelp(Ei2, "A fast, accurate approximation of the exponential integral for real arguments\nusing a method due to Donald Wadsworth, ``Improved Asymptotic Expansion for the Exponential Integral with Positive Argument'' (~1964)");


predicate(ff, lim)={
	my(v=[]);
	for(i=1,lim,
		if(ff(i), v=concat(v,i))
	);
	v
};
addhelp(predicate, "predicate(ff, lim): Returns a vector containing those k in 1..lim for which ff holds.");


abc.pr={
	my(v=[]);
	if (type(abc) != "t_VEC",
		error("Not a vector!")
	);
	for(i=1,#abc,if(isprime(abc[i]),v=concat(v,abc[i])));
	v
};


abc.prat={
	my(v=[]);
	if (type(abc) != "t_VEC",
		error("Not a vector!")
	);
	for(i=1,#abc,if(isprime(abc[i]),v=concat(v,i)));
	v
};


quad(a,b,c)={
	local(d,s1,s2,g);
	d=b^2-4*a*c;
	if (d == 0,
		b=-b/2/a;
		print("Solution: " b);return()
	);
	if (d < 0,
		print("Note: complex answer!")
	);
	a=2*a;
	s1=-b/a;
	g=sqrtint(largestSquareFactor(abs(d)));
	d/=g^2;
	s2=gcd(g,a);
	a/=s2;
	g/=s2;
	if(g==1,
		if(a==1,
			print("Solutions: "s1" +/- sqrt("d")");
		,
			print("Solutions: "s1" +/- sqrt("d")/"a);
		)
	,
		if(a==1,
			print("Solutions: "s1" +/- "g"sqrt("d")");
		,
			print("Solutions: "s1" +/- "g"sqrt("d")/"a);
		);
	);
};
addhelp(quad, "quad(a,b,c): Solves the quadratic equation ax^2 + bx + c in a simple fashion appropriate for high school algebra.");


default(timer, timervalue);
