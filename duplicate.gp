timervalue = default(timer);
default(timer, 0);

\\ This file is a collection of functions that are already implemented in
\\ auto.gp.c.  Occasionally it is useful to see a script for these functions,
\\ for example giving code to others who cannot load the .run files.


dsumTable=[1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,11,3,4,5,6,7,8,9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18];
dsum(n:int)={
	my(s=0,tmp);
	while(n>9,
		tmp = n%100;
		s += if (tmp < 10, tmp, dsumTable[tmp]);
		n \= 100
	);
	s+n
};
addhelp(dsum, "dsum(n): Digit sum of n. Sloane's A007953.");


\\ Really bad performance here: >= 32 bits per number in sieve, about 100 times
\\ the size of a basic wheel mod 6
forbigprime(from,to,ff)={
	my(bits, sq, lim = default(primelimit));
	
	\\ Small enough for forprime
	if (to <= lim,
		forprime(p=from,to,ff(p));
		return()
	);
	
	\\ Too big for segmented sieve, go one-at-a-time
	if (to > lim^2,
		if (default(debug) >= 2,
			print("forbigprime: testing primes individually");
		);
		from -= 1;
		while (from < to,
			from = nextprime(from+1);
			ff(from)
		)
	);
	
	\\ Remove all small primes and set variables
	if (from <= lim,
		forprime(p=from,lim,ff(p));
		from = nextprime(lim + 1)
	);
	if (bigprimechunk == 'bigprimechunk, bigprimechunk = 10^6);
	sq = precprime(sqrt(to));
	from = ceil(from - 1);
	to = floor(to);

	\\ Segment loop
	while(from < to,
		bits = vector(min(bigprimechunk, to - from), n, 1);
		forprime(p=2,sq,
			forstep(i = p - (from%p), #bits, p,
				bits[i] = 0;
			)
		);
		for(i=1,#bits,
			if(bits[i], ff(from + i))
		);
		from += bigprimechunk;
	)
};
addhelp(forbigprime, "forbigprime(X=a,b,seq): the sequence is evaluated, X running over the primes between a and b.");


graeffe(f:pol)={
	my(d=poldegree(f),g=vector(d\2+1,i,polcoeff(f,2*i-2)),h=vector((d+1)\2,i,polcoeff(f,2*i-1)));
	Polrev(g)^2-x*Polrev(h)^2
};
addhelp(graeffe, "graeffe(f): Assumes that f is a nonzero t_POL with t_INT coefficients.  Returns a polynomial with roots that are precisely the squares of the roots of f.");


poliscyclo(f:pol)={
	my(f1=graeffe(f),f2,fn);
	if(f==f1,return(1));
	fn=subst(f,x,-x);
	if(f1==fn&iscyclo(fn),return(1));
	issquare(f1,&f2)&iscyclo(f2)
};
addhelp(poliscyclo, "poliscyclo(f): Is f a cyclotomic polynomial?  Uses the Bradford-Davenport algorithm.");


\\addhelp(poliscycloproduct, "poliscycloproduct(f, {flag=0}): Is f a product of distinct cyclotomic polynomials?  If flag is 1, return the least n such that f | x^n-1.");


istotient(n:int)={
	my(k:int,p:int,d:int);
	if(n<2,return(n==1));
	if(n%2,return(0));
	k=n;
	while(1,
		if(totientHelper(k,2), return(1));
		if(k%2,break);
		k\=2
	);
	fordiv(n>>1,dd,
    	d=dd<<1;
		if(!isprime(p=d+1),next);
		k=n\d;
		while(1,
			if(totientHelper(k,p), return(1));
			if(k%p,break);
			k\=p
		)
	);
	0
};
addhelp(istotient, "istotient(n): Does there exist some m such that eulerphi(m) = n?");


totientHelper(n:int,m:int=1)={
	my(k:int,p:int,d:int);
	if(n==1,return(1));
	if(n%2,return(0));
	fordiv(n>>1,dd,
    	d=dd<<1;
		if(d<m|!isprime(p=d+1),next);
		k=n\d;
		while(1,
			if(totientHelper(k,p), return(1));
			if(k%p,break);
			k\=p
		)
	);
	0
};
addhelp(totientHelper, "totientHelper(n,m): Helper function for totient.");


Eng(n:int)={
	my(tmp, s="");
	if (n < 1000,
		if (n > 0, return(Eng(n)));
		if (n == 0, return("zero"));
		error("Can't do negatives")
	);
	tmp = n \ 1000000000000000;
	if (tmp,
		if (tmp >= 1000, error("Too big"));
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(tmp), " quadrillion");
		n -= tmp * 1000000000000000
	);
	tmp = n \ 1000000000000;
	if (tmp,
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(tmp), " trillion");
		n -= tmp * 1000000000000
	);
	tmp = n \ 1000000000;
	if (tmp,
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(tmp), " billion");
		n -= tmp * 1000000000
	);
	tmp = n \ 1000000;
	if (tmp,
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(tmp), " million");
		n -= tmp * 1000000
	);
	tmp = n \ 1000;
	if (tmp,
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(tmp), " thousand");
		n -= tmp * 1000
	);
	if (n,
		if (#s, s = Str(s, ", "));
		s = Str(s, Eng(n))
	);
	s
};
addhelp(Eng, "Eng(n): English name of the number n.");


fibmod(n:int, m:int)={
	((Mod([1,1;1,0],m))^n)[1,2]
	/*
	my(f,t);
	if (m < 6,
		if (m < 0, m = -m);
		if (m == 0, error("division by 0"));
		if (m == 1, return(0));
		if (m == 2, return(n%3>0));
		if (m == 3, return(fibonacci(n%8)%3));	\\ 0 1 1 2 0 2 2 1
		if (m == 4, return(fibonacci(n%6)%4));	\\ 0 1 1 2 3 1
		if (m == 5, return(fibonacci(n%20)%5));	\\ 0 1 1 2 3 0 3 3 1 4 0 4 4 3 2 0 2 2 4 1
	);
	if (n < 1000,
		return(fibonacci(n)%m);
	);
	f = factor(m);
	t = Mod(fibonacci(n%Pisano(f[1,1], f[1,2])),f[1,1]^f[1,2]);
	for(i=2,#f[,1],
		t = chinese(t,Mod(fibonacci(n%Pisano(f[i,1], f[i,2])),f[i,1]^f[i,2]));
	);
	lift(t)
	*/
};
addhelp(fibmod,"fibmod(n,m): Returns the nth Fibonacci number mod m. Same as finonacci(n)%m, but faster for large n.");


\\ Helper function for fibmod; returns a multiple of the period of Fibonacci numbers mod p^e.
\\ Assumes p is prime and e > 0.
Pisano(p:int,e:small)={
	my(t);
	if(p==5, return(4 * 5^e));	\\ powuu? powis?
	t = p%5;
	if (t == 2 || t == 3,
		t = 2 * (p + 1)
	,
		t = p - 1
	);
	t * p^(e-1)
};
addhelp(Pisano, "Pisano(p, e): Returns a multiple of the period of the Fibonacci numbers mod p^e. Assumes p is prime and e is a (small) positive integer.");


\\ Version 2... worse, I think?
/*
Pisano(p:int)={
	if (!isprime(p), return(0));
	my(m=p%5,m1,f);
	if(!m, return(20));
	m = if (m == 2 || m == 3,2 * (p + 1),p - 1);
	f = factor(m);
	for (i = 1, #f[,1],
		m1 = m / f[i,1];
		while (fibonacci(m1)%p == 0 && fibonacci(m1 - 1)%p == 1,
			m = m1;
			m1 /= f[i,1];
			if (denominator(m1) > 1, break)
		)
	);
	m
};*/

\\ ***************************************************************************************************
\\ *				Prime-related arithmetic functions				*
\\ ***************************************************************************************************

rad(n:int)={
	if(n < 0, n = -n);
	if (issquarefree(n), return(n));
	vecprod(factor(n)[,1])
};
addhelp(rad, "rad(n): Radical of n, the largest squarefree number dividing n. Sloane's A007947.");


isprimepower(n:int)={
	ispower(n,,&n);
	isprime(n)
};
addhelp(isprimepower, "isprimepower(n): Is n a prime power? Characteristic sequence of Sloane's A000961.");


isPowerful(n:int)={
	if (bitand(n,3)==2,return(0));
	if (n%3 == 0 && n%9 > 0, return(0));
	if (n%5 == 0 && n%25 > 0, return(0));
	vecmin(factor(n)[,2])>1
};
addhelp(isPowerful, "isPowerful(n): Is n powerful (min exponent 2)? Sloane's A112526; characteristic function of Sloane's A001694.");


prp(n:int,b:int=2)={
	Mod(b,n)^(n-1)==1
};
addhelp(prp, "prp(n,b=2): Is n a b-probable prime?");


sprp(n:int,b:int=2)={
	my(d = n, s = 0);
	until(bitand(d,1), d >>= 1; s++);
	d = Mod(b, n)^d;
	if (d == 1, return(1));
	for(i=1,s-1,
		if (d == -1, return(1));
		d = d^2;
	);
	d == -1
};
addhelp(sprp, "sprp(n,b=2): Is n a b-strong probable prime?");


sopfr(n:int)={
	my(f=factor(n));
	sum(i=1,#f[,1],f[i,1]*f[i,2])
};
addhelp(sopfr, "sopfr(n): Sum of prime factors of n (with multiplicity). Sloane's A001414.");


sopf(n:int)={
	my(f=factor(n)[,1]);
	sum(i=1,#f,f[i])
};
addhelp(sopf, "sopf(n): Sum of distinct prime factors of n. Sloane's A008472.");


primorial(n)={
	my(t1=1,t2=1,t3=1,t4=1);
	forprime(p=2,n>>3,t1=t1*p);
	forprime(p=(n>>3)+1,n>>2,t2=t2*p);
	t1=t1*t2;t2=1;
	forprime(p=(n>>2)+1,(3*n)>>3,t2=t2*p);
	forprime(p=((3*n)>>3)+1,n>>1,t3=t3*p);
	t1=t1*(t2*t3);t2=1;t3=1;
	forprime(p=(n>>1)+1,(5*n)>>3,t2=t2*p);
	forprime(p=((5*n)>>3)+1,(3*n)>>2,t3=t3*p);
	t2=t2*t3;t3=1;
	forprime(p=((3*n)>>2)+1,(7*n)>>3,t3=t3*p);
	forprime(p=((7*n)>>3)+1,n,t4=t4*p);
	t1*(t2*(t3*t4))
};
addhelp(primorial, "Returns the product of each prime less than or equal to n. Sloane's A034386.");


gpf(n:int)={
	my(f=factor(n)[,1]);
	if (n <= 1,
		if (n >= -1, return (1));	\\ Convention
		n = -n
	);
	f[#f]
};
addhelp(gpf, "The greatest prime factor of a number. Sloane's A006530.");


lpf(n:int)={
	if (n <= 1,
		if (n >= -1, return (1));	\\ Convention
		n = -n
	);
	forprime(p=2,2e4,
		if(n%p==0,return(p));
	);
	factor(n)[1,1]
};
addhelp(lpf,"The least prime factor of a number. Sloane's A020639.");


\\ ***************************************************************************************************
\\ *				Other arithmetic functions					*
\\ ***************************************************************************************************

isTriangular(n:int)={
	issquare(n<<3+1)
};
addhelp(isTriangular, "isTriangular(n): Is n a triangular number? Sloane's A010054.");


isFibonacci(n:int)={
	my(k=n^2);
	k+=((k + 1) << 2);
	issquare(k) || (n > 0 && issquare(k-8))
};
addhelp(isFibonacci, "isFibonacci(n): Is n a Fibonacci number? Sloane's A010056; characteristic function of Sloane's A000045.");


largestSquareFactor(n:int)={
	my(f=factor(n));
	prod(i=1,#f[,1],f[i,1]^(f[i,2]>>1))^2
};
addhelp(largestSquareFactor, "largestSquareFactor(n): Largest square dividing n. Sloane's A008833.");


hamming(n:int)={
	my(v=binary(n));
	sum(i=1,#v,v[i])
};
addhelp(hamming, "hamming(n): Hamming weight of n (considered as a binary number). Sloane's A000120.");


istwo(n:int)={
	my(f);
	if (n < 3,
		return (n >= 0);
	);
	f = factor(oddres(n));
	for (i=1,#f[,1],
		if(bitand(f[i,2],1)==1 && bitand(f[i, 1], 3)==3, return(0))
	);
	1
};
addhelp(istwo, "Is the number a sum of two squares? Characteristic function of Sloane's A001481.");


ways2(n:int)={
	\\sum(k=0,sqrtint(n>>1),issquare(n-k^2))
	my(f=factor(oddres(n)),t=1);
	for(i=1,#f[,1],
		if (f[i,1]%4==3,
			if (f[i,2]%2, return(0))
		,
			t *= f[i,2] + 1;
		)
	);
	(t + 1) >> 1
};
addhelp(ways2, "Number of ways that n can be represented as a sum of two squares. Sloane's A000161.");


isthree(n:int)={
	my(tmp=valuation(n, 2));
	bitand(tmp, 1) || bitand(n >> tmp, 7) != 7
};
addhelp(isthree, "isthree(n): Is the number the sum of three squares? Sloane's A071374; characteristic function of Sloane's A000378.");


ways3(n:int)={
	my(t);
	sum(k=0,sqrtint(n\3),
		t=n-k^2;
		sum(m=k,sqrtint(t>>1),
			issquare(t-m^2)
		)
	)
};
addhelp(ways3, "Number of ways that n can be represented as a sum of three squares. Sloane's A000164.");


msb(n:int)={
	my(k:int=0);
	n=floor(n);
	while(n>15,
		n>>=4;
		k+=4
	);
	if(n>3,
		n>>=2;
		k+=2
	);
	min(n,2)<<k
};
addhelp(msb, "msb(n): Most significant bit of n: returns the greatest power of 2 <= the number. Sloane's A053644.");


\\ ***************************************************************************************************
\\ *					Other number theory						*
\\ ***************************************************************************************************

Faulhaber(e:small,a='x)={
	substpol(sum(i=0,e,binomial(e+1,i)*bernfrac(i)*'x^(e+1-i))/(e+1) + 'x^e, 'x, a)
};
addhelp(Faulhaber, "Faulhaber(e,{a='x}): Returns the polynomial for the sum 1^e + 2^e + ... + x^e, evaluated at a.");


rp(b:small)={
	b=1<<(b-1);
	nextprime(b+random(b))
};
addhelp(rp, "rp(b): Returns a random b-bit prime.");


countPowerful(lim)={
	sum(k=1,(lim+0.5)^(1/3),
		if(issquarefree(k), sqrtint(lim\k^3))
	)
};
addhelp(countPowerful, "countPowerful(lim): Number of powerful numbers up to lim.");


countSquarefree(lim)={
	my(b=floor(sqrt(lim/2)));
	sum(k=1,b,moebius(k)*(lim\k^2))+
	sum(k=b+1,sqrt(lim),moebius(k))
};
addhelp(countSquarefree, "countSquarefree(lim): Counts the number of squarefree numbers up to lim.");


\\ ***************************************************************************************************
\\ *						Factoring						*
\\ ***************************************************************************************************

Mfactor(p:int,lim,start=2)={
	my(v=[],i,k);
	if (p < 1, error("p must be positive"));
	if (!isprime(p), error("p must be prime"));
	if (bitand(p, 3) != 3, print("  ###   warning: p must be a Mersenne exponent equal to 3 mod 4"));

	\\ Really, only k in [0, 5, 8, 9] mod 12 need to be checked (at last for Mersenne prime exponents?).
	\\ So this should use a quick check, then a full loop with step [10p, 6p, 2p, 6p].
	\\ Check for a factor of 3 also before starting loop.
	k = ceil((start-1)/(2*p));
	forstep(q=2*p*k+1,lim,2*p,
		if (bitand(q,7) != 1 && bitand(q,7) != 7, next);
		if (Mod(2,q)^(p%(q-1))!=1, next);
		
		v = concat(v, q);
		i = 2;
		while (Mod(2,q^i)^(p%((q-1)*q^(i-1)))==1,
			v = concat(v, q);
			i++
		)
	);
	v
};
addhelp(Mfactor, "Mfactor(p,lim,start=2): Returns factors of the Mersenne number 2^p-1 up to lim, starting at start, provided p is a prime = 3 mod 4. Same as bigfactor(2,p,1,lim,start) but faster because it checks only factors of the form 2kp+1 that are +/- 1 mod 8.");


bigfactor(a:int,b:int,c:int,lim,start=2)={
	my(v=[],i);
	if (b < 0,
	
		\\ These two should have their formats changed to match that given below.
		if (a == 1, return(factor(1 - c)));
		if (a == -1, return(factor(1-2*bitand(b, 1)-c)));
		
		\\ a^b not in Z
		error("not an integer argument in an arithmetic function")
	);
	forprime(p=start,min(a,lim),
		if (Mod(a,p)^if (gcd(a, p) > 1, b, b%(p-1))!=c, next);
		v = concat(v, p);
		i = 2;
		while (Mod(a,p^i)^if (gcd(a, p) > 1, b, b%((p-1)*p^(i-1)))==c,
			v = concat(v, p);
			i++
		)
	);
	forprime(p=precprime(max(start, min(a,lim)))+1,lim,
		if (Mod(a,p)^(b%(p-1))!=c, next);
		v = concat(v, p);
		i = 2;
		while (Mod(a,p^i)^(b%((p-1)*p^(i-1)))==c,
			v = concat(v, p);
			i++
		)
	);
	v
};
addhelp(bigfactor, "bigfactor(a,b,c,lim,start): Find small prime factors of a^b - c (up to lim). Optional parameter start gives a starting point below which primes are not checked.");


bigdiv(a:int,b:int,c:int,d:int)={
	if (b < 0,
		if (a == 1, return((1-c)%d == 0));
		if (a == -1, return((1-2*bitand(b, 1)-c)%d == 0));
		
		\\ a^b not in Z
		error("not an integer argument in an arithmetic function")
	);
	if (d <= 2,
		if (d == 0, error("division by zero"));
		if (d < 0, d = -d);
		if (d == 1, return(1));
		
		\\ d == 2
		if (b == 0,
			return(bitand(1-c,1) == 0)
		);
		return(bitand(a-c, 1) == 0)
	);
	if (gcd(a, d) > 1,
		Mod(a,d)^b==c
	,
		Mod(a,d)^(b%eulerphi(d))==c
	)
};
addhelp(bigdiv, "bigdiv(a,b,c,d): Does d divide a^b - c?  Same as (a^b-c)%d == 0, but faster for large b. Example: bigdiv(2,p,1,d) checks if d divides the p-th Mersenne number.");


\\ ***************************************************************************************************
\\ *				Real and complex functions					*
\\ ***************************************************************************************************

Bell(n:small)={
	my(pr:small=default(realprecision),B:real,Br:int,sz:small);
	sz = n * lg(n/exp(1)) \ lg(10) + 5;	\\ estimate of # of bits needed
	default(realprecision, sz);
	my(f:real=1.0,t:real,k:small=1);
	B=0.;
	while(1,
		f /= k;
		t = k^n * f;
		k++;
		if(t < 1e-9, break, B+=t)
	);
	B*=exp(-1);
	Br=round(B);
	if(abs(B-Br)>1e-6,error("not enough precision"));
	default(realprecision,pr);
	Br
};
addhelp(Bell, "Bell(n): Returns the n-th Bell or exponential number, Sloane's A000110.");


deBruijnXi(x)={
	my(left, right, m);
	if (x < 1, error ("deBruijnXi: Can't find a xi given x < 1."));
	if (x > 1, left = log(x), left = eps());
	right = 1.35 * log(x) + 1;		\\ Heuristic

	\\Bisection
	while (right - left > left * eps(),
		m = (left + right) / 2;
		if (exp(m) - 1 > x * m, right = m, left = m)
	);
	(left + right) / 2
};
addhelp(deBruijnXi, "deBruijnXi(x): Helper function for rhoest.  Finds a xi such that e^xi - 1 = x * xi.");


rhoest(x)={
	my(xi = deBruijnXi(x));
	\\exp(Euler) / sqrt(2 * Pi * x) * exp(1 - exp(xi) + intnum(s = eps(), xi, (exp(s) - 1) / s))
	exp(-eint1(-xi) - x * xi) / sqrt(2 * Pi * x) / xi
};
addhelp(rhoest, "de Bruijn's asymptotic approximation for rho(x), rewritten as in van de Lune and Wattel 1969.  Curiously, their paper shows values for this estimate that differ from those calculated by this function, often as soon as the second decimal place -- but as the difference is in the direction of the true value, I have not looked further into this.");


rhoTable = [1, 3.068528194e-1, 4.860838829e-2, 4.910925648e-3, 3.547247005e-4, 1.964969635e-5, 8.745669953e-7, 3.232069304e-8, 1.016248283e-9, 2.770171838e-11, 6.644809070e-13, 1.419713165e-14, 2.729189030e-16, 4.760639989e-18, 7.589908004e-20];
DickmanRho(x)={
	local(left, right, scale);
	if (x <= 2, return (1 - log(max(x, 1))));
	if (x <= 3, return(
		1 - (1 - log(x - 1))*log(x) + real(dilog(1 - x)) + Pi^2 / 12
	));

	\\ Asymptotic estimate (scaled for continuity)
	if (x > #rhoTable,
		scale = rhoTable[#rhoTable] / rhoest(#rhoTable);
		
		\\ Let the scale factor dwindle away, since the estimate is (presumably)
		\\ better in the long run than any scaled version of it.  The exponent
		\\ of 0.25 has been chosen to give the best results for 10 < x < 100
		\\ with a table size of 10.
		scale = (scale - 1) * (#rhoTable / x)^.25 + 1;
		return (precision(rhoest(x) * scale, 9))
	);

	\\ Scaling factors: the factor by which the true value of rho differs from
	\\ the estimates at the endpoints.
	left = rhoTable[floor(x)] / rhoest(floor(x));
	right = rhoTable[ceil(x)] / rhoest(ceil(x));
	
	\\ Linear interpolation on the scale factors.
	scale = left + (right - left) * (x - floor(x));
	
	\\ Return a result based on the scale factor and the asymptotic formula.
	precision(rhoest(x) * scale, 9)
};
addhelp(DickmanRho, "Estimates the value of the Dickman rho function. For x <= 3 the exact values are used, up to rounding; up to "#rhoTable" the value is interpolated using known values and rhoest; after "#rhoTable" rhoest is used, along with a correction factor based on the last value in rhoTable.");


contfracback(v,terms=-1)={
	my(x);
	if(terms==-1, terms=#v-1);
	x = v[terms+1];
	forstep(i=terms,1,-1,
		x = v[i] + x^-1;
	);
	x
};
addhelp(contfracback, "contfracback(v, terms): Given a continued fraction v, gives the real number back. If terms is given, use only that many terms.");


LambertW(x)={
	local(e,t,w,ep);
	if(x <= 0,
		if (x == 0, return (0));
		if (x == -exp(-1), return (-1));
		if (x < -exp(-1), error("LambertW: "x" is out of range, exiting."))
	);

	\\ Initial approximation for iteration
	if (x < 1,
		ep = sqrt(2*exp(1)*x+2);	\\ Using ep as a temporary variable
		w = ep * (1 - ep * (1/3 + 11/72*ep)) - 1
	,	w = log(x));
	if (x>3, w = w - log(w));

	t = 1;
	ep = eps() * (1 + abs(w));	\\ ep = epsilon

	while (abs(t) > ep, \\ Halley loop
		e = exp(w);
		t = w*e - x;
		t = t/(e*(w+1) - .5*(w+2)*t/(w+1));
		w = w - t
	);
	w
};
addhelp(LambertW,"Primary branch of Lambert's W function. Finds an L >= -1 such that L * exp(L) = x, where x >= -1/e.");


lg(x)={
	log(x)/log(2);
};
addhelp(lg, "Base-2 logarithm.");


\\ ***************************************************************************************************
\\ *					Convenience						*
\\ ***************************************************************************************************

vecsum(v)={
	sum(i=1,#v,v[i])
};
addhelp(vecsum, "vecsum(v): Sum of the elements of v.");


vecprod(v)={
	prod(i=1,#v,v[i])
};
addhelp(vecprod, "vecprod(v): Product of the elements of v.");


vecgcd(v)={
	local(g);
	g = 0;
	for (i=1,#v, g=gcd(g, v[i]));
	g
};
addhelp(vecgcd, "Vector gcd: returns the gcd of all elements in the vector.");


veclcm(v)={
	local(l);
	l = 1;
	for (i=1,#v, l=lcm(l, v[i]));
	l
};
addhelp(veclcm, "Vector lcm: returns the lcm of all elements in the vector.");


oddres(n)={
	n>>valuation(n,2)
};
addhelp(oddres, "oddres(n): Returns the greatest odd number dividing n.");


toC(n)={
	my(words,t);
	if (type(n) != "t_INT",
		error("can't handle non-integers");
	);
	if (n < 3,
		if (n == 2, print ("gen_2"));
		if (n == 1, print ("gen_1"));
		if (n == 0, print ("gen_0"));
		if (n == -1, print ("gen_m1"));
		if (n >= -1, return);
		error("can't handle negatives yet")
	);
	words=ceil(log(n)/log(1<<32));
	if (n>>(words*32) > 0, words++);
	if (n>>((words-1)*32) == 0, words--);
	if (words == 1,
		print("stoi("n")");
		return
	);
	if (words == 2,
		print1("uutoi(");
		print1(n>>32);
		print(", "bitand(n, 2^32-1)")");
		return
	);

	\\ Large numbers
	print1("mkintn("words);
	forstep(i=words-1,0,-1,
		t=n>>(i*32);
		print1(", "t);
		n -= t<<(i*32)
	);
	print(")");
};
addhelp(toC, "toC(n): Format n for use with the Pari library (e.g., with gp2c programs).");


digits(x)={
	my(s=sizedigit(x)-1);
	if(x<10^s,s,s+1)
};
addhelp(digits, "digits(n): Number of decimal digits in n. Sloane's A055642.");


\\ 2 ^ -(decimal_precision * log(10)/log (2) / 32) * 32 - 1)
\\ precision(2. >> (32 * ceil(precision(1.) * 0.103810252965230073370947482171543443)), 9)
eps()={
	precision(2. >> (32 * ceil(default(realprecision) * 38539962 / 371253907)), 9)
};
addhelp(eps,"Returns machine epsilon for the current precision.");


\\ ***************************************************************************************************
\\ *						I/O							*
\\ ***************************************************************************************************

bfile(name, v:vec, offset=1)={
	my(cur=offset-1,Anum);
	if (type(name) == "t_INT",
		name = Vec(Str(name));
		while(#name < 6, name = concat(["0"], name));
		Anum = concat(name);
		name = Str("b"Anum".txt");
	,
		if (type(name) != "t_STR", error("name must be an integer (A-number) or filename (\"b000040.txt\")."));
		Anum = concat(vecextract(Vec(name), 126))
	);
	for(i=1,#v,
		if (digits(v[i]) > 1000,
			print("Next term has "digits(v[i])" digits; exiting.");
			break
		);
		write(name, cur++" "v[i]);
	);
	print("%H A"Anum" Author, <a href=\"b"Anum".txt\">Table of n, a(n) for n = "offset".."cur"</a>");
};
addhelp(bfile, "bfile(name, v, offset=1): Creates a b-file with the values of v for A-number name (given as a number or a filename).");


fnice(n)={
	my(f,s="",s1);
	if(n < 0,
		s = "-";
		n = -n
	);
	if (n < 2, return(concat(s, n)));
	
	f = factor(n);
	s = Str(s, f[1,1]);
	if (f[1, 2] != 1, s=Str(s, "^", f[1,2]));
	for(i=2,#f[,1],
		s1 = Str(" * ", f[i, 1]);
		if (f[i, 2] != 1, s1 = Str(s1, "^", f[i, 2]));
		s = Str(s, s1)
	);
	s
};
addhelp(fnice, "fnice(n): Returns a string with a 'nice' factorization of n.");


nice(o)={
	my(s,t,t1);
	if (type(o) == "t_POL",	\\ Can handle only single-variable polynomials.
		t = content(o);
		if (denominator(t) == 1,
			my(v=variable(o));
			t = poldegree(o);
			t1 = polcoeff(o, t, v);
			s = if(t1 < 0,
				Str("-", monomialnice(-t1, t, v))
			,
				monomialnice(t1, t, v)
			);
			o -= t1 * v^t;
			t = poldegree(o,v);
			while (t > 0,
				t1 = polcoeff(o, t, v);
				s = Str(s, if(t1 > 0, " + ", " - "), monomialnice(abs(t1), t, v));
				o -= t1 * v^t;
				t = poldegree(o,v);
			);
			o = simplify(o);
			if (type(o) != "t_INT", error("nice: cannot display multivariate polynomials!"));
			o = polcoeff(o, 0);
			if (o > 0, s = Str(s, " + ", o));
			if (o < 0, s = Str(s, " - ", -o));
			return(s)
		);
		o = nice(denominator(t) * o);
		if (numerator(t) == 1, return(Str("(", o, ")/", 1/t)));
		return(Str(numerator(t), "(", o, ")/", denominator(t)))
	);
	Str(o)
};
addhelp(nice, "nice(o): Reformats the object o 'nicely' when possible. Currently chokes on multivariable polynomials.");


\\ Degree assumed to be positive
monomialnice(coeff:int, degree:int, v='x)={
	if (coeff == 1,
		return(if(degree == 1, Str(v), Str(v, "^", degree)))
	);
	if(degree == 1, Str(coeff, v), Str(coeff, v, "^", degree))
};


\\ ***************************************************************************************************
\\ *						Set stuff						*
\\ ***************************************************************************************************

sumset(a,b)={
	local(c);
	c=vector(#a * #b);
	for(i=1,#a,
		for (j=1,#b,
			c[(i - 1)*#b + j] = a[i] + b[j];
		)
	);
	vecsort(c,,8)
};
addhelp(sumset,"sumset(A, B) is the set of all numbers of the form a+b, a in A, b in B.");


diffset(a,b)={
	local(c);
	c=vector(#a * #b);
	for(i=1,#a,
		for (j=1,#b,
			c[(i - 1)*#b + j] = a[i] - b[j];
		)
	);
	vecsort(eval(Set(c)))
};
addhelp(diffset,"diffset(A, B) is the set of all numbers of the form a-b, a in A, b in B.");


\\ ***************************************************************************************************
\\ *					Statistics						*
\\ ***************************************************************************************************


normd(a,b)={
	\\ Infinities
	if(a==[-1],
		return(if (b==[+1], 1, -erfc(b/sqrt(2))/2))
	);
	if(b==[1], return(erfc(a/sqrt(2))/2));
	
	(erfc(a/sqrt(2))-erfc(b/sqrt(2)))/2
};
addhelp(normd, "normd(a,b): Amount of the normal distribution between a and b standard deviations. Plus/minus infinity coded as [+1]/[-1].");


rnormal()={
	my(pr=32*ceil(default(realprecision)*log(10)/log(4294967296)),u1=random(2^pr)*1.>>pr,u2=random(2^pr)*1.>>pr);
	sqrt(-2*log(u1))*cos(2*Pi*u1)
	\\ Could easily be extended with a second normal at very little cost.
};
addhelp(rnormal, "rnormal(): Returns a random normal variable with mean 0 and standard deviation 1 at the current precision.");


default(timer, timervalue);
