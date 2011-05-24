timervalue = default(timer);
default(timer, 0);


\\ ***************************************************************************************************
\\ *					Working space						*
\\ ***************************************************************************************************


cons(x:real)={
	my(t=precision(x),s=log(x)\log(10),v=vector(t));
	x /= 10^s;
	for(i=1,t,
		v[i]=floor(x);
		x = (x - floor(x)) * 10
	);
	v
};
addhelp(cons, "cons(x): Converts a constant into a sequence of decimal digits.");


\\ Output format (work-in-progress):
\\ [x1, y1, x2, y2, ..., xk, yk, c]: represents density c * x^x1 / log(x)^y1 * x^x2 / log(x)^y2 * ... * x^xk / log(x)^yk
\\ [x1, y1, x2, y2, ..., xk, yk]: same, but with unknown constant.
\\ Examples:
\\ [k]			Exactly k primes
\\ []			O(1) primes
\\ [1, 1, 1]	(1 + o(1)) x/log(x) primes (the maximum)
\\ [1/2, 1]		Theta(sqrt(x)/log(x)) primes -- not sure about which symbol to use, maybe O is better
primesin(P:pol)={
	my(t=type(P));
	if (t=="t_INT",
		return(if (isprime(P),
			print("Contains 1 prime.");
			[1]
		,
			print("Contains 0 primes.");
			[0]
		));
	);
	if (t=="t_FRAC" || t=="t_REAL",
		print("Contains 0 primes.");
		return([0])
	);
	if (t!="t_POL", error("bad type"));
	
	if (poldegree(P) == 1,
		my(a=polcoeff(P,1),b=polcoeff(P,0),g,t);
		if (type(a) == "t_INT" && type(b) == "t_INT",
			g=gcd(a,b);
			if (g == 1,
				t=eulerphi(a);
				print("Infinitely many primes (Dirichlet), density n/log n * ",1/t);
				return([1,1,1/t])
			);
			return(if (isprime(g) && (1 - b/g) % (a/g) == 0,
				print("Contains 1 prime.");
				[1]
			,
				print("Contains 0 primes.");
				[0]
			));
		);
	);
	if (poldegree(P) > 2, error("Cannot be determined"));
};
addhelp(primesin, "primesin(P): Attempts to give the asymptotic number of primes in the polynomial P, either unconditionally or under standard conjectures. The output format is a vector, which represents:\n[x1, y1, x2, y2, ..., xk, yk, c]: density c * x^x1 / log(x)^y1 * x^x2 / log(x)^y2 * ... * x^xk / log(x)^yk\n[x1, y1, x2, y2, ..., xk, yk]: same, but with unknown leading term.");


ConjectureFEmp(a:int,b:int,c:int,lim=1e6)={
	my(s=0,n=0,k);
	for(n=1,lim,
		k=(a*n+b)*n+c;
		if(k>lim,break);
		s+=isprime(k)
	);
	2*s/Li(sqrt(lim))
};
addhelp(ConjectureFEmp, "ConjectureFEmp(a,b,c,{lim}): Helper function for ConjectureF, given integers a, b, and c as coefficients of ax^2 + bx + c.");


ConjectureF(a:int,b:int,c:int,lim=1e6)={
	my(D=b^2-4*a*c,g=gcd(a,b),f);
print("Delta = "D);
	if (a <= 0 | issquare(D) | gcd(g,c) > 1, return(0));
	if ((a+b)%2 ==0 & c%2 == 0, return(0));

	f=factor(g)[,1];
	if(a+b%2, 1, 2)/sqrt(a) * prod(i=1,#f,f[i]/(f[i]-1)) * ConjectureF_C(D,lim)
};
addhelp(ConjectureF, "ConjectureF(a,b,c): Returns a number C such that there are ~ C*sqrt(x)/log(x) primes of the form an^2 + bn + c up to x, or 0 if there are o(sqrt(x)/log(x)) or Omega(sqrt(x)) such primes. Relies on the Hardy-Littlewood Conjecture F and possible the ERH.");


\\ prodinf(i=1,p=prime(i);1-kronecker(D,p)/(p-1))
ConjectureF_C(D,lim=1e6)={
	my(cf, pr=1., l=log(lim));
	if (D < -4,
		print("Fung & Williams 1990");
		my(f=factor(oddres(-D))[,1], c = if(delta%8 == 1, 5/2, if (delta%8 == 5, 1/2, 15/16)));
		\\ Corrects for the 1 - 1/p^3 (since kronecker(D,q) == 1 about half the time)
		forprime(q=3,lim,if(kronecker(D,q)==1,pr*=1-2/q/(q-1)^2));
		cf = -2*eint1(2*l) - (1-2*l)/2/lim^2/l^2;
		pr *= 1 - cf;
print("Correction: ", cf);
		return(c * Pi^3 * sqrt(-D) / 90 / qfbclassno(D) / Lquad(D, 2) * prod(i=1,#f,1-1/f[i]^4) * pr)
	);
	if (D > 0,
		print("Jacobson 1995");
		\\ Possibly needed:
		\\ QFBclassno(D) = qfbclassno(D) * if (D < 0 || norm(quadunit(D)) < 0, 1, 2)
		my(f=factor(oddres(D))[,1], c = if(delta%8 == 1, 5/2, if (delta%8 == 5, 1/2, 15/16)));
		forprime(q=3,lim,if(kronecker(D,q)==1,pr*=1-2/q/(q-1)^2));
		\\ Corrects for the 1 - 1/p^3 (since kronecker(D,q) == 1 about half the time)
		cf = -2*eint1(2*l) - (1-2*l)/2/lim^2/l^2;
		pr *= 1 - cf;
print("Correction: ", cf);
		return(c * Pi^4 * sqrt(D) / 180 / quadregulator(D) / qfbclassno(D) / Lquad(D, 2) * prod(i=1,#f,1-1/f[i]^4) * pr)
	);
	if(1,
		print("KConrad");
		my(t);
		\\ Corrects for the 1 - 2/p^3
		forprime(p=3,lim,t=kronecker(D,p);pr*=(1-t/(p-1))/(1-t/p));
		cf = -4*eint1(2*l) - (1-2*l)/lim^2/l^2;
		pr *= 1 - cf;
print("Correction: ", cf);
		return(pr/Lquad(D, 1));
	);
	print("Naive method");
	\\ cf: the correction factor that estimates suminf(p=lim,-1/p^2) over the primes
	\\ Corrects for the average of (1 - 1/p) and (1 + 1/p); not very useful.
	cf=1/lim/l-eint1(l);
	forprime(p=3,lim,pr*=1-kronecker(D,p)/(p-1));
print("Correction: ", cf);
	pr * (1 - cf)
};
addhelp(ConjectureF_C, "ConjectureF_C(D,{lim}): Helper function for ConjectureF, finding the infinite product in the Hardy-Littlewood Conjecture F.");


hurwitz(s,x)={
	my(a,res,tes,in,sig,t,m,pr,lim,nn,in2,s1,s2);
	sig=real(x);
	if (sig>1.5,
		m=floor(sig-0.5);
		return (hurwitz(s,x-m)-sum(i=1,m,(x-i)^(-s)))
	);
	if (sig<=0,
		m=ceil(-sig+0.5);
		return (hurwitz(s,x+m)+sum(i=0,m-1,(x+i)^(-s)))
	);
	pr=precision(1.);
	sig=real(s);
	t=imag(s);
	default(realprecision,9);
	res=s-1.;
	res=if(abs(res)<0.1,-1,log(res));
	lim=(pr*log(10)-real((s-.5)*res)+(1.*sig)*log(2.*Pi))/2;
	lim=max(2,ceil(max(lim,abs(s*1.)/2)));
	nn=ceil(sqrt((lim+sig/2-.25)^2+(t*1.)^2/4)/Pi);
	default(realprecision,pr+5);
	a=(x+nn+0.)^(-s);
	res=sum(n=0,nn-1,(x+n)^(-s),a/2);
	in=x+nn;
	in2=1./(in*in);
	s1=2*s-1;
	s2=s*(s-1);
	tes=bernreal(2*lim);
	forstep (k=2*lim-2,2,-2,
		tes=bernreal(k)+in2*(k*k+s1*k+s2)*tes/((k+1)*(k+2))
	);
	tes=in*(1+in2*s2*tes/2);
	res+=tes*a/(s-1);
	res=precision(res,pr);
	default(realprecision,pr);
	res
};


\\ (VIII.2) Complex L function, vector form. Chivec is a vector of complex
\\ values, assumed to be the values from 1 to m of a periodic function. In
\\ addition, with zero sum, such as a nontrivial character, if s=1.
\\ simple implementation, for small m.
Lsimp(chivec,s)={
	my(m=#chivec);
	if (s==1,
		-sum(r=1,m,chivec[r]*psi(r/m))/m
	,
		sum(r=1,m,chivec[r]*hurwitz(s,r/m))/m^s
	)
};


\\ (VIII.3) Complex L function of quadratic character, simple implementation
\\ for small |D|.
Lquad(D, s)={
	my(v=vector(abs(D),i,kronecker(D,i)));
	Lsimp(v,s)
};









forpseudo(ff)={
	if(default(parisize) < 2e8,
		print("Not enough memory!  Resizing, please rerun command.");
		allocatemem(200<<20)	\\ Need ~200 MB of stack space!
	);
	for(i=0,89,
		read(concat("psp/psp-chunk", i));
		trap(user,
			print("User error, ending loop...");
			pspChunk=0;
			return()
		,
			for(j=1,#pspChunk,
				ff(pspChunk[j])
			)
		);
		pspChunk=0
	)
};
addhelp(forpseudo, "forpseudo(ff): Runs the command (closure) ff on each 2-pseudoprime up to 2^64.");


\\ Characteristic function of A
\\ output sequence b_n = 1 if n if exist k : a_k == n
\\ otherwise b_n = 0
\\ 2nd parameter is length of output sequence
trv_char(A, n=-1)={
	local(B,v);
	if(n<0, n = A[#A]);
	B = vector(n);
	for(i=1, #A,
		v=A[i];
		if(type(v)=="t_INT" && v>=0 && v<n, B[v+1]=1);
	);
	B
};
addhelp(trv_char, "trv_char(A, {n}): Characteristic function of vector A. Gives the first n terms, or all terms that can be calculated if omitted.");


\\ Euler transform
trv_euler(A)={
	my(B=vector(#A-1,n,1/n),C);
	B = dirmul(vecextract(A,"2.."),B);
	C = exp(Ser(concat(0,B)));
	for(i=1, #A, A[i] = polcoeff(C,i-1));
	A
};
addhelp(trv_euler, "trv_euler(A): Euler transform of vector A.");


ways(v, n)={
	my(w=vector(n));
	for(i=1,n,
		for(j=1,#v,
			if (i > v[j],
				w[i] += w[i - v[j]];
			,
				if (i == v[j], w[i]++);
				break
			)
		)
	);
	w
};
addhelp(ways, "ways(v, n): Given an increasing vector v of positive integers, return an n-element vector w so that w[i] is the number of ways (with ordering) to sum elements of v to get i.");


/*

time(ff,lim,sz=0)={
	my(tEmpty,tNull,tFunc,timeRange,appx,gg);
	sz=floor(sz)+1;
	(gg=n->0);
	gettime();
	for(n=sz,lim+sz,gg(n));
	tEmpty=gettime();
	for(n=sz,lim+sz,ff(n));
	tFunc=gettime();
	for(n=sz,lim+sz,);
	tNull=gettime();
	timeRange=[tFunc - tEmpty, tFunc - tNull] / lim;
	appx = (timeRange[1] + timeRange[2])/2;
	if (appx < .001, print("About "round(appx*1000000)" ns"),
	if (appx < 1, print("About "round(appx*1000)" µs"),
	if (appx < 1000, print("About "round(appx)" ms"),
	print("About "round(appx/1000)" s")
	)));
	precision(timeRange/1000,9)	\\ Convert to seconds
};


find(m,mn,startAt)={
	my(w,startAt);
	if(startAt==0,
		startAt = m
	,
		startAt = (startAt \ (m + m)) * (m + m) + m
	);
	forstep(n=startAt,9e999,m+m,
		w=ways2(n);
		if(w>mn,
			print(w" "n" = "fnice(n));
			return([w,n])
		)
	)
};

intersect(v1,v2)={
	my(i=1,j=1,v=[]);
	while(i<=#v1&&j<=#v2,
		if(v1[i]==v2[j],
			v=concat(v,v1[i]);
			i++;
			j++
		,
			if(v1[i]>v2[j],j++,i++)
		)
	);
	v
};
addhelp(intersect, "intersect(v1,v2): Intersection of the sorted sets v1 and v2.");


doesintersect(v1,v2)={
	my(i=1,j=1,v=[]);
	while(i<=#v1&&j<=#v2,
		if(v1[i]==v2[j],
			return(1)
		,
			if(v1[i]>v2[j],j++,i++)
		)
	);
	0
};


ord(a,p)={
	my(f=factor(p-1),T=1,q,b);
	for(i=1,#f[,1],
		q=f[i,1];
		b=Mod(a,p)^((p-1)/q^f[i,2]);
		while(b!=1,
			T*=q;
			b=b^q
		)
	);
	T
};


isperf(n)={
	my(f=factor(n),spf,nd,sig);
\\print([spf,nd,sig]);
	sig=prod(i=1,#f[,1],
		if(f[i,2]==1,f[i,1]+1,(f[i,1]^(f[i,2]+1)-1)/(f[i,1]^f[i,2]-1))
	);
	if(issquare(sig,&sig),
		nd=prod(i=1,#f[,1],f[i,2]+1);
		spf=sum(i=1,#f[,1],f[i,1]);
		sig==nd*spf
	,
		0
	)
};


English(n:int)={
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


vdiff(v1,v2)={
	my(i=1,j=1);
	while(i<=#v1 && j<=#v2,
		if (v1[i] == v2[j],
			i++;
			j++
		,
			if (v1[i] < v2[j],
				print("Second vector missing "v1[i]);
				i++
			,
				print("First vector missing "v2[j]);
				j++
			)
		)
	);
};


fibmod(n:int, m:int)={
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

\\ Version 2... worse, I think?
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
}; */


findrec(v:vec, verbose:bool=1)={
	my(c,d = (#v - 1) >> 1, firstNonzero = 0);
	if (#v < 3,
		warning("findrec: Not enough information to determine if the sequence is a recurrence relation: matrix is underdetermined. Give more terms and try again.");
		return
	);
	forstep(i=d,1,-1,
		if(v[i] != 0, firstNonzero = i)
	);
	if(firstNonzero == 0,
		if(verbose,
			print1("Recurrence relation is a(n) = 0.");
		);
		return([0]~);
	);
	for (i=firstNonzero,d,
		c = findrecd(v,i,verbose);
		if(c, return(c))
	);
	if(verbose,print("Cannot be described by a homogeneous linear recurrence relation with "d" or fewer coefficients."));
	0
};
addhelp(findrec, "findrec(v, {verbose=1}): Tries to find a homogeneous linear recurrence relation with constant coefficients to fit the vector v.");


findrecd(v:vec, d:int, verbose:bool=1)={
	my(M,c);
	if (#v < 2*d,
		warning("findrec: Not enough information to determine if the sequence is a "d"-coefficient recurrence relation; matrix is underdetermined. Give more terms or reduce d and try again.");
		return
	);
	M = matrix(d,d,x,y,v[d+x-y]);
	\\print(M" * c = "vector(d,i,v[d+i])~);
	if(matdet(M) == 0, return(0));	\\ Non-invertible matrix, no solutions for this d
	c = matsolve(M,vector(d,i,v[d+i])~);
	for(n=2*d+1,#v,
		if(v[n] != sum(i=1,d,v[n-i] * c[i]),return(0))
	);
	if(verbose,
		my(init=1,s);
		print1("Recurrence relation is a(n) = ");
		for(i=1,#c,
			if(c[i] == 0, next);
			if(init,
				s = initial(c[i], Str("a(n-", i, ")"));
				init = 0
			,
				s = Str(s, medial(c[i], Str("a(n-", i, ")")))
			)
		);
		print(s".");
		if((vecmax(c) == 1 && vecmin(c) == 0 && vecsum(c) == 1) || c == [1]~,
			print("Sequence has period "#c".");		
		,
			my(g=0);
			for(i=1,#c,
				if(c[i] != 0, g = gcd(g, i))
			);
			if (g > 1,
				my(gvec = vector(#c/g, i, c[i*g]),s,init=1);
				for(i=1,#gvec,
					if(gvec[i] == 0, next);
					if(init,
						s = initial(gvec[i], Str("a(n - ", i, ")"));
						init = 0
					,
						s = Str(s, medial(gvec[i], Str("a(n - ", i, ")")))
					)
				);
				print("Can be though of as "g" interlocking sequences, each of the form a(n) = "s".")
			)
		);
		print1("<a href=\"/Sindx_Rea.html#recLCC\">Index to sequences with linear recurrences with constant coefficients</a>, signature ("c[1]);
		for(i=2,#c,print1(","c[i]));
		print(").");
		print((#v-2*d)" d.f.")
	);
	c
};
addhelp(findrecd, "findrecd(v, d, {verbose=1}): Helper function for findrec. Tries to find a d-coefficient homogeneous linear recurrence relation with constant coefficients to fit the vector v.");


dotproduct(a:vec, b:vec)={
	if(#a != #b,
		error("dotproduct: vectors are not the same size")
	);
	sum(i=1,#a,a[i]*b[i])
};
addhelp(dotproduct, "dotproduct(a, b): Returns the dot product of vectors a and b.");


/*
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
*/

\\ Varies as n*log(log(n))/log(n)
countsemi(n:int)={
	my(s=0,B=sqrtint(n));
	forprime(p=2,B,
		s+=primepi(n\p)
	);
	B=primepi(B);
	s-B*(B-1)/2
};
addhelp(countsemi, "countsemi(n): Number of semiprimes up to n. Sloane's A072000.");


/*
\\ TODO: This generates too few even semiprimes.  For example, 19% of the
\\ 20-bit semiprimes are even, but only 16% of the numbers generated by rsp(20)
\\ are even:
\\ shouldbe(b)=(primepi(2^(b-1)-1)-primepi(2^(b-2)-1))/(countsemi(2^b-1)-countsemi(2^(b-1)-1))
\\ is(b)=my(M=.261497212847642783755);(log(log(3))+M)/(log(b*log(2)/2)+M)
\\
\\ This is because the first prime 2 takes essentially all the error in the
\\ approximation log log x + M on itself.
\\
\\ Other issues include: (1) primes with big prime gaps following take too much
\\ weight, (2) primes around the square root may be over- or under-represented,
\\ and (3) semiprimes are larger than 2^b on occasion. (Also under 2^(b-1), but
\\ this is much less common and not really a concern.) Problem #1 can be
\\ handled with the problem with 2 by handling small primes separately -- maybe
\\ through 127. This might require handling the semiprimes under 15 bits as
\\ special cases.
rsp(b:small,flag:small=0)={
	if (bitand(flag,2),
		if(b<3, error("No semiprimes with fewer than 3 bits"));
		my(p,q,lower,upper,t);
		if(b%2,
			lower=sqrtint(1<<(b-2))+1;
			upper=1<<((b+1)/2);
		,
			lower=1<<((b-2)/2);
			upper=sqrtint(1<<(b+1));
		);
		while(p = precprime(random(upper - lower) + lower);
			q = nextprime(random(upper - lower) + lower);
			if (p > q,
				t = p;
				p = q;
				q = t;
			);
			p * 2 < q || #binary(p * q) != b
		,);
		return(if(bitand(flag,1),
			if(p==q,
				matrix(1,2,i,j,if(j==1,p,2))
			,
				[p,1;q,1]
			)
		,
			p*q
		))
	);
	if(b<26,	\\ Cases below 15 must be handled here so that primes up to
				\\ 127 can be handled separately. Cases below 26 are just
				\\ faster this way, depending on precision.
		my(n);
		return(if(b<8,
			if(b<3, error("No semiprimes with fewer than 3 bits"));
			if(b==3,n=[4,6]);
			if(b==4,n=[9,10,14,15]);
			if(b==5,n=[21,22,25,26]);
			if(b==6,n=[33,34,35,38,39,46,49,51,55,57,58,62]);
			if(b==7,n=[65,69,74,77,82,85,86,87,91,93,94,95,106,111,115,118,119,121,122,123]);
			n=n[random(#n)+1];
			if(flag,factor(n),n)
		,
			b=2^(b-1);
			while(bigomega(n=random(b)+b)!=2,);	\\ Brute force
			n
		))
	);

	\\ Choose the small prime as 2 with probability 1/2, 3 with probability 1/3, etc.
	my(M=.261497212847642783755,weight=log(b*log(2)/2)+M,r=random(round(weight << 64))*1.>>64,p=precprime(exp(exp(r-M))),q);
	q=2^(b-1)\p;
	q=nextprime(q+random(q));
	if(p > q,
		b = p;
		p = q;
		q = b;
	);
	if(bitand(flag,1),
		if(p==q,
			matrix(1,2,i,j,if(j==1,p,2))
		,
			[p,1;q,1]
		)
	,
		p*q
	)
};
addhelp(rsp, "rsp(b, {flag=0}): Generates a random b-bit semiprime. The bits of flag mean: 1: returns the factorization of the number instead of the number, 2: return worst-case semiprimes rather than average semiprimes.");


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


bestappr2(x, n)={
	my(b=bestappr(x, n), cf=contfrac(x, #contfrac(b)+1), lower=ceil(cf[#cf]/2), upper=cf[#cf]);
	cf[#cf] = lower;
	contfracback(cf)
};
addhelp(bestappr2, "bestappr2(x, n): Work-in-progress.  Should determine the best rational approximation of x with denominator up to n.");


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
*/

primezeta(s)={
	suminf(k=1,moebius(k)/k*log(abs(zeta(k*s))))
};
addhelp(primezeta, "primeZeta(s): Returns the prime zeta function of s, the sum of p^-s over all primes p.");


/*
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
*/


initial(n:int, s:genstr)={
	if (n == 0, return(""));
	if (n == -1, return(Str("-", s)));
	if (n == 1, return(s));
	Str(n, s)
};


medial(n:int, s:genstr)={
	if (n == 0, return(""));
	if (n == -1, return(Str(" - ", s)));
	if (n == 1, return(Str(" + ", s)));
	Str(if (n < 0, " - ", " + "), abs(n), s)
};

/*
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
*/

binomialIntervalExact(n:int,k:int,conf=.05)={
	my(a,v,flip,result);
	if (k==0, return([0, solve(p=0,1,(1-p)^n-conf)]));
	if (k==n, return([solve(p=0,1,p^n-conf), 1]));

	conf = conf/2;
	if (k + k > n,
		flip = 1;
		k = n - k;
	);
	a = if (k*8 < n - 10, .5, 1);
	
	v = vector(k+1,i,binomial(n,i-1));
	result = [solve(p=0,.5,1-conf-sum(i=0,k-1,v[i+1]*p^i*(1-p)^(n-i))),
	solve(p=0,a,sum(i=0,k,v[i+1]*p^i*(1-p)^(n-i))-conf)];
	
	if(flip,[1 - result[2], 1 - result[1]], result)
};
addhelp(binomialIntervalExact, "binomialIntervalExact(n,k,conf=.05): Gives a confidence interval for the probability of an event which happens k times in n independent trials with fixed probability. Uses a direct inversion of the exact combinatorial probability; seems to be the Clopper-Pearson interval.");


binomialInterval(n,k,conf=.05)={
	my(kappa=ierfc(conf/2)*sqrt(2), kest=k+kappa^2/2, nest=n+kappa^2, pest=kest/nest, radius=kappa*sqrt(pest*(1-pest)/nest));
	[max(0, pest - radius), min(1, pest + radius)]
};
addhelp(binomialInterval, "binomialInterval(x,k,conf=.05): Computes the Agresti-Coull interval. Conservative: tends to give slightly wider ranges than binomialIntervalExact.");


\\ ierfc(1 - x) = ierf(x)
ierfc(x)={
	if (x <= 1.999999999998 && x >= 2e-12, return(solve(t=-5,5,erfc(t)-x)));
	if (x <= 0, error("ierfc: Out of range/overflow"));
	if (x >= 2, error("ierfc: Out of range/overflow"));
	solve(t=-999,999,erfc(t)-x)
};
addhelp(ierfc, "Inverse complementary error function.");


\\ ***************************************************************************************************
\\ *					Verbose monstrosities					*
\\ ***************************************************************************************************

\\ These first three aren't 'verbose monstrosities', but are needed by them.

\\ suminf(n=1,moebius(n)/n*Li(x^(1/n)))
R(x)={
	local(l2,u);
	l2 = li(2);
	u = -log(x);
	-suminf(n=1,
		moebius(n)/n * (eint1(u/n) + l2)
	)
};
addhelp(R, "Riemann's R function, an estimate for primepi.");

li(x)=-eint1(-log(x));
addhelp(li, "Logarithmic integral.");

Li(x)=eint1(-2) - eint1(-log(x));
addhelp(Li, "Offset logarithmic integral, an estimate for primepi. Crandall and Pomerance call this li_0.");

piBounds(x,verbose:bool=0)={
	local(lower,upper,lx,t,lowerRH,upperRH,roundFlag=default(realprecision)>sizedigit(x));
	
	if(roundFlag, x = floor(x));
	if (x < default(primelimit),
		lower=primepi(x);
		print("primepi("x") = "lower);
		return([lower,lower]);
	);

	lx = log(x);
	lower = x/lx * (1 + 1/lx + 1.8/lx^2);	\\ Dusart, x >= 32299
	upper = x/(lx - 1.1);					\\ Dusart, x >= 60184
	
	if (x >= 355991,
		t = x/lx * (1 + 1/lx + 2.51/lx^2);	\\ Dusart, x >= 355991
		if (upper > t, upper = t)
	);

	if (x >= 1.332e10,
		t = x/lx * (1 + 1.0992/lx);			\\ Dusart, x >= 1.332e10
		if (upper > t, upper = t)
	);
	
	if (roundFlag,
		lower = ceil(lower);
		upper = floor(upper);
	);
	
	t = sqrt(x)/8/Pi*log(x);
	lowerRH = li(x)-t;					\\ Schoenfield, x >= 2657
	upperRH = li(x)+t;					\\ Schoenfield, x >= 2657
	if (roundFlag,
		lowerRH = ceil(lowerRH);
		upperRH = floor(upperRH);
	);

	print("For primepi("x"):");
	print(lower" (lower bound)");
	if (lowerRH > lower,
		print(lowerRH" (lower bound under the RH)");
	);
	if(roundFlag,
		print(round(R(x))" (Riemann R, approximate)");
		print(round(Li(x))" (logarithmic integral, apx)");
	,
		print(R(x)" (Riemann R, approximate)");
		print(Li(x)" (logarithmic integral, apx)");
	);
	if (upperRH < upper,
		print(upperRH" (upper bound under the RH)");
	);
	print(upper" (upper bound)");
	
	if (verbose,
		print("\nPierre Dusart, 'Autour de la fonction qui compte le nombre de nombres");
		print("premiers', doctoral thesis for l'Université de Limoges (1998).");
		
		if (lowerRH > lower || upperRH < upper,
			print("Lowell Schoenfeld, 'Sharper Bounds for the Chebyshev Functions theta(x)");
			print("and psi(x). II'. Mathematics of Computation, Vol 30, No 134 (Apr 1976),");
			print("pp. 337-360.");
		);
		print();
	);
	[lower, upper]
};
addhelp(piBounds, "piBounds(x, verbose=0): Bounds on primepi(x). Set verbose=1 to get a list of sources for the results.");


pBounds(n, verbose:bool=0)={
	my(lower,upper,appx,l,ll);
	if (n < 6548,
		return(if (n < 1,
			print("There are no negative primes.");
			[0,0]
		,
			n = floor(n);
			lower=prime(n);
			print("p_"n" = "lower" (exactly)");
			[lower,lower]
		))
	);
	n = floor(n);
	l = log(n);
	ll = log(l);
	
	lower = n * (l + ll - 1);						\\ Dusart, n >= 2
	if (n > 13196,
		lower = n * (l + ll - 1 + ll/l - 2.25/l)	\\ Dusart, n >= 2
	);
	appx  = n * (l + ll - 1 + ll/l - 2/l - ll^2/2/l^2 + 3*ll/l^2 + 11/2/l^2);	\\ + O(ll^3/l^3)
	upper = n * (l + ll);							\\ ?, n >= 6
	if (n >= 27076,
		upper = n * (l + ll - 1 + ll/l - 1.8/l)		\\ Dusart, n >= 27076
	);
	if (n >= 39017,
		upper = min(upper, n * (l + ll - .9484))	\\ Dusart, n >= 39017
	);
	
	lower = ceil(lower);
	upper = floor(upper);
	
	print(lower" (lower bound)");
	if (lower<appx && appx<upper,
		print(appx" (approximate)")
	);
	print(upper" (upper bound)");
	
	if (verbose,
		print("\nPierre Dusart, 'Autour de la fonction qui compte le nombre de nombres");
		print("premiers', doctoral thesis for l'Université de Limoges (1998).");
		if (lower<appx && appx<upper,
			print("Ernest Cesàro (1894). \"Sur une formule empirique de M. Pervouchine\". Comptes");
			print("rendus hebdomadaires des séances de l'Académie des sciences 119, pp. 848-849.")
		);
	);
	
	\\t = sqrt(x)/8/Pi*log(x);
	\\lowerRH = ceil(li(x)-t);					\\ Schoenfield, x >= 2657
	\\upperRH = floor(li(x)+t);					\\ Schoenfield, x >= 2657
	[lower, upper]
};
addhelp(pBounds, "pBounds(n, verbose=0): Estimates the nth prime. Set verbose=1 to get a list of sources for the results.");

\\ ***************************************************************************************************
\\ *					Works-in-progress						*
\\ ***************************************************************************************************

/*
\\ Really bad here: >= 32 bits per number in sieve, about 100 times the size of a basic wheel mod 6
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

*/
Psi(x, B)={
	local(log3x, log2, tmp, u, result);
	if (x < 1, return (0));
	if (B < 7,
		if (B < 2, return (B >= 1));
		log2 = log(2);
		if (B < 3, return (floor(log(x) / log2) + 1));
		log3x = log(x) / log(3) + 1;	\\ The + 1 accounts for the pure powers of 2.
		tmp = log2 / log(3);
		if (B < 5,
			\\ Snould rewrite with loop up to log_3 (x) for improved speed
			return (sum(n=0, log(x) / log2, floor(log3x - tmp * n)))
		);
		result = 0;
		while (x >= 1,
			\\ Rewrite as above
			log3x = log(x) / log(3) + 1;
			result = result + sum(n=0, lg(x), floor(log3x - tmp * n));
			x = floor (x / 5);
		);
		return (result);
	);
	if (x <= B,
		if(x < 0, return(0));
		return(if (x * eps() < 2,
			floor (x)
		,
			x
		))
	);
	if (B < 100,
		\\ Buchstab identity
		return(sum(i=2,primepi(B),Psi(x\prime(i),prime(i)))+Psi(x,2))
	);

	if (B >= sqrt(x),
		trap(accurer,,
			return (floor(x) - primepi(x) + primepi(B))
		);
		if (x - B <= 10000,
			return (floor(x) - sum(n=floor(B)+1,x,ispseudoprime(n)))
		);

		result = x - li(x) + li(B);		\\ Schoenfeld's famous error bound, applied to both x and B.
		tmp = (sqrt(x)*log(x) + sqrt(B)*log(B))/(8*Pi);
		if (tmp < 1e9,
				u = 0;
				if (x - B < 1e9 && x * eps() < 2,
					u = floor(x - B);
					u = (u - u%210)*8/35 + (u%210)/3 + 1;	\\ Crude upper bound on primes on the interval.  Assumes B > 7.
					if (result > x - u,
						u = ceil(max(x-result, result-x+u));
						if (u < 10 * tmp, print("Estimate.  The error is (unconditionally) less than "u".")),
						u = 0
					);
				);
				if (u > tmp, print("Estimate.  On the RH, the error is less than "ceil(tmp)"."))
			,
			print("Estimate.  On the RH, the error is less than: "precision(tmp,9))
		);
		if (result * eps() < 2 && result < 1e19,
			return (floor(result)),
			return (result)
		)
	);

	\\ B^3 > x section?

	if (x <= 1e6,
		x = floor(x);
		
		return (sum(n = 1, x, gpf(n) <= B));
	);

	\\ At this point, x >= 1e6 and 7 <= B < sqrt(x).
	u = log(x) / log(B);	\\ Psi (x, B) = Psi (x, x ^ (1/u))
	tmp = x/log(x)^u;
	if (tmp > 500,
		if (tmp < 1e7 && tmp * eps() < 2,
			print("Konyagin and Pomerance 1997 lower bound: "floor(tmp)),
			print("Konyagin and Pomerance 1997 lower bound: "precision(tmp, 9))
		);
	);

	if (B < log(x),
		tmp = log(x)^primepi(B)/primepi(B)!;
		forprime(p = 2, B,
			tmp = tmp / log(p)
		);
		print ("Asymptotic estimate for small B.");
		return (precision(tmp, 9))
	);
	print ("Asymptotic estimate using Dickman's rho.");
	x * DickmanRho(u)
};
addhelp(Psi, "Calculates or estimates the count of smooth number up to the given bound. Psi(1000, 10) is the number of 10-smooth numbers up to 1000.");


default(timer, timervalue);

