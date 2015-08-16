timervalue = default(timer);
default(timer, 0);


\\ ***************************************************************************************************
\\ *					Working space
\\ ***************************************************************************************************

isRicePoly(P:pol, verbose:bool=1)={
	if(type(P)!="t_POL", error("Must be a polynomial"));
	if(poldegree(P)<2, error("Degree must be >= 2"));
	for(i=0,poldegree(P),
		if(type(polcoeff(P,i))!="t_INT",error("Must be an integer polynomial"))
	);
	if(pollead(P)!=1, return(0));	\\ Must be monic
	if(polcoeff(P,1), return(0));	\\ Linear coefficient must be 0
	my(c=polcoeff(P,0));
	if(abs(c)>=2, return(1));
	if(c==0, return(0));		\\ P(0) cannot be 0
	if(subst(P,variable(P),c)==0, return(0));	\\ P(P(0)) cannot be 0
	if(verbose,
		print("The sequence starting with "c" and iterating");
		print(variable(P)," |--> ",nice(P));
		print("is a superrigid divisibility sequence.")
	);
	1
};
addhelp(isRicePoly, "isRicePoly(P): Checks if a polynomial P fits the conditions for Proposition 3.2 in Brian Rice, 'Primitive Prime Divisors in Polynomial Arithmetic Dynamics', Integers 7:1 (2007).");

\\ ***************************************************************************************************
\\ *					Uncategorized
\\ ***************************************************************************************************

yafu(n)={
	if(type(n) != "t_INT", error("Bad type in yafu."));
	system(Str("~/mth/yafu 'factor(",n,")'"))
};
addhelp(yafu, "yafu(n): Factor a number with yafu.");


src(funcName:genstr)={
	my(cmd=Str("egrep --color=auto --after-context=10 --before-context=1 -R '^", funcName, "'"),C,t);
	C=[Str(cmd, " basemath"), Str(cmd, " kernel/gmp"), Str(cmd, " ."), Str("egrep --color=auto --after-context=5 -R '#define.*", funcName, "' .")];
	for(i=1,#C,
		t=Str("cd ~/mth/pari/src;", C[i]);
		if (#externstr(t),
			system(t);
			return()
		)
	)
};
addhelp(src, "src(funcName): Search the PARI source for a function (or #define) called funcName and display with context.");


polygonArea(v:vec)={
       if(type(v)!="t_VEC", error("Not a vector!"));
       if(#v<6 || #v%2, error("Need an even number of coordinates representing at least 3 points."));
       abs(sum(i=1,#v/2-1,v[2*i+2]*v[2*i-1]-v[2*i]*v[2*i+1],v[#v-1]*v[2]-v[#v]*v[1]))/2
};
addhelp(polygonArea, "polygonArea(v): Find the area of the polygon with points (v[1], v[2]), (v[3], v[4]), ..., (v[#v-1], v[#v]).");


chVec(u:vec,v:vec,umod:int=0,vmod:int=0)={
	if(umod,u=apply(n->Mod(n,umod),u));
	if(vmod,v=apply(n->Mod(n,vmod),v));
	vector(#u*#v,i,chinese(u[(i-1)%#u+1],v[(i-1)\#u+1]))
};
addhelp(chVec, "chVec(v, u, umod, vmod): Given a vector v mod vmod and a vector u mod umod, return a vector with one representative mod lcm(umod, vmod) which is u[i] mod umod and v[j] mod vmod.");


isWardC(v:vec)={
	my(aux=vector(#v),t);
	if(#v==0,return(1));
	aux[1]=v[1];
	for(n=2,#v,
		t=v[n];
		fordiv(n,d,if(d<n,t/=aux[d]));
		if(denominator(t)>1,
			print("Not a Property C sequence at n = "n", value should be a multiple of "v[n]/t);
			return(0);
		);
		aux[n]=t
	);
	aux
};
addhelp(isWardC, "isWardC(v): Does the vector v have Ward's (1939) property C? If not, return 0; if so, return the auxiliary sequence.");


plt(mn,mx,ff,flags:small=0,n:small=0)={
	my(o=mn\1-1,k);
	if(mx-mn < 2*n || (!n && mx-mn < 500),
		my(v=vector(ceil(mx)-o,i,ff(i+o)));
		ploth(x=max(mn,mn\1+eps()),min(mx,ceil(mx)-eps()),
			k=x\1;v[k-o]*(k+1-x)+v[k-o+1]*(x-k), flags, n
		)
	,
		ploth(x=max(mn,mn\1+eps()),min(mx,ceil(mx)-eps()),
			k=x\1;ff(k-o)*(k+1-x)+ff(k-o+1)*(x-k), flags, n
		)
	)
};
addhelp(plt, "plt(mn, mx, ff, flags, n): Make a high-resolution plot of the int");


shortestPath(G, startAt)={
	my(n=#G[,1],dist=vector(n,i,9e99),prev=dist,Q=2^n-1);
	dist[startAt]=0;
	while(Q,
		my(t=vecmin(vecextract(dist,Q)),u);
		if(t==9e99, break);
		for(i=1,#v,if(dist[i]==t && bittest(Q,i-1), u=i; break));
		Q-=1<<(u-1);
		for(i=1,n,
			if(!G[u,i],next);
			my(alt=dist[u]+G[u,i]);
			if (alt < dist[i],
				dist[i]=alt;
				prev[i]=u;
			)
		)
	);
	dist
};


gammainv(x:real)={
	if(x<5040,return(solve(y=1,8,gamma(y)-x)));
	my(L=log(x));
	solve(y=L/log(L),L,lngamma(y)-L)
};
addhelp(gammainv, "gammainv(x): Inverse gamma function.");


rand(b:small)={
	if(b<30,b--;return(factor(random(1<<b)+1<<b)));
	while(1,
		my(t=1,v=List());
		forprime(p=2,127,	\\ prime factors with up to 7 bits
			while(!random(p),t*=p; listput(v,p))
		);
		if (#binary(t)>b, next);
		if(#binary(t)==b,return(list2fac(v)));
		for(n=8,b,
			my(lambda=log(n/(n-1)),x=random(1.)/exp(-lambda),k=0);
			while(1,
				x-=lambda^k/k!;
				if(x<0,break);
				listput(v,rp(n));
				t*=v[#v];
				if(#binary(t)>b,next(3));
				if(#binary(t)==b,return(list2fac(v)));
			)
		);
			
	)
};
addhelp(rand, "rand(b): Gives a random factored n-bit integer.");


coin(v:vec)={
	if(type(v)!="t_VEC", error("Must be a vector"));
	for(i=1,#v,
		if(type(v[i])!="t_INT" || v[i]<1,error("Must be a vector of positive integers"))
	);
	v=vecsort(v,,8);
	my(g=vecgcd(v));
	if(g>1,
		print("All coins are a multiple of "g"; returning largest nonrepresentable multiple.");
		return(g*coin(v/g))
	);
	if(#v==1,return(0));
	if(#v==2,return(v[1]*v[2]-v[1]-v[2]));
	if(#v==3,
		my(a:int=v[1],b:int=v[2],c:int=v[3],A=gcd(b,c),B=gcd(a,c),C=gcd(a,b));
		if(A>1||B>1||C>1,
			\\ Johnson reduction
			\\ not quite right -- G(a,a,c) != G(a,c)
			error("fail");
			my(a1=a/B/C,b1=b/A/C,c1=c/A/B);
			return(A*B*C*(coin([a1,b1,c1])+a1+b1+c1)-a-b-c)
		);
		error("not quite implemented")
	);
	error("not implemented")
};
addhelp(coin, "coin(v): What is the largest number such that change cannot be made with zero or more coins of denomination v[1], v[2], .., v[#v]?  Usually called the coin problem or the Frobenius problem.");


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
	if (a <= 0 || issquare(D) || gcd(g,c) > 1, return(0));
	if ((a+b)%2 ==0 & c%2 == 0, return(0));

	f=factor(g)[,1];
	if(a+b%2, 1, 2)/sqrt(a) * prod(i=1,#f,f[i]/(f[i]-1)) * ConjectureF_C(D,lim)
};
addhelp(ConjectureF, "ConjectureF(a,b,c): Returns a number C such that there are ~ C*sqrt(x)/log(x) primes of the form an^2 + bn + c up to x, or 0 if there are o(sqrt(x)/log(x)) or Omega(sqrt(x)) such primes. Relies on the Hardy-Littlewood Conjecture F and possibly the ERH.");


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
		print("Keith Conrad"); \\ http://mathoverflow.net/questions/31150/calculating-the-infinite-product-from-the-hardy-littlewood-conjecture-f
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


forgap(ff)={
	if(default(parisize) < 1e8,
		print("Not enough memory!  Resizing, please rerun command.");
		allocatemem(100<<20)	\\ Need ~200 MB of stack space!
	);
	for(i=0,13,
		read(Str("gap200-chunk", if(i<10,"0",""), i));
		trap(e_USER,
			print("User error, ending loop...");
			pg=0;
			return()
		,
			for(j=1,#pg,
				ff(pg[j])
			)
		);
		pg=0
	)
};
addhelp(forgap, "forgap(ff): Runs the command (closure) ff on each prime 2 < p < 10^12 starting a prime gap of length 200 or greater.");


forpseudo(ff)={
	if(default(parisize) < 2e8,
		print("Not enough memory!  Resizing, please rerun command.");
		allocatemem(200<<20)	\\ Need ~200 MB of stack space!
	);
	for(i=0,89,
		read(concat("psp/psp-chunk", i));
		trap(e_USER,
			print("User error, ending loop...");
			pspChunk=0;
			return()
		,
			for(j=1,#pspChunk,
				my(t=ff(pspChunk[j]));
				if(t,return(t))
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
*/


dotproduct(a:vec, b:vec)={
	if(#a != #b,
		error("dotproduct: vectors are not the same size")
	);
	sum(i=1,#a,a[i]*b[i])
};
addhelp(dotproduct, "dotproduct(a, b): Returns the dot product of vectors a and b.");


\\ ***************************************************************************************************
\\ *	                                    Polynomials                                              *
\\ ***************************************************************************************************

finiteOrbit(f:pol, x:mp)={
	my(a=pollead(f), y=variable(f), v, d, b);
	if(a < 0, error("Not implemented: negative leading coefficient"));
	if(type(f) != "t_POL", error("Type error in finiteOrbit: should be polynomial"));
	d=poldegree(f);
	b=vecsum(abs(Vec(f)[2..d+1]))/a;
	if(b < 1, b = 1);
	v=List([x]);
	while(1,
		x=subst(f,y,x);
		if(x > b, return(0));
		for(i=1,#v,	\\ Slow loop -- need a better data structure.
			if(x==v[i], return(1));
		);
		listput(v, x)
	)
};
addhelp(finiteOrbit, "finiteOrbit(f, x): Is the forward orbit of x (under the polynomial f) finite?");


isCriticallyFinite(f:pol)={
	my(cp=nfroots(,f'), d=poldegree(f), y=variable(f));
	if(#cp < d-2 && sum(i=1,#roots,valuation(f,y-cp[i])) < d-1,
		error("Not implemented: polynomials with irrational critical points")
	);
	for(i=1,#cp,
		if(!finiteOrbit(f, cp[i]), return(0))
	);
	1
};
addhelp(isCriticallyFinite, "isCriticallyFinite(f): ");


RgX_check_ZX(P:pol, s:genstr)={
	for(i=0,poldegree(P),
		if(type(polcoeff(P,i)) != "t_INT", error("Type error: must be an integer polynomial in ",s))
	)
};
addhelp(RgX_check_ZX, "RgX_check_ZX(P, s): Given a polynomial P, checks if its coefficients are integers and returns an error if not. s is a string with the name of the calling function (for the error message). Mimics the library function of the same name.");


hasIntegerRoot(P:pol)={
	#select(x->type(x)=="t_INT", nfroots(,P)) > 0
};
addhelp(hasIntegerRoot, "hasIntegerRoot(P): Returns 1 if the polynomial P has at least one integer root and 0 otherwise.");


Eisenstein(P:pol)={
	my(d=poldegree(P),x=variable(P),g,t,times);
	P/=content(P);
	while(1,
		g=0;
		for(i=0,d-1,
			g=gcd(g,polcoeff(P,i));
			if(g==1,break)
		);
		while((t=gcd(g,polcoeff(P,d)))>1,
			g/=t
		);
		g=factor(g)[,1];
		for(i=1,#g,if(polcoeff(P,0)%g[i]^2,return([g[i],times,P])));
		times++;
		P=subst(P,x,x+1)
	)
};
addhelp(Eisenstein, "Eisenstein(P): Given an irreducible polynomial P, searches for some k such that P(x+k) can be proven irreducible by Eisenstein's criterion using a prime q. Output is [q, k, P(x+k)].");

\\ ***************************************************************************************************
\\ *	                                  Number theory                                              *
\\ ***************************************************************************************************

sqrtformal(n:int)={
	my([d,f]=core(n,1));
	if(d>1,quadgen(4*d),1)*f
};
addhelp(sqrtformal, "sqrtformal(n): Returns a number (t_QUAD or t_INT) representing the square root of n.");


factordb(n)={
	my(N=n,res,pr);
	if(type(n)=="t_STR",N=eval(n));
	if(type(N)!="t_INT", error("Bad input"));
	my(t=factor(N,2^24));	\\ Check if number is already factored (or tiny)
	if(#select(k->k>2^24 || ispseudoprime(k), t[,1]) == #t[,1],
		return(t)
	);
	\\t=alarm(1,factor(N));	\\ Pretest to remove easy composites
	\\if(type(t)!="t_ERROR",return(t));
	if(N<1e1000 && ispseudoprime(N),
		if(default(factor_proven) && !isprime(N),
			warning("Input appears to be a BPSW pseudoprime, PLEASE REPORT")
		,
			return(Mat([N,1]))
		)
	);
	res = extern(concat("~/mth/fdb_factor.pl ",n));
	if(!res, error("Could not find factors -- did you hit your hourly limit?"));
	if(type(res)=="t_VEC",res=Mat(res));
	if(type(res)!="t_MAT", error("fdb_factor.pl returned a non-matrix: ", res));
	pr=prod(i=1,#res~, res[i,1]^res[i,2]);
	if(N!=pr,
		if(N%pr, error("factordb error on input "n));
		res=concat(res,[N/pr,1])
	);
	for(i=1,#res~,
		my(p=res[i,1]);
		if(if(default(factor_proven),isprime(p),ispseudoprime(p)),
			if(p>2^24 && default(factor_add_primes), addprimes(p))
		,
			warning("composite factor C",#Str(p)," = ",precision(p*1.,9))
		)
	);
	res
};
addhelp(factordb, "factordb(n): Look up the factorization of n on factordb.com.");


\\ Varies as n*log(log(n))/log(n)
countsemi(n:int)={
	my(B=sqrtint(n),s=-binomial(primepi(B),2));
	forprime(p=2,B,
		s+=primepi(n\p)
	);
	s
};
addhelp(countsemi, "countsemi(n): Number of semiprimes up to n. Sloane's A072000.");


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
		while(1,
			p = precprime(random(upper - lower) + lower);
			q = nextprime(random(upper - lower) + lower);
			if (p > q,
				t = p;
				p = q;
				q = t;
			);
			if (2*p > q && #binary(p * q) == b, break)
		);
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


/*
bestappr2(x, n)={
	my(b=bestappr(x, n), cf=contfrac(x, #contfrac(b)+1), lower=ceil(cf[#cf]/2), upper=cf[#cf]);
	cf[#cf] = lower;
	contfracback(cf)
};
addhelp(bestappr2, "bestappr2(x, n): Work-in-progress.  Should determine the best rational approximation of x with denominator up to n.");
*/


\\ ***************************************************************************************************
\\ *                               Operations on generic lists                                       *
\\ ***************************************************************************************************

ternarySearch(f, left, right)={
	my(e=2*eps());\\left+right));
	left *= 1.;	\\ convert to t_REAL
	while(right-left>e,
		my(L=(2*left+right)/3,R=(left+2*right)/3);
		if(f(L) < f(R), left=L, right=R)
	);
	(left+right)/2
};
addhelp(ternarySearch, "ternarySearch(f, left, right): Assuming f is unimodal, find left <= x <= right such that f(x) is maximal.");


chop(v:vec,n:int)={
	if(#v<n,v,vector(n,i,v[i]))
	\\ vecextract?
};
addhelp(chop, "chop(v, n): Returns the first n terms of v, or all the terms if n > #v.");


Aitken(v:vec)={
	vector(#v-2, i, (v[i+2]*v[i]-v[i+1]^2)/(v[i+2]-2*v[i+1]+v[i]))
};
addhelp(Aitken, "Aitken(v): Applies Aitken's delta-squared process to v. The precision of v needs to be at least twice the desired precision to get meaningful results.");


LIP(X=0,Y)={
	my(n=#Y);
	if(type(X)!="t_VEC",X=vector(n,i,i));	\\ Default to X = [1,2,3,...]
	if (#X!=n, error("Vectors must be the same size!"));
	sum(j=1,n,
		Y[j]*prod(k=1,n,if(j==k,1,(x-X[k])/(X[j]-X[k])))
	)
};
addhelp(LIP, "LIP(X,Y): Gives the Lagrange interpolating polynomial for X=[x1,x2,...] and Y=[y1,y2,...]. If X is omitted, use [1,2,...].");


isperiodic(v:vec)={
	\\ k is the period
	for(k=1,(#v+1)\2,
		for(i=k+1,#v,if(v[i]!=v[i-k], next(2)));
		return(k)
	);
	print("Not periodic with period <= ",(#v+1)\2);
	0
};


isrec(v:vec,P:pol)={
	my(d=poldegree(P),c=vector(d,i,-polcoeff(P,d-i)));
	if (#v < 2*d,
		warning(Str("isrec: Not enough information to determine if the sequence is a ",d,"-coefficient recurrence relation; matrix is underdetermined. Give more terms or reduce d and try again."));
		return
	);
	for(n=d+1,#v,
		if(v[n] != sum(i=1,d,v[n-i] * c[i]),return(0))
	);
	1
};
addhelp(isrec, "isrec(v, P): Is the vector v a recurrence relation with characteristic polynomial P?");


minrec(v:vec,P:pol)={
	if (!isrec(v,P),error("Sequence is not generated by the given recurrence."));
	my(f=factor(P));
	for(i=1,#f[,1],
		for(j=1,f[i,2],
			if(!isrec(v,P/f[i,1]),break);
			P/=f[i,1]
		)
	);
	P
};
addhelp(minrec, "minrec(v, P): Given a vector v which is a recurrence relation with characteristic polynomial P, return the lowest-degree polynomial which represents a recurrence relation for v.");


ggf(v)={
	my( p, q, B=#v\2 ); if( B<4, "Need at least 8 values!",
 if( !#q=qflll(matrix(B,B,x,y,v[x-y+B+1]),4)[1], "QFLLL returned empty result.",
  polcoeff(p=Ser(q[,1]),0)<0 && p=-p; /* adjust sign */
  q=Pol(Ser(v)*p);
  if( Ser(v)==q/=Pol(p), q,
   Str("I guessed ",q,", but it doesn't reproduce the last value!")
)))
};
pgf(f)={
	my(p=P->concat( vecextract( Vec(Str( P+O(x^99) )), "..-11")));
	Str("(",p(numerator(f)),")/(",p(denominator(f)),")")  /* prettyprint rational g.f. */
};


\\ Helper function.
tidy(x)={
	if(type(x) == "t_COMPLEX",
		if(abs(imag(x))<=eps(), x=real(x));
		if(abs(real(x))<=eps(), x=imag(x)*I);
		return(x)
	);
	if(type(x) == "t_VEC",
		vector(#x,i,tidy(x[i]))
	,
		x
	)
};


\\ Only works in the case of distinct roots -- doesn't know about polynomials
findBinet(v:vec,offset:int=1)={
	my(c=findrec(v,0),poly=Pol(concat(-1,c)),t=polroots(poly),roots=tidy(t),coeff);
	coeff=matsolve(matrix(#c,#c,i,j,roots[i]^(j+offset-1)),vectorv(#c,i,v[i]));
	coeff=tidy(coeff);
	for(i=1,#c,
		print("("coeff[i]") * ("roots[i]")^n")
	);
	coeff
};
addhelp(findBinet, "findBinet(v, {offset}): Find a Binet-style formula for a linear recurrence given by the terms in the vector v. Assumes the roots are distinct -- no polynomial multipliers.");


findrec(v:vec, verbose:bool=1)={
	my(c,d = (#v - 1) >> 1, firstNonzero = 0);
	if (#v < 3,
		warning(Str("findrec: Not enough information to determine if the sequence is a recurrence relation: matrix is underdetermined. Give more terms and try again."));
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


rec(v:vec, verbose:bool=1)={
	my(c=findrec(v,verbose),gf=-Pol(Ser(c)*'x-1),poly=Pol(concat(-1,c)),f=factor(poly),roots=polroots(poly));
	if(c == 0, return(0));
	print("G.f. denominator: "gf);
	print("Roots have magnitude up to "precision(vecmax(abs(roots)),9)".");
	print("Roots: "roots);
	if(#f[,1]>1 || f[1,2]>1,
		print("Polynomial factors as: "nice(f));
	);
	c
};


findrecd(v:vec, d:int, verbose:bool=1)={
	my(M,c);
	if (#v < 2*d,
		warning(Str("findrec: Not enough information to determine if the sequence is a "d"-coefficient recurrence relation; matrix is underdetermined. Give more terms or reduce d and try again."));
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
			/*
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
			)*/
		);
		print1("<a href=\"/index/Rec#order_"if(d<10,"0","")d"\">");
		print1("Index to sequences with linear recurrences with constant coefficients</a>, signature ("c[1]);
		for(i=2,#c,print1(","c[i]));
		print(").");
		print((#v-2*d)" d.f.")
	);
	c
};
addhelp(findrecd, "findrecd(v, d, {verbose=1}): Helper function for findrec. Tries to find a d-coefficient homogeneous linear recurrence relation with constant coefficients to fit the vector v.");


initial(n:int, s:str)={	\\ Helper function
	if (n == 0, return(""));
	if (n == -1, return(Str("-", s)));
	if (n == 1, return(s));
	Str(n, s)
};


medial(n:int, s:genstr)={	\\ Helper function
	if (n == 0, return(""));
	if (n == -1, return(Str(" - ", s)));
	if (n == 1, return(Str(" + ", s)));
	Str(if (n < 0, " - ", " + "), abs(n), s)
};


diff(v:vec)={
	vector(#v-1,i,v[i+1]-v[i])
};
addhelp(diff, "diff(v): First difference of the vector v.");


complement(v:vec)={
	my(u=List(),n=1,i=1);
	while(i<=#v,
		if(v[i]==n,
			i++;
			n++
		,
			listput(u,n);
			n++
		)
	);
	Vec(u)
};
addhelp(complement, "complement(v): Returns to complement of v with respect to {1, 2, ..., v[#v]}. Assumes v is sorted.");


findCycle(ff:closure,startAt)={
	my(power=1,len=1,tortoise=startAt,hare=ff(startAt));
	while (tortoise != hare,
		if (power == len,
			tortoise = hare;
			power <<= 1;
			len = 0
		);
		hare = ff(hare);
		len++
	);
	
	my (mu=0);
	tortoise=hare=startAt;
	for (i=1,len,
		hare = ff(hare)
	);
	while(tortoise != hare,
		tortoise=ff(tortoise);
		hare=ff(hare);
		mu++
	);
print("mu = "mu);
	len
};
addhelp(findCycle, "findCycle(ff, startAt): Finds the length of the first cycle that startAt, ff(startAt), ff(ff(startAt)), ... enters into. Prints the prefix length (steps taken before the cycle begins).");


LinearRecurrence(sig:vec, initial:vec, terms:small)={
	if(terms<=#initial, return(vector(terms,i,initial[i])));
	if(#sig>#initial, error("Not enough terms to uniquely determine sequence"));
	my(v=vector(terms));
	for(i=1,#initial,v[i]=initial[i]);
	for(i=#initial+1,terms,
		v[i]=sum(j=1,#sig,sig[j]*v[i-j])
	);
	v
};
addhelp(LinearRecurrence, "LinearRecurrence(sig, initial, terms): Given a signature sig and a vector of initial terms, gives the first terms terms of that linear recurrence.");


\\ ***************************************************************************************************
\\ *					Statistics						*
\\ ***************************************************************************************************


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
	-real(suminf(n=1,
		moebius(n)/n * (eint1(u/n) + l2)
	))
};
addhelp(R, "Riemann's R function, an estimate for primepi.");

li(x)=-real(eint1(-log(x)));
addhelp(li, "Logarithmic integral.");

Li(x)=real(eint1(-log(2)) - eint1(-log(x)));
addhelp(Li, "Offset logarithmic integral, an estimate for primepi. Crandall and Pomerance call this li_0.");

piBounds(x,verbose=0)={
	my(lower,upper,t,srcLower,srcUpper,lowerRH,upperRH,lowerC,upperC,roundFlag);
	roundFlag=default(realprecision)>sizedigit(x);
	
	if(roundFlag, x = floor(x));
	if (x < default(primelimit) || x < 65557,
		lower=primepi(x);
		print("primepi("x") = "lower);
		return([lower,lower]);
	);
	my(lx=log(x),lix=li(x),Rx=R(x));	\\ Precompute these, we'll use them several times.

	\\ "PARI currently guarantees that the first 6547 primes, up to and
	\\ including 65557, are precomputed, even if primelimit is 1."
	\\ --User's Guide
	lower = x/lx * (1 + 1/lx + 1.8/lx^2);	\\ Dusart 1998, x >= 32299
	upper = x/(lx - 1.1);					\\ Dusart 1998, x >= 60184
	srcLower = srcUpper = "Dusart1998";
	
	if (x >= 88783,
		lower = x/lx * (1 + 1/lx + 2/lx^2);		\\ Dusart 2010, x >= 88783
		srcLower = "Dusart2010"
	);
	
	if (x >= 1332450001,
		lower = x/lx * (1 + 1/lx + 2/lx^2 + 5.65/lx^3 + 23.65/lx^4 + 118.25/lx^5 + 709.5/lx^6 + 4966.5/lx^7);		\\ Axler 2014, x >= 1332450001
		srcLower = "Axler2014"
	);
	
	if (x >= 355991,
		t = x/lx * (1 + 1/lx + 2.51/lx^2);	\\ Dusart 1998, x >= 355991
		if (upper > t,
			upper = t;
			srcUpper = "Dusart1998"
		)
	);
	
	if (x >= 2953652302,
		t = x/lx * (1 + 1/lx + 2.334/lx^2);	\\ Dusart 2010, x >= 2 953 652 302
		if (upper > t,
			upper = t;
			srcUpper = "Dusart2010"
		)
	);
	
	if (x >= 13041027276,
		\\ This is the best in [13041027276, 16526414332]
		t = x/lx * (1 + 1.0992/lx);			\\ Dusart 1998, x >= 1.332e10
		if (upper > t,
			upper = t;
			if(x < 1.332e10,
				srcUpper = "Dusart1998+finite"
			,
				srcUpper = "Dusart1998"
			)
		)
	);
	
	if (x >= 2,
		t = x/lx * (1 + 1/lx + 2/lx^2 + 6.35/lx^3 + 24.35/lx^4 + 121.75/lx^5 + 730.5/lx^6 + 6801.4/lx^7);	\\ Axler 2014, x >= 2
		if (upper > t,
			upper = t;
			srcUpper = "Axler2014"
		)
	);
	
	if(x >= 2 && x <= 1e18,
		if(upper > lix,
			upper = lix;
			srcUpper = "StollDemichael2011"
		)
	);

	t = sqrt(x)/8/Pi*log(x);
	lowerRH = lix-t;					\\ Schoenfield, x >= 2657
	upperRH = lix+t;					\\ Schoenfield, x >= 2657
	lowerC = lix-sqrt(x);				\\ Kotnik, x >= 2
	upperC = Rx+sqrt(x);				\\ Kotnik, x >= 2
	
	\\ Note: Stoll & Demichael have even tighter, even more conjectural bounds:
	\\ lix +- (lix-Rx)*(1+(log(log(lx))+1)/exp(1)). But this is so close --
	\\ +- 1/e vs. Omega_{+-}(1) -- to Littlewood's bound that it scares me a bit.
	
	\\ primepi(10770325941) = 488450930
	\\ 10770325941 being the Dusart-Axler crossover
	
	if (roundFlag,
		lowerRH = ceil(lowerRH);
		upperRH = floor(upperRH);
		lower = ceil(lower);
		upper = floor(upper);
		lowerC = ceil(lowerC);
		upperC = floor(upperC);
	);

	if (x <= 14316504791907170806555192,	\\ Büthe gives x <= 1.4e25, but this is the exact number based on Gourdon's verification up to height 2445999556030.
		if (lowerRH > lower,
			lower = lowerRH;
			srcLower = "Buthe2014"
		);
		if (upperRH < upper,
			upper = upperRH;
			srcUpper = "Buthe2014"
		);
	);

	\\print("For primepi("x"):");
	print(lower" (lower bound)");
	if (lowerRH > lower,
		print(lowerRH" (lower bound under the RH)");
	);
	if (lowerC > max(lowerRH, lower),
		print(lowerC" (conjectural lower bound)");
	);
	if(roundFlag,
		print(round(Rx)" (Riemann R, approximate)");
		print(round(lix)" (logarithmic integral, apx)");
	,
		print(Rx" (Riemann R, approximate)");
		print(lix" (logarithmic integral, apx)");
	);
	if (upperC < min(upperRH, upper),
		print(upperC" (conjectural upper bound)");
	);
	if (upperRH < upper,
		print(upperRH" (upper bound under the RH)");
	);
	print(upper" (upper bound)");
	
	if (verbose,
		if (srcLower == "Dusart1998" || srcUpper == "Dusart1998" || srcUpper == "Dusart1998+finite",
			if(srcLower == "Dusart1998" && (srcUpper == "Dusart1998" || srcUpper == "Dusart1998+finite"),
				print("\nUpper and lower bounds:");
			,
				print(if(srcLower == "Dusart1998", "\nLower bound:", "\nUpper bound:"))
			);
			print("Pierre Dusart, 'Autour de la fonction qui compte le nombre de nombres");
			print("premiers', doctoral thesis for l'Universite de Limoges (1998).");
			if(srcUpper == "Dusart1998+finite",
				print("  plus finite checking (crg4)")
			)
		);
		
		if (srcLower == "Dusart2010" || srcUpper == "Dusart2010",
			if(srcLower == "Dusart2010" && srcUpper == "Dusart2010",
				print("\nUpper and lower bounds:")
			,
				print(if(srcLower="Dusart2010","\nLower bound:","\nUpper bound:"))
			);
			print("Pierre Dusart, 'Estimates of some functions over primes without R.H.',");
			print("preprint (2010), arXiv:1002.0442.");
		);
		
		if (srcLower == "Axler2014" || srcUpper == "Axler2014",
			if(srcLower == "Axler2014" && srcUpper == "Axler2014",
				print("\nUpper and lower bounds:")
			,
				print(if(srcLower="Axler2014","\nLower bound:","\nUpper bound:"))
			);
			print("Christian Axler, New bounds for the prime counting function pi(x), preprint");
			print("(2014), arXiv:1409.1780.");
		);
		
		if(srcUpper == "StollDemichael2011",
			print("\nUpper bound:");
			print("Douglas A. Stoll and Patrick Demichel, The impact of zeta(s) complex zeros on");
			print("pi(x) for x < 10^10^13, Math. Comp. 80:276 (2011), pp. 2381-2394.")
		);
		
		t=(srcUpper == "Buthe2014")+(srcLower == "Buthe2014");
		if(t,
			if(t==2,
				print("\nUpper and lower bounds:")
			, srcUpper == "Buthe2014",
				print("\nUpper bound:")
			,
				print("\nLower bound:")
			);
			print("Jan Buthe, Estimating pi(x) and related functions under partial RH assumptions,");
			print("preprint (2014), arXiv:1410.7015.");
		);
		
		t=(lowerRH > lower) + (upperRH < upper);
		if (t,
			if(t==1,
				print("\nBound on the RH:");
			,
				print("\nBounds on the RH:");
			);
			print("Lowell Schoenfeld, 'Sharper Bounds for the Chebyshev Functions theta(x) and");
			print("psi(x). II'. Mathematics of Computation, Vol 30, No 134 (Apr 1976),");
			print("pp. 337-360.");
		);

		t=(lowerC > max(lowerRH, lower)) + (upperC < min(upperRH, upper));
		if(t,
			if(t==1,
				print("\nConjectural bound:")
			,
				print("\nConjectural bounds:")
			);
			print("Tadej Kotnik, The prime-counting function and its analytic approximations,");
			print("Adv. Comput. Math. 29 (2008), pp. 55-70.");
		);

		print();
	);
	[lower, upper]
};
addhelp(piBounds, "piBounds(x, verbose=0): Bounds on primepi(x). Set verbose=1 to get a list of sources for the results.");


pBounds(n, verbose:bool=0)={
	my(lower,upper,appx,l,ll,ali=x->solve(y=x*log(x),x*log(x)*log(log(x)),x+eint1(-log(y))));
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
	\\appx  = n * (l + ll - 1 + ll/l - 2/l - ll^2/2/l^2 + 3*ll/l^2 + 11/2/l^2);	\\ + O(ll^3/l^3)
	appx = ali(n);
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

Psi(x, B)={
	local(log3x, log2, tmp, u, result);
	if (x < 1, return (0));
	if (B < 7,
		if (B < 2, return (B >= 1 && x >= 1));
		log2 = log(2);
		if (B < 3, return (floor(log(x\1+.5) / log2) + 1));
		log3x = log(x) / log(3) + 1;	\\ The + 1 accounts for the pure powers of 2.
		tmp = log2 / log(3);
		if (B < 5,
			\\ Should rewrite with loop up to log_3 (x) for improved speed
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
	if (B >= x,
		if(x < 0, return(0));
		return(if (x * eps() < 2,
			floor (x)
		,
			x
		))
	);

	if (B^2 >= x,
		trap(/*accurer*/,,
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

	if (B < 257,
		\\ Buchstab identity
		return(sum(i=2,primepi(B),Psi(x\prime(i),prime(i)))+Psi(x,2))
	);
	error("B too large: "B);
	
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
addhelp(Psi, "Psi(x, B): Calculates or estimates the count of B-smooth number up to x. Psi(1000, 10) counts the number of numbers up to 1000 which have no prime factor greater than 10.");


default(timer, timervalue);

