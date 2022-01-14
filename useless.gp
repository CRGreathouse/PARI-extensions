\\ ***************************************************************************************************
\\ *                                    One-time use functions                                       *
\\ ***************************************************************************************************


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


A024916(z)=
{
	my(s=z^2,p,n);
	for(d=1, sqrtint(z),
		n = z\d - z\(d+1); 
		if(n<2, p=d; break, s -= (z%d2+(n-1)*d)*n/2);
	);
	for(d=2, z\if(p, p, sqrtint(z)+1),
		s -= z%d
	);
	s
};

pack(v)={
	sum(i=1,#v,1<<v[i])>>1
};
unpack(n)={
	my(t=1,v=vector(hammingweight(n)),i);
	while(n,
		if(bitand(n,1), v[i++]=t);
		t++;
		n>>=1
	);
	v
};
step(V,n)= {
	my(N=2^(n-1),u=List(V),v,t);
	for(i=1,#V,
		v=unpack(V[i]);
		forprime(p=2,n,
			if(#vecsort(v%p,,8)==p,next(2))
		);
		listput(u,V[i]+N)
	);
	Vec(u)
};
steps()={
	my(V=[0],i);
	for(n=1,99,
		if(gcd(n,10)>1,next);
		i++;
		V=step(V,n);
		my(MB=sizebyte(V)/2.^20);
		print1(n" "#V" ("round(MB)" MB; est. ");
		print(round(#V^(40/i-1)*MB)," MB to finish)")
	)
};



\\ Find 2-element members of J_Z[i] such that
\\ * Neither element is 0
\\ * The elements are not the negatives of each other
findGaussianPair(lim)={
	lim\=1;
	my(ok=v->isJ(v)&&v[1]!=-v[2]&&v[1]&&v[2]&&real(v[1])>1&&real(v[2])>1);
	for(a=-lim,lim,
		for(b=-lim,lim,
			my(m=3*(b^2-a^2)+6*a+1,n=6*b*(1-a),x=m+n*I);
			if(n,
				fordiv(n/2,d,
					my(e=n/(2*d));
					if((d+e*I)^2==x,
						my(v);
						v=[a+b*I,(a+b*I+1+d+e*I)/2];
						if(ok(v),print(v));
						v=[a+b*I,(a+b*I+1-d-e*I)/2];
						if(ok(v),print(v));
					)
				)
			,
				my(k);
				if(issquare(m,&k),
					my(v);
					v=[a+b*I,(a+b*I+1+k)/2];
					if(ok(v),print(v));
					v=[a+b*I,(a+b*I+1-k)/2];
					if(ok(v),print(v));
				)
			)
		)
	)
};
isJ(v)={
	sum(i=1,#v,v[i]^3)==sum(i=1,#v,v[i])^2
};


\\ A015704
\\ best(p, cur) gives an upper bound on the value of sigma(n)/n in [0..L]
\\ for numbers divisible by cur but not divisible by p.
best(avoid,cur=1,L=336280120525440)={
	my(q=2,r,rAt,t);
	if(q==avoid,q=nextprime(q+1));
	while(cur<L,
		r=0;
		forprime(p=2,q,
			if(p==avoid,next);
			t=valuation(cur,p);
			t=sigma(p^(t+1),-1)/sigma(p^t,-1);
			if(t>r,r=t;rAt=p)
		);
		cur*=rAt;
		if(rAt==q,
			q=nextprime(q+1);
			if(q==avoid,q=nextprime(q+1))
		)
	);
	sigma(cur,-1)
};


dosomething(lim,startAt)={
	my(q=precprime(startAt-2),p=precprime(q-1),k,pr=primorial(200)/6);
	forbigprime(r=q+2,lim,
		if(r-p!=6 || r-q!=2 || (k=r+2)%715,
			p=q;
			q=r;
			next
		,
			p=q;
			q=r;
		);
		\\k=q+2;
		if(gcd((k-32)*(k-64)*(k-128)*(k-256)*(k-512)*(k-1024)*(k-2048)*(k-4096),pr)>1,next);
		for(x=4,13,if(!isprime(k-2^x),next(2)));
		return(k)
	)
};


\\ Right now this generates A187060, though it could easily be modified to generate similar sequences.
makeList(lim,startAt=2)={
	my(v=vector(14),u=List(),p);
	v[2]=nextprime(startAt);
	for(i=3,#v,v[i]=nextprime(v[i-1]+1);lim=nextprime(lim+1));
	lim=nextprime(nextprime(lim+1)+1);
	forbigprime(q=v[#v]+1,lim,
		v=concat(vecextract(v,16382),q);
		p=v[1];
		if(v[3]!=p+6 || v[2]!=p+2, next);
		for(x=3,7,if(!setsearch(v,p+x+x^2),next(2)));
		listput(u,p)
	);
	Vec(u)
};


\\ Functions for displaying + and * tables for finite fields
randpol(p,n)={
	my(P=Pol(vector(n+1,i,random(p)))*Mod(1,p),x='x);
	if(!polcoeff(P,0),P+=1);
	if(!polcoeff(P,n),P+=x^n*Mod(1,p));
	while(!polisirreducible(P) || poldegree(P)<n,
		P+=x^random(n+1)*Mod(1,p)
	);
	P
};
tables(P)={	\\ tables(randpol(5,2))
	my(p=polcoeff(P,0).mod,n=poldegree(P),v=vector(p^n,i,Polrev(vector(n,j,(i-1)\p^(j-1)%p))),S,t,disp);
	disp=if(p^n-p<=26,
		t->my(n=sum(i=0,n-1,polcoeff(t,i)*p^i));if(n<p,n,Strchr([97+n-p]))
	,
		t->my(n=sum(i=0,n-1,polcoeff(t,i)*p^i));if(n<p,n,Str("e_{",n-p,"}"))
	);
	
	print1("\\begin{tabular}{");
	for(i=0,p^n,print1("c"));
	print1("}\n$+$&");
	for(i=0,p^n-2,
		print1(disp(i)"&")
	);
	print(disp(p^n-1),"\\\\ ");
	for(i=1,p^n,
		print1(disp(i-1));
		for(j=1,p^n,
			t=disp(lift(v[i]+v[j]*Mod(1,p)));
			print1("&"t)
		);
		print("\\\\ ")
	);
	print("\\end{tabular}");

	print("\\vspace{1em}");

	print1("\\begin{tabular}{");
	for(i=0,p^n,print1("c"));
	print1("}\n$\\times$&");
	for(i=0,p^n-2,
		print1(disp(i)"&")
	);
	print(disp(p^n-1),"\\\\ ");
	for(i=1,p^n,
		print1(disp(i-1));
		for(j=1,p^n,
			t=disp(lift(v[i]*v[j]*Mod(1,p)));
			print1("&"t)
		);
		print("\\\\ ")
	);
	print("\\end{tabular}");
};

rootsmod(f:pol,m:int)={
	my(v=List(),x=variable(f));
	for(n=0,m-1,
		if(subst(f,x,Mod(n,m))==0,listput(v,n))
	);
	Vec(v)
};


smallestNeededToAdd(v:vec,n:int)={
	v=vecsort(v,,8);
	if(v[1] <= 0,
		if(v[1]==0,v=vecextract(v,"2.."));
		if(v[1]<0,error("can't handle negatives"))
	);
	if(v[1]==1,
		snta1(v,n,#v)
	,
		error("not implemented")
	)
};
snta1(v:vec,n:int,a:small)={
	if(a==1,return(n));
	my(best=snta1(v,n,a-1),t);
	for(k=1,n\v[a],
		t=k+snta1(v,n-v[a]*k,a-1);
		if(t<best,best=t)
	);
	best
};


Witness(n)={
       my(b=2);
       while(ispower(b) || sprp(n,b), b++);
       b
};
Wspecial(n)={
       \\ For numbers known to be {3, 7, 11, 13, 17}-strong pseudoprimes
       if(!sprp(n,2) || !sprp(n,5), return(0));
       my(v=[6,10,12,14,15,18],b=19);
       for(i=1,#v,if(!sprp(n,v[i]),return(v[i])));
       while(ispower(b) || sprp(n,b), b++);
       b
};
acceptable(w)={
       if(w<18, return(0));
       w != 19 && w != 22 && w != 23 && w != 26 && w != 34 && w != 37
};
\\ Construct pseudoprimes with Arnault's method.
construct(startAt=0)={
       \\ 2 and 7 to 17: Mod(54607, 204204), Mod(8910, 204204)
       my(a=54607,b=8910,m=204204,p1,p2,N);
       a += startAt * m;
       b += startAt * m;
       until(ispseudoprime(p2=a^2+b^2) && isprime(p1=(p2+1)/2) && isprime(p2) && acceptable(Wspecial(N=p1*p2)),
               a += m;
               b += m
       );
       print("Least witness "Witness(N));
       print("Coefficient of m: ",a\m);
       N
};
construct1(startAt=0)={
       \\ 2 and 7 to 19: Mod(2913463, 3879876), Mod(2255154, 3879876)
       my(a=2913463,b=2255154,m=3879876,p1,p2,N);
       a += startAt * m;
       b += startAt * m;
       until(ispseudoprime(p2=a^2+b^2) && isprime(p1=(p2+1)/2) && isprime(p2=a^2+b^2) && acceptable(Wspecial(N=p1*p2)),
               a += m;
               b += m
       );
       print("Least witness "Witness(N));
       print("Coefficient of m: ",a\m);
       N
};


\\ A058445
\\ ok(n)=n=eval(Vec(Str(n)));for(i=1,#n,if(n[i]&&n[i]!=5&&n[i]!=6,return(0)));1
\\ do(m,u=0)=if(m==1,return([4,5,6]));my(v=List(),M=10^m,n);if(!u,u=do(m-1));forstep(d=0,M-1,M/10,for(i=1,#u,n=d+u[i];if(ok(n^2%M),listput(v,n))));Vec(v)
\\ v=do(1);
\\ v=do(digits(v[#v])+1,v);#v
searchA058445()={
	my(sz=digits(v[#v]),start,stop,T);
	stop=setsearch(v,ceil(sqrt(20/3)*10^(sz-2)),1);
	print("Searching for tiny solutions (through "floor(sqrt(5)*10^(sz-1))")...");
	for(i=1,stop,
		if(ok(v[i]^2),print("** "v[i]))
	);
	
	start=setsearch(v,floor(sqrt(5)*10^(sz-1)),1);
	stop=setsearch(v,ceil(sqrt(20/3)*10^(sz-1)),1);
	print("Searching for small solutions (through "floor(sqrt(5)*10^sz)")...");
	for(i=start,stop,
		if(ok(v[i]^2),print("** "v[i]))
	);
	
	start=setsearch(v,floor((sqrt(5)-2)*10^sz),1);
	stop=setsearch(v,ceil((sqrt(20/3)-2)*10^sz),1);
	T=2*10^sz;
	print("Searching for moderate solutions (through "floor(sqrt(5)*10^(sz+1))")...");
	for(i=start,stop,
		if(ok((T+v[i])^2),print("** ",T+v[i]))
	);
	
	T=10^sz;
	print("Searching for large solutions (through "floor(sqrt(5)*10^(sz+2))")...");
	forstep(d=22*T,25*T,T,
		for(i=1,#v,
			if(ok((d+v[i])^2),print("** ",d+v[i]))
		)
	);

	T=[223,224,225,226,234,235,236,237,238,239,244,245,246,247,254,255,256,257,258]*10^sz;
	print("Searching for huge solutions (through "floor(sqrt(5)*10^(sz+3))")...");
	for(i=1,#T,
		for(j=1,#v,
			if(ok((T[i]+v[j])^2),print("** ",T[i]+v[j]))
		)
	)
};


\\ Find terms of A145746.
okA145746(n,m)={
	if(m%2,
		if(n%2,
			!issquare(n)
		,
			issquare(n) || issquare(n/2)
		)
	,
		\\ This is the 'hot' path, taken 1 - (5/9)^d of the time, or 99.5% for 9-digit numbers.
		if(n%2,
			\\ About all the benefit is here.
			issquare(n)
		,
			1
			\\ Don't pretest here at all -- the effort saved is minimal and the cost of testing
			\\ is relatively high.
			\\!issquare(n)&&!issquare(n>>1)
		)
	) && sigma(n)-n==m
};
P=vector(1000,i,if(i<112,0,dprod(i-1)));
check1()={
	my(P0);
	for(t=1,999,
		P0=dprod(t);
		if(sigma(t)-t==P0, print1(t", "))
	)
};
check2()={
	my(PA,A,PB,t);
	for(a=1,999,
		PA=dprod(a);
		if(!PA,next);
		A=1000*a;
		for(b=111,999,
			PB=P[b+1];
			if(!PB,next);
			PB*=PA;
			t=A+b;
			if(sigma(t)-t==PB, print1(t", "))
		)
	)
};
check3()={
	my(PA,A,PB,B,PC,t);
	for(a=1,999,
		PA=dprod(a);
		if(!PA,next);
		A=1000*a;
		for(b=111,999,
			PB=P[b+1];
			if(!PB,next);
			PB*=PA;
			B=1000*(A+b);
			for(c=111,999,
				PC=P[c+1];
				if(!PC,next);
				PC*=PB;
				t=B+c;
				if(okA145746(t,PC), print1(t", "))
			)
		)
	)
};
check4()={
	my(PA,A,PB,B,PC,C,PD,t);
	for(a=1,999,
		PA=dprod(a);
		if(!PA,next);
		A=1000*a;
		for(b=111,999,
			PB=P[b+1];
			if(!PB,next);
			PB*=PA;
			B=1000*(A+b);
			for(c=111,999,
				PC=P[c+1];
				if(!PC,next);
				PC*=PB;
				C=1000*(B+c);
				for(d=111,999,
					PD=P[d+1];
					if(!PD,next);
					PD*=PC;
					t=C+d;
					if(okA145746(t,PD), print1(t", "))
				)
			)
		)
	)
};


\\ A176494 (without term 341)
\\[2,3,1,1,2,1,2,1,2,4,1,3,2,1,2,4,4,1,3,2,1,3,2,4,3,2,1,2,1,2,47,2,6,1,8,1,3,5,2,4,4,1,6,1,2,1,5,5,2,1,2,4,1,8,4,6,8,1,3,2,1,4,7,2,1,2,9,791,4,1,2,8,3,9,5,2,4,3,2,3,8,1,6,1,3,2,4,3,2,1,2,4,3,2,3,2,9,6,1,5,7,4,4,8,1,3,4,4,16,1,3,9,2,1,5,6,1,2,8,4,1,5,2,6,3,14,3,8,3,5,9,2,3,11,2,3,2,7,6,636,1,6,1,2,1,4,5,2,1,2,11,2,1,2,291,2,3,8,3,2,6,4,7,2,10,4,3,11,5,2,10,1,6,1,3,4,1,6,1,3,35,2,1,2,4,4,3,5,5,6,1,8,3,6,4,6,3,5,2,8,4,1,3,5,12,6,1,2,8,1,3,2,1,2,4,1,3,6,6,8,3,5,8,9,2,1,2,4,3,2,1,3,5,10,1,2,1,2,4,6,6,3,5,11,2,4,3,2,3,2,587,2,6,1,2,12,1,3,4,39,9,2,1,9,2,1,4,6,1,6,3,7,5,17,11,17,2,1,9,6,6,3,4,7,11,2,1,2,1,4,10,10,8,6,1,4,1,14,8,3,9,2,1,2,3,7,4,1,8,11,6,4,6,1,2,1,4,10,1,4,1,3,2,1,4,3,9,93,2,22,3,4,1,2,3,4,1,2,3,11,"?",2,4,1,6,8,1,3,2,4,5,15,2,1,3,2,4,6,8,10,7,2,6,8,3,5,2,13,7,6,3,2,8,1,3,20,10,1,4,3,2,4,11,6,1,2,3,7,17,2,1,2,4,3,2,1,3,4,1,6,3,2,11,6,12,1,3,2,1,4,5,2,22,3,7,2,12,3,6,4,12,3,8,8,9,2,8,4,1,12,1,10,3,2,7,2,3,9,7,11,2,8,11,12,1,4,101,2,1,13,15,2,4,8,3,2,3,6,1,2,1,11,12,1,3,9,7,2,4,1,57,2,4,1,6,1,4,1,4,10,13,4,3,2,1,2,1,10,3,9,7,4,1,2,12,1,3,4,1,5,6,1,6,11,5,2,4,3,15,13,4,6,1,2,5,487,2,3,8,3,5,19,6,8,1,12,5,2,6,10,1,6,7,2,1,4,19,2,3,7,2,6,1,2,8,1,5,2,23,8,18,1,2,12,4,1,3,8,1,3,4,8,6,1,3,5,4,1,2,12,21,2,1,7,7,10,4,1,20,1,6,1,10,4,1,14,1,6,4,3,58091,6,1,8,3,9,4,6,10,8,6,1,5,11,2,4,4,10,1,10,11,5,2,1,2,6,1,9,15,4,3,11,9,7,4,1,2,16,1,3,13,4,4,130,6,6,1,3,2,10,3,5,2,1,2,4,1,381,2,5,7,10,4,16,5,5,4,10,1,2,10,3,11,4,1,2,11,9,5,2,4,1,6,1,4,6,3,15,2,6,144,1,3,5,5,12,10,7,2,3,8,3,7,7,2,1,2,6,6,8,1,4,22,12,349,8,11,2,10,6,12,12,3,5,2,1,5,7,2,1,3,6,6,1,2,6,1,2,4,1,3,2,15,5,7,2,3,6,5,4,1,3,2,1,2,1,4,5,5,2,3,6,4,15,2,1,3,7,2,3,6,4,6,3,11,5,2,12,1,3,2,4,1,8,1,4,6,13,2,4,9,10,6,7,2,12,3,5,2,4,7,45,5,4,1,6,8,3,6,1,6,3,5,14,6,1,2,3,7,2,3,6511,10,4,16464,1,3,4,8,1,6,4,8,8,3,5,4,4,1,3,9,9,4,3,5,7,56,1,11,2,3,6,34,3,35,2,1,4,4,1,3,2,8,3,7,103,4,24,1,29,7,4,1,6,1,3,4,5,2,10,1,4,1,34,1,8,5,2,1,2,3,4,8,108,1,4641,4,3,2,8,6,1,8,1,3,2,4,4,3,7,2,4,8,3,7,2,11,11,8,6,16,6,8,1,18,3,71,10,4,7,2,1,3,4,3,5,2,4,6,15,6,1,7,38,1,4,1,9,41,7,19,4,10,6,1,7,2,12,1,8,3,14,4,12,3,2,14,1,12,1,5,2,30,4,1,25,2,5,5,2,8,7,2,3,7,2,3,2,7,5,2,9,13,2,1,5,2,6,6,8,14,6,19,5,2,1,2,26,10,12];


/*
\\ L: Prime at the beginning of a large (>= 400) prime gap
\\ U: Corresponding upper prime.  There are no primes p with L[i] < p < U[i].
\r pg.gp
L=vector(#U,i,precprime(U[i]-1));for(i=1,#L,if(!isprime(L[i]),print("Bad news at index #"i)))
for(n=401,1e9,t=ff2(n);if(type(t)!="t_INT",return(n-1));print1(t","))
*/

\\ Find A110835(n). Designed for n <= 400.  First 400:
\\ 8, 4, 8, 6, 18, 15, 17, 25, 13, 20, 29, 44, 87, 81, 35, 83, 79, 74, 70, 67, 118, 330, 58, 223, 172, 229, 179, 471, 292, 360, 506, 367, 586, 577, 645, 545, 424, 743, 503, 637, 766, 467, 937, 579, 698, 683, 542, 1443, 641, 628, 616, 604, 2026, 1661, 571, 1834, 551, 2989, 1820, 2242, 2842, 2515, 2475, 2938, 2399, 2849, 4960, 2293, 5227, 5290, 5215, 4695, 2136, 5004, 2079, 4872, 2025, 1999, 4687, 7210, 1925, 12852, 4461, 4408, 4243, 15273, 4256, 10544, 12809, 5468, 4069, 11944, 3878, 3939, 13826, 3857, 13992, 13849, 14589, 16718, 3666, 13306, 16231, 16075, 19150, 12804, 17662, 18618, 17338, 23980, 12158, 12118, 41935, 16626, 13582, 17334, 17186, 32424, 16897, 13016, 31620, 84785, 47876, 30855, 16086, 47891, 60010, 15709, 105630, 55618, 46063, 28985, 132722, 85978, 34462, 61921, 61469, 98741, 144060, 121798, 143866, 124310, 32534, 123864, 117598, 275658, 75605, 257360, 114441, 135234, 137956, 238629, 257451, 235530, 254129, 252500, 293037, 131844, 296359, 130196, 576404, 235119, 290348, 565860, 278829, 125490, 443629, 123996, 413263, 350806, 528878, 654653, 457617, 404577, 622617, 610992, 510950, 117030, 264395, 689854, 1119070, 671235, 258616, 663939, 255820, 578143, 253084, 110805, 962318, 249088, 639660, 471032, 632978, 107378, 242701, 555908, 1164043, 622961, 454463, 959564, 449941, 604776, 1931012, 940749, 936160, 1034749, 927115, 1041676, 918243, 1751793, 1834579, 905249, 890590, 1329399, 1007761, 3014698, 1972560, 880334, 2358631, 862253, 1773735, 854485, 2343606, 1642306, 1937410, 2285576, 1885663, 2849036, 3450284, 1236919, 5602894, 3054872, 2382584, 1654257, 1856481, 1847107, 6070801, 1833080, 4372938, 5343597, 6009329, 4968769, 5980116, 6175469, 5234544, 1906470, 5192159, 5859549, 3791939, 5177074, 9167102, 5089140, 5914788, 5049068, 6193360, 7161197, 8953084, 5576666, 9804705, 15479668, 6484011, 4894898, 4921173, 6410329, 7000178, 6959606, 9510931, 6220845, 10419040, 1615826, 8367364, 5342530, 8428361, 9267951, 15733435, 5265102, 15436508, 10225681, 5208488, 8098413, 19972759, 15162811, 11245931, 14171527, 20677825, 20605525, 19555210, 14846919, 7961739, 8391828, 14693858, 13159626, 19154762, 18675031, 29996960, 14445651, 35087063, 8166544, 29920233, 14341358, 13015240, 14246382, 14111923, 7568890, 28213972, 20559987, 38079851, 35800634, 34714902, 33615670, 27669651, 41949273, 31934314, 37278592, 41497004, 33945901, 13572263, 13529583, 36694288, 13445023, 40721359, 59414803, 51988611, 32163141, 36016855, 34636202, 45847407, 51246260, 39731174, 48974757, 11609096, 45156934, 32212927, 48388233, 32020611, 68645124, 41739096, 59859254, 49583402, 47534323, 86284290, 58680755, 111128157, 72837448, 72626325, 58117850, 30913270, 67022439, 107967669, 57453646, 30560982, 30474161, 63362847, 65154694, 70580513, 87038962, 30047352, 62477891, 64247247, 175502367, 140695500, 138973543, 61617314, 56807825, 174838837, 61112254, 120615817, 128788491, 137685707, 60451581, 114422833, 135237695, 67174483, 59805040, 66816219, 169723871, 134855982, 168825861, 130096671, 95191395, 65763995, 65591838, 65420580, 94199818, 165756300, 93711736, 219786991, 122149909, 164051865, 92750590, 163212725, 108807700, 102460630, 63594117, 101941842, 91345278, 63113557, 125414410, 106898793, 106631546
\\ Next several:
\\ 316919127,316130771,62173901,105575788,203681025,61714488,324216795,61411966,346812469,321844477,396569020,308457694,103275105,280083082,87163205,307635420,153036392,355200739,152305908,363925867,432842984,587648645,150865663,100595798,100359102,698359333,99889036,742400756,424771320,297619383,550053070,305454249,674913703,555669641,303347668,464520255,697024047,695432668,691051151,96937769,1320493795,969040287,940425950,288234988,826441680,934100215,286300525,1007176777,662585915,1641850951,329210441,671175786,1356484445,1646323897,1350521876,1068611381,663832506,1586503974,648150492,1606158539,395286109,521992693,655229925,995022651,921109262,517512069,1951172971,890189521,2024453776,1307420114,1304644275,642736134,509853328,1028031202,1293657797,1237600167,1947917667,1838182582,635907116,867934783,2525614752,3297039415,1822606623,1526514314,613404280,626747960,1517110735,2219954014,2573693081,3048598334,2811841042,1501692943,2730820103,2529424186,1241388795,1489582516,2977996441,609179629,834887166,3936377113,4402180384,2806166322,4219755279,1219221138,2789496027,912431838,821713404,4096868979,5912718456,3859193248,2701788555,2751358386,899981501,8445255552,4962751468,1190867158,5057914334,2599023766,2714249506,4081801741,4627622940,2579107875,4773702121,4990346776,6422015302,5414029916,4961938730,2549799831,2797852989,6361430252,6369357446,3912047822,2525880508,3313115872,2633075689,5050473671,6298191441,4860486451,4851468851,7714135917,4833533661,11838904881,2594282677,3617993670,8940613046,3604740946,9361483090,2456741443,3585042908,4754439474,13107805402,4737213244,9259911845,4720111391,10913263380,7769343669,3533552166,12653686592,18748243789,4669538769,12793516242,13795952741,2391286520,8639422181,3483519569,4620038358,3471232022,3465120698,4595679632,12591513354,7295329939,13554764756,17081555860,9581514843,2341381410,7536217628,8874750867,13858930276,7194530907,17432644005,19074688582,12331894522,18273609420,9417447808,15972808608,4462357868,12226852831,4447179780,14739060524,7060395585,8664519882,4417131268,4409682480,16477662471,12062458171,15911880607,17845082566,15858663615,13514565138
\\ 
\\vv=[8, 4, 8, 6, 18, 15, 17, 25, 13, 20, 29, 44, 87, 81, 35, 83, 79, 74, 70, 67, 118, 330, 58, 223, 172, 229, 179, 471, 292, 360, 506, 367, 586, 577, 645, 545, 424, 743, 503, 637, 766, 467, 937, 579, 698, 683, 542, 1443, 641, 628, 616, 604, 2026, 1661, 571, 1834, 551, 2989, 1820, 2242, 2842, 2515, 2475, 2938, 2399, 2849, 4960, 2293, 5227, 5290, 5215, 4695, 2136, 5004, 2079, 4872, 2025, 1999, 4687, 7210, 1925, 12852, 4461, 4408, 4243, 15273, 4256, 10544, 12809, 5468, 4069, 11944, 3878, 3939, 13826, 3857, 13992, 13849, 14589, 16718, 3666, 13306, 16231, 16075, 19150, 12804, 17662, 18618, 17338, 23980, 12158, 12118, 41935, 16626, 13582, 17334, 17186, 32424, 16897, 13016, 31620, 84785, 47876, 30855, 16086, 47891, 60010, 15709, 105630, 55618, 46063, 28985, 132722, 85978, 34462, 61921, 61469, 98741, 144060, 121798, 143866, 124310, 32534, 123864, 117598, 275658, 75605, 257360, 114441, 135234, 137956, 238629, 257451, 235530, 254129, 252500, 293037, 131844, 296359, 130196, 576404, 235119, 290348, 565860, 278829, 125490, 443629, 123996, 413263, 350806, 528878, 654653, 457617, 404577, 622617, 610992, 510950, 117030, 264395, 689854, 1119070, 671235, 258616, 663939, 255820, 578143, 253084, 110805, 962318, 249088, 639660, 471032, 632978, 107378, 242701, 555908, 1164043, 622961, 454463, 959564, 449941, 604776, 1931012, 940749, 936160, 1034749, 927115, 1041676, 918243, 1751793, 1834579, 905249, 890590, 1329399, 1007761, 3014698, 1972560, 880334, 2358631, 862253, 1773735, 854485, 2343606, 1642306, 1937410, 2285576, 1885663, 2849036, 3450284, 1236919, 5602894, 3054872, 2382584, 1654257, 1856481, 1847107, 6070801, 1833080, 4372938, 5343597, 6009329, 4968769, 5980116, 6175469, 5234544, 1906470, 5192159, 5859549, 3791939, 5177074, 9167102, 5089140, 5914788, 5049068, 6193360, 7161197, 8953084, 5576666, 9804705, 15479668, 6484011, 4894898, 4921173, 6410329, 7000178, 6959606, 9510931, 6220845, 10419040, 1615826, 8367364, 5342530, 8428361, 9267951, 15733435, 5265102, 15436508, 10225681, 5208488, 8098413, 19972759, 15162811, 11245931, 14171527, 20677825, 20605525, 19555210, 14846919, 7961739, 8391828, 14693858, 13159626, 19154762, 18675031, 29996960, 14445651, 35087063, 8166544, 29920233, 14341358, 13015240, 14246382, 14111923, 7568890, 28213972, 20559987, 38079851, 35800634, 34714902, 33615670, 27669651, 41949273, 31934314, 37278592, 41497004, 33945901, 13572263, 13529583, 36694288, 13445023, 40721359, 59414803, 51988611, 32163141, 36016855, 34636202, 45847407, 51246260, 39731174, 48974757, 11609096, 45156934, 32212927, 48388233, 32020611, 68645124, 41739096, 59859254, 49583402, 47534323, 86284290, 58680755, 111128157, 72837448, 72626325, 58117850, 30913270, 67022439, 107967669, 57453646, 30560982, 30474161, 63362847, 65154694, 70580513, 87038962, 30047352, 62477891, 64247247, 175502367, 140695500, 138973543, 61617314, 56807825, 174838837, 61112254, 120615817, 128788491, 137685707, 60451581, 114422833, 135237695, 67174483, 59805040, 66816219, 169723871, 134855982, 168825861, 130096671, 95191395, 65763995, 65591838, 65420580, 94199818, 165756300, 93711736, 219786991, 122149909, 164051865, 92750590, 163212725, 108807700, 102460630, 63594117, 101941842, 91345278, 63113557, 125414410, 106898793, 106631546,316919127,316130771,62173901,105575788,203681025,61714488,324216795,61411966,346812469,321844477,396569020,308457694,103275105,280083082,87163205,307635420,153036392,355200739,152305908,363925867,432842984,587648645,150865663,100595798,100359102,698359333,99889036,742400756,424771320,297619383,550053070,305454249,674913703,555669641,303347668,464520255,697024047,695432668,691051151,96937769,1320493795,969040287,940425950,288234988,826441680,934100215,286300525,1007176777,662585915,1641850951,329210441,671175786,1356484445,1646323897,1350521876,1068611381,663832506,1586503974,648150492,1606158539,395286109,521992693,655229925,995022651,921109262,517512069,1951172971,890189521,2024453776,1307420114,1304644275,642736134,509853328,1028031202,1293657797,1237600167,1947917667,1838182582,635907116,867934783,2525614752,3297039415,1822606623,1526514314,613404280,626747960,1517110735,2219954014,2573693081,3048598334,2811841042,1501692943,2730820103,2529424186,1241388795,1489582516,2977996441,609179629,834887166,3936377113,4402180384,2806166322,4219755279,1219221138,2789496027,912431838,821713404,4096868979,5912718456,3859193248,2701788555,2751358386,899981501,8445255552,4962751468,1190867158,5057914334,2599023766,2714249506,4081801741,4627622940,2579107875,4773702121,4990346776,6422015302,5414029916,4961938730,2549799831,2797852989,6361430252,6369357446,3912047822,2525880508,3313115872,2633075689,5050473671,6298191441,4860486451,4851468851,7714135917,4833533661,11838904881,2594282677,3617993670,8940613046,3604740946,9361483090,2456741443,3585042908,4754439474,13107805402,4737213244,9259911845,4720111391,10913263380,7769343669,3533552166,12653686592,18748243789,4669538769,12793516242,13795952741,2391286520,8639422181,3483519569,4620038358,3471232022,3465120698,4595679632,12591513354,7295329939,13554764756,17081555860,9581514843,2341381410,7536217628,8874750867,13858930276,7194530907,17432644005,19074688582,12331894522,18273609420,9417447808,15972808608,4462357868,12226852831,4447179780,14739060524,7060395585,8664519882,4417131268,4409682480,16477662471,12062458171,15911880607,17845082566,15858663615,13514565138];
ff1(n)={
	my(s=n^2);
	for(i=1,#u1,if(u1[i]>n,s=max(s,(u[i]\n)*n);break));forstep(k=s,9e99,n,if(nextprime(k)-k>n,return(k/n)))
};

ff2(n)={
	if(n < 399, error("Too small, use ff1"));
	my(t);
	for(i=1,#L,
		t = L[i]\n + 2;
		if (n * t < U[i], return(t-1))
	);
	"?"
};



\\ ??
searchUnknown(lim)={
	my(t,A,B,C);
	for(a=0,lim,
		if(a%100==0,print1(a" "));	\\ show progress
		A = a^2;
		if(a%3,
			forstep(c=a+a%3,lim,if(a%3==1,[1,2],[2,1]),
				C = c^2;
				forstep(b=1,lim,[1,2],
					if(c!=b & a!=b & (t=try(A,b^2,C)),print("\n**"t))
				)
			)
		,
			forstep(b=0,lim,3,
				if (a == b, next);
				B = b^2;
				forstep(c=a+3,lim,3,
					if(b!=c & (t=try(A,B,c^2)),print("\n**"t))
				)
			)
		)
	)
};
try(a,b,c)={
	my(d=(b+4*c-2*a)/3,e,f,g,h,i);
	if(a!=d&b!=d&c!=d&issquare(d)&issquare(e=(a+b+c)/3)&a!=e&b!=e&c!=e&d!=e&issquare(f=(4*a+b-2*c)/3)&a!=f&b!=f&c!=f&d!=f&e!=f&issquare(g=(2*a+2*b-c)/3)&a!=g&b!=g&c!=g&d!=g&e!=g&f!=g&issquare(h=(2*a-b+2*c)/3)&a!=h&b!=h&c!=h&d!=h&e!=h&f!=h&g!=h&issquare(i=(2*b+2*c-a)/3)&a!=i&b!=i&c!=i&d!=i&e!=i&f!=i&g!=i&h!=i&issquare(a)&issquare(b)&issquare(c)&a!=b&a!=c&b!=c
	,
		[a,b,c;d,e,f;g,h,i]
	)
};


\\ All through findCE related to http://mymathforum.com/viewtopic.php?f=40&t=19010

\\ Find the number of primes that the claim says can be found
limitkiller(p)={
	my(pr=1.,pp=primepi(p),ret);
	forprime(q=3,p,pr*=1-2/q);
	ret=floor(1-pp+p^2/2*pr);
	print("Claim: There are at least "ret" primes p with "p" < p < "p^2"\nwhich are not congruent to (x_2, x_3, ..., x_"pp") mod (3, 5, ..., "p").");
	ret
};

\\ maximum number of occurances of a number up to n in vector v
mx(v,n)={
	my(ct=vector(n));
	for(i=1,#v,ct[v[i]]++);
	vecmax(ct)
};

\\ "List primes"
LP(p)={
	my(P=primes(primepi(p^2)),pp=primepi(p));
	vector(#P-pp,i,P[i+pp])
};

\\ Find residue classes that leave few primes
heuristicSearch(p)={
	my(v=LP(p), remaining=vector(primepi(p)-1,i,prime(i+1)),best,bestAt,t,q,residue);
	while(#remaining,
		\\ Find the prime that removes the largest number beyond the expected number
		best=0;
		for(i=1,#remaining,
			t=mx(v%remaining[i],remaining[i])-ceil(#v/(remaining[i]-1));
			if (t >= best,
				best = t;
				bestAt = i
			)
		);
		q = remaining[bestAt];
		best=mx(v%q,q);
		\\print1("  Removing "best" with");
		
		\\ Find what residue class to remove
		for(i=1,q-1,
			if(sum(j=1,#v,v[j]%q==i)==best,
				t = i;
				break
			)
		);
		print("  "t" mod "q);
		
		\\ Actually remove those residue classes from v and the modulus from remaining
		v=select(n->n%q!=t,v);
\\print1("  "remaining" -> ");
		remaining=concat(vector(bestAt-1,i,remaining[i]),vector(#remaining-bestAt,i,remaining[i+bestAt]))
\\;print(remaining)
	);
	v
};

\\ Search for counterexamples
findCE(lim)={
	my(need,have);
	forprime(p=3,lim,
		need = limitkiller(p);
		have = heuristicSearch(p);
		if (#have < need,
			error("This looks like a counterexample mod "p": needed "need", have only ",#have)
		);
		print(#have" >= "need", ok")
	);
};


A066888(lim,sz)={
	my(v=vector(sz),n=1,t=1,s=0);
	forprime(p=2,lim,
		if(p>t,
			print(s" primes between T(",n-1,") = ",n*(n-1)/2," and "t);	\\ show progress
			if(s>0&s<=sz&!v[s],
				v[s]=n-1
			);
			n++;
			t=n*(n+1)/2;
			s=1
		,
			s++
		)
	);
	v
};


\\ A056757 etc.
builder(p=2,n=1,x=0)={	\\ called by user; also handles p <= 7
\\print("builder("p","n","x")");
	my(v=[],u,y,q=nextprime(p+1),f=if(p==7,builderSmall,builder));
	for(e=0,if(p==2,10,7),
		y=x+3*log(e+1)-e*log(p);
		u=f(q,n*p^e,y);
		if(#u, v=concat(v,u));
		if(y>0, v=concat(v,n*p^e))
	);
	if(p==2,vecsort(v),v)
};
builderSmall(p,n,x)={	\\ p > 7
\\print("builderSmall("p","n","x")");
\\if(p<5,error("builderSmall called with inappropriate argument"));
	if(x+3*log(2)<log(p), return([]));	\\ Not enough mojo to continue with more primes
	my(v=[],u,y,N,q=nextprime(p+2));
	v=builderSmall(q,n,x);
	for(e=1,5,
		y=x+3*log(e+1)-e*log(p);
		if(y>0,
			u=builderSmall(q,N=n*p^e,y);
			v=if(#u, concat(v,concat(u,N)), concat(v,N))
		,
			break
		)
	);
	v
};

doit()={
	\\ A125827, etc.
	my(V=vector(20,i,[]),total=vector(20),t,n=0);
	forbigprime(p=2,1e14,
		n++;
		t=1;
		for(i=1,20,
			t*=p;
			total[i] += t;
			if (total[i]%n == 0,
				V[i]=concat(V[i],[n]);
				print("Sequence #"i": ",V[i])
			)
		)
	)
};


\\ Related to A082794, etc.
intminus(a,b)={
	my(v);
	for(i=1,#V,
		v=V[i];
		if(v[2]<=a|v[1]>=b,next);
		if(a<=v[1]&b>=v[2],
			V[i]=0
		,
			if(a<=v[1],
				V[i]=[b,v[2]]
			,
				if(b>=v[2],
					V[i]=[v[1],a]
				,
					V[i]=[v[1],a];
					V=concat(V,[[b,v[2]]])
				)
			)
		)
	);
	V=vecsort(V,lex,8);
	if(#V & V[1]==0,V=vector(#V-1,i,V[i+1]));
	V=vecsort(V,1);
};
go(N)=local(V=[[1,10]]);for(n=1,99,oV=V;intminus(N/n,(N+1)/n);intminus(10*N/n,10*(N+1)/n);intminus(100*N/n,100*(N+1)/n);if(V!=oV,print1(n"n, ")));V*1.;


A037053(n)={
	my(N=10^(n+1),NN);
	forstep(a=N,9*N,N,
		if(ispseudoprime(a+1), return(a+1));
		if(ispseudoprime(a+3), return(a+3));
		if(ispseudoprime(a+7), return(a+7));
		if(ispseudoprime(a+9), return(a+9))
	);
	N *= 10;
	forstep(a=N,9*N,N,
		NN=1;
		while((NN*=10) < N,
			forstep(b=a+NN,a+9*NN,NN,
				if(ispseudoprime(b+1), return(b+1));
				if(ispseudoprime(b+3), return(b+3));
				if(ispseudoprime(b+7), return(b+7));
				if(ispseudoprime(b+9), return(b+9))
			)
		)
	);
	"?"
};


coverage(ff, n, conf=.05)={
	intnum(p=0,1,
		cover(ff, n, p, conf)
	)
};
cover(ff, n, p, conf=.05)={
	sum(k=0,n,if(inInterval(p,ff(n,k,conf)), p^k*(1-p)^(n-k)*binomial(n,k)))
};
inInterval(p,interval)={
	p>=interval[1] && p <= interval[2]
};


\\ See Li 2003.
L(X,n)={
	my(t,half=ceil(X/2));
	sum(k=ceil(sqrt(X)/2),sqrtint(X),
		t = n - k^2;
		half <= t && t < X
	)
};
addhelp(L, "L(X,n): Li's (2003) function L.");


roulette(n, m, c, p)={
	my(mx,tmp);
	if (n == 0,
		return (c >= p)
	);
	
	n--;
	mx = roulette(n, m, c - 1, p) * 19/37 + roulette(n, m, c + 1, p) * 18/37;
	for (b = 1, m,
		tmp = P(n, m, c - b, p);
		mx = max(max(mx, tmp * 19/37 + roulette(n, m, c + b, p) * 18/37), max(tmp * 25/37 + roulette(n, m, c + 2 * b, p) * 12/37, tmp * 36/37 + roulette(n, m, c + 35 * b, p) * 1/37))
	);
	mx	
};
addhelp(roulette, "roulette(n, m, c, p) be the probability of winning (a simplified version of) European roulette with an optimal strategy and n rounds, a max bet of m, c chips, and a target profit of p");


zander(product)={
  my(x,y,sqrty);
  for(n = 1, 10000,
    x = sqrtint(product * n) + 1;
    y = (x^2) % product; \\ or y = lift(Mod(x, product)^2);
    if (issquare(y, &sqrty), return(rt(x + sqrty, product)))
  );
};


findA141768(a,b,cur,writeToFile)={
	my(t,i=0);
	forstep(n=bitor(a,1),b,2,
		if(!isprime(n) && Bpsp(n) > cur,
			t = Bspsp(n);
			if(t > cur,
				print(t" "n);
				cur = t;
				if(writeToFile,
					write("b141768.txt", i++" "n)
				)
			)
		)
	)
};
addhelp(findA141768, "findA141768(a,b,cur,writeToFile): Find terms of Sloane's A141768 in [a, b].");


\\ Throwaway functions for testing a conjecture of Z-W Sun.
\\ http://listserv.nodak.edu/cgi-bin/wa.exe?A2=ind0812&L=nmbrthry&T=0&F=&S=&P=2140
A059389(n)={ \\ For n > 1
	my(t=trinv(n-2));
	fibonacci(t+2)+fibonacci(n-t*(t-1)/2)
};
addhelp(A059389, "A059389(n): Sums of two nonzero Fibonacci numbers. Sloane's A059389.");
sunconj(n,s=100000)={
	if(fib2=='fib2,
		print1("Initializing vector...");
		fib2=concat(2,vector(99999,n,A059389(n+1)));
		print("done.")
	);
	forstep(i=s,1,-1,
		if(isprime(n-fib2[i]),return(0))
	);
	1
};
testrange(a,b)={
	my(m=9);
	while(fib2[m+1]<a,m++);
	forstep(n=floor(a),b,3,
		if(fib2[m+1]<n,m++);
		if(sunconj(n,m),print(n);write("log.txt",n));
		if(sunconj(n+1,m),print(n+1);write("log.txt",n+1));
		if(sunconj(n+2,m),print(n+2);write("log.txt",n+2));
	);
	write("log.txt", "testrange: Range "floor(a)" to "floor(b)" finished.")
};


bad(p,d)={
	forstep(n=p,p*d,p,for(b=1,d,
		if(isprime(n+b)||isprime(n-b),return(0))
	));
	1
};
addhelp(bad, "bad(p,d): Is p d-nice? A problem proposed by davar55; see http://mersenneforum.org/showthread.php?p=152030");
first(d,start=2)={
	forprime(p=start,default(primelimit),
		if(bad(p,d),return(p))
	);
	-1
};
firstbig(d,mx=1e9,p=0)={
	p=max(p,default(primelimit))-1;
	while(p<mx,
		p=nextprime(p+1);
		if(bad(p,d),return(p))
	);
	-1
};


Pierce(a,b)={
	my(k=0);
	while(b>0,b=a%b;k++);
	k
};
Pmax(a)={
	my(k=-1);
	for(b=2,a-1,
		k=max(k,Pierce(a,b))
	);
	k
};


least(n,m)={
	(m-1)-(n-1)%m
};
addhelp(least, "least(n, m): Returns the least number >= n that is 0 mod m.");
p43()={
	my(s=0,kk);
	forstep(k=012,987,3,
		kk = k\10;
		if(bitand(kk,1),next);
		if(k%10 == kk%10 || k%10 == kk\10 || kk%10 == kk\10, next);
		for(n=01234+least(1234,17)+k*10^5,98765+k*10^5,
			if((n\10000)%5==0&&((n\1000)%1000)%7==0&&((n\100)%1000)%11==0&&((n\10)%1000)%13==0&&(n%1000)%17==0&&#vecsort([n%10,(n\10)%10,(n\100)%10,(n\1000)%10,(n\10000)%10,(n\100000)%10,(n\1000000)%10,n\10000000],,8)==8,
				print1(n" (k = "k")\t");
				s+=hundred(n);
			)
		)
	);
	s
};
hundred(n)={
	my(s=0,ss=0);
	forstep(k=10^9+n,98*10^8+n,10^8,
		if(pan10(k),s+=k;ss++)
	);
	print(ss);
	s
};

Sokol(n,limit=1e7)={
	if (isA129912(n-1), return(-1));
	if (isA129912(n+1), return(1));
	forprime(p=3,n,if(isA129912(n-p),return(-p)));
	forprime(p=3,limit,if(isA129912(n+p),return(p)));
	0
};
addhelp(Sokol, "Does the Sokol conjecture hold for n?  Returns p such that n + p is in A129912 and |p| is either 1 or an odd prime, or 0 if it cannot find such a p.");


walk(n)={
	local(vv,v,cur,nxt);
	n=floor(abs(n/2));
	vv=vector(n);
	v=1.;
	for(step=1,n,
		cur = vv[1];
		nxt = vv[2];
		vv[1] = vv[1] + v/2;
		v = v / 2;
		vv[2] = vv[2] + cur/4;
		v = v + cur/4;
		for(i=2,n-1,
			cur=nxt;
			nxt = vv[i+1];
			vv[i-1] = vv[i-1] + cur/4;
			vv[i+1] = vv[i+1] + cur/4;
			vv[i] = cur/2;
		);
		vv[n] = vv[n] + cur/4;
	);
	print(v);
	print(vv);
	sum(k=1,n,vv[k]*k)*2/sqrt(2*n)
};
addhelp(walk,"Simulates a 1-D random walk.");


polynum(a,b,n)={
	local(tmp);
	if(a>b,
		tmp = a;
		a = b;
		b = tmp;
	);
	
	if (a == 2,
		forstep(i=floor(n^(1/b)),0,-1,
			if(issquare(n-i^b),return(1))
		),
		forstep(i=floor(n^(1/b)),0,-1,
			tmp = ispower(n-i^b);
			if(tmp&&(tmp%a == 0),return(1))
		)
	);
	0
};
addhelp(polynum,"polynum(a,b,n): can n be expressed in the form x^a + y^b with x,y non-negative integers?");


LegendreC(start,stop,b)={
	my(n);
	for (i=ceil(start),stop,
		n=nextprime(i^2)-i^2;
		if (n>b,
			b=n;
			print(n"\t"i);
			if(!isprime(i^2 + n),print("**The above line does not meet isprime! **"))
		)
	)
};
addhelp(LegendreC, "LegendreC(start,stop,b): Related to the Legendre conjecture.");


\\ Smallest prime in each E-S class.
A005113=[2, 13, 37, 73, 1021, 2917, 15013, 49681, 532801, 1065601, 8524807, 68198461, 545587687, 1704961513, 23869461181, 288310406533];


\\ Gives the lower bound on E-S classes -- see A005113.
lbES(n)={
	if (n <= #A005113,
		A005113[n],
		((A005113[#A005113]+1)<<(n-#A005113))-1
	)
};


\\ Macro to do lots of calculation on Erdos-Selfridge prime classes.
doES()={
	lst4=writeBaseES(4,"4.txt");
	print("List 4 size: "#lst4);
	lst5=buildES(lst4,"5.txt");
	print("List 5 size: "#lst5);
	lst6=buildES(lst5,"6.txt");
	print("List 6 size: "#lst6);
	lst7=buildES(lst6,"7.txt");
	print("List 7 size: "#lst7);
	lst8=buildES(lst7,"8.txt");
	print("List 8 size: "#lst8);
	lst9=buildES(lst8,"9.txt");
	print("List 9 size: "#lst9);
	lst10=buildES(lst9,"10.txt");
	print("List 10 size: "#lst10);
	lst11=buildES(lst10,"11.txt");
	print("List 11 size: "#lst11);
	lst12=buildES(lst11,"12.txt");
	print("List 12 size: "#lst12);
};


\\ Build as many E-S primes of the appropriate class as possible from those one class smaller.
buildFancyES(lastlst,filename="b.txt")={
	local(limit,lst);
	limit = lastlst[#lastlst]*2-1;
	lst=[];
	for(i=1,#lastlst,
		forstep(p=2*lastlst[i]-1,limit,2*lastlst[i],
			if(isprime(p), lst=concat(lst,p))
		)
	);
	lst=vecsort(lst);
	for(i=1,#lst,
		write("c:\\"filename,i" "lst[i])
	);
	#lst
};


\\ Build as many E-S primes of the appropriate class as possible from those one class smaller.
buildES(lastlst,filename="b.txt")={
	local(limit,lst);
	limit = lastlst[#lastlst]*2-1;
	lst=[];
	for(i=1,#lastlst,
		forstep(p=2*lastlst[i]-1,limit,2*lastlst[i],
			if(isprime(p), lst=concat(lst,p))
		)
	);
	lst=vecsort(lst);
	for(i=1,#lst,
		write("c:\\"filename,lst[i])
	);
	lst
};


writeFancyES(class,filename="b.txt")={
	local(n);
	n=0;
	for(i=1,#ES,if(ES[i]==class,
		n=n+1;
		write("c:\\"filename,n" "prime(i));
		if(n >= 5000, return(n))
	));
	n
};


\\ This writes the precalculated (with buildES) E-S primes of the given class
\\ to a file and then returns the result as a vector.
writeBaseES(class,filename="p.txt")={
	for(i=1,#ES,if(ES[i]==class,
		write("c:\\"filename,prime(i))
	));
	readvec("c:\\"filename)
};


\\ Lazy evaluation of Erdos-Selfridge+ prime class.  This is a good method for finding small
\\ ES+ prime classes; large ones should be constructed out of lists of smaller ones using buildES.
ESClass(p)={
	local(f,m,class,start);
	if (ES == 'ES, ES = [1,1]);
	start=#ES+1;
	
	\\ Allocate space to hold all new entries.
	ES=concat(ES,vector(primepi(p)-#ES));
	for(n=start,#ES,
		f=factor(prime(n)+1);
		m=matsize(f)[1];
		class=1;
		for(i=2,m,
			if(f[i,1] > 3, 
				class=max(class,ES[primepi(f[i,1])]+1)
			)
		);
		ES[n]=class
	);
	if(isprime(p), return(ES[primepi(p)]));

	\\ At this point the E-S class is undefined, but it is at times
	\\ convenient to define this as the maximum of the E-S classes of
	\\ its prime factors.
	0
};


\\ Returns the ES+ number.
singleES(p)={
	local(f,m,class);
	if (p < #ES,
		if (p  < 13, return (1));
		return (ES[primepi(p)])
	);
	f = factor(p+1);
	m = matsize(f)[1];
	class = 1;
	for (i = 2, m,
		if (f[i, 1] > 3, class = max(class, singleES(f[i,1]) + 1))
	);
	class
};


\\ Estimates log(dag(t)) based on R. P. Stanley (2006)
estlogdag(t)={
	0.5544947694708433699946880820231264067578639900861156+t*(t-1)*log(2)/2+lngamma(t+1)-
	t*0.3974857210390778070092136788410687008793342476338963345239587
};


\\ Estimates the parameter alpha ~= 1.48807... in the estimation of dag.
estalpha(n)={
	(1.0 << (n-1)) * n * dag(n-1) / dag(n)
};


\\ Estimates the parameter C ~= 1.74106... in the estimation of dag.
estC(n, alpha)={
	dag(n)/(n!) * 2^(-n*(n-1)/2) * alpha ^ n
};


\\ Rough number of correct decimal places in b, given true value a.
precdec(a, b)={
	if(a*b <= 0, return(0));
	if(a==b, return(precision(b)));
	-floor(log(abs(a-b)) / log(10))
};


Stanley(n)={
	local(a,a1,c,c1);
	a = estalpha(n);
	a1 = estalpha(n-1);
	c = estC(n, a);
	c1 = estC(floor(.9*n)-2, a);
	print("Parameters from R. P. Stanley, 'Acyclic orientation of graphs' (2006)");
	print("alpha: "a);
	print("C:     "c);
	print("Decimal places: "precdec(a, a1)" for alpha and "precdec(c, c1)" for C.");
};


\\ ***************************************************************************************************
\\ *                                         Momo stuff                                              *
\\ ***************************************************************************************************

	\\ For S = {1, 2, ..., a+b}, phi(a, b) is the number of coprimes of a in S, plus the number of coprimes to b in S.
\\ See http://mymathforum.com/viewtopic.php?f=40&t=12181
phi(a,b)={
	sum(n=1,b%a,gcd(a,n)==1)+sum(n=1,a,gcd(b,n)==1)+eulerphi(a)*(b\a+1)+eulerphi(b)
};

\\ sumphi(n) is ∑ φ(i,2n-i) with i in {1, 2, ..., n}.
\\ momo writes this as φφ(2n)
\\ See http://mymathforum.com/viewtopic.php?f=40&t=12181
sumphi(n)={
	sum(i=1,n,phi(i,2*n-i))
};

\\ Phi(n) calculates φ(1) + φ(2) + ... + φ(n)
Phi(n)={
	my(t,s=0);for(d=1,n,t=n\d;s+=moebius(d)*t*(t+1));s/2
};

nPhi(n)={
	n*Phi(n)
};

phisetest(v)={
	my(S=sum(i=1,#v,v[i]));
	#v*S-sum(i=1,#v,ceil(S/v[i]*(v[i]-eulerphi(v[i]))))
};

testPhi(lim,startAt=1)={
	my(diff=0,b,t);
	for(s=startAt,lim,
		for(a=1,s>>1,
			b=s-a;
			t=phisetest([a,b])-phiset([a,b]);
			if(abs(t)>abs(diff),
				diff=t;
				print(t," ",[a,b])
			)
		)
	);
	diff
};
testPhi1(lim,startAt=1)={
	my(diff=0,b,t);
	for(s=startAt,lim,
		for(a=1,s>>1,
			b=s-a;
			t=phisetest([a,b])-phiset([a,b]);
			if(t>diff,
				diff=t;
				print(t," ",[a,b])
			)
		)
	)
};

\\ phiset(v) calculates φ(v1, v2, ..., vk)
\\ phiset generalizes phi: phiset([a, b]) = phi(a,b)
\\ momo calls this the "uple totient"
\\ See http://mymathforum.com/viewtopic.php?f=40&t=12181&start=42
phiset(v)={
	my(s=sum(i=1,#v,v[i]),t);
	sum(i=1,#v,
		t=v[i];
		eulerphi(t) * (s\t) +
		sum(j=1,s%t,
			gcd(t,j)==1
		)
	)
};

partial(n,k)={
	sum(i=1,k,gcd(n,i)==1)
};

phiset1(v)={
	my(n=#v,S=sum(i=1,n,v[i]));
	n*(S-sum(i=1,n,
		if(isprime(v[i]),
			S\v[i]
		,
		if(S%v[i],
			S/v[i]*(v[i]-eulerphi(v[i]))
		,
			(S\v[i])*(v[i]-eulerphi(v[i]))+
			partial(v[i],S%v[i])
		))
	))
};

findphiprime(v)={
	v=vector(#v+1,i,if(i<=#v,v[i],1));
	for(i=1,9e9,
		v[#v]=i;
		if(isprime(phiset(v)),return(i))
	)
};

\\ Finds 1, 1, 4, 15, 64, 301, 1688, 10097, 70868, 519073, 4349551, 40628112, 387002087, 4105064694, 44541193708, ...
\\ See http://mymathforum.com/viewtopic.php?f=40&t=12307
phiprimes(lim)={
	my(v=[1]);
	for(i=2,lim,
		v=concat(v,findphiprime(v))
	);
	v
};

parallelPrime2(lim=25)={
	my(syze=100000,arr,a=1,lim2=997,minn,v=[1]);
	arr=List(vector(syze, i, i));

	for(ctr=1,lim,
		ptr=a;
		move=2;

		forprime(move=2,lim2,
			ptr += move;      \\eg 1+2, then 1+3, then 1+5,1+7,1+11,...
			if(ptr>syze, error(ptr" = ptr > syze = "syze,k5));
			arr[ptr]=0;      \\strike-out sieve elements
		);

		minn=a*10;    \\Feb 6 tweak
		for(l7=1,syze,
			elem=arr[l7];
			if(elem<minn && elem>a,minn=elem)  \\want lowest nonzero value NOT already setaside
		);  \\end FOR l7
		a=minn;
		\\write("c:\\Pari\\momo3.txt",minn," , ")
		v=concat(v,minn)

	);  \\end while
	v
};

