\\ Makes programs for sequences like A179643.
generate(v)={
	my(s);
	if(type(v) != "t_VEC", error("Input should be a vector"));
	if (#v == 1,
		s=Str("a(n)=",gPow("prime(n)", v[1]));
		return(gH(s))
	);
	v=vecsort(v,,4);
	if(#v>3,
		my(pVar=if(#v>6,vector(#v,i,Str("p",i)),["p","q","r","s","t","u"][1..#v]),tVar=vector(#v-1,i,Str("t",i)),t,m);
		t=gRoot(Str("lim\\",tVar[#v-1]), v[#v]);
		s=Str("forprime(",pVar[#v],"=2,",t,", if(",gOr(pVar[1..#v-1], pVar[#v]),", next); listput(v, ",tVar[#v-1],"*",gPow(pVar[#v],v[#v]),"))");
		forstep(i=#v-1,2,-1,
			m=factorback(primes(#v-i),v[i+1..#v]);
			t=gRoot(Str("lim\\(",m,"*",tVar[i-1],")"), v[i]);
			s=Str("forprime(",pVar[i],"=2,",t,", if(",gOr(pVar[1..i-1], pVar[i]),", next); ",tVar[i],"=",gPow(pVar[i],v[i]),"*",tVar[i-1],"; ",s,")");
		);
		m=factorback(primes(#v-1),v[2..#v]);
		t=gRoot(Str("lim\\",m), v[1]);
		s=Str("forprime(",pVar[1],"=2,",t,", ",tVar[1],"=",gPow(pVar[1],v[1]),"; ",s,")");
		s=Str("list(lim)=my(v=List()",concat(apply(x->Str(",",x), tVar)),"); ",s);
	);
	if (#v == 2,
		s="list(lim)=my(v=List(),t);forprime(p=2, ";
		s=if(v[1]==v[2],
			Str(s, gPow("lim", 1/(2*v[1])), ", ")
		,
			Str(s, gPow(Str("lim\\", 2^v[2]), 1/v[1]), ", ")
		);
		s = Str(s, "t=", gPow("p", v[1]), ";");	\\ Not needed for v[1] == 1, but that's not too important -- only happens for [1,1].
		s=if(v[1]==v[2],
			Str(s, "forprime(q=p+1, ", gPow("lim\\t", 1/v[2]), ", ")
		,
			Str(s, "forprime(q=2, ", gPow("lim\\t", 1/v[2]), ", if(p==q, next);")
		);
		s = Str(s, "listput(v,t*", gPow("q", v[2]), ")))");
	);
	if(#v == 3,
		s="list(lim)=my(v=List(),t1,t2);forprime(p=2, ";	\\ We don't always need the temp variables... but it's a pain to remove them automatically.  For now you can remove them by hand or just ignore the issue.
		s=Str(s,if(v[1]==v[2],
			if(v[1]==v[3],
				gPow("lim", 1/(3*v[1]))
			,
				gPow(Str("lim\\",2^v[3]), 1/(2*v[1]))
			)
		,
			gPow(Str("lim\\", 2^v[2]*3^v[3]), 1/v[1])
		),", t1=", gPow("p", v[1]), ";");
		s=Str(if(v[1]==v[2],
			if(v[1]==v[3],
				Str(s, "forprime(q=p+1, ", gPow("lim\\t1", 1/(2*v[2])), ", ")
			,
				Str(s, "forprime(q=p+1, ", gPow("lim\\t1", 1/v[2]), ", ")
			)
		,
			Str(s, "forprime(q=2, ", gPow("lim\\t1", 1/v[2]), ", if(p==q, next);")
		),"t2=t1*", gPow("q", v[2]), ";");

		s=Str(if(v[2]==v[3],
			if(v[1]==v[3],
				Str(s, "forprime(r=q+1, ", gPow("lim\\t2", 1/v[3]), ", ")
			,
				Str(s, "forprime(r=q+1, ", gPow("lim\\t2", 1/v[3]), ", if(p==r,next);")
			)
		,
			Str(s, "forprime(r=2, ", gPow("lim\\t2", 1/v[3]), ", if(p==r||q==r, next);")
		));
		s = Str(s, "listput(v,t2*", gPow("r", v[3]), "))))");
	);
	s = Str(s, "; Set(v)");
	return(gH(s))
};
gH(s)={	\\ Helper function for printing or returning values.
	my(v,lm=625);
	print("(PARI) "s" \\\\ ~~~~");
	eval(s);
	while(#(v=list(lm))<30, lm*=2);
	if(v!=Set(v),
		if(v==vecsort(v),
			warning("repeats")
		,
			warning("not in order")
		)
	);
	v
};

gPow(s, n)={ \\ n-th power of s.
	if(n == 1, return(s));
	if(type(n)=="t_INT" || type(n)=="t_REAL",
		return(Str(if(needsParens(s),Str("("s")"),s), "^", n))
	);
	if (type(n) != "t_FRAC", error("bad type"));
	if (n == 1/2,
		return(Str("sqrt",if(hasParens(s),s,Str("("s")"))))
	);
	if(needsParens(s),
		Str("("s")^("n")")
	,
		Str(s"^("n")")
	)
};

\\ n-th root of s.
gRoot(s, n)=
{
	if(n == 1, return(s));
	if(n == 2, return(Str("sqrtint(", s, ")")));
	Str("sqrtnint(", s, ", ", n, ")");
}

gOr(v:vec, x)=
{
	if (#v==1, return(Str(x,"==",v[1])));
	Str(x,"==",v[1],concat(apply(s->Str(" || ",x,"==",s), v[2..#v])));
}

needsParens(s)={
	my(v=Vec(s),start);
	if(alpha(v[1]),
		while(start < #v && varchar(v[start+1]), start++);
		if (start == #v, return(0))
	);
	v[start+1] != "(" || v[#v] != ")"
};	\\ variable, function(foo), (expression) are ok; a+b is not since (e.g.) a+b^2 is not the same as (a+b)^2.  Doesn't yet handle pure numbers.
hasParens(s)={
	my(v=Vec(s));
	v[1] == "(" && v[#v] == ")"
};	\\ variable bad, (expression) good since you can't do sqrtvariable, but sqrt(expression) is fine.
alpha(s)={
	(s <= "z" && s >= "a") || (s <= "Z" && s >= "A")
};	\\ Is s a letter?
varchar(s)={
	alpha(s) || s == "_" || (s >= "0" && s <= "9")
};	\\ Is s a valid variable character?  Letters, 0-9, and _.
