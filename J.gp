time(x)={
	if(x < 1,
		return(Strprintf("%d.2 ms", x*1000))
	);
	if(x < 120,
		return(Strprintf("%d.2 s", x))
	);
	if(x < 7200,
		return(Strprintf("%d.1 min", x/60))
	);
	if(x < 48*3600,
		return(Strprintf("%d.1 hr", x/3600))
	);
	if(x < 770*24*3600,
		return(Strprintf("%d.1 d", x/24/3600))
	);
	if(x < 36525*24*3600,
		return(Strprintf("%d.1 yr", x/365.25/24/3600))
	);
	Strprintf("%g yr", x/365.25/24/3600)
};

fac(v)={
   if(#v<4 || isprime(#v),
      \\ Length 0: []
      \\ Length 1: [1]   (identity)
      \\ Length 2: [1, 2] (prime)
      \\ Length 2: [2, 2] (prime)
      \\ Length 3: [1, 2, 3] (prime)
      \\ Length 3: [3, 3, 3] (prime)
      return(if(#v==1,[],[v]))
   );
   
   \\ Check for factors of the form [d, ..., d] with d elements.
   my(g=gcd(gcd(v),#v));
   fordiv(g,d,
      if(d>1 && d<#v,
         for(i=2,#v,if(v[i]!=v[(i-1)\d+1],next(2)));
         return(concat(fac(vector(d,i,d)),fac(vector(#v/d,i,v[d*i-1]/d))))
      )
   );
   
   \\ Check if [1, 2] divides v.
   if(#v%2==0 && sum(i=1,#v,v[i]%2)<=#v/2,
      my(need=[],has=[]);
      forstep(i=#v,1,-1,
         if(#need && need[1]==v[i],
            need=if(#need==1,[],need[2..#need]);
            has=concat(v[i],has);
            next
         );
         if(v[i]%2,has=0;break);   \\ Not a multiple of [1, 2]
         need=concat(v[i]/2,need);
      );
      if(#has==#v/2,
         my(f=fac(has));   \\ Factor, if possible, otherwise return partial factorization
         return(concat([[1,2]],if(f==-1,has,f)))
      );
   );
   
   \\ Check if [1, 2, 3] divides v.
   if(#v%3==0 && v[#v]%3==0 && sum(i=1,#v,v[i]%3==0)>=#v/3 && sum(i=1,#v,gcd(v[i],6)>1)>=#v*2/3 && sum(i=1,#v,v[i]%2==0)>=#v/3,
      \\ ...
      if(#v/gpf(#v)>2, return(-1));
   ,
      if(#v/gpf(#v)>3, return(-1));
   );
   
   
   
   [v]
};
Jmult(u,v)={
   my(V=vector(#u*#v,i,u[(i-1)\#v+1]*v[(i-1)%#v+1]));
   V=lift(select(n->n,V))
   \\;vecsort(V)
};
isJ(v)={
   sum(i=1,#v,v[i]^3)==sum(i=1,#v,v[i])^2
};
genmod(n,m)={
   my(v=vector(n,i,1));
   while(v[#v]<m-1,
      v[1]++;
      if(v[1]>=m,
         for(i=2,n,
            if(v[i]<m-1,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               break
            )
         )
      );
      if(sum(i=1,n,v[i]^3)==Mod(sum(i=1,n,v[i]),m)^2,
         print(Vecrev(v))
      );
   )
};
gen(n,verbose=1,mx=biggest(n))={
   my(v=vector(n,i,1),V=List(),tried);
   while(v[#v]<mx,
      v[1]++;
      if(v[1]>mx,
         for(i=2,n,
            if(v[i]<mx,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               break
            )
         )
      );
      tried++;
      if(sum(i=1,n,v[i]^3)==sum(i=1,n,v[i])^2,
         my(u=Vecrev(v),f);
         if(verbose,
            f=fac(u);
            if(f==-1,
               print(u" unknown")
            ,#f>1,
               print(u" = "f)
            ,
               print(u" prime")
            )
         );
         listput(V,u)
      );
   );
   print("Tried "tried" multisets");
   Vec(V)
};


biggest(n,start=n^(4/3)\1)={
   if(n<5,return(n));
   forstep(k=start,n,-1,
      my(P=k^3+(n-1)*'x^3-(k+(n-1)*'x)^2);
      if(!issquarefree(P) || polsturm(P,1),
         return(k)
      )
   )
};
addhelp(biggest, "biggest(n): Gives an upper bound on the size of the largest element of a J-multiset of cardinality n.");


count(n,mx=biggest(n))={
   if(n<2,return(1));
   my(v=vectorsmall(n,i,1),tried,found,s=n,s3=n,best,bestAt,S3,n2=n^2);
   S3=vector(mx,i,i^3-(i-1)^3);
   print("Effort required: "binomial(mx+#v-1,#v)" (about "time(binomial(mx+#v-1,#v)/5e5)")");
   while(v[#v]<mx,
      v[1]++;
      if(v[1]>mx,
         for(i=2,n,   \\ carry
            if(v[i]<mx,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               s-=(i-1)*(mx-v[i])-1;
               /*
               if(s>n2,
                  print(v"?");
                  s-=i*v[1]-1;
                  for(j=i+1,n,
                     s-=v[j]++;
                     if(j*v[j]+s<=n2,
                        print1(s" ");
                        for(k=1,j-1,v[k]=v[j]);
                        s+=v[j]*j;
                        print(v" "s);
                        s3=sum(k=1,n,v[k]^3);
                        break
                     )
                  );
                  break(2)
               ,
                  s3-=(i-1)*(mx^3-v[i]^3)+(v[i]-1)^3-v[i]^3;
               );*/
                  s3-=(i-1)*(mx^3-v[i]^3)+(v[i]-1)^3-v[i]^3;
               break
            )
         ) \\ end fixup
      ,
         s++;
         s3+=S3[v[1]]
      );
      
      \\print(sum(i=1,n,v[i]^3)"=="s3"\t"sum(i=1,n,v[i])"=="s);
      \\print(v);
      
      tried++;
      if(s3==s^2,
      \\if(sum(i=1,n,v[i]^3)==sum(i=1,n,v[i])^2,
         found++;
         if(vecmax(v)>best,
            best=vecmax(v);
            bestAt=v
         )
      );
   );
   print("Tried "tried" multisets");
   print("Best: "vecsort(Vec(bestAt)));
   found
};
count1(n,mx=biggest(n),pre=[])={
   if(n<2,return(1));
   my(v=vectorsmall(n-#pre,i,1),tried,found,s=n+sum(i=1,#pre,pre[i]),s3=n+sum(i=1,#pre,pre[i]^3),best,bestAt,S3,n2=n^2);
   S3=vector(mx,i,i^3-(i-1)^3);
   print("Effort required: "binomial(mx+#v-1,#v)" (about "time(binomial(mx+#v-1,#v)/5e5)")");
   while(v[#v]<mx,
      v[1]++;
      if(v[1]>mx,
         for(i=2,#v,   \\ carry
            if(v[i]<mx,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               s-=(i-1)*(mx-v[i])-1;
                  s3-=(i-1)*(mx^3-v[i]^3)+(v[i]-1)^3-v[i]^3;
               break
            )
         ) \\ end fixup
      ,
         s++;
         s3+=S3[v[1]]
      );
      
      tried++;
      if(s3==s^2,
         found++;
         if(vecmax(v)>best,
            best=vecmax(v);
            bestAt=v
         )
      );
   );
   print("Tried "tried" multisets");
   if(best,print("Best: "vecsort(concat(pre,Vec(bestAt)))));
   found
};

find1(n,mx=biggest(n),pre=[])={
   \\error("broken -- sometimes gives solutions which don't work");
   if(n<2,return(1));
   my(v=vectorsmall(n-#pre,i,1),s=n+sum(i=1,#pre,pre[i]),s3=n+sum(i=1,#pre,pre[i]^3),S3,n2=n^2);
   S3=vector(mx,i,i^3-(i-1)^3);
   print("Effort required: "binomial(mx+#v-1,#v)" (up to about "time(binomial(mx+#v-1,#v)/5e5)")");
   while(v[#v]<mx,
      v[1]++;
      if(v[1]>mx,
         for(i=2,#v,   \\ carry
            if(v[i]<mx,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               s-=(i-1)*(mx-v[i])-1;
                  s3-=(i-1)*(mx^3-v[i]^3)+(v[i]-1)^3-v[i]^3;
               break
            )
         ) \\ end fixup
      ,
         s++;
         s3+=S3[v[1]]
      );
      
      if(s3==s^2,
         return(vecsort(concat(pre,Vec(v))))
      );
   )
};


findSpecial(n,v1=gen(n-1),v2=gen(n))={
   my(s1=vector(#v1,i,sum(j=1,#v1[i],v1[i][j])),s2=vector(#v2,i,sum(j=1,#v2[i],v2[i][j])));
   for(i=1,#v1,
      for(j=1,#v2,
         my(s=s2[j]-s1[i]);
         if(!vecsearch(v1[i],s), next);
         my(needed=multisetminus(v2[j],v1[i]),left);
         if(#needed!=2,next);
         left=multisetminus(v1[i],v2[j]);
         if(#left!=1,next);
         print(v1[i]" -> "v2[j])
      )
   )
};


multisetminus(u:vec,v:vec)={   \\ Assumes u and v are sorted
   my(V=List(),i=1,j=1);
   while(i<=#u && j<=#v,
      if(u[i]==v[j],
         i++;
         j++
      ,
         if(u[i]<v[j],
            listput(V,u[i]);
            i++
         ,
            j++
         )
      )
   );
   while(i<=#u,
      listput(V,u[i]);
      i++
   );
   Vec(V)
};


generate(n)={
   my(v=vector(n,i,1),s=n,s3=n);
   while(v[#v]<n,
      v[1]++;
      if(v[1]>n,
         for(i=2,n,
            if(v[i]<n,
               v[i]++;
               for(j=1,i-1,v[j]=v[i]);
               s-=(i-1)*(n-v[i])-1;
               s3-=(i-1)*(n^3-v[i]^3)+(v[i]-1)^3-v[i]^3;
               break
            )
         )
      ,
         s++;
         s3+=v[1]^3-(v[1]-1)^3
      );
      \\print(v"\t"sum(i=1,n,v[i]^3)"=="s3"\t"sum(i=1,n,v[i])"=="s);
      \\if(sum(i=1,n,v[i]^3)==sum(i=1,n,v[i])^2,
      if(s3==s^2,
         my(u=Vecrev(v));
         if(u[1]==n,next);   \\ (n+n+...+n)^3 = (n+n+...+n)^2
         if(u==vector(n,i,i),next);   \\ Ignore this infinite family
         
         \\ n can be assumed to be >= 4 at this point
         my(ok=0);   \\ ok = 1 means that u is not of the form vecsort(apply(numdiv, divisors(N))) for any N.
         if(u[1]>1, ok=1);   \\ 1 has 1 divisor
         if(u[2]!=2, ok=1);   \\ Must be divisible by at least one prime
         if(u[3]!=2, ok=1); \\ Must be divisible by at least two distinct primes, or else the only possibility is [1, 2, ..., n] which is excluded above.
         if(v[1]!=n, ok=1);   \\ The last element must have d(N) = n divisors...
         if(v[2]==n, ok=1);   \\ Can have only one element with the maximal number of divisors
         if(!vecsearch(u,4), ok=1);   \\ Must be divisible by the product of the primes above
         print(u,if(ok,"","*"))
      );
   )
};


check2(n,mx=n)={
   for(b=1,mx,for(a=1,b,
      if(n-2+a^3+b^3==(n+a+b-2)^2,
         print("{1, ..., 1, "a", "b"}")
      )
   ))
};
