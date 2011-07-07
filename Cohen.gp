/* The following file contains a large number of useful GP scripts, 
concatenated from scripts written by Henri Cohen in 1998-2005. */

/*********************************************************************/
/*  (Z/nZ)^* and Dirichlet characters 				     */
/*  A znstruct is a new technical structure useful for computations: */
/*  znstruct[1] contains the group structure in the usual way        */
/*  [h,cyc,gen], znstruct[2] and znstruct[3] contains technical data */
/*  essential for computing discrete logs, znstruct[4] contains      */
/*  factor(n), znstruct[5] contains n itself                         */
/*********************************************************************/

/* (0) Computing lambda(n) */

lambda(n)=
{
 my(fa,lipr,liex,lfa,res,p,ex);
 fa=factor(n);lipr=fa[,1];liex=fa[,2];lfa=#lipr;res=1;
 for (i=1,lfa,
  p=lipr[i]; ex=liex[i];
  if (p==2, 
   if (ex<=2,res=ex,res=2^(ex-2)),
   res=lcm(res,p^(ex-1)*(p-1))
  )
 );
 return (res);
}

/* (I) Computing znstructs */

/* (I.1) Case of prime powers pa = p^a */

znstar_prp(pa)=
{
 if (pa%2, return (znstar(pa)));
 if (pa==2, return ([1,[],[]]));
 if (pa==4, return ([2,[2],[Mod(-1,4)]]));
 return ([pa/2,[pa/4,2],[Mod(5,pa),Mod(-1,pa)]]);
}

/* (I.2) The main program: znstar_2(n) returns the znstruct of (Z/nZ)^* */

znstar_2(n)=
{
 my(fa,lipr,liex,lfa,lipa,vpa,cln,w,wgen,wpro,gen,uvd,d,uinv,wgenmods,lgsnf);
 if (n==0,error("n = 0"));
 n=abs(n); fa=factor(n);
 if (n<=2,return([[1,[],[]],[],[;],fa,n]));
 lipr=fa[,1]; liex=fa[,2]; lfa=#lipr; lipa=vector(lfa,i,lipr[i]^liex[i]);
 vpa=vector (lfa,i,znstar_prp(lipa[i]));
 cln=prod (i=1,lfa,vpa[i][1]);
 w=[]; for (i=1,lfa,w=concat(w,vpa[i][2]));
 lgsnf=#w;
 wgen=[];
 for (i=1,lfa,
  wpro=vpa[i][3];
  for (j=1,#wpro,
   gen=chinese(wpro[j],Mod(1,n/lipa[i]));
   wgen=concat(wgen,gen)
  )
 );
 uvd=matsnf(matdiagonal(w),1);
 u=uvd[1]; d=uvd[3]; d=vector(lgsnf,i,d[i,i]); uinv=u^(-1);
 wgenmods=vector(lgsnf);
 for (i=1,lgsnf,
  wgenmods[i]=prod(j=1,lgsnf,wgen[j]^uinv[j,i])
 );
 return ([[cln,d,wgenmods],vpa,u,fa,n]);
}

/* given zns1 and zns2 corr. to n1 and n2 with n2\mid n1, find map from
   zns1 to zns2 */
znstar_map(zns1,zns2)=
{
 my(n1,n2,clgp1,clgp2,cyc1,cyc2,gen1,gen2,lg1,lg2,m);
 n1=zns1[5]; n2=zns2[5]; if (n1%n2,error("n2 nmid n1 in znstar_map"));
 clgp1=zns1[1]; cyc1=clgp1[2]; gen1=clgp1[3]; lg1=#cyc1;
 clgp2=zns2[1]; cyc2=clgp2[2]; gen2=clgp2[3]; lg2=#cyc2;
 m=matrix(lg1,0);
 for (j=1,lg2,
  m=concat(m,Mat(znlog_2(zns1,lift(gen2[j]))))
 );
 return (m);
}

/*------------------------------------------------------------------*/

/* (II) Compute discrete logs in (Z/nZ)^* */

/* (II.1) Case of odd prime power; x is an int or an intmod a multiple of p */

znlog_prpodd(zprpstruct,x)=
{
 my(p,prpa,ex,gen,u,b,xk,qk);
 gen=zprpstruct[3][1];
 prpa=component(gen,1);
 p=prpa/(prpa-zprpstruct[1]);
 ex=valuation(prpa,p);
 xk=znlog(Mod(lift(x),p),Mod(lift(gen),p));
 if (ex==1,return([xk]~));
 u=lift(Mod(p/lift(gen^(p-1)-1),p)); 
 b=1/Mod(lift(x),prpa);
 qk=p-1;
 for (k=2,ex,
  qk*=p;
  xk=(xk-(p-1)*u*((lift(gen^xk*b)-1)/p))%qk
 );
 return ([xk]~);
}

/* (II.2) Case of power of 2. */

znlog_prp2(zprpstruct,x)=
{
 my(prpa,ex,gen,b,xk,qk,eps);
 prpa=zprpstruct[1];
 b=1/Mod(lift(x),2*prpa);
 if (prpa==1,return([]~));
 if (lift(x)%4==1,eps=0,eps=1);
 if (prpa==2,return([eps]~));
 prpa*=2;
 b*=(-1)^eps;
 ex=valuation(prpa,2);
 if (lift(b)%8==1,xk=0,xk=1);
 qk=2; gen=Mod(5,prpa);
 for (k=4,ex,
  qk*=2;
  xk=(xk-(lift(gen^xk*b)-1)/4)%qk
 );
 return ([xk,eps]~);
}

/* (II.3) Case of general prime power. */

znlog_prp(zprpstruct,x)=
{
 if (zprpstruct[1]==1,return([]~));
 if (component(zprpstruct[3][1],1)%2,
  return(znlog_prpodd(zprpstruct,x)),
  return(znlog_prp2(zprpstruct,x))
 );
 return (1);
}

/* (II.4) The main program: discrete log of x on znstruct. */

znlog_2(znstruct,x)=
{
 my(res1,res2,lfa,u,d,vlog,lfb);
 res1=znstruct[1]; res2=znstruct[2]; u=znstruct[3]; lfa=#res2;
 vlog=[]~;for(i=1,lfa,vlog=concat(vlog,znlog_prp(res2[i],x)));
 d=res1[2];
 lfb=#vlog; if (#d!=lfb,error("lengths"));
 vlog=u*vlog;
 vlog=vectorv(lfb,i,vlog[i]%d[i]);
 return(vlog);
}

/*---------------------------------------------------------------------*/

/* (III) Dirichlet characters */

/* given a znstruct [[cl, [cyc], [gen]], aux1, aux2, aux3] modulo m and    */
/* lcyc=#cyc, a character chi is given by one of the following structures, */
/* where the first component is a codeword:                                */
/* [0,[chi(gen[1]),...,chi(gen[lcyc])]] with chi(gen[i]) polmods,          */
/* [1,[chi(gen[1]),...,chi(gen[lcyc])]] with chi(gen[i]) complex,          */
/* [2,[...w[i]=v[i]*cyc[1]/cyc[i]...]] with chi(gen[i])=z_{d_1}^w[i],      */
/* [8,[chi(1),...,chi(m)]] with chi(i) polmods,                            */
/* [9,[chi(1),...,chi(m)]] with chi(i) complex,                            */
/* [64,[chi(gen[1]),...,chi(gen[lcyc])],[chi(1),..,chi(m)]] with polmods   */
/* [65,[chi(gen[1]),...,chi(gen[lcyc])],[chi(1),..,chi(m)]] with complex   */

/* (III.01) Product of two chars of the same code and modulus. */

prod_char_s(char1, char2)=
{
 my(cod,lcod,res,ch1,ch2);
 cod=char1[1]; 
 if (char2[1]!=cod,error("not the same code in prod_char_s"));
 lcod=#char1;
 res=vector(lcod); res[1]=cod;
 for (j=2,lcod,
  ch1=char1[j]; ch2=char2[j];
  res[j]=vector(length(ch1),i,ch1[i]*ch2[i])
 );
 return (res);
}

/* (III.01) Product of two chars of the same code different modulus. */

prod_char(zns1,char1,zns2,char2)=
{
 my(n1,n2,n3,zns3);
 if (char2[1]!=char1[1],error("not the same code in prod_char"));
 n1=zns1[5]; n2=zns2[5];
 if (n1==n2,return (prod_char_s(char1,char2)));
 n3=lcm(n1,n2); zns3=znstar_2(n3);
 return (prod_char_s(char_map(zns3,zns1,char1),char_map(zns3,zns2,char2)));
}

/* (III.02) Power of a char. */

pow_char(chara, n)=
{
 my(lcod,res,ch,code);
 lcod=#chara;
 res=vector(lcod); code=chara[1]; res[1]=code;
 if (code==2||code==3,
  ch=chara[2]; res[2]=vector(length(ch),i,n*ch[i]),
  for (j=2,lcod,
   ch=chara[j];
   res[j]=vector(length(ch),i,ch[i]^n)
  )
 );
 return (res);
}

/* (III.021) Conjugate (or inverse) of a character */

inv_char(chara)=return(pow_char(chara,-1));

/* (III.025) Given chara=[chi(gen[1]),...,chi(gen[lcyc])], transform it into */
/* character with code, assuming values already polmods or complex          */

char_transform(znstruct,chara,code)=
{
 my(res); 
 if (code<8,return([code,chara]));
 res=gentoper_char(znstruct,chara);
 if (code<64,return([code,res]),return([code,chara,res]));
 return (1);
}

/* given a character on zns2, lift it to a character on zns1, assuming
  n2 \mid n1 */

char_map(zns1,zns2,chara)=
{
 my(clgp1,cyc1,gen1,lg1,n1,n2,code,v,w);
 clgp1=zns1[1]; cyc1=clgp1[2]; gen1=clgp1[3]; lg1=#cyc1;
 n1=zns1[5]; n2=zns2[5]; code=chara[1];
 if (code<=3 || code >=64,
  v=vector(lg1,i,eval_char(zns2,chara,lift(gen1[i])))
 );
 if (code<=3,return([code,v]));
 if (code<=9,w1=chi[2],w1=chi[3]);
 w=vector(n1); for (i=1,n1,w[i]=w1[(i-1)%n2+1]);
 if (code<=9,return ([code,w]),return ([code,v,w]));
}

/* (III.03) List all characters of (Z/nZ)^*, primitive or not */
/* output characters using code as above.                     */

all_charscom(znstruct,code=0)=
{
 my(clgp,cyc,gen,z,lchar,chi,ct);
 if (type(znstruct)=="t_INT",znstruct=znstar_2(znstruct));
 if (znstruct[5]<=2,
  if (code%2,res=1,res=Mod(1,x-1));
  return (char_transform(znstruct,[res],code))
 );
 clgp=znstruct[1]; cyc=clgp[2]; gen=clgp[3]; lcyc=length(cyc);
 if (code%2,z=exp(2*I*Pi/cyc[1]),z=Mod(x,polcyclo(cyc[1])));
 lchar=vector(clgp[1]); ct=0;
 vfor=vector(lcyc,i,[0,cyc[i]-1]);
 forvec (v=vfor,
  chi=char_transform(znstruct,vector(lcyc,i,z^(v[i]*cyc[1]/cyc[i])),code);
  ct++; lchar[ct]=chi
 );
 return (lchar);
}

all_chars(znstruct)=all_charscom(znstruct,0);
all_chars_c(znstruct)=all_charscom(znstruct,1);

/* (III.1) Evaluate a character on an int or intmod x */

eval_char(znstruct,chara,x)=
{
 my(vlog,lfb,n,code);
 n=znstruct[5]; if (n==1,return(1));
 if (gcd(n,lift(x))>1,return(0));
 code=chara[1];
 if (code<8,
  chara=chara[2]; vlog=znlog_2(znstruct,x); lfb=#vlog;
  if (code<=1,
   return (prod(i=1,lfb,chara[i]^vlog[i])),
   clgp=znstruct[1]; cyc=clgp[2];
   return (TAB[(sum(i=1,lfb,vlog[i]*chara[i])-1)%cyc[1]+1])
  ),
  if (code<64,chara=chara[2],chara=chara[3]);
  return (chara[lift(x)%n])
 );
}

/* (III.15) Evaluate a vector of characters on an int or intmod x */

eval_vecchar(znstruct,vecchi,x)=
{
 my(vlog,lfb,n,code,lv,chara);
 lv=#vecchi; if(lv==0,return([]));
 n=znstruct[5]; if (n==1,return(vector(lv,i,1)));
 if (gcd(n,lift(x))>1,return(vector(lv)));
 chara=vecchi[1]; code=chara[1];
 if (code<8,
  vlog=znlog_2(znstruct,x); lfb=#vlog;
  if (code<=1,
   return (vector(lv,j,prod(i=1,lfb,vecchi[j][2][i]^vlog[i]))),
   clgp=znstruct[1]; cyc=clgp[2];
   return (vector(lv,j,TAB[(sum(i=1,lfb,vlog[i]*vecchi[j][2][i])-1)%cyc[1]+1]))
  ),
  if (code<64,
   return (vector(lv,j,vecchi[j][2][lift(x)%n])),
   return (vector(lv,j,vecchi[j][3][lift(x)%n]))
  )
 );
}

/* (III.2) Is q equal to 1. */

isone_char(q,n)=
{
 my(tq,p);
 tq=type(q);
 if (tq=="t_REAL" || tq=="t_COMPLEX", return (abs(q-1)<n^2*10^(3-precision(1.))));
 if (tq=="t_PADIC", 
  p = component(q,1);
  if (valuation (q,p), return(0));
  return (q+O(p^(padicprec(q,p)-1))==1)
 );
 return (q==1);
}

/* (III.3) Is character char even ? */

iseven_char(znstruct,chara)=
{
 my(n);
 n=znstruct[5]; 
 return (isone_char(eval_char(znstruct,chara,Mod(n-1,n)),n));
}

all_charseven(znstruct)=
{
 my(li,lli,v,c);
 li=all_chars(znstruct); lli=#li;
 v=vector(lli); c=0;
 for (i=1,lli,if(iseven_char(znstruct,li[i]),c++;v[c]=li[i]));
 vector(c,i,v[i]);
}

/* (III.4) Is character char defined modulo d dividing n ? */

isdefinedmod_char(znstruct,chara,d)=
{
 my(n,b);
 n=znstruct[5]; 
 if (n%d, error("d does not divide n"));
 forstep (a=d+1,n,d,
  b=a; while (gcd(b,n)>1, b+=d);
  if (!isone_char(eval_char(znstruct,chara,b),n),return (0))
 );
 return (1);
}

/* (III.4bis) Is character trivial ? */

istrivial_char(znstruct,chara)=
{
 return (isdefinedmod_char(znstruct,chara,1));
}

isreal_char(znstruct,chara)=
{
 return (istrivial_char(znstruct,pow_char(chara,2)));
}

/* (III.5) Compute the conductor of char. */

conductor_char(znstruct,chara)=
{
 my(fa,lipr,liex,lfa,n,m);
 fa=znstruct[4]; n=znstruct[5]; lipr=fa[,1]; liex=fa[,2]; lfa=#lipr; m=n;
 for (i=1,lfa,
  for (j=1,liex[i],
   if (isdefinedmod_char(znstruct,chara,m/lipr[i]), m/=lipr[i])
  )
 );
 return (m);
}

/* (III.52) Compute the primitive character equivalent to char. Only for
complex-valued (not checked). Returns [znstar_2(f),chi_f] */

getprim_char(znstruct,chara)=
{
 my(ff,n,zns,clgp,cyc,gg,res);
 n=znstruct[5]; ff=conductor_char(znstruct,chara);
 if (ff==n,return ([znstruct,chara]));
 zns=znstar_2(ff);
 clgp=zns[1]; cyc=clgp[2]; gen=clgp[3]; lcyc=length(cyc);
 res=vector(lcyc);
 for (i=1,lcyc,
  gg=lift(gen[i]);
  while (gcd(gg,n)>1,gg+=ff);
  res[i]=eval_char(znstruct,chara,gg)
 );
 return ([zns,res]);
}

/* (III.55) Compute the order of char. */

order_char(znstruct,chara)=
{
 my(clgp,h,fa,p);
 clgp=znstruct[1]; h=clgp[1]; fa=factor(h)[,1];
 for (i=1,length(fa),
  p=fa[i];
  while (h%p==0,
   if (istrivial_char(znstruct,pow_char(chara,h/p)),
    h/=p,
    break()
   )
  )
 );
 return (h);
}

/* (III.6) Is char a primitive character ? */

isprimitive_char(znstruct,chara)=
{
 my(fa,lipr,lfa,n);
 fa=znstruct[4]; n=znstruct[5]; lipr=fa[,1]; lfa=#lipr;
 for (i=1,lfa,
  if (isdefinedmod_char(znstruct,chara,n/lipr[i]),return(0))
 );
 return (1);
}

numprimchar(n)=
{
 my(fa,lipr,liex,P,p);
 fa=factor(n); lipr=fa[,1]; liex=fa[,2]; P=n;
 for (i=1,#lipr,
  p=lipr[i];
  if (liex[i]==1,P*=(1-2/p),P*=(1-1/p)^2)
 );
 return (P);
}

/* (III.65) List all primitive characters of (Z/nZ)^* */

allprim_charscom(znstruct,code=0,half=0)=
{
 my(n,clgp,cyc,gen,z,lchar,chi,sd,ct,w,ik,ctreal,llchar);
 if (type(znstruct)=="t_INT",znstruct=znstar_2(znstruct));
 clgp=znstruct[1]; cyc=clgp[2]; gen=clgp[3]; lcyc=length(cyc);
 n=znstruct[5];
 if (n==1,return([Mod(1,x-1)]));
 if (n%4==2,return([]));
 if ((code%2)==0,z=Mod(x,polcyclo(cyc[1])),z=exp(2*I*Pi/cyc[1]));
 if (code>1,
  TAB=vector(cyc[1]); TAB[1]=z; for(j=2,cyc[1],TAB[j]=z*TAB[j-1])
 );
 lchar=vector(clgp[1]); ct=0; ctreal=0;
 vfor=vector(lcyc,i,[0,cyc[i]-1]);
 forvec (v=vfor,
  if (half,
   w=vector(lcyc,i,if(v[i],cyc[i]-v[i],0));
   ik=0;for(i=1,lcyc,if(w[i]!=v[i],ik=i;break()));
  );
  if ((half&&((ik==0)||v[ik]<w[ik]))|| (half==0),
   if (code<=1,
    chi=char_transform(znstruct,vector(lcyc,i,z^(v[i]*cyc[1]/cyc[i])),code),
    chi=[code,vector(lcyc,i,v[i]*cyc[1]/cyc[i])]
   );
   if (isprimitive_char(znstruct,chi), 
    ct++; lchar[ct]=chi; if(half&&(ik==0),ctreal++)
   )
  )
 );
 lchar=vector(ct,i,lchar[i]);
 if (half,llchar=2*ct-ctreal,llchar=ct);
 if (llchar != numprimchar(n),
  print("lchar = ",llchar,", numchar = ",numprimchar(n));
  error("primitive characters not correct")
 );
 return (lchar);
}

allprim_chars(znstruct)=allprim_charscom(znstruct,0,0);
allprim_chars_c(znstruct)=allprim_charscom(znstruct,1,0);

/* list of all wild characters of conductor p^v */

allwild_chars(p,v)=
{
 my(zns,li,w);
 if (p==2,return([]));
 zns=znstar_2(p^v);
 li=allprim_chars(zns);
 w=[];
 for (i=1,#li,
  if (order_char(zns,li[i])==p^(v-1),
   if (!iseven_char(zns,li[i]),error("something wrong in allwildchars"));
   w=concat(w,[li[i]])
  )
 );
 return (w);
}

/* (III.7) Given a character modulo m as a vector of length m representing */
/* its values on 1..m, transform it into a char */

pertogen_char(znstruct,per)=
{
 my(gen,chara);
 gen=znstruct[1][3];
 if (#per!=component(gen[1],1),error("period not equal to n"));
 chara=vector(#gen,i,per[lift(gen[i])]);
 return (chara);
}

/* (III.8) Same in the other direction. */

gentoper_char(znstruct,chara)=
{
 my(n,clgp,cyc,gen,lcyc,res,vfor,ind);
 n=znstruct[5]; clgp=znstruct[1]; cyc=clgp[2]; gen=clgp[3]; lcyc=#cyc;
 res=vector(n);
 vfor=vector(lcyc,i,[0,cyc[i]-1]);
 forvec (v=vfor,
  ind=lift(prod(i=1,lcyc,gen[i]^v[i]));
  res[ind]=prod(i=1,lcyc,chara[i]^v[i])
 );
 return (res);
}

/* (III.9) Transform the kronecker character (D/n) into a char */

krotogen_char(znstruct,D)=
{
 my(gen,chara);
 gen=znstruct[1][3];
 if (abs(D)!=component(gen[1],1),error("period not equal to |D|"));
 chara=vector(#gen,i,kronecker(D,lift(gen[i])));
 return (chara);
}

/* temporary */
ll(znstruct, chara, m, p,k) = my(s, ap); s=0;for(a=1,m-1,if(a%p,ap=a+O(p^20);s+=eval_char(znstruct,chara,a)*log(ap/teichmuller(ap))^k/k!));s/m;
doll(m, p,k) = my(zns); zns=znstar_2(m);allp=allprim_chars(zns);valuation(vector(length(allp),i,ll(zns,allp[i],m,p,k)),p);

ispowk(m)=my(fa);fa=factor(m)[,1]; if (length(fa)==1,return(fa[1]),return(0));
gcdinf(a,binf)=my(fa);fa=factor(gcd(binf,a))[,1];prod(i=1,length(fa),fa[i]^valuation(a,fa[i]));

sninit(ff,n)=
{
 my(zns,li,lli,p,vord,orc);
 if (ff==4 && n%2,return(0));
 zns=znstar_2(ff);
 li=allprim_chars(zns); lli=#li;
 p=ispowk(ff);
 if (p==0 || ff%2==0,
  vord=vector(lli,i,li[i]),
  vord=[];
  for (i=1,lli,
   orc=order_char(zns,li[i]);
   if ((p-1!=gcd(p-1,n)*gcd(p-1,orc)) || (orc%p==0 && n%p==0),
    vord=concat(vord,[li[i]])
   )
  )
 );
 return([zns,vord]);
}

sncomp(ff,n)=
{
 my(zz,zns,vord,w,chara,lv);
 zz=sninit(ff,n); 
 if (zz==0,return([]));
 zns=zz[1]; vord=zz[2]; lv=#vord;
 w=vector(lv);
 for (i=1,lv,
  chara=vord[i];
  w[i]=content(lift(sum(r=1,ff-1,eval_char(zns,chara,r)*r^n)))
 );
 w;
}

/*********************************************************************/
/*  p-adic gamma functions and related functions                     */
/*  padic precision should be given, but if not defaults to          */
/*  DFPADP 							     */
/*********************************************************************/

global(DFPADP);
global(TAB);

DFPADP=20;

/* (IV) Auxilliary p-adic functions */

/* (IV.1) recognize a p-adic number */

bestappr_p(x,k,p)=
{
 my(m,res);
 m=[1,0;truncate(x+O(p^k)),p^k];
 res=(m*qflll(m,1))[,1];
 if (res[1]<0,res=-res);
 return (res[2]/res[1]);
}

/* (IV.2) linear dependence of p-adics */

lindep_p(v,k,p)=
{
 my(lv,m,res,w,wmin,wval);
 lv=#v;
 if (v[lv]==0,return(vectorv(lv,i,i==lv)));
 m=matid(lv);
 w=vector(lv-1,i,-v[i]/v[lv]); wval=vector(lv-1,i,valuation(w[i],p));
 wmin=p^vecmin(wval);
 for (i=1,lv-1,m[lv,i]=truncate(w[i]/wmin+O(p^k)));
 m[lv,lv]=p^k;
 res=(m*qflll(m,1))[,1]; res[lv]*=wmin;
 res=res/content(res);
 if (res[1]<0,res=-res);
 return (res);
}

/* (IV.3) algebraic dependence of p-adics */

algdep_p(a,n,k,p)=
{
 my(v,res);
 v=vector(n+1,i,a^(n+1-i));
 res=lindep_p(v,k,p);
 res=sum(j=0,n,res[j+1]*('x)^(n-j));
 return (res/x^valuation(res,x));
}

/* (IV.4) Compute sump(a=1,qp,1/a^n) */

Har_p(n,p,padp=DFPADP)=
{
 my(qp);
 if (p==2,qp=4,qp=p);
 return (sum (a=1,qp-1,if(a%p,1/a^n,0),O(p^padp)));
}

/* (IV.5) Computes B_{-k,p}, so k>0 here */

bernneg(k,p,padp=DFPADP)=
{
 my(qp,lim);
 if (p==2,qp=4,qp=p);
 lim=ceil((padp+1)/valuation(qp,p));
 return (sum(j=0,lim\2,binomial(k+2*j-1,k-1)*qp^(2*j-1)*bernfrac(2*j)*Har_p(k+2*j,p,padp),(k/2)*Har_p(k+1,p,padp)+O(p^padp)));
}

/* (IV.6) Computes p-adic Euler's constant */

Euler_p(p,padp=DFPADP)=
{
 my(qp,lim);
 if (p==2,qp=4,qp=p);
 lim=ceil((padp+1)/valuation(qp,p));
 return (-log(-(qp-1)!+O(p^padp))/qp+sum(j=1,lim\2,qp^(2*j-1)*bernfrac(2*j)/(2*j)*Har_p(2*j,p,padp),Har_p(1,p,padp)/2+O(p^padp)));
}

/* (IV.7) Compute p^kk!u(pk+r) */

ukpkrdwork(p,k,r,padp=DFPADP)=
{
 my(s,pr);
 s=O(p^padp); pr=1;
 for (m=0,k,
  s+=pr/(p*m+r)!;
  pr*=p*(k-m)
 );
 s;
}

/* (IV.8) Compute u(k). Not used, but why not ? */

ukdwork(p,k,padp=DFPADP)=
{
 my(s,pr,lim);
 s=O(p^padp); pr=1; lim=k\p;
 for (j=0,lim,
  s+=1/(pr*((k-p*j)!));
  pr*=p*(j+1)
 );
 s;
}

/*-------------------------------------------------------------*/
 
/* (V) Morita's p-adic gamma function */


/* (V.1) Morita's Gamma_p(z), where z is a p-adic integer, which */
/* can have its own p-adic prec, of course given priority        */
/* Version 1: Rodriguez-Villegas, using u(pk+r)                  */

gamma_p(z,p,padp=DFPADP)=
{
 my(D,rs,r,x,pr,s,lim);
 D=padicprec(z,p);
 if (D>=2^30,z=z+O(p^padp); D=padp);
 rs=z%p;
 if (rs==0,r=0;x=z/p,r=p-rs;x=(z+r)/p);
 pr=1; s=0;
 lim=1;
 for (i=1,3,
  lim=ceil((D+r*(2+1/(p-1))/p^2+log(lim)/log(p))*p/(p-1));
 );
 for(k=0,lim,
  s+=pr*ukpkrdwork(p,k,r,D);
  pr*=(x+k)/(k+1)
 );
 s;
}

/* (V.2) Morita's Gamma_p(z).                                      */
/* Version 2: using power series. Slower, but similar to psi, only */
/* for comparison                                                  */

gamma_pslow(s,p,padp=DFPADP)=
{
 my(res,u,r,fa);
 s=s-1;
 r=s%p; fa=(-1)^r*prod(i=1,r,s+1-i);s=s-r;
 res=-(Euler_p(p,padp)*s+sum(k=1,ceil(padp/valuation(s,p)),bernneg(2*k,p,padp)/((2*k)*(2*k+1))*s^(2*k+1),O(p^padp)));
 res=-exp(res);
 if (p==2,
  u=s*(s-2)/8;
  if (u%2==0,res=-res)
 );
 return (res*fa);
}


/* (V.3) Morita's psi function */

psi_p(s,p,padp=DFPADP)=
{
 my(r,fa);
 s=s-1;
 r=s%p;  fa=sum(i=1,r,1/(s+1-i),O(p^padp)); s=s-r;
 return(fa-(Euler_p(p,padp)+sum(k=1,ceil(padp/valuation(s,p)),bernneg(2*k,p,padp)/(2*k)*s^(2*k),O(p^padp))));
}

/* (V.4) Morita's psi' function */

psi1_p(s,p,padp=DFADP)=
{
 my(r,fa);
 s=s-1;
 r=s%p;  fa=-sum(i=1,r,1/(s+1-i)^2,O(p^padp)); s=s-r;
 return(fa-sum(k=1,ceil(padp/valuation(s,p)),bernneg(2*k,p,padp)*s^(2*k-1),O(p^padp)));
}

/* (V.5) Morita's psi^(m) function */

psim_p(s,p,m,padp=DFADP)=
{
 my(r,fa);
 if (m<0,error("negative m"));
 if (m==0,return(psi_p(s,p,padp)));
 if (m==1,return(psi1_p(s,p,padp)));
 s=s-1;
 r=s%p;  fa=(-1)^m*sum(i=1,r,1/(s+1-i)^(m+1),O(p^padp)); s=s-r;
 return(fa-sum(k=ceil(m/2),ceil(padp/valuation(s,p)),bernneg(2*k,p,padp)*s^(2*k-m)*prod(i=1,m-1,2*k-i),O(p^padp)));
}

/*********************************************************************/
/*  p-adic L-functions and related functions                         */
/*********************************************************************/

/* (VI) p-adic Hurwitz zeta function and related functions */

/* (VI.1) p-adic Hurwitz zeta function, function of x. */

hurwitz_p(s,x,p,padp=DFPADP)=
{
 my(S,vx);
 vx=valuation(x,p);
 if (vx>=0,error("x not of strictly negative valuation"));
 S=sum(j=0,(1+ceil((padp+1)/(-vx)))\2,binomial(1-s,2*j)*x^(-2*j)*bernfrac(2*j),(s-1)/(2*x)+O(p^padp));
 return (exp((1-s)*log(x+O(p^padp)))*S/(s-1));
}

/* (VI.1bis) same at s=1, minus 1/(s-1) if s=1 */

hurwitz0_p(s,x,p,padp=DFPADP)=
{
 if (s==1,return(-psiL_p(x,p,padp)),return(hurwitz_p(s,x,p,padp)));
}

/* (VI.2) p-adic Hurwitz zeta function, homogeneous, function of a, m. */

hurwitzhom_p(s,a,m,p,padp=DFPADP)=
{
 my(S);
 if (m%p, error("m not divisible by p"));
 S=hurwitz_p(s,a/m,p,padp);
 return (exp((1-s)*log(m+O(p^padp)))*S/m);
}

/* (VI.2bis) same at s=1, minus 1/(m(s-1)) if s=1 */

hurwitzhom0_p(s,a,m,p,padp=DFPADP)=
{
 if (s==1,
  return(-(log(m+O(p^(padp+1)))+psiL_p(a/m+O(p^padp),p,padp))/m),
  return(hurwitzhom_p(s,a,m,p,padp))
 );
}

/* (VI.3) Diamond's p-adic log gamma function LG */

LG_p(x,p,padp=DFPADP)=
{
 my(vx,S);
 vx=valuation(x,p);
 if (vx>=0,error("x not of strictly negative valuation"));
 S=sum(j=1,1+ceil((padp+1)/(-2*vx)),x^(1-2*j)*bernfrac(2*j)/(2*j*(2*j-1)),O(p^padp));
 return ((x-1/2)*log(x*(1+O(p^padp)))-x+S);
}

/* (VI.4) Diamond's p-adic psi function */

psiL_p(x,p,padp=DFPADP)=
{
 my(vx,S);
 vx=valuation(x,p);
 if (vx>=0,error("x not of strictly negative valuation"));
 S=sum(j=1,1+ceil((padp+1)/(-2*vx)),x^(-2*j)*bernfrac(2*j)/(2*j),O(p^padp));
 return (log(x*(1+O(p^padp)))-1/(2*x)-S);
}

/* (VI.5) Diamond's p-adic psi' function */

psiL1_p(x,p,padp=DFPADP)=
{
 my(vx,S);
 vx=valuation(x,p);
 if (vx>=0,error("x not of strictly negative valuation"));
 S=sum(j=0,1+ceil((padp+1)/(-2*vx)),x^(-2*j-1)*bernfrac(2*j),O(p^padp));
 return (1/(2*x^2)+S);
}

/* (VI.6) Diamond's p-adic psi^(m) function */

psiLm_p(x,p,m,padp=DFPADP)=
{
 my(vx,S);
 if (m<0,error("negative m"));
 if (m==0,return (psiL_p(x,p,padp)));
 if (m==1,return (psiL1_p(x,p,padp)));
 vx=valuation(x,p);
 if (vx>=0,error("x not of strictly negative valuation"));
 S=sum(j=0,1+ceil((padp+1)/(-2*vx)),x^(-2*j-m)*bernfrac(2*j)*prod(i=1,m-1,2*j+i),O(p^padp));
 return ((-1)^(m-1)*(S+m!/(2*x^(m+1))));
}

/* (VI.7) special Teichmuller as in book */

teichspec_p(x,p,padp=DFPADP)=
{
 my(v);
 v=valuation(x,p);
 x=(x/p^v)*(1+O(p^padp));
 return (p^v*teichmuller(x));
}

/* (VI.8) special <x> */

teichangle_p(x,p,padp=DFPADP)=
{
 my(v);
 v=valuation(x,p);
 x=(x/p^v)*(1+O(p^padp));
 return (x/teichmuller(x));
}

/*----------------------------------------------------------------*/

/* (VII) Kubota-Leopoldt p-adic L functions */

/* (VII.1) KL p-adic zeta function */

zeta_p(s,p,padp=DFPADP)=
{
 if (p>2,
  return (2*sum(a=1,(p-1)/2,hurwitzhom_p(s,a,p,p,padp))),
  return (2*hurwitzhom_p(s,1,4,2,padp))
 );
}

/* (VII.2) KL p-adic L function of general character */

L_p(znstruct,chara,s,p,padp=DFPADP)=
{
 my(S);
 if (!iseven_char(znstruct,chara),return(0));
 if (istrivial_char(znstruct,chara) && s==1,
  error("padic L function of trivial character at s = 1")
 );
 n=znstruct[5];
 if (p==2,m=lcm(n,4),m=lcm(n,p));
 S=O(p^padp);
 for(a=1,floor(m-1)/2,
  if(gcd(a,p*n)==1,
   S+=2*eval_char(znstruct,chara,Mod(a,n))*hurwitzhom0_p(s,a,m,p,padp)
  )
 );
 if (m%2==0 && gcd(m/2,p*n)==1,
  S+=eval_char(znstruct,chara,Mod(m/2,n))*hurwitzhom0_p(s,m/2,m,p,padp)
 );
 return (S);
}

/* (VII.3) KL p-adic L function for quadratic character mod D>0 */

Lquad_p(D,s,p,padp=DFPADP)=
{
 my(znstruct,chara);
 if (D==1,return(zeta_p(s,p,padp)));
 znstruct=znstar_2(D);
 chara=krotogen_char(znstruct,D);
 return (L_p(znstruct,chara,s,p,padp));
}

/* (VII.4) KL p-adic L function for Teichmuller^k */

Lteich_p(k,s,p,padp=DFPADP)=
{
 my(qp,gen,chara);
 if (p==2,qp=4,qp=p);
 znstruct=znstar_2(qp);
 gen=znstruct[1][3];
 chara=vector(#gen,i,teichmuller(lift(gen[i])*(1+O(p^padp)))^k);
 return (L_p(znstruct,chara,s,p,padp));
}

/*********************************************************************/
/*  L(chi,1-k) for k>=1 in Z                                         */
/*********************************************************************/

/* Bernoulli polynomial B_n(x) */

B(n,x)=sum(j=0,n,binomial(n,j)*bernfrac(j)*x^(n-j));

/* chi-Bernoulli polynomial B_n(chara,x) */

Bchipol(znstruct,chara,n,x)=
{
 my(m);
 m=znstruct[5];
 m^(n-1)*sum(r=0,m-1,eval_char(znstruct,chara,r)*B(n,(r+x)/m));
}

/* chi-Bernoulli number B_n(chara) */

Bchi(znstruct,chara,n)=
{
 e=1-iseven_char(znstruct,chara);
 if ((e+n)%2 && n>1, return (0));
 return (Bchipol(znstruct,chara,n,0));
}

/* L(chi,1-n) */

Lchineg(znstruct,chara,n)=
{
 if (n<=0, error("n must be >= 1 in Lchineg"));
 if (znstruct[5]==1,
  if (n==1,return (-1/2),return (-bernfrac(n)/n)),
  return (-Bchi(znstruct,chara,n)/n)
 );
}

/* L(leg{D}{.},1-n) */

LchinegD(D,n)=
{
 my(zns,chara);
 if (!isfundamental(D), error("D must be fundamental in LchinegD"));
 zns=znstar_2(abs(D));chara=krotogen_char(zns,D);
 Lchineg(zns,chara,n);
}

/* compute vector of all L(chi,1-n) for primitive chi mod m */

all_lchinegs(m,n)=
{
 my(znstruct,allpc);
 znstruct=znstar_2(m);
 allpc=allprim_chars(znstruct);
 return (vector(length(allpc),i,Lchineg(znstruct,allpc[i],n)));
}

/* Naive method 1 */

gausssum(zns,chi)=
{
  my(m,q,S);
  m=zns[5]; if (m==1,return(1)); q=exp(2*I*Pi/m); S=0.;
  forstep (n=m-1,1,-1,S=eval_char(zns,chi,n)+q*S); 
  return (q*S);
}

jacobisum(zns,chi1,chi2)=
{
  my(m,q,S);
  m=zns[5]; if (m<=2,return(0)); 
  return (sum(n=2,m-1,eval_char(zns,chi1,n)*eval_char(zns,chi2,m+1-n)));
}

jac3spec(zns,chi)=
{
 o=order_char(zns,chi); if (o==1,error("jac3spec"));
 if (o==3,return(-eval_char(zns,chi,-1)*jacobisum(zns,chi,chi)));
 if (o==2,return(zns[5]*eval_char(zns,chi,-1)));
 return (jacobisum(zns,chi,chi)*jacobisum(zns,chi,prod_char_s(chi,chi)));
}


rootno(zns,chi)=
{
  my(m,res);
  m=zns[5]; res=gausssum(zns,chi)/sqrt(m);
  if (!iseven_char(zns,chi),res=res/I);
  return (res);
}

allrootno(zns)=
{
 my(li);
 li=allprim_chars_c(zns);
 return (vector(#li,i,rootno(zns,li[i])));
}

allroot1(N)=
{
 my(res);
 res=allrootno(znstar_2(N));
 return (vector(#res,i,if(abs(res[i]^N-1)<10^(-20),1,0)));
}

allrootnonotreal(zns)=
{
 my(li,lc);
 li=allprim_chars_c(zns);
 lc=[]; 
 for (i=1,#li,if(!isreal_char(zns,li[i]),lc=concat(lc,rootno(zns,li[i]))));
 return (lc);
}

allroot1notreal(N)=
{
 my(res);
 res=allrootnonotreal(znstar_2(N));
 return (vector(#res,i,if(abs(res[i]^N-1)<10^(-20),1,0)));
}
 
isall1(v)=for(i=1,#v,if(v[i]==0,return(0)));return(1);

clc(m,n)=content(lift(all_lchinegs(m,n)));

thchieven(zns,chi)=
{
 my(S,q,q1,eps);
 S=0.;q=exp(-Pi/zns[5]);q1=q;n=1;eps=10.^(-precision(1.));
 while(q1>eps,S+=eval_char(zns,chi,n)*q1;n++;q1=q^(n^2));
 return(S);
}

thveceven(zns,vecchi)=
{
 my(lv,S,q,q1,eps);
 lv=#vecchi;S=vector(lv,i,0.);
 q=exp(-Pi/zns[5]);
 q1=q;q2=q^2;q3=q;n=1;
 eps=10.^(-precision(1.));
/* while(q1>eps, */
  while(abs(q1)>eps,
  S+=eval_vecchar(zns,vecchi,n)*q1;
  n++;q3*=q2;q1*=q3
 );
 return(S);
}

thchiodd(zns,chi)=
{
 my(S,q,q1,eps);
 S=0.;q=exp(-Pi/zns[5]);q1=q;n=1;eps=10.^(-precision(1.));
 while(q1>eps,S+=n*eval_char(zns,chi,n)*q1;n++;q1=q^(n^2));
 return(S);
}

thvecodd(zns,vecchi)=
{
 my(lv,S,q,q1,eps);
 lv=#vecchi;S=vector(lv,i,0.);q=exp(-Pi/zns[5]);q1=q;q2=q^2;q3=q;n=1;
 eps=10.^(-precision(1.));
 while(q1>eps,
  S+=n*eval_vecchar(zns,vecchi,n)*q1;
  n++;q3*=q2;q1*=q3
 );
 return(S);
}

dotheven(N)=
{
 my(zns,li,lie,ct);
 zns=znstar_2(N);li=allprim_chars_c(zns);lie=vector(#li); ct=0;
 for(i=1,#li,
  if(iseven_char(zns,li[i]),ct++; lie[ct]=thchieven(zns,li[i]))
 );
 return(vector(ct,i,lie[i]));
}

dothvecevenold(N)=
{
 my(zns,li,lie,ct,cu,res,lif);
 zns=znstar_2(N);li=allprim_charscom(zns,1,1);lie=vector(#li); ct=0;
 for(i=1,#li,if(iseven_char(zns,li[i]),ct++; lie[ct]=li[i]));
 lie=vector(ct,i,lie[i]);
 res=thveceven(zns,lie);
 lif=vector(2*ct); cu=0;
 for(i=1,ct,
  if(isreal_char(zns,lie[i]),
   cu++;lif[cu]=res[i],
   cu+=2;lif[cu-1]=res[i];lif[cu]=conj(res[i])
  )
 );
 return(vector(cu,i,lif[i]));
}

dothveceven(N)=
{
 my(zns,li,lie,ct,cu,res,lif);
 if(isprime(N),return(dothvecevenprime(N)));
 zns=znstar_2(N);li=allprim_charscom(zns,3,1);lie=vector(#li); ct=0;
 for(i=1,#li,if(iseven_char(zns,li[i]),ct++; lie[ct]=li[i]));
 lie=vector(ct,i,lie[i]);
 res=thveceven(zns,lie);
 lif=vector(2*ct); cu=0;
 for(i=1,ct,
  if(isreal_char(zns,lie[i]),
   cu++;lif[cu]=res[i],
   cu+=2;lif[cu-1]=res[i];lif[cu]=conj(res[i])
  )
 );
 return(vector(cu,i,lif[i]));
}

dothvecevenprime(N)=
{
 my(g,ze,NS4,S,q,q1,q2,q3,n,eps,vpro,vlog,lif,NS2);
 g=znprimroot(N); ze=exp(4*I*Pi/(N-1)); NS4=(N-1)\4; NS2=(N-1)\2;
 S=vector(NS4,j,0.);q=exp(-Pi/N);q1=q;q2=q^2;q3=q;n=1;
 eps=10.^(-precision(1.));
 vpro=vector(NS2); vpro[1]=ze;
 for(j=2,NS2,vpro[j]=ze*vpro[j-1]);
 while(q1>eps,
  if (gcd(n,N)==1,
   vlog=lift(znlog(n,g));
   S+=vector(NS4,j,vpro[((vlog*j-1)%NS2)+1])*q1
  );
  n++;q3*=q2;q1*=q3
 );
 fl=(N%4)==1;
 lif=vector(2*NS4-fl);
 for(j=1,NS4-fl,lif[2*j-1]=S[j];lif[2*j]=conj(S[j]));
 if(fl,lif[2*NS4-1]=S[NS4]);
 return(lif);
}

dothodd(N)=
{
 my(zns,li,lie,ct);
 zns=znstar_2(N);li=allprim_chars_c(zns);lie=vector(#li); ct=0;
 for(i=1,#li,
  if(iseven_char(zns,li[i])==0,ct++; lie[ct]=thchiodd(zns,li[i]))
 );
 return(vector(ct,i,lie[i]));
}

dothvecoddold(N)=
{
 my(zns,li,lie,ct,cu,res,lif);
 zns=znstar_2(N);li=allprim_charscom(zns,1,1);lie=vector(#li); ct=0;
 for(i=1,#li,if(iseven_char(zns,li[i])==0,ct++; lie[ct]=li[i]));
 lie=vector(ct,i,lie[i]);
 res=thvecodd(zns,lie);
 lif=vector(2*ct); cu=0;
 for(i=1,ct,
  if(isreal_char(zns,lie[i]),
   cu++;lif[cu]=res[i],
   cu+=2;lif[cu-1]=res[i];lif[cu]=conj(res[i])
  )
 );
 return(vector(cu,i,lif[i]));
}

dothvecodd(N)=
{
 my(zns,li,lie,ct,cu,res,lif);
 if(isprime(N),return(dothvecoddprime(N)));
 zns=znstar_2(N);li=allprim_charscom(zns,3,1);lie=vector(#li); ct=0;
 for(i=1,#li,if(iseven_char(zns,li[i])==0,ct++; lie[ct]=li[i]));
 lie=vector(ct,i,lie[i]);
 res=thvecodd(zns,lie);
 lif=vector(2*ct); cu=0;
 for(i=1,ct,
  if(isreal_char(zns,lie[i]),
   cu++;lif[cu]=res[i],
   cu+=2;lif[cu-1]=res[i];lif[cu]=conj(res[i])
  )
 );
 return(vector(cu,i,lif[i]));
}

dothvecoddprime(N)=
{
 my(g,ze,NS4,S,q,q1,q2,q3,n,eps,vpro,vlog,lif,NM1);
 g=znprimroot(N); NM1=N-1; ze=exp(2*I*Pi/NM1); NS4=(N+1)\4; 
 S=vector(NS4,j,0.);q=exp(-Pi/N);q1=q;q2=q^2;q3=q;n=1;
 eps=10.^(-precision(1.));
 vpro=vector(NM1); vpro[1]=ze;
 for(j=2,N-1,vpro[j]=ze*vpro[j-1]);
 while(q1>eps,
  if (gcd(n,N)==1,
   vlog=lift(znlog(n,g));
   S+=n*vector(NS4,j,vpro[((vlog*(2*j-1)-1)%NM1)+1])*q1
  );
  n++;q3*=q2;q1*=q3
 );
 fl=(N%4)==3;
 lif=vector(2*NS4-fl);
 for(j=1,NS4-fl,lif[2*j-1]=S[j];lif[2*j]=conj(S[j]));
 if(fl,lif[2*NS4-1]=S[NS4]);
 return(lif);
}

dothall(N)=
{
 my(zns,li,lie);
 zns=znstar_2(N);li=allprim_chars_c(zns);lie=vector(#li);
 for(i=1,#li,
  if(iseven_char(zns,li[i]),
   lie[i]=thchieven(zns,li[i]),
   lie[i]=thchiodd(zns,li[i])
  )
 );
 return(lie);
}

myvm(v)=if(#v,return(vecmin(v)),return(1));

alleven(lim1,lim2)=
{
 my(res);
 for(p=lim1,lim2,
  if(p%4!=2,
   res=myvm(norm(dotheven(p)));
   if(res<10^(-11),print("\n",p,": ",res),print1("."))
  )
 );
}

allveceven(lim1,lim2)=
{
 my(res,tmp);
 for(p=lim1,lim2,
  if(p%4!=2,
   tmp=dothveceven(p);
   res=myvm(norm(tmp));
   if(res<10^(-11),print("\n",p,": ",res),print1("."));
   if(p%100==0,print1("(",p,")"))
  )
 );
}

/* done up to 10600 */

allodd(lim1,lim2)=
{
 my(res);
 for(p=lim1,lim2,
  if(p%4!=2,
   res=myvm(norm(dothodd(p)));
   if(res<10^(-11),print("\n",p,": ",res),print1("."))
  )
 );
}

allvecodd(lim1,lim2)=
{
 my(res);
 for(p=lim1,lim2,
  if(p%4!=2,
   res=myvm(norm(dothvecodd(p)));
   if(res<10^(-11),print("\n",p,": ",res),print1("."))
  )
 );
}

allall(lim1,lim2)=for(p=lim1,lim2,res=myvm(norm(dothall(p)));if(res<10^(-11),print("\n",p,": ",res)));

do300(LIM)=zns=znstar_2(300);li=allprim_chars(zns);chi=li[7];v=vector(300);for(i=0,19,for(j=0,1,for(k=0,1,res=(277^i*151^j*101^k)%300;v[res]=[i,j,k])));z=exp(4*I*Pi/5);w=vector(300,n,if(#v[n]==3,z^(v[n][1])*(-1)^(v[n][2]+v[n][3]),0));sum(n=1,LIM,w[(n-1)%300+1]*exp(-Pi*n^2/300));

putargbin(li)=
{
 for(i=1,#li,
  res=ceil(((arg(li[i])+Pi)/(2*Pi))*100); if (res==0,res=100);
  if (res<1||res>100,print(li[i]);error("arg out of bounds"));
  BIN[res]++
 );
}

putnormbin(li)=
{
 for(i=1,#li,
  res=ceil(((arg(li[i])+Pi)/(2*Pi))*100); if (res==0,res=100);
  if (res<1||res>100,print(li[i]);error("arg out of bounds"));
  BIN[res]++
 );
}

/*********************************************************************/
/*  Complex zeta and L function                                      */
/*********************************************************************/


/* (VIII.1) hurwitz zeta function zeta(s,x). Valid for all reasonable */
/* values of x and s, with s not equal to 1.                          */

hurwitz_c(s,x)=
{
 my(a,res,tes,in,sig,t,m,pr,lim,nn,in2,s1,s2);
 sig=real(x);
 if (sig>1.5,
  m=floor(sig-0.5);
  return (hurwitz_c(s,x-m)-sum(i=1,m,(x-i)^(-s)))
 );
 if (sig<=0,
  m=ceil(-sig+0.5);
  return (hurwitz_c(s,x+m)+sum(i=0,m-1,(x+i)^(-s)))
 );
 pr=precision(1.); sig=real(s); t=imag(s);
 default(realprecision,9);
 res=s-1.;if(abs(res)<0.1,res=-1,res=log(res));
 lim=(pr*log(10)-real((s-.5)*res)+(1.*sig)*log(2.*Pi))/2;
 lim=max(2,ceil(max(lim,abs(s*1.)/2)));
 nn=ceil(sqrt((lim+sig/2-.25)^2+(t*1.)^2/4)/Pi);
 default(realprecision,pr+5);
 a=(x+nn+0.)^(-s);
 res=sum(n=0,nn-1,(x+n)^(-s),a/2);
 in=x+nn; in2=1./(in*in);
 s1=2*s-1; s2=s*(s-1);
 tes=bernreal(2*lim);
 forstep (k=2*lim-2,2,-2,
  tes=bernreal(k)+in2*(k*k+s1*k+s2)*tes/((k+1)*(k+2))
 );
 tes=in*(1+in2*s2*tes/2);
 res+=tes*a/(s-1);
 res=precision(res,pr); default(realprecision,pr);
 return(res);
}

/* (VIII.2) Complex L function, vector form. Chivec is a vector of complex */
/* values, assumed to be the values from 1 to m of a periodic function. In */
/* addition, with zero sum, such as a nontrivial character, if s=1.        */
/* simple implementation, for small m                                      */

Lsimp_c(chivec,s)=
{
 my(m);
 m=length(chivec);
 if (s==1,
  return(-sum(r=1,m,chivec[r]*psi(r/m))/m),
  return(sum(r=1,m,chivec[r]*hurwitz_c(s,r/m))/m^s)
 );
}

/* (VIII.3) Complex L function of quadratic character, simple implementation */
/* for small |D|                                                             */

Lsimpquad_c(D,s)=
{
 my(v);
 v=vector(abs(D),i,kronecker(D,i));
 return (Lsimp_c(v,s));
}

/* (VIII.4) Complex L function, character form, using functional equation. */
/* a ecrire */

/* L_c(znstruct,chara,s)= */

/* (VIII.5) Complex L function of quadratic character, using functional */
/* equation                                                             */
/* a ecrire */

/***********************************************************************/
/* Sums of series of rational functions                                */
/***********************************************************************/

/* (IX.1) partial fraction decomposition of rational function F. The answer  */
/* is a vector of 4-component vectors [al,v,rr,sal], where al ranges through */
/* the distinct poles of F, and only one among 2 conjugates if F is real, v  */
/* is the order of the pole, rr=2 if F is real and al nonreal, and sal is    */
/* the power series giving the polar decomposition around al, so that if     */
/* sal=sum(k=0,v-1,ck*x^k)+O(x^v) then around al we have                     */
/*   F=sum(k=0,v-1,ck/(x-al)^(v-k))+O(1).                                    */

ratdec(F)=
{
 my(vx,D,N,d,fl,ropro,ro,vfl,al,ct,rr,Dal,sal);
 vx=variable(F);
 D=denominator(F); N=numerator(F);
 d=poldegree(D);
 fl=(imag(F)==0);
 ropro=polroots(D);
 ro=[];
 vfl=vector(d,i,1);
 for (j=1,d,
  if (vfl[j],
   al=ropro[j]; ct=1; vfl[j]=0;
   for (i=j+1,d,
    if (vfl[i],
     if (ropro[i]==al,
      ct++;vfl[i]=0,
      if (fl && ropro[i]==conj(al),
       vfl[i]=0
      )
     )
    )
   );
   rr=1;
   if (imag(al)==0,al=real(al),if (fl,rr=2));
   Dal=D\((vx-al)^ct);
   sal=subst(N/Dal,vx,vx+al)+O(vx^ct);
   ro=concat(ro,[[al,ct,rr,sal]]);
  )
 );
 ro;
}

/* (IX.2) Computes the integral from N to infinity of rational function F. */

{
intinfrat(F,N,flerr=0)=
 my(s,vx,pr,lim,G,ro,r);
 s=0.;
 if (poldegree(F)>=-1,error("infinite integral in intinfrat"));
 vx=variable(F);
 G=subst(F,vx,1/vx);
 pr=precision(1.);
 default(realprecision,18);
 ro=polroots(denominator(G));
 r=1;for(k=1,length(ro),r=min(r,norml2(ro[k])));
 r=1/r;
 if (N<r+1,
  if (flerr,return([r]),
  error("N not large enough in intinflograt"))
 );
 lim=ceil(pr*log(10)/log(N/sqrt(r)))+1;
 default(realprecision,pr+5);
 G=1.*G+O(vx^(lim+2));
 s=0.;
 forstep (k=lim,1,-1,s=(polcoeff(G,k+1)/k+s)/N);
 s=precision(s,pr);
 default(realprecision,pr);
 s;
}

/* (IX.3) F rational fraction. Computes the sum from in to infinity of F(n). */
/* By default in=0. Optional parameter la>=1 (default 5) is only used for    */
/* speeding up the computation.                                              */

sumrat(F,in=0,la=5)=
{
 my(vx,pr,lim,nn,res,sal);
 vx=variable(F);
 if (in,F=subst(F,vx,vx+in));
 if (poldegree(F)>=-1,error("infinite sum in sumrat"));
 pr=precision(1.);
 default(realprecision,9);
 la=max(la+0.,1);
 lim=ceil((pr*log(10)+2)/(2*(1+log(la))));
 nn=ceil((lim+.25)/Pi*la);
 default(realprecision,pr+5);
 res=intinfrat(F,nn,1);
 if (type(res)=="t_VEC",
  nn=max(nn,ceil(2*res[1]));
  default(realprecision,9);
  lim=ceil(nn*Pi/la+.75);
  default(realprecision,pr+5);
  res=intinfrat(F,nn);
 );
 res+=sum(ijkl=0,nn-1,subst(F,vx,ijkl),subst(F,vx,nn)/2.);
 sal=1.*subst(F,vx,vx+nn)+O(vx^(2*lim));
 bernreal(2*lim);
 for (k=1,lim,
  res-=bernreal(2*k)/(2*k)*polcoeff(sal,2*k-1)
 );
 res=precision(res,pr);
 default(realprecision,pr);
 res;
}

/* (IX.4) F rational fraction. Computes the sum from in to infinity of   */
/* (-1)^n*F(n). By default in=0. Optional parameter la is here set to 3, */
/* but only controls the speed.                                          */

sumratalt(F,in=0,la=3)=
{
 my(vx,res);
 if (type(in)!="t_INT",error("not an integer starting value in sumratalt"));
 vx=variable(F);
 if (in,F=subst(F,vx,vx+in));
 res=sumrat(2*subst(F,vx,2*vx)-F,0,la);
 if (poldegree(F)==-1,
  res-=pollead(numerator(F))/pollead(denominator(F))*log(2)
 );
 if (in%2,res=-res);
 res;
}

/* (IX.5) F rational fraction. Computes the sum from in to infinity of      */
/* (-1)^n*F(n). By default in=0. Uses another reduction. Optional parameter */
/* la is here set to 3, but only controls the speed. Usually slower than    */
/* above. */

sumratalt2(F,in=0,la=3)=
{
 my(vx,res);
 if (type(in)!="t_INT",error("not an integer starting value in sumratalt2"));
 vx=variable(F);
 if (in,F=subst(F,vx,vx+in));
 res=sumrat(subst(F,vx,2*vx)-subst(F,vx,2*vx+1),0,la);
 if (in%2,res=-res);
 res;
}

/* (IX.6) F rational function. Computes the integral from N to infty of    */
/* log(F).                                                                 */

intinflograt(F,N,flerr=0)=
{
 my(s,vx,pr,lim,G,ro,r);
 s=0.;
 if (poldegree(F-1)>=-1,error("infinite integral in intinflograt"));
 vx=variable(F);
 G=subst(F,vx,1/vx);
 G=G'/G;
 pr=precision(1.);
 default(realprecision,18);
 ro=polroots(denominator(G));
 r=1;for(k=1,length(ro),r=min(r,norml2(ro[k])));
 r=1/r;
 if (N<r+1,
  if (flerr,return([r]),
  error("N not large enough in intinflograt"))
 );
 lim=ceil(pr*log(10)/log(N/sqrt(r)))+1;
 default(realprecision,pr+5);
 G=1.*G+O(vx^(lim+1));
 s=0.;
 forstep (k=lim,1,-1,s=(polcoeff(G,k)/(k*(k+1))+s)/N);
 s=precision(s,pr);
 default(realprecision,pr);
 s;
}

/* (IX.7) F rational function. Computes the product from in to infinity of */
/* F(n), default in=0. We must have F(x)=1+O(1/x^2) as x tends to infty.   */

{
prodrat(F,in=0,la=3)=
 my(vx,pr,lim,nn,res,tes,sal);
 vx=variable(F);
 if (in,F=subst(F,vx,vx+in));
 if (poldegree(F-1)>=-1,error("infinite product in prodrat"));
 pr=precision(1.);
 default(realprecision,9);
 la=max(la+0.,1);
 lim=ceil((pr*log(10)+2)/(2*(1+log(la))));
 nn=ceil((lim+.25)/Pi*la);
 default(realprecision,pr+5);
 tes=intinflograt(F,nn,1);
 if (type(tes)=="t_VEC",
  nn=max(nn,ceil(2*tes[1]));
  default(realprecision,9);
  lim=ceil(nn*Pi/la+.75);
  default(realprecision,pr+5);
  tes=intinflograt(F,nn);
 );
 res=prod(ijkl=0,nn-1,subst(F,vx,ijkl),sqrt(subst(F,vx,nn)));
 sal=1.*subst(F'/F,vx,vx+nn)+O(vx^(2*lim-1));
 bernreal(2*lim);
 for (k=1,lim,
  tes-=bernreal(2*k)*polcoeff(sal,2*k-2)/((2*k)*(2*k-1))
 );
 res*=exp(tes);
 res=precision(res,pr);
 default(realprecision,pr);
 res;
}

/* (IX.8) clear */

logzetaa(N,A)=log(zeta(N)*prodeuler(p=2,A,1-1/p^N));

/* (IX.9) F rational function. Sum of F(p), p over all prime numbers */

sumeulerrat(F,pin,A=30)=
{
 my(vx,pr,mro,lim,sal,res);
 vx=variable(F);
 if (poldegree(F)>=-1,error("infinite sum in sumeulerrat"));
 pr=precision(1.);
 default(realprecision,9);
 mro=max(1,vecmax(abs(polroots(denominator(F)))));
 A=ceil(max(A,3*mro));
 A=max(A,pin);
 lim=ceil(pr*log(10)/log(A/mro))+1;
 default(realprecision,pr+5);
 sal=1.*subst(F,vx,1/x)+O(x^(lim+1)); 
 res=sum(N=2,lim,logzetaa(N,A)*sumdiv(N,k,moebius(k)*polcoeff(sal,N/k)/k));
 forprime (p=pin,A,res+=subst(F,vx,p));
 res=precision(res,pr);
 default(realprecision,pr);
 res;
}

/* (IX.10) F rational function. Product of F(p), p over all prime numbers */

prodeulerrat(F,pin,A=30)=
{
 my(vx,pr,mro,lim,sal,res);
 vx=variable(F);
 if (poldegree(F-1)>=-1,error("infinite product in prodrat"));
 pr=precision(1.);
 default(realprecision,9);
 mro=max(1,vecmax(abs(polroots(denominator(F)))));
 mro=max(mro,vecmax(abs(polroots(numerator(F)))));
 A=ceil(max(A,3*mro));
 A=max(A,pin);
 lim=ceil(pr*log(10)/log(A/mro))+1;
 default(realprecision,pr+5);
 sal=log(1+1.*(subst(F,vx,1/x)-1)+O(x^(lim+1))); 
 res=exp(sum(N=2,lim,logzetaa(N,A)*sumdiv(N,k,moebius(k)*polcoeff(sal,N/k)/k)));
 forprime (p=pin,A,res*=subst(F,vx,p));
 res=precision(res,pr);
 default(realprecision,pr);
 res;
}

/*******************************************************/
/* p-adic log conjectures                              */
/*******************************************************/

/* Compute vector(ff,k,log_p(1-z_m^(u*p^k))-p*log_p(1-z_m^(u*p^(k-1)))) */
/* when m is not a power of p and m \nmid u p^\infty                   */

logpmaux(p,m,u,ff,padp=DFPADP)=
{
 my(pol,q,lim,polp,ratp,s,qs);
 pol=polcyclo(m); 
 q=(1+O(p^padp))*vector(ff,k,1/Mod(x^(u*p^(k-1))-1,pol)); 
 lim=padp*(p-1);
 polp=(x+1)^p-x^p;
 ratp=-(polp'+O(x^(lim+1)))/polp;
 s=vector(ff,k,O(p^padp)); qs=vector(ff,k,1+O(p^padp));
 for (i=1,lim, 
  qs=vector(ff,k,q[k]*qs[k]);
  s=vector(ff,k,s[k]+qs[k]*polcoeff(ratp,i-1)/i)
 );
 return (-s);
}

/* Compute log_p(1-z_m^(u*p))-p*log_p(1-z_m^u) under same cond. */
 

logpm1(p,m,u,padp=DFPADP)=logpmaux(p,m,u,1,padp)[1];
 

/* Compute log_p(1-z_m^u) when p \nmid m               */

logpm2(p,m,u,padp=DFPADP)=
{
 my(ff,res);
 if (m%p==0,error("m divisible by p"));
 ff=znorder(Mod(p,m));
 res=logpmaux(p,m,u,ff,padp);
 return (sum(k=1,ff,res[k]/p^k)/(1/p^ff-1));
}

/* Compute log_p(1-z_m^u) when m is not a power of p */

logpm3(p,m,u,padp=DFPADP)=
{
 my(v,res1,res2);
 v=valuation(m,p); m1=m/p^v;
 if (v==0,return(logpm2(p,m,u,padp)));
 res1=logpmaux(p,m,u,v,padp);
 res2=Mod(subst(lift(logpm2(p,m1,u,padp)),x,x^(p^v)),polcyclo(m));
 return (res2/p^v-sum(k=1,v,res1[k]/p^k));
}

/* Compute vector(v-1,k,log_p(1-z_m^(u*p^k))-p*log_p(1-z_m^(u*p^(k-1)))) */
/* when m = p^v                                                        */

logpmaux2(p,u,v,padp=DFPADP)=
{
 my(pol,q,lim,polp,ratp,s,qs);
 if (v<=1,return([]));
 m=p^v;
 pol=polcyclo(m); 
 q=(1+O(p^(7*padp)))*vector(v-1,k,1/Mod(x^(u*p^(k-1))-1,pol)); 
 lim=ceil(padp*(p-1)*m/(m-p));
 polp=(x+1)^p-x^p;
 ratp=-(polp'+O(x^(lim+1)))/polp;
 s=vector(v-1,k,O(p^(7*padp))); qs=vector(v-1,k,1+O(p^(7*padp)));
 for (i=1,lim, 
  qs=vector(v-1,k,q[k]*qs[k]);
  s=vector(v-1,k,s[k]+qs[k]*polcoeff(ratp,i-1)/i);
/* print(i,": ",polcoeff(lift(s[1]),0)) */
 );
 return (-s*(1+O(p^padp)));
}

/* Compute log_p(1-z_p^u) when p \nmid u */

logp(p,u,padp=DFPADP)=
{
 my(lim,z,t);
 lim=(p-1)*padp;
 z=Mod(x,polcyclo(p));
 t=(1+(1-z^u)^(p-1)/p);
 return (-sum(i=1,lim,t^i/i*(1+O(p^padp)))/(p-1));
}

/* Compute log_p(1-z_m^u) when m \nmid u and m=p^v */

logpv(p,u,v,padp=DFPADP)=
{
 local (res,w);
 w=valuation(u,p); if (w>=v,error("logpv"));
 v-=w; u/=p^w;
 res=logpmaux2(p,u,v,padp);
 res2=lift(logp(p,u,padp));
 res2=Mod(subst(res2,x,x^(p^(v-1))),polcyclo(p^v));
 return (res2/p^(v-1)-sum(k=1,v-1,res[k]/p^k));
}

/* Compute log_p(1-z_m^u) for m \nmid u */

logpall(p,m,u,padp=DFPADP)=
{
 my(d,res);
 d=gcd(m,u); m/=d; u/=d;
 if (m/p^valuation(m,p)!=1,
  res=logpm3(p,m,u,padp),
  if (m==p,
   res=logp(p,u,padp),
   res=logpv(p,u,valuation(m,p),padp)
  )
 );
 if (d>1,
  res=Mod(1,polcyclo(m*d))*subst(lift(res),x,x^d)
 );
 return (res);
}

/* check conjecture 1 */

conj1(r,ff,p)=psi_p(r/ff,p)-polcoeff(lift(sum(a=1,ff-1,Mod(x,polcyclo(ff))^(-a*r)*(logpall(p,ff,a)-(1/p)*logpall(p,ff,a*p)))),0)+Euler_p(p)+(1-1/p)*log(ff+O(p^DFPADP));

/* check conjecture 2 */

conj2(r,ff,p)=psiL_p(r/ff,p)-polcoeff(lift(sum(a=1,ff-1,Mod(x,polcyclo(ff))^(-a*r)*logpall(p,ff,a))),0)+p/(p-1)*Euler_p(p)+log(ff+O(p^DFPADP));

/* must be zero when p \mid ff and p \nmid r */

A1(m)=my(fa);fa=factor(m)[,1];sum(i=1,#fa,log(fa[i])/(fa[i]-1));

L1(zns,chara,z)=my(m);m=zns[5];if(conductor_char(zns,chara)>1,(-1/m)*sum(r=1,m,subst(lift(eval_char(zns,chara,r)),x,z)*psi(r/m)),Euler*eulerphi(m)/m);

av1all(m)=my(zns,li,z);zns=znstar_2(m);z=exp(2*Pi/eulerphi(m));li=all_chars(zns);sum(i=1,#li,L1(zns,li[i],z))/eulerphi(m);

av1prim(m)=my(zns,li,z);if(m%4==2,error("no primitive characters"));zns=znstar_2(m);z=exp(2*Pi/eulerphi(m));li=allprim_chars(zns);sum(i=1,#li,L1(zns,li[i],z))/#li;

test(f,k,eps)=iee=if((-1)^k*eps==1,1,0);zns=znstar_2(f);li=allprim_chars(zns);for(i=1,#li,ch=li[i];if(iseven_char(zns,ch)==iee,print(lift((sum(r=1,f/2-1,eval_char(zns,ch,r)*r^k)+2*Bchi(zns,ch,k+1)/(k+1))/f))));
