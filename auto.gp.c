/*-*- compile-command: "/usr/bin/gcc -c -o auto.gp.o -O3 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -I"/usr/local/include" auto.gp.c && /usr/bin/gcc -o auto.gp.so -shared -O3 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari"; -*-*/
#include <pari/pari.h>
#include <pari/paripriv.h>
#include "float.h"
/*
GP;install("log_2","Gp","lg","./auto.gp.so");
GP;addhelp(lg, "lg(x): Binary logarithm of x.");
GP;install("rad","G","rad","./auto.gp.so");
GP;addhelp(rad, "rad(n): Radical of n, the largest squarefree number dividing n. Sloane's A007947.");
GP;install("isprimepower","lG","isprimepower","./auto.gp.so");
GP;addhelp(isprimepower, "isprimepower(n): Is n a prime power? Sloane's A010055; characteristic function of Sloane's A000961.");
GP;install("isPowerful","lG","isPowerful","./auto.gp.so");
GP;install("isPowerfulCorrect","lG","isPowerfulC","./auto.gp.so");
GP;addhelp(isPowerful, "isPowerful(n): Is n powerful (min exponent 2)? Sloane's A112526; characteristic function of Sloane's A001694.");
GP;install("prp","lGDG","prp","./auto.gp.so");
GP;addhelp(prp, "prp(n,b=2): Is n a b-probable prime?");
GP;install("sprp","lGDG","sprp","./auto.gp.so");
GP;addhelp(sprp, "sprp(n,b=2): Is n a b-strong probable prime?");
GP;install("sopf","G","sopf","./auto.gp.so");
GP;addhelp(sopf, "sopf(n): Sum of prime factors of n. Sloane's A008472.");
GP;install("sopfr","G","sopfr","./auto.gp.so");
GP;addhelp(sopfr, "sopfr(n): Sum of prime factors of n (with multiplicity). Sloane's A001414.");
GP;install("primorial","G","primorial","./auto.gp.so");
GP;addhelp(primorial, "Returns the product of each prime less than or equal to n. Sloane's A034386.");
GP;install("gpf","G","gpf","./auto.gp.so");
GP;addhelp(gpf, "The greatest prime factor of a number. Sloane's A006530.");
GP;install("lpf","G","lpf","./auto.gp.so");
GP;addhelp(lpf, "The least prime factor of a number. Sloane's A020639.");
GP;install("isTriangular","lG","isTriangular","./auto.gp.so");
GP;addhelp(isTriangular, "isTriangular(n): Is n a triangular number? Sloane's A010054; characteristic function of Sloane's A000217.");
GP;install("isHexagonal","lG","isHexagonal","./auto.gp.so");
GP;addhelp(isHexagonal, "isHexagonal(n): Is n a hexagonal number? Characteristic function of Sloane's A000384.");
GP;install("isFibonacci","lG","isFibonacci","./auto.gp.so");
GP;addhelp(isFibonacci, "isFibonacci(n): Is n a Fibonacci number? Sloane's A010056; characteristic function of Sloane's A000045.");
GP;install("largestSquareFactor","G","largestSquareFactor","./auto.gp.so");
GP;addhelp(largestSquareFactor, "largestSquareFactor(n): Largest square dividing n. Sloane's A008833.");
GP;install("hamming","lG","hamming","./auto.gp.so");
GP;addhelp(hamming, "hamming(n): Hamming weight of n (considered as a binary number). Sloane's A000120.");
GP;install("ispow2","lG","ispow2","./auto.gp.so");
GP;addhelp(ispow2, "ispow2(n): Is n a power of two? Characteristic function of Sloane's A000079.");
GP;install("issemiprime","lG","issemi","./auto.gp.so");
GP;addhelp(issemi, "issemi(n): Is n a semiprime? Sloane's A064911; characteristic function of Sloane's A001358.");
GP;install("istwo","lG","istwo","./auto.gp.so");
GP;addhelp(istwo, "Is the number a sum of two squares? Characteristic function of Sloane's A001481.");
coGP;install("ways2","G","ways2","./auto.gp.so");
GP;addhelp(ways2, "Number of ways that n can be represented as a sum of two squares. Sloane's A000161.");
GP;install("isthree","lG","isthree","./auto.gp.so");
GP;addhelp(isthree, "isthree(n): Is the number the sum of three squares? Sloane's A071374; characteristic function of Sloane's A000378.");
GP;install("ways3","G","ways3","./auto.gp.so");
GP;addhelp(ways3, "Number of ways that n can be represented as a sum of three squares. Sloane's A000164.");
GP;install("msb","G","msb","./auto.gp.so");
GP;addhelp(msb, "msb(n): Most significant bit of n: returns the greatest power of 2 <= the number. Sloane's A053644.");
GP;install("Faulhaber","LDG","Faulhaber","./auto.gp.so");
GP;addhelp(Faulhaber, "Faulhaber(e,{a='x}): Returns the polynomial for the sum 1^e + 2^e + ... + x^e, evaluated at a.");
GP;install("rp","L","rp","./auto.gp.so");
GP;addhelp(rp, "rp(b): Returns a random b-bit prime.");
GP;install("countPowerful","D0,G,","countPowerful","./auto.gp.so");
GP;addhelp(countPowerful, "countPowerful(lim): Number of powerful numbers up to lim. Partial sum of characteristic function of of Sloane's A001694.");
GP;install("countSquarefree","D0,G,","countSquarefree","./auto.gp.so");
GP;addhelp(countSquarefree, "countSquarefree(lim): Counts the number of squarefree numbers up to lim.");
GP;install("Mfactor","GD0,G,DG","Mfactor","./auto.gp.so");
GP;addhelp(Mfactor, "Mfactor(p,lim,{start=2}): Returns factors of the Mersenne number 2^p-1 up to lim, starting at start, provided p is a prime = 3 mod 4. Same as bigfactor(2,p,1,lim,start) but faster because it checks only factors of the form 2kp+1 that are +/- 1 mod 8.");
GP;install("bigfactor","GGGD0,G,DG","bigfactor","./auto.gp.so");
GP;addhelp(bigfactor, "bigfactor(a,b,c,lim,{start=2}): Find small prime factors of a^b - c (up to lim). Optional parameter start gives a starting point below which primes are not checked.");
GP;install("bigdiv","lGGGG","bigdiv","./auto.gp.so");
GP;addhelp(bigdiv, "bigdiv(a,b,c,d): Does d divide a^b - c? Same as (a^b-c)%d == 0, but faster for large b. Example: bigdiv(2,p,1,d) checks if d divides the p-th Mersenne number.");
GP;install("contfracback","D0,G,DG","contfracback","./auto.gp.so");
GP;addhelp(contfracback, "contfracback(v, terms): Given a continued fraction v, gives the real number back. If terms is given, use only that many terms.");
GP;install("W","D0,G,p","W","./auto.gp.so");
GP;addhelp(W, "W(x): Primary branch of Lambert's W function. Finds an L >= -1 such that L * exp(L) = x, where x >= -1/e.");
GP;install("vecsum","G","vecsum","./auto.gp.so");
GP;addhelp(vecsum, "vecsum(v): Sum of the elements of v.");
GP;install("vecprod","G","vecprod","./auto.gp.so");
GP;addhelp(vecprod, "vecprod(v): Product of the elements of v.");
GP;install("vecgcd","G","vecgcd","./auto.gp.so");
GP;addhelp(vecgcd, "Vector gcd: returns the gcd of all elements in the vector.");
GP;install("veclcm","G","veclcm","./auto.gp.so");
GP;addhelp(veclcm, "Vector lcm: returns the lcm of all elements in the vector.");
GP;install("oddres","G","oddres","./auto.gp.so");
GP;addhelp(oddres, "oddres(n): Returns the greatest odd number dividing n.");
GP;install("toC","vG","toC","./auto.gp.so");
GP;addhelp(toC, "toC(n): Format n for use with the Pari library (e.g., with gp2c programs).");
GP;install("digits","D0,G,","digits","./auto.gp.so");
GP;addhelp(digits, "digits(n): Number of decimal digits in n. Sloane's A055642.");
GP;install("eps","p","eps","./auto.gp.so");
GP;addhelp(eps, "Returns machine epsilon for the current precision.");
GP;install("fnice","G","fnice","./auto.gp.so");
GP;addhelp(fnice, "fnice(n): Returns a string with a 'nice' factorization of n.");
GP;install("tonice","D0,G,p","nice","./auto.gp.so");
GP;addhelp(nice, "nice(o): Reformats the object o 'nicely' when possible. Currently chokes on multivariable polynomials.");
GP;install("sumset","D0,G,D0,G,","sumset","./auto.gp.so");
GP;addhelp(sumset, "sumset(A, B) is the set of all numbers of the form a+b, a in A, b in B.");
GP;install("diffset","D0,G,D0,G,","diffset","./auto.gp.so");
GP;addhelp(diffset, "diffset(A, B) is the set of all numbers of the form a-b, a in A, b in B.");
GP;install("normd","D0,G,D0,G,p","normd","./auto.gp.so");
GP;addhelp(normd, "normd(a,b): Amount of the normal distribution between a and b standard deviations. Plus/minus infinity coded as [+1]/[-1].");
GP;install("rnormal","p","rnormal","./auto.gp.so");
GP;addhelp(rnormal, "rnormal(): Returns a random normal variable with mean 0 and standard deviation 1 at the current precision.");
GP;install("checkVDW","D0,G,DG","checkVDW","./auto.gp.so");
GP;addhelp(checkVDW, "checkVDW(vv): Given a partition vv = [p1, p2, ...] with union(p1, p2, ...) = [1, 2, ..., n], finds a lower-bound proof for van der Waerden numbers based on the partition. Returns 0 if vv is not a partition of any initial segment, and k if vv proves that W(#vv, k) > n.");
GP;install("longestProgression","D0,G,","longestProgression","./auto.gp.so");
GP;addhelp(longestProgression, "longestProgression(v): Finds the longest arithmetic progression in v. Assumes that v is a vector of integers sorted from smallest to largest. Uses a space-efficient naive algorithm.");
GP;install("longestProgression1","D0,G,","longestProgression1","./auto.gp.so");
GP;addhelp(longestProgression1, "longestProgression1(v): Uses a quadratic algorithm of Jeff Erickson, which is worst-case optimal; better algorithms are available when there are long progressions (> lg #v lg lg #v).");
GP;install("fusc","D0,G,","fusc","./auto.gp.so");
GP;addhelp(fusc, "fusc(n): Stern's diatomic series, which has many interpretations. Sloane's A002487.");
GP;install("forodd","vV=GGI","forodd","./auto.gp.so");
GP;addhelp(forodd, "forodd(X=a,b,seq): the sequence is evaluated, X running over the odds between a and b.");
GP;install("fortwin","vV=GGI","fortwin","./auto.gp.so");
GP;addhelp(fortwin, "fortwin(X=a,b,seq): the sequence is evaluated, X running over the twin primes between a and b.");
GP;install("forthinprime","vV=GGI","forthinprime","./auto.gp.so");
GP;addhelp(forthinprime, "forthinprime(X=a,b,seq): the sequence is evaluated, X running over the primes between a and b, even if b > primelimit. EXPERIMENTAL!");
GP;install("forbigprime","vV=GGI","forbigprime","./auto.gp.so");
GP;addhelp(forbigprime, "forbigprime(X=a,b,seq): the sequence is evaluated, X running over the primes between a and b. EXPERIMENTAL!");
GP;install("sumformal","V=GGG","sumformal","./auto.gp.so");
GP;addhelp(sumformal, "sumformal(X=start,end,expr): Formal version of sum(X=start,end,expr). Start and end can be expressions instead of numbers.");
GP;install("fibmod","GG","fibmod","./auto.gp.so");
GP;addhelp(fibmod, "fibmod(n,m): Returns the nth Fibonacci number mod m. Same as finonacci(n)%m, but faster for large n.");
GP;install("bfile","GDGDG","bfile","./auto.gp.so");
GP;addhelp(bfile, "bfile(name, v, offset=1): If v is given, creates a b-file with the values of v, using name as the A-number (given as a number or a filename). If v is not given, open the b-file name (again, as a filename or number) and return a vector of its values.");
GP;install("dsum","G","dsum","./auto.gp.so");
GP;addhelp(dsum, "dsum(n): Digit sum of n. Sloane's A007953.");
GP;install("checkmult","lGD1,L,","checkmult","./auto.gp.so");
GP;addhelp(checkmult, "checkmult(v,{verbose=1}): Is the sequence v multiplicative?");
GP;install("checkcmult","lGD1,L,","checkcmult","./auto.gp.so");
GP;addhelp(checkcmult, "checkcmult(v,{verbose=1}): Is the sequence v completely multiplicative?");
GP;install("checkdiv","lGD1,L,","checkdiv","./auto.gp.so");
GP;addhelp(checkdiv, "checkdiv(v,{verbose=1}): Is v a divisibility sequence?");
GP;install("solvePell","Gp","solvePell","./auto.gp.so");
GP;addhelp(solvePell, "solvePell(n): Returns a solution to the equation x^2 - ny^2 = 1.");
GP;install("Engel","Gp","Engel","./auto.gp.so");
GP;addhelp(Engel, "Engel(x): Engel expansion of x.");
GP;install("Eng","G","Eng","./auto.gp.so");
GP;addhelp(Eng, "Eng(n): English name of the number n.");
GP;install("composite","L","composite","./auto.gp.so");
GP;addhelp(composite, "composite(n): Returns the n-th composite. Sloane's A002808.");
GP;install("deBruijnXi","G","deBruijnXi","./auto.gp.so");
GP;addhelp(deBruijnXi, "deBruijnXi(x): Helper function for rhoest.  Finds a xi such that e^xi - 1 = x * xi.");
GP;install("rhoest","Gp","rhoest","./auto.gp.so");
GP;addhelp(rhoest, "de Bruijn's asymptotic approximation for rho(x), rewritten as in van de Lune and Wattel 1969.  Curiously, their paper shows values for this estimate that differ from those calculated by this function, often as soon as the second decimal place -- but as the difference is in the direction of the true value, I have not looked further into this.");
GP;install("DickmanRho","Gp","DickmanRho","./auto.gp.so");
GP;addhelp(DickmanRho, "Estimates the value of the Dickman rho function. For x <= 3 the exact values are used, up to rounding; up to 15 the value is interpolated using known values and rhoest; after 15 rhoest is used, along with a correction factor based on the last value in rhoTable.");
// New \/
GP;install("tetrMod","GGG","tetrMod","./auto.gp.so");
GP;addhelp(tetrMod, "tetrMod(a,b,M): Returns a^^b mod M.");
GP;install("iscyclo","lG","iscyclo","./auto.gp.so");
GP;addhelp(iscyclo, "iscyclo(f): Is f a cyclotomic polynomial?  Uses the Bradford-Davenport algorithm.");
GP;install("istotient","lG","istotient","./auto.gp.so");
GP;addhelp(istotient, "istotient(n): Does there exist some m such that eulerphi(m) = n?");
* // New /\
// No associated help
GP;install("init_auto","v","init_auto","./auto.gp.so");
GP;install("consistency","l","consistency","./auto.gp.so");
GP;install("ucountPowerfuli","lD0,G,","cP","./auto.gp.so");
GP;install("ucountSquarefree","lL","cS","./auto.gp.so");
GP;install("graeffe","G","graeffe","./auto.gp.so");
GP;install("totientHelper","lGDG","totientHelper","./auto.gp.so");
*/

// Removed: //GP;install("pBounds","vD0,G,DGp","pBounds","./auto.gp.so");
// Removed: //GP;addhelp(pBounds, "pBounds(n, verbose=0): Estimates the nth prime. Set verbose=1 to get a list of sources for the results.");



//////////////////////////////////////////////////////////// New
GEN graeffe(GEN f);
long iscyclo(GEN f);
long istotient(GEN n);
long totientHelper(GEN n, GEN m);
//////////////////////////////////////////////////////////// New


long checkmult(GEN v, long verbose);
long checkcmult(GEN v, long verbose);
long checkdiv(GEN v, long verbose);
GEN solvePell(GEN n, long prec);
GEN tetrMod(GEN a, GEN b, GEN M);
GEN Engel(GEN x, long prec);
GEN Eng(GEN n);
GEN Edigit(GEN n);
GEN composite(long n);
GEN deBruijnXi(GEN x);
GEN rhoest(GEN x, long prec);
GEN DickmanRho(GEN x, long prec);
GEN dsum(GEN n);
ulong dsum_small(ulong n);
long issemiprime(GEN n);
long uissemiprime(ulong n);
GEN phiset(GEN v);
GEN rad(GEN n);
long isprimepower(GEN n);
long uisprimepower(ulong n);
INLINE long valu(ulong n);
long isPowerful_small(ulong n);
long isPowerful(GEN n);
long prp(GEN n, GEN b);
long sprp(GEN n, GEN b);
GEN sopf(GEN n);
GEN sopfr(GEN n);
GEN primorial(GEN n);
GEN prodtree(GEN A, long start, long stop);
GEN gpf(GEN n);
GEN lpf(GEN n);
long isTriangular(GEN n);
long isHexagonal(GEN n);
long isFibonacci(GEN n);
long isSmallFib(long n);
GEN fibmod(GEN n, GEN m);
long Pisano(long p, long e);
GEN largestSquareFactor(GEN n);
INLINE long hamming_word(ulong w);
long hamming(GEN n);
long ispow2(GEN n);
long istwo(GEN n);
GEN ways2(GEN n);
long isthree(GEN n);
INLINE long uissquare(ulong n);
long sways3s(ulong n);
GEN ways3(GEN n);
GEN msb(GEN n);
GEN sumformal(GEN expr, GEN start, GEN end);
GEN Faulhaber(long e, GEN a);
GEN rp(long b);
GEN countPowerful(GEN lim);
GEN countSquarefree(GEN lim);
ulong ucountSquarefree(ulong lim);
GEN Mfactor(GEN p, GEN lim, GEN start);
GEN bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start);
long bigdiv(GEN a, GEN b, GEN c, GEN d);
GEN contfracback(GEN v, GEN terms);
GEN W(GEN x, long prec);
GEN vecsum(GEN v);
GEN vecprod(GEN v);
GEN vecgcd(GEN v);
GEN veclcm(GEN v);
GEN oddres(GEN n);
void toC(GEN n);
GEN digits(GEN x);
GEN eps(long prec);
GEN fnice(GEN n);
GEN tonice(GEN o, long prec);
GEN initial(GEN n, char *s);
GEN medial(GEN n, char *s);
GEN monomialnice(GEN coeff, GEN degree, GEN v);
GEN sumset(GEN a, GEN b);
GEN diffset(GEN a, GEN b);
GEN normd(GEN a, GEN b, long prec);
GEN rnormal(long prec);
void pBounds(GEN n, GEN verbose, long prec);
GEN checkVDW(GEN vv, GEN verbose);
GEN longestProgression(GEN v);
GEN longestProgression1(GEN v);
GEN fusc(GEN n);
GEN fusc_large(GEN n);
ulong fusc_small(GEN n);
ulong cuberoot(ulong n);
ulong ucountPowerfulu(ulong lim);
long issquarefree_small(ulong n);
void forodd(GEN a, GEN b, GEN code);
char* getBValue(char*);
GEN bfile(GEN name, GEN v, GEN offset);
GEN bfilein(char* name);
void fortwin(GEN ga, GEN gb, GEN code);
void forbigprime(GEN ga, GEN gb, GEN code);
void forbigprime_sieve(ulong a, ulong b, GEN code);
void forthinprime(ulong a, ulong b, GEN code);
INLINE GEN gtor(GEN x, const char* funcName, long prec);
void init_auto(void);
/*End of prototype*/

#include "othergpincludes.h"


#define FAKE_PREC 0		// Used when a precision is required but will not be used
#define NEVER_USED 0	// Used to initialize values so the compiler doesn't complain
GEN rnormal_cached;
void
init_auto(void)
{
	pari_sp ltop = avma;
	rnormal_cached = 0;
	avma = ltop;
}


/*
 * TODO: Split file into several parts, comment dependencies,
 * and upload to website
*/

/*
default(debug,4)
consistency()
*/


long
consistency()
{
	pari_sp ltop = avma;
	ulong n;
	long bad = 0;
	long failures = 0;
	long verbose = DEBUGLEVEL > 3;
	
	pari_printf("Cube roots...");
	pari_flush();
#ifdef LONG_IS_64BIT
	for (n = 1; n <= 2642245; n++) {
#else
	for (n = 1; n <= 1625; n++) {
#endif
		ulong c = cuberoot(n * n * n);
		ulong c1 = cuberoot(n * n * n - 1);
		if (c != n) {
			if (verbose)
				pari_printf("\n	cuberoot(%lu) should be %lu, was %lu", n * n * n, n, c);
			bad = 1;
			break;
		}
		if (c1 != n - 1) {
			if (verbose)
				pari_printf("\n	cuberoot(%lu) should be %lu, was %lu", n * n * n - 1, n - 1, c1);
			bad = 1;
			break;
		}
	}
	if (bad)
		pari_printf(verbose ? "; failed.\n" : "failed.\n");
	else
		pari_printf(" ok.\n");
	failures += bad;
	bad = 0;


	pari_printf("issquarefree_small...");
	pari_flush();
	if (issquarefree_small(4293001441)) {
		if (verbose)
			pari_printf("\n	issquarefree_small(4293001441) should have been false");
		bad = 1;
	}
	if (bad)
		pari_printf(verbose ? "; failed.\n" : "failed.\n");
	else
		pari_printf(" ok.\n");
	failures += bad;
	bad = 0;


	pari_printf("ucountPowerfulu...");
	pari_flush();
	if (cmpiu(countPowerful(gen_0), ucountPowerfulu(0)) != 0) {
		if (verbose)
			pari_printf("\n	ucountPowerfulu(0) != countPowerful(0): %lu != %Ps", ucountPowerfulu(n), countPowerful(utoi(n)));
		bad = 1;
	}
	for (n = 100; n%100 == 0; n *= 10) {
		if (cmpiu(countPowerful(utoi(n)), ucountPowerfulu(n)) != 0) {
			if (verbose)
				pari_printf("\n	ucountPowerfulu(%lu) != countPowerful(%Ps): %lu != %Ps", n, utoi(n), ucountPowerfulu(n), countPowerful(utoi(n)));
			bad = 1;
		}
		avma = ltop;
	}
	if (bad)
		pari_printf(verbose ? "; failed.\n" : "failed.\n");
	else
		pari_printf(" ok.\n");
	failures += bad;
	bad = 0;

	return failures;
}


INLINE GEN
gtor(GEN x, const char* funcName, long prec)
{
	switch (typ(x)) {
		case t_REAL:
			return x;	// x, not a copy of x
		case t_INT:
		case t_FRAC:
			return cxcompotor(x, prec);
		default:
			pari_err(typeer, funcName);
	}
	return NEVER_USED;
}


/******************************************************************************/
/**                     Prime-related arithmetic functions                   **/
/******************************************************************************/

// TODO: Optimize multi-word inputs
long
issemiprime(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "issemiprime");
	if (signe(n) <= 0)
		return 0;
	
	ulong nn = itou_or_0(n);
	if (nn)
		return uissemiprime(nn);
		
	pari_sp ltop = avma;
	if (!mpodd(n)) {
		long ret = mod4(n) && isprime(shifti(n, -1));
		avma = ltop;
		return ret;
	}
	
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > 997)	// What is a good breakpoint here?
			break;
		if (dvdis(n, p))	// "1 if p divides n and 0 otherwise"
		{
			long ret = isprime(diviuexact(n, p));
			avma = ltop;
			return ret;
		}
	}
	
	if (isprime(n))
		return 0;
	
	GEN fac = Z_factor_until(n, shifti(n, -1));	// Find a nontrivial factor -- returns just the factored part
	GEN expo = gel(fac, 2);
	GEN pr = gel(fac, 1);
	long len = glength(expo);
	if (len > 2) {
		avma = ltop;
		return 0;
	}
	if (len == 2) {
		if (cmpis(gel(expo, 1), 1) > 0 || cmpis(gel(expo, 2), 1) > 0) {
			avma = ltop;
			return 0;
		}
		GEN p = gel(pr, 1);
		GEN q = gel(pr, 2);
		long ret = !cmpii(mulii(p, q), n) && isprime(p) && isprime(q);
		avma = ltop;
		return ret;
	}
	if (len == 1) {
		long e = itos(gel(expo, 1));
		if (e != 2) {
			avma = ltop;
			return 0;
		}
		GEN p = gel(pr, 1);
		long ret = !cmpii(sqri(p), n) && isprime(p);
		avma = ltop;
		return ret;
	}
	
	pari_err(talker, "Z_factor_until returned an unexpected value %Ps at n = %Ps, exiting...", fac, n);
	avma = ltop;
	return NEVER_USED;
}


void
dostuff(GEN lm) {
	byteptr d = diffptr, d1 = diffptr;
	ulong p = 0, q = 0;
	ulong lim = maxprime();
	if (lim < 1000000)
		pari_err(primer1, "<-- you really should have a billion primes precalculated for this function, but you need at least a million or so.");
	lim -= 1000000;
	
	lim = minuu(itou_or_0(lm), lim);

	ulong sum = 0;	// Is a ulong big enough?  Needs to hold approximately p^2 log p.
	for(;;)
	{
		ulong pdiff = *d + 1;	// Increase number of primes by g_n: lose one at the beginning and add g_n + 1 at the end
		NEXT_PRIME_VIADIFF(p,d);
		while(pdiff) {
			NEXT_PRIME_VIADIFF(q,d1);
			sum += q;
			pdiff--;
		}
		sum -= p;
		
		if (uisprime(sum))
			printf("%llu, ", (long long unsigned int)p);

		if (q > lim) break;
	}
	
	pari_printf("\n");
}

// FIXME: Handle the case of small primelimit
/*
P=primorial(661)/2;
v=vector(10^4);i=0;forstep(n=2^64-1,1,-2,if(1==gcd(n,P),v[i++]=n;if(i==#v,return(#v))))
#
sum(i=1,#v,issemi(v[i]))
v=vector(10^5);i=0;forstep(n=2^64-1,1,-2,if(1==gcd(n,P),v[i++]=n;if(i==#v,return(#v))))
sum(i=1,#v,issemi(v[i]))

sum(i=1,#v,print(v[i]);issemi(v[i]))

No rho, crossover  661: 27.68s for 10k, 4:38.21s for 100k
   Rho, crossover  661: 28.06s for 10k, 4:45.50s for 100k
   Rho, crossover 1000: 28.05s for 10k, 4:43.82s for 100k
   Rho, crossover 3000: 27.98s for 10k, 4:44.21s for 100k

With prime test fronting:
   Rho, crossover  661: 13.60s for 10k, 2:12.36s for 100k
   Rho, crossover 1000: 13.72s for 10k, 2:12.53s for 100k
   Rho, crossover 1500: 13.36s for 10k, 2:12.98s for 100k
   Rho, crossover 2000: 13.65s for 10k, 2:12.80s for 100k
   Rho, crossover 2500: 13.57s for 10k, 2:12.98s for 100k
   Rho, crossover 3000: 13.50s for 10k, 2:13.30s for 100k
   Rho, crossover 4000: 13.48s for 10k, 2:12.98s for 100k
   Rho, crossover 9000: 13.54s for 10k, 2:14.17s for 100k

With prime test fronting, before even rho:
   Rho, crossover  661: 12.43s for 10k, 2:01.02s for 100k
   Rho, crossover  750: 12.50s for 10k, 2:01.03s for 100k
   Rho, crossover 1000: 12.38s for 10k, 2:00.59s for 100k
   Rho, crossover 1250: 12.40s for 10k, 2:01.12s for 100k
   Rho, crossover 1500: 12.38s for 10k, 2:00.77s for 100k
   Rho, crossover 2000: 12.58s for 10k, 2:01.75s for 100k
   Rho, crossover 3000: 12.42s for 10k, 2:01.08s for 100k
*/
long
uissemiprime(ulong n)
{
#define CUTOFF 1000ULL
#ifdef LONG_IS_64BIT
	#if CUTOFF <= 2642245ULL
		#define CUTOFF_CUBE CUTOFF * CUTOFF * CUTOFF
	#else
		#define CUTOFF_CUBE 18446744073709551615ULL
	#endif
#else
	#if CUTOFF <= 1625ULL
		#define CUTOFF_CUBE CUTOFF * CUTOFF * CUTOFF
	#else
		#define CUTOFF_CUBE 4294967295ULL
	#endif
#endif
#if CUTOFF < 661
	#error uissemiprime misconfigured, needs more primes to use uisprime_nosmalldiv.
#endif

	// Remove even numbers. Half of random inputs are caught here.
	if (!(n&1))
		return uisprime(n >> 1);

	// If n is small, simply test up to the cube root; no need for fancy stuff
	ulong lim;
	if (n <= CUTOFF_CUBE) {
		if (n < 9)
			return 0;
		lim = cuberoot(n);

		long p = 0;
		byteptr primepointer = diffptr;
		NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
		for (;;)
		{
			NEXT_PRIME_VIADIFF(p, primepointer);
			if (p > lim)
				break;
			if (n%p == 0)
				return uisprime(n / p);
		}
		
		return !uisprime(n);
	}
	
	// First trial division loop, catches 'easy' numbers.
	// 83% of 'random' odd numbers trapped by this loop for CUTOFF = 661,
	// or more for larger values.
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > CUTOFF)
			break;
		if (n%p == 0)
			return uisprime(n / p);
	}

	// Test for primality. About 27% of 661-rough numbers are caught here.
	if (uisprime_nosmalldiv(n))
		return 0;
	
	// Check for a small prime factor with rho. Catches about 70% of remaining
	// composites, based on testing with CUTOFF = 1000.
	pari_sp ltop = avma;
	GEN fac = pollardbrent(utoipos(n));
	if (fac == NULL) {
		avma = ltop;
	} else if (typ(fac) == t_INT) {
		ulong f = itou(fac);
		avma = ltop;
		return uisprime_nosmalldiv(f) && uisprime_nosmalldiv(n / f);
	} else if (typ(fac) == t_VEC) {
		// TODO: Slight speedup possible by paying attention to format instead
		// of just taking first factor:
		//   "a vector of t_INTs, each triple of successive entries containing
		//   a factor, an exponent (equal to one),  and a factor class (NULL
		//   for unknown or zero for known composite)"
		ulong f = itou(gel(fac, 1));
		avma = ltop;
		return uisprime_nosmalldiv(f) && uisprime_nosmalldiv(n / f);
	}
	
	// Second part of trial division loop: avoids the cube root calculation
	// for numbers with a tiny prime divisor, and allows the use of
	// uisprime_nosmalldiv instead of uisprime.  Neither really matter for
	// hard numbers, but for 'average' numbers the first, at least, is
	// worthwhile.
	lim = cuberoot(n);
	for (;;)
	{
		if (p > lim)
			break;
		if (n%p == 0)
			return uisprime_nosmalldiv(n / p);
		NEXT_PRIME_VIADIFF(p, primepointer);
	}
	
	return 0;	//!uisprime_nosmalldiv(n);
#undef CUTOFF_CUBE
#undef CUTOFF
}


GEN
rad(GEN n)
{
	pari_sp ltop = avma;
	GEN ret;
	if (typ(n) != t_INT)
		pari_err(arither1, "rad");
	if (signe(n) < 0)
		n = negi(n);
	ret = vecprod(gel(Z_factor(n), 1));
	ret = gerepileupto(ltop, ret);
	return ret;
}


// TODO: This is slow and I'm not sure why.
long
isprimepower(GEN n)
{
	pari_sp ltop = avma;
	long ret;
	if (typ(n) != t_INT)
		pari_err(arither1, "isprimepower");
	else if (signe(n) < 1)
		return 0;
	else if (!mod2(n))
		return ispow2(n);
		
	// Deal with small odd values
	ulong nn = itou_or_0(n);
	if (nn)
		return uisprimepower(nn);
	
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p >= 103)
			break;
		if (!smodis(n, p))
		{
			ret = !cmpii(n, powuu(p, Z_lval(n, p)));
			avma = ltop;
			return ret;
		}
	}

	Z_isanypower(n, &n);	// Expensive test!
	ret = isprime(n);
	avma = ltop;
	return ret;
}


// n is assumed to be odd and positive.
long
uisprimepower(ulong n)
{
	// CUTOFF should be at least 89 for best performance. Tests suggest that
	// 200-300 is the best range for 64-bit platforms.
#define CUTOFF 199ULL
	long ret;
	//if (n < 6)
	//	return n > 0;
	//if (!(n&1))
	//	return hamming_word(n) == 1;

	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p >= CUTOFF)
			break;
		if (!(n%p))
		{
			pari_sp ltop = avma;
			ret = !cmpui(n, powuu(p, u_lval(n, p)));
			avma = ltop;
			return ret;
		}
	}
	
	if (n < CUTOFF * CUTOFF * CUTOFF)
	{
#if CUTOFF < 102
		if (n < CUTOFF * CUTOFF || uisprime(n))
#else
		if (n < CUTOFF * CUTOFF || u_IsLucasPsP(n))
#endif
			return 1;
		return uissquare(n);
	}

// These are the primes preceeding the appropriate root of ULONG_MAX.
#ifdef LONG_IS_64BIT
	#define ROOT11 53
	#define ROOT9 137
	#define ROOT8 251
	#define ROOT7 563
	#define ROOT5 7129
	#define ROOT4 65521
#else
	#define ROOT11 7
	#define ROOT9 11
	#define ROOT8 13
	#define ROOT7 23
	#define ROOT5 83
	#define ROOT4 251
#endif

	pari_sp ltop = avma;
#if CUTOFF <= ROOT4
	if (uissquareall(n, &n))
#endif
#if CUTOFF <= ROOT8
	if (uissquareall(n, &n))
#endif
	uissquareall(n, &n);
	
	GEN nn = utoi(n);
#if CUTOFF <= ROOT11
	pari_warn(warner, "isprimepower: cutoff marginal, performace suffers");
	Z_isanypower(nn, &nn);
#endif
#if CUTOFF <= ROOT7
	ulong mask = 7;
	if (n < CUTOFF * CUTOFF * CUTOFF * CUTOFF * CUTOFF * CUTOFF * CUTOFF) {
		if (n < CUTOFF * CUTOFF * CUTOFF * CUTOFF * CUTOFF)
			mask = 1;
		else
			mask = 3;
	}
#elif CUTOFF <= ROOT5
	ulong mask = 3;
	if (n < CUTOFF * CUTOFF * CUTOFF * CUTOFF * CUTOFF)
		mask = 1;
#else
	ulong mask = 1;
#endif
#if CUTOFF <= ROOT9
	if (is_357_power(nn, &nn, &mask))
#endif
	is_357_power(nn, &nn, &mask);
	
	ret = isprime(nn);
	avma = ltop;
	return ret;
#undef CUTOFF
#undef ROOT11
#undef ROOT9
#undef ROOT8
#undef ROOT7
#undef ROOT5
#undef ROOT4
}


// 2-adic valuation of n
INLINE long
valu(ulong n)
{
#if 1
	return n ? __builtin_ctzll(n) : -1;
#else
	if (n == 0)
		return -1;
	long count = 0;
	while (!(n & 1)) {
		n >>= 1;
		count++;
	}
	return count;
#endif
}


long
isPowerful_small(ulong n)
{
	long u = valu(n);
	if (u == 1)
		return 0;
	n >>= u;
	//if (n < 3)
	//	return 1;
	
	long p = 0, lim = usqrtsafe(usqrtsafe(n));
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);
	for (;p <= lim;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (n % p == 0)
		{
			n /= p;
			if (n % p)
				return 0;
			do {
				n /= p;
			} while (n % p == 0);
			lim = usqrtsafe(usqrtsafe(n));
		}
	}
	
	if (n == 1)
		return 1;
	
	// if the input was a powerful number, n the square of a prime or the cube of a prime.
	pari_sp btop = avma;
	long ret = n == 1 || Z_isanypower(utoipos(n), NULL);	// FIXME: Check for correctness and speed, something was wrong here.
	avma = btop;
	return ret;
}


long
isPowerful(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isPowerful");
	ulong nn = itou_or_0(n);
	if (nn)
		return isPowerful_small(nn);

	if (!signe(n) || is_pm1(n))
		return 1;
	if (mod4(n) == 2)
		return 0;
	if (!smodis(n, 3) && smodis(n, 9))
		return 0;
	if (!smodis(n, 5) && smodis(n, 25))
		return 0;
	if (!smodis(n, 7) && smodis(n, 49))
		return 0;
	if (!smodis(n, 11) && smodis(n, 121))
		return 0;
	if (!smodis(n, 13) && smodis(n, 169))
		return 0;
	if (!smodis(n, 17) && smodis(n, 289))
		return 0;

	pari_sp ltop = avma;
	n = shifti(n, -vali(n));
	GEN fac, expo, pr;
	
	while (!is_pm1(n) && !isprime(n)) {
		fac = Z_factor_until(n, shifti(n, -1));	// Find a nontrivial factor -- returns just the factored part
		expo = gel(fac, 2);
		pr = gel(fac, 1);
		long len = glength(expo);
		int i;
		
		for (i = 1; i <= len; i++) {
			GEN q = gel(pr, i);	// The found factor, probably prime
			long e = itos(gel(expo, i));
			GEN nn = n = diviiexact(n, powis(q, e));
			GEN g = gcdii(nn, q);
			if (is_pm1(g)) {
				if (e > 1 || isPowerful(q)) {
					n = nn;
				} else {
					avma = ltop;
					return 0;
				}
			} else {
				// This really shouldn't happen!
				GEN nnn;
				long ee = Z_pvalrem(n, q, &nnn);
				if (ee > 1 || isPowerful(nnn)) {
					n = nnn;
				} else {
					avma = top;
					return 0;
				}
			}
		}
	}
	
	avma = ltop;	// slight abuse of ltop: n is used below
	return is_pm1(n);
}


///////////////////////////////
long
isPowerfulCorrect(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isPowerful");

	if (!signe(n) || is_pm1(n))
		return 1;
	if (mod4(n) == 2)
		return 0;
	if (!smodis(n, 3) && smodis(n, 9))
		return 0;
	if (!smodis(n, 5) && smodis(n, 25))
		return 0;
	if (!smodis(n, 7) && smodis(n, 49))
		return 0;
	if (!smodis(n, 11) && smodis(n, 121))
		return 0;
	if (!smodis(n, 13) && smodis(n, 169))
		return 0;
	if (!smodis(n, 17) && smodis(n, 289))
		return 0;

	pari_sp ltop = avma;
	n = shifti(n, -vali(n));
	GEN fac, expo, pr;
	
	while (!is_pm1(n) && !isprime(n)) {
		fac = Z_factor_until(n, shifti(n, -1));	// Find a nontrivial factor -- returns just the factored part
		expo = gel(fac, 2);
		pr = gel(fac, 1);
		long len = glength(expo);
		int i;
		
		for (i = 1; i <= len; i++) {
			GEN q = gel(pr, i);	// The found factor, probably prime
			long e = itos(gel(expo, i));
			GEN nn = n = diviiexact(n, powis(q, e));
			GEN g = gcdii(nn, q);
			if (is_pm1(g)) {
				if (e > 1 || isPowerful(q)) {
					n = nn;
				} else {
					avma = ltop;
					return 0;
				}
			} else {
				// This really shouldn't happen!
				GEN nnn;
				long ee = Z_pvalrem(n, q, &nnn);
				if (ee > 1 || isPowerful(nnn)) {
					n = nnn;
				} else {
					avma = top;
					return 0;
				}
			}
		}
	}
	
	avma = ltop;	// slight abuse of ltop: n is used below
	return is_pm1(n);
}
///////////////////////////////


long
prp(GEN n, GEN b)
{
	pari_sp ltop = avma;
	if (typ(n) != t_INT)
		pari_err(arither1, "prp");
	if (!b)
		b = gen_2;
	else if (typ(b) != t_INT)
		pari_err(arither1, "prp");
	long ret = gequal1(powgi(gmodulo(b, n), subis(n, 1)));
	avma = ltop;
	return ret;
}


long
sprp(GEN n, GEN b)
{
	pari_sp ltop = avma;
	GEN d;
	long l2;
	if (typ(n) != t_INT)
		pari_err(arither1, "sprp");
	else if (cmpis(n, 3) < 0)
		return cmpis(n, 1) > 0;		// Doesn't like even primes
	if (!b)
		b = gen_2;
	else if (typ(b) != t_INT)
		pari_err(arither1, "sprp");

	d = shifti(n, -1);	// At least 1
	long s = vali(d);
	d = shifti(d, -s);
	s++;
	d = gpow(gmodulo(b, n), d, FAKE_PREC);
	if (gequal1(d))
	{
		avma = ltop;
		return 1;
	}

	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 1; i <= s; i++)
	{
		if (gequalm1(d))
		{
			avma = ltop;
			return 1;
		}
		d = gsqr(d);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &d);
	}
	l2 = gequalm1(d);
	avma = ltop;
	return l2;
}


GEN
sopf(GEN n)
{
	pari_sp ltop = avma;
	GEN f, ret = gen_0;
	long l1;
	if (typ(n) != t_INT)
		pari_err(arither1, "sopf");
	f = Z_factor(n);
	l1 = glength(gel(f, 1));
	long i;
	for (i = 1; i <= l1; ++i)
	{
		ret = addii(ret, gcoeff(f, i, 1));
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
sopfr(GEN n)
{
	pari_sp ltop = avma;
	GEN f, ret = gen_0;
	long l1;
	if (typ(n) != t_INT)
		pari_err(arither1, "sopfr");
	f = Z_factor(n);
	l1 = glength(gel(f, 1));
	long i;
	for (i = 1; i <= l1; ++i)
	{
		ret = addii(ret, mulii(gcoeff(f, i, 1), gcoeff(f, i, 2)));
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


static long smallpr[] = {
	1, 1, 2, 6, 6, 30, 30, 210, 210, 210, 210, 2310, 2310, 30030, 30030, 30030,
	30030, 510510, 510510, 9699690, 9699690, 9699690, 9699690, 223092870,
	223092870, 223092870, 223092870, 223092870, 223092870
#ifdef LONG_IS_64BIT
	, 6469693230, 6469693230, 200560490130, 200560490130, 200560490130,
	200560490130, 200560490130, 200560490130
#endif
};


GEN
gpf(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "gpf");
	if (is_pm1(n))
		return gen_1;
	if (!signe(n))
		return gen_0;
		// My choice of convention: gpf(0) = 0
	
	pari_sp ltop = avma;
	GEN f, ret;
	f = gel(Z_factor(n), 1);
	ret = gel(f, glength(f));
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
prodtree(GEN A, long start, long stop)
{
	pari_sp ltop = avma;
	//pari_sp st_lim = stack_lim(ltop, 1);
	long diff = stop - start;
	if (diff >= 8) {
		diff >>= 1;
		GEN leftprod = prodtree(A, start, start + diff);
		//if (low_stack(st_lim, stack_lim(ltop, 1)))
			leftprod = gerepileupto(ltop, leftprod);
		GEN rightprod = prodtree(A, start + diff + 1, stop);
		//if (low_stack(st_lim, stack_lim(ltop, 1)))
			gerepileall(ltop, 2, &leftprod, &rightprod);
		GEN ret = mulii(leftprod, rightprod);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	
	GEN ret = NEVER_USED, a, b, c, d;
	switch (diff) {
		case 7:
			a = mulss(A[start], A[start+7]);
			b = mulss(A[start+1], A[start+6]);
			c = mulss(A[start+2], A[start+5]);
			d = mulss(A[start+3], A[start+4]);
			ret = mulii(mulii(a, b), mulii(c, d));
			break;
		case 6:
			a = mulss(A[start], A[start+3]);
			b = mulss(A[start+1], A[start+2]);
			c = mulss(A[start+4], A[start+6]);
			ret = mulii(mulii(a, b), mulis(c, A[start+5]));
			break;
		case 5:
			a = mulss(A[start], A[start+5]);
			b = mulss(A[start+1], A[start+4]);
			ret = mulii(mulis(a, A[start+2]), mulis(b, A[start+3]));
			break;
		case 4:
			a = mulss(A[start], A[start+2]);
			b = mulss(A[start+3], A[start+4]);
			ret = mulii(mulis(a, A[start+1]), b);
			break;
		case 3:
			a = mulss(A[start], A[start+3]);
			b = mulss(A[start+1], A[start+2]);
			ret = mulii(a, b);
			break;
		default:
			pari_err(talker, "prodtree passed small argument");
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
primorial(GEN n)
{
	pari_sp ltop = avma;
	long nn = NEVER_USED;
	GEN ret;
	if (typ(n) == t_REAL) {
		nn = itos_or_0(floorr(n));
		avma = ltop;
	} else if (typ(n) == t_INT) {
		nn = itos_or_0(n);
	} else {
		pari_err(arither1, "primorial");
	}
	
	if (signe(n) <= 0)
		return gen_1;
	if (nn > maxprime() || nn == 0)	// nn == 0 if n didn't fit into a word
		pari_err(primer1, n);
	if (nn < 37) {
#ifdef LONG_IS_64BIT
		ret = stoi(smallpr[nn]);
#else
		avma = ltop;
		if (nn < 29)
			ret = stoi(smallpr[nn]);
		else if (nn < 31)
			ret = uu32toi(1, 2174725934);
		else
			ret = uu32toi(46, 2991994514);
#endif		
		return ret;
	}

	ulong primeCount = uprimepi(nn);
	GEN pr = primes_zv(primeCount);
	ret = prodtree(pr, 1, primeCount);
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
lpf(GEN n)
{
	pari_sp ltop = avma;
	GEN res;
	if (typ(n) != t_INT)
		pari_err(arither1, "lpf");
	if (!signe(n))
		return gen_0;	// My choice of convention: lpf(0) = 0
	if (!mod2(n))
		return gen_2;
	if (cmpis(n, 2) < 0)
	{
		if (cmpis(n, -1) >= 0)
			return gen_1;
		n = negi(n);
	}
	pari_sp btop = avma;
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);
	for (;p < 9999; avma = btop)	// TODO: Find appropriate cutoff here
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (!smodis(n, p))
		{
			avma = ltop;
			return stoi(p);
		}
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (!smodis(n, p))
		{
			avma = ltop;
			return stoi(p);
		}
	}
	res = gcoeff(Z_factor(n), 1, 1);	// TODO: Partial factorization?  Tricky to do right...
	res = gerepileupto(ltop, res);
	return res;
}

/******************************************************************************/
/**											 Other arithmetic functions													*/
/******************************************************************************/

long
isTriangular(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isTriangular");
	pari_sp ltop = avma;
	long ret = Z_issquare(addis(shifti(n, 3), 1));
	avma = ltop;
	return ret;
}


long isHexagonal(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isHexagonal");
	if (signe(n) < 1)
		return signe(n) == 0;	// 0 is hexagonal
	pari_sp ltop = avma;
	GEN root;
	long ret = Z_issquareall(addis(shifti(n, 3), 1), &root);
	if (ret)
		ret = mod4(root) == 3;
	avma = ltop;
	return ret;
}


// Approximate speed on a Phenom II at various input sizes:
// size	 cycles		time
// 1e00			110	 40 ns
// 1e10			150	 50 ns
// 1e20			268	100 ns
// 1e30			282	100 ns
// 1e50			297	110 ns
// 1e75			316	110 ns
// 1e100		 322	120 ns
// 1e200		 385	140 ns
// 1e500		 569	200 ns
// 1e1000		941	300 ns
// 1e2000	 1710	600 ns
// 1e5000	 3910	1.4 µs
// 1e10000	7540		3 µs
// 1e20000 15200		5 µs
// 1e50000 37000	 13 µs
long
isFibonacci(GEN n)		/* bool */
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isFibonacci");
	if (!is_bigint(n))
		return isSmallFib(itos(n));
	pari_sp ltop = avma;

	// Good residue classes: 55, 76, 144, 199, 377, 521, 987, 1364, 2584, 3571, 6765, 9349, 17711, 24476, 46368, 64079, 121393, 167761, 317811, 439204, 832040, 1149851, 2178309, 3010349, 5702887, 7881196, 14930352, 20633239, 39088169, ...
	long rem = smodis(n, 17711);
	if (rem & 1) {
		if (rem > 5 && rem != 13 && rem != 21 && rem != 55 && rem != 89 && rem != 233 && rem != 377 && rem != 987 && rem != 1597 && rem != 4181 && rem != 6765 && rem != 15127 && rem != 17567 && rem != 17703)
			return 0;
	} else {
		if (rem > 2 && rem < 17708 && rem != 8 && rem != 34 && rem != 144 && rem != 610 && rem != 2584 && rem != 10946 && rem != 16724 && rem != 17334 && rem != 17656 && rem != 17690)
			return 0;
	}
	GEN k = sqri(n);
	k = addis(mulis(k, 5), 4);	// Multiplication via shifts is slower here (by ~4 cycles)
	long res = Z_issquare(k) || (signe(n) > 0 && Z_issquare(subis(k, 8)));
	avma = ltop;
	return res;
}


static long smallfib[] = {
#ifdef LONG_IS_64BIT
	-7540113804746346429,-2880067194370816120,-1100087778366101931,
	-420196140727489673,-160500643816367088,-61305790721611591,
	-23416728348467685,-8944394323791464,-3416454622906707,
	-1304969544928657,-498454011879264,-190392490709135,-72723460248141,
	-27777890035288,-10610209857723,-4052739537881,-1548008755920,
	-591286729879,-225851433717,-86267571272,-32951280099,-12586269025,
	-4807526976,-1836311903,-701408733,-267914296,-102334155,-39088169,
	-14930352,-5702887,-2178309,-832040,-317811,-121393,-46368,-17711,
	-6765,-2584,-987,-377,-144,-55,-21,-8,-3,-1,0,1,1,2,3,5,8,13,21,34,
	55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657,
	46368,75025,121393,196418,317811,514229,832040,1346269,2178309,
	3524578,5702887,9227465,14930352,24157817,39088169,63245986,
	102334155,165580141,267914296,433494437,701408733,1134903170,
	1836311903,2971215073,4807526976,7778742049,12586269025,20365011074,
	32951280099,53316291173,86267571272,139583862445,225851433717,
	365435296162,591286729879,956722026041,1548008755920,2504730781961,
	4052739537881,6557470319842,10610209857723,17167680177565,
	27777890035288,44945570212853,72723460248141,117669030460994,
	190392490709135,308061521170129,498454011879264,806515533049393,
	1304969544928657,2111485077978050,3416454622906707,5527939700884757,
	8944394323791464,14472334024676221,23416728348467685,
	37889062373143906,61305790721611591,99194853094755497,
	160500643816367088,259695496911122585,420196140727489673,
	679891637638612258,1100087778366101931,1779979416004714189,
	2880067194370816120,4660046610375530309,7540113804746346429
#else
	-1836311903,-701408733,-267914296,-102334155,-39088169,-14930352,
	-5702887,-2178309,-832040,-317811,-121393,-46368,-17711,-6765,-2584,
	-987,-377,-144,-55,-21,-8,-3,-1,0,1,1,2,3,5,8,13,21,34,55,89,144,
	233,377,610,987,1597,2584,4181,6765,10946,17711,28657,46368,75025,
	121393,196418,317811,514229,832040,1346269,2178309,3524578,5702887,
	9227465,14930352,24157817,39088169,63245986,102334155,165580141,
	267914296,433494437,701408733,1134903170,1836311903
#endif
};


long
isSmallFib(long n)
{
	int left = 0;
#ifdef LONG_IS_64BIT
	int right = 139;
#else
	int right = 70;
#endif
	while (left + 1 < right) {
		int mid = (left + right) >> 1;
		if (n < smallfib[mid])
			right = mid;
		else
			left = mid;
	}
	return n == smallfib[left] || n == smallfib[right];
}


static void
lucasmod(ulong n, GEN m, GEN *a, GEN *b)
{
	GEN z, t, zt;
	if (!n) { *a = gen_2; *b = gen_1; return; }
	lucasmod(n >> 1, m, &z, &t);
	zt = mulii(z, t);
	switch(n & 3) {
		case  0: *a = addsi(-2,sqri(z)); *b = addsi(-1,zt); break;
		case  1: *a = addsi(-1,zt);      *b = addsi(2,sqri(t)); break;
		case  2: *a = addsi(2,sqri(z));  *b = addsi(1,zt); break;
		case  3: *a = addsi(1,zt);       *b = addsi(-2,sqri(t));
	}
	*a = modii(*a, m);
	*b = modii(*b, m);
}


GEN
fibomod(long n, GEN m)
{
	if (typ(m) != t_INT)
		pari_err(arither1, "fibomod");
	pari_sp av = avma;
	GEN a, b;
	if (!n) return gen_0;
	lucasmod((ulong)(labs(n)-1), mulis(m, 5), &a, &b);
	a = modii(diviuexact(addii(shifti(a,1),b), 5), m);
	if (n < 0 && !odd(n)) setsigne(a, -1);
	return gerepileuptoint(av, a);
}


static void
lucasmod_tiny(ulong n, ulong m, ulong *a, ulong *b)
{
	ulong z, t;
	if (!n) { *a = 2; *b = 1; return; }
	lucasmod_tiny(n >> 1, m, &z, &t);
	switch(n & 3) {
		case  0: *a = z*z - 2; *b = z*t - 1; break;
		case  1: *a = z*t - 1; *b = t*t + 2; break;
		case  2: *a = z*z + 2; *b = z*t + 1; break;
		case  3: *a = z*t + 1; *b = t*t - 2;
	}
	*a = *a % m;
	*b = *b % m;
}


ulong
fibomod_tiny(long n, ulong m)
{
	ulong a, b;
	if (!n) return 0;
	lucasmod_tiny((ulong)(labs(n)-1), (ulong)(5*m), &a, &b);
	a = (((a << 1) + b) / 5) % m;
	
	return n < 0 && !odd(n) ? -a : a;
}


GEN
fibmod(GEN n, GEN m)
{
	if (typ(n) != t_INT || typ(m) != t_INT)
		pari_err(arither1, "fibmod");
	long nn = itos_or_0(n);
	if (nn) {
		ulong mm = itou_or_0(m);
		if (mm) {
#ifdef LONG_IS_64BIT
			if (mm <= 858993459UL)
#else
			if (mm <= 13107UL)
#endif
				return utoi(fibomod_tiny(nn, mm));
		}
		return fibomod(nn, m);
	}
	if (!signe(n))
		return gen_0;

	GEN f, t, res;
	long l;
	pari_sp ltop = avma;
	if (cmpis(m, 6) < 0)
	{
		if (signe(m) < 0)
			m = negi(m);
		if (!signe(m))
			pari_err(user, "division by 0");
		if (equali1(m))
		{
			avma = ltop;
			return gen_0;
		}
		if (cmpis(m, 2) == 0)
		{
			l = smodis(n, 3) > 0;
			avma = ltop;
			return stoi(l);
		}
		if (cmpis(m, 3) == 0)
		{
			/* 0 1 1 2 0 2 2 1 */
			res = modis(fibo(mod8(n)), 3);
			res = gerepileupto(ltop, res);
			return res;
		}
		if (cmpis(m, 4) == 0)
		{
			/* 0 1 1 2 3 1 */
			res = remi2n(fibo(smodis(n, 6)), 2);
			res = gerepileupto(ltop, res);
			return res;
		}
		if (cmpis(m, 5) == 0)
		{
			/* 0 1 1 2 3 0 3 3 1 4 0 4 4 3 2 0 2 2 4 1 */
			res = modis(fibo(smodis(n, 20)), 5);
			res = gerepileupto(ltop, res);
			return res;
		}
	}
	// Rough estimate of the break-even point.	Yes, it's quite large.
	/*
	if (cmpis(n, 1000) < 0)
	{
		res = modii(fibo(itos(n)), m);
		res = gerepileupto(ltop, res);
		return res;
	}*/
	
	f = Z_factor(m);
	long e;
	t = gmodulo(gen_0, gen_1);
	l = glength(gel(f, 1));
	pari_sp btop = avma;
	long i;
	for (i = 1; i <= l; ++i)
	{
		e = itos(gcoeff(f, i, 2));
		t = chinese(t, gmodulo(fibo(smodis(n, Pisano(itos(gcoeff(f, i, 1)), e))), powis(gcoeff(f, i, 1), e)));
		t = gerepileupto(btop, t);
	}
	res = lift(t);
	res = gerepileupto(ltop, res);
	return res;
}


// Pisano periods of odd primes (with 0s elsewhere)
static long PisanoArr[] = {
	0,8,20,16,0,10,28,0,36,18,0,48,0,0,14,30,0,0,76,0,40,88,0,32,0,0,108,0,0,58,
	60,0,0,136,0,70,148,0,0,78,0,168,0,0,44,0,0,0,196,0,50,208,0,72,108,0,76,0,
	0,0,0,0,0,256,0,130,0,0,276,46,0,0,0,0,148,50,0,0,316,0,0,328,0,336,0,0,348,
	0,0,178,90,0,0,0,0,190,388,0,396,22,0,0,0,0,0,42,0,0,0,0,0,448,0,456,114,0,
	52,0,0,238,240,0,0,0,0,250,0,0,516,0,0,176,0,0,268,270,0,0,556,0,56,568,0,0,
	0,0,588,0,0,0,0,0,0,88,0,310,628,0,636,0,0,0,0,0,0,110,0,0,676,0,0,0,0,232,
	174,0,236,0,0,358,0,0,0,736,0,0,748,0,0,378,0,768,0,0,388,0,0,0,796,0,200,0,
	0,0,408,0,0,0,0,418,84,0,0,0,0,430,868,0,0,438,0,888,0,0,448,0,0,0,916,0,46,
	928,0,936,0,0,0,0,0,478,0,0,0,976,0,490,0,0,0,498,0,1008,0,0,254,0,0,0,0,0,
	26,1048,0,0,0,0,0,0,0,0,90,0,0,1096,0,0,0,0,124,0,0,376,0,0,568,570,0,0,
	1156,0,0,0,0,1176,0,0,1188,0,0,598,600,0,0,1216,0,0,1228,0,1236,206,0,0,0,0,
	0,630,0,0,0,0,640,1288,0,1296,0,0,1308,0,0,658,220,0,0,0,0,0,1348,0,452,0,0,
	1368,0,0,0,138,0,0,0,0,700,0,0,0,118,0,0,0,0,718,0,0,0,1456,0,0,1468,0,0,
	738,0,496,0,0,0,750,0,0,1516,0,380,0,0,0,192,0,1548,0,0,0,0,0,0,1576,0,0,0,
	0,228,0,0,0,0,0,202,270,0,0,0,0,820,1648,0,1656,276,0,0,0,0,838,0,0,0,0,0,0,
	1708,0,1716,78,0,1728,0,0,0,0,0,0,1756,0,176,1768,0,1776,0,0,0,0,0,0,0,0,0,
	1816,0,70,0,0,0,102,0,0,0,0,928,0,0,0,1876,0,470,0,0,1896,0,0,212,0,0,0,0,0,
	0,176,0,970,0,0,652,0,0,1968,0,0,0,198,0,0,1996
};


/* Helper function for fibmod; returns a multiple of the period of Fibonacci numbers mod p^e. */
/* Assumes p is prime and e > 0. */
long
Pisano(long p, long e)
{
	long p1;
	if (p < 999) {
		if (p == 5)
			return (long)upowuu(5, (ulong)e) << 2;
		e--;
		if (p == 2)
			return 3 << e;
		p1 = PisanoArr[p>>1];
		return e ? p1 * (long)upowuu(p, e) : p1;
	}
	long t = p%5;
	if (t == 2 || t == 3)
		p1 = (p + 1) << 1;
	else
		p1 = p - 1;
	return e == 1 ? p1 : p1 * (long)upowuu(p, e - 1);
}


GEN
largestSquareFactor(GEN n)
{
	pari_sp ltop = avma;
	GEN f, res;
	long l1, e, i;
	if (typ(n) != t_INT)
		pari_err(arither1, "largestSquareFactor");
	if (!signe(n))
		return gen_0;

	f = Z_factor(n);
	l1 = glength(gel(f, 1));
	res = gen_1;
	for (i = 1; i <= l1; ++i)
	{
		e = itos(gcoeff(f, i, 2));
		if (e > 1)
			res = mulii(res, e >= 4 ? powis(gcoeff(f, i, 1), e >> 1) : gcoeff(f, i, 1));
	}
	res = gsqr(res);	// Remove this line to instead calculate A000188.
	res = gerepileupto(ltop, res);
	return res;
}


INLINE long
hamming_word(ulong w)
{
	return __builtin_popcountll(w);
#if 0
static long byte_weight[] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
	long sum = 0;
	while (w) {
		//sum += w & 1;
		//w >>= 1;

		//sum++;
		//w &= w - 1;

		sum += byte_weight[w & 255];
		w >>= 8;
	}
	return sum;
#endif
}


long
hamming(GEN n)
{
	pari_sp ltop = avma;
	if (typ(n) != t_INT)
		pari_err(arither1, "hamming");
	if (!signe(n))
		return 0;

	GEN xp = int_MSW(n);
	long lx = lgefint(n);
	ulong u = *xp;
	long sum = 0;
	long i = 3;
	for (; i < lx; ++i)
	{
		sum += hamming_word(u);
		xp=int_precW(xp);
		u = *xp;
	}
	avma = ltop;
	return sum + hamming_word(u);
}


static ulong DS[] ={
	0,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,11,3,4,5,6,7,8,
	9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,
	12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,
	12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,9,10,11,3,4,5,6,7,8,
	9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,
	12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,
	12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,2,3,4,5,6,7,8,9,10,11,3,
	4,5,6,7,8,9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,
	9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,
	9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,
	16,17,18,19,20,3,4,5,6,7,8,9,10,11,12,4,5,6,7,8,9,10,11,12,13,5,6,7,8,9,10,
	11,12,13,14,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,
	12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,
	19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,21,4,5,6,7,8,9,
	10,11,12,13,5,6,7,8,9,10,11,12,13,14,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,
	12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,
	11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,
	18,19,20,21,13,14,15,16,17,18,19,20,21,22,5,6,7,8,9,10,11,12,13,14,6,7,8,9,
	10,11,12,13,14,15,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,
	10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,
	17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,20,21,22,14,
	15,16,17,18,19,20,21,22,23,6,7,8,9,10,11,12,13,14,15,7,8,9,10,11,12,13,14,
	15,16,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,12,13,
	14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,
	21,13,14,15,16,17,18,19,20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,17,18,
	19,20,21,22,23,24,7,8,9,10,11,12,13,14,15,16,8,9,10,11,12,13,14,15,16,17,9,
	10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,13,14,15,16,
	17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,20,21,22,14,
	15,16,17,18,19,20,21,22,23,15,16,17,18,19,20,21,22,23,24,16,17,18,19,20,21,
	22,23,24,25,8,9,10,11,12,13,14,15,16,17,9,10,11,12,13,14,15,16,17,18,10,11,
	12,13,14,15,16,17,18,19,11,12,13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,
	19,20,21,13,14,15,16,17,18,19,20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,
	17,18,19,20,21,22,23,24,16,17,18,19,20,21,22,23,24,25,17,18,19,20,21,22,23,
	24,25,26,9,10,11,12,13,14,15,16,17,18,10,11,12,13,14,15,16,17,18,19,11,12,
	13,14,15,16,17,18,19,20,12,13,14,15,16,17,18,19,20,21,13,14,15,16,17,18,19,
	20,21,22,14,15,16,17,18,19,20,21,22,23,15,16,17,18,19,20,21,22,23,24,16,17,
	18,19,20,21,22,23,24,25,17,18,19,20,21,22,23,24,25,26,18,19,20,21,22,23,24,
	25,26,27
};


// TODO: Binary splitting first, then by thousands.
GEN
dsum(GEN n)	  /* int */
{
	if (typ(n) != t_INT)
		pari_err(typeer, "dsum");
	long nn = itou_or_0(n);
	if (nn)
		return utoipos(dsum_small(nn));
#ifdef LONG_IS_64BIT
	if (lgefint(n) > 4003199668773774)
#else
	if (lgefint(n) > 49540182)
#endif
		pari_err(overflower, "freaking giant number in dsum");
		// TODO: Handle very large numbers that overflow ulong?
	
	pari_sp ltop = avma;
	ulong s = 0;
	GEN t, ret;
	GEN thou = stoi(1000);
	
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while (signe(n) > 0)
	{
		t = divrem(n, thou, -1);
		s += DS[itos(gel(t, 2))];
		n = gel(t, 1);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &n);
	}
	ret = stoi(s);
	ret = gerepileuptoint(ltop, ret);
	return ret;
}


ulong
dsum_small(ulong n)
{
	ulong s = 0;
	while (n)
	{
		s += DS[n % 1000];
		n /= 1000;
	}
	return s;
}


long
istwo(GEN n)
{
	// TODO: Serious improvements available with incremental factoring
	pari_sp ltop = avma;
	GEN f = gen_0;
	long ret, l;
	if (typ(n) != t_INT)
		pari_err(arither1, "istwo");
	
	// Remove factors of 2
	long tmp = vali(n);
	if (tmp)
		n = shifti(n, -tmp);
	
	// Remove negatives
	if (cmpis(n, 3) < 0)
	{
		ret = signe(n) >= 0;
		avma = ltop;
		return ret;
	}

	// End early if possible: numbers that are 3 mod 4 always have some prime factor 3 mod 4 raised to an odd power.
	if (mod4(n) == 3)
		return 0;
	
	f = Z_factor(n);
	l = glength(gel(f, 1));
	long i;
	for (i = 1; i <= l; ++i)
	{
		if (mod4(gcoeff(f, i, 1)) == 3 && mod2(gcoeff(f, i, 2)))
		{
			avma = ltop;
			return 0;
		}
	}
	avma = ltop;
	return 1;
}


GEN
ways2(GEN n)
{
	// TODO: Serious improvements available with incremental factoring
	if (typ(n) != t_INT)
		pari_err(arither1, "ways2");

	pari_sp ltop = avma;
	GEN res = gen_1, f = factor(oddres(n));
	long l1 = glength(gel(f, 1));
	pari_sp btop = avma;
	long i;
	for (i = 1; i <= l1; ++i)
	{
		if (mod4(gcoeff(f, i, 1)) == 3)
		{
			if (mod2(gcoeff(f, i, 2)))
			{
				avma = ltop;
				return gen_0;
			}
		} else {
			res = mulii(res, addis(gcoeff(f, i, 2), 1));
		}
		res = gerepileupto(btop, res);
	}
	res = shifti(addis(res, 1), -1);
	res = gerepileuptoleaf(ltop, res);
	return res;
}


long
isthree(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "isthree");
	if (signe(n) < 0)
		return 0;

	long tmp = vali(n);
	if (tmp & 1)
		return 1;
	pari_sp ltop = avma;
	long l1 = mod8(shifti(n, -tmp)) != 7;
	avma = ltop;
	return l1;
}


INLINE long
uissquare(ulong n)
{
	ulong ignore;
	return uissquareall(n, &ignore);
}


long
sways3s(ulong n)
{
	pari_sp ltop = avma;
	long p1 = (long)usqrtsafe(n);
	long res = 0;
	long k, t, m, p2, tmod4;
	for (k = 0; k <= p1; k++)
	{
		tmod4 = (n - (k&1))&3;
		if (tmod4 == 3)
			continue;
		t = n - k * k;
		if (!istwo(stoi(t)))
			continue;
		avma = ltop;
		p2 = (long)usqrtsafe(t>>1);
		switch (tmod4) {
			case 0:
				for (m = ((k + 1)|1) - 1; m <= p2; m += 2)
					if (uissquare(t - m * m))
						res++;
				break;
			case 1:
				for (m = k; m <= p2; m++)
					if (uissquare(t - m * m))
						res++;
				break;
			case 2:
				for (m = k|1; m <= p2; m += 2)
					if (uissquare(t - m * m))
						res++;
				break;
		}
	}
	return res;
}


GEN
ways3(GEN n)
{
	pari_sp ltop = avma;
	GEN p1, t, res;
	if (typ(n) != t_INT)
		pari_err(typeer, "ways3");
	if (signe(n) < 1)
		return signe(n) ? gen_0 : gen_1;
	
	// ways3(4n) = ways3(n); ways3(4^k(8n + 7)) = 0
	if (!mod4(n)) {
		n = shifti(n, -(vali(n) >> 1) << 1);
	}
	if (mod8(n) == 7) {
		avma = ltop;
		return gen_0;
	}
	ulong small = itou_or_0(n);
	if (small) {
		avma = ltop;
		return stoi(sways3s(small));
	}
	
	p1 = sqrtint(truedivis(n, 3));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k, p3, p4;
	res = gen_0;
	for (k = gen_0; cmpii(k, p1) <= 0; k = addis(k, 1))
	{
		t = subii(n, sqri(k));
		p3 = sqrtint(shifti(t, -1));
{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN m;
		p4 = gen_0;
		for (m = k; cmpii(m, p3) <= 0; m = addis(m, 1))
		{
			if (Z_issquare(subii(t, gsqr(m))))
				p4 = addis(p4, 1);
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 2, &p4, &m);
		}
}
		res = addii(res, p4);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 4, &res, &k, &p3, &p4);
			//gerepileall(btop, 2, &res, &k);	// This is significantly (15%+) slower -- I don't know why.
	}
	res = gerepileupto(ltop, res);
	return res;
}


// Looking for A000523?  Try expi.  A029837 and A070939 are similar.
// See also __builtin_ffsll
GEN
msb(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "msb");
	if (signe(n) < 1) {
		if (signe(n))
			pari_err(talker, "msb: negative argument");	// TODO: What error type to use?
		return gen_0;	// Convention from A053644
	}

	return int2n(expi(n));
}


GEN
fusc(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "fusc");
	if (signe(n) < 1)
		return gen_0;
	pari_sp ltop = avma;
	GEN ret;
	
	// Note: Different constants depending on word size!
#ifdef LONG_IS_64BIT
	long l = lgefint(n);
	if (l < 4 || (l == 4 && *int_MSW(n) <= 13153337344ULL))
#else
	if (cmpii(n, u2toi(9386, 2863311530ULL)) <= 0)
#endif
		ret = utoi(fusc_small(n));
	else
		ret = fusc_large(n);
	ret = gerepileupto(ltop, ret);
	return ret;
}


ulong
fusc_small(GEN n)
{
	ulong a = 1, b = 0;
	pari_sp ltop = avma;
	
	GEN xp = int_LSW(n);
	long lx = lgefint(n);
	ulong u = *xp;
	
	// If n has two words, handle the least-significant one (otherwise it has
	// only one word).
	if (lx > 3)
	{
		if (u & 0x1)
			b += a;
		else
			a += b;
		if (u & 0x2)
			b += a;
		else
			a += b;
		if (u & 0x4)
			b += a;
		else
			a += b;
		if (u & 0x8)
			b += a;
		else
			a += b;
		if (u & 0x10)
			b += a;
		else
			a += b;
		if (u & 0x20)
			b += a;
		else
			a += b;
		if (u & 0x40)
			b += a;
		else
			a += b;
		if (u & 0x80)
			b += a;
		else
			a += b;
		if (u & 0x100)
			b += a;
		else
			a += b;
		if (u & 0x200)
			b += a;
		else
			a += b;
		if (u & 0x400)
			b += a;
		else
			a += b;
		if (u & 0x800)
			b += a;
		else
			a += b;
		if (u & 0x1000)
			b += a;
		else
			a += b;
		if (u & 0x2000)
			b += a;
		else
			a += b;
		if (u & 0x4000)
			b += a;
		else
			a += b;
		if (u & 0x8000)
			b += a;
		else
			a += b;
		if (u & 0x10000)
			b += a;
		else
			a += b;
		if (u & 0x20000)
			b += a;
		else
			a += b;
		if (u & 0x40000)
			b += a;
		else
			a += b;
		if (u & 0x80000)
			b += a;
		else
			a += b;
		if (u & 0x100000)
			b += a;
		else
			a += b;
		if (u & 0x200000)
			b += a;
		else
			a += b;
		if (u & 0x400000)
			b += a;
		else
			a += b;
		if (u & 0x800000)
			b += a;
		else
			a += b;
		if (u & 0x1000000)
			b += a;
		else
			a += b;
		if (u & 0x2000000)
			b += a;
		else
			a += b;
		if (u & 0x4000000)
			b += a;
		else
			a += b;
		if (u & 0x8000000)
			b += a;
		else
			a += b;
		if (u & 0x10000000)
			b += a;
		else
			a += b;
		if (u & 0x20000000)
			b += a;
		else
			a += b;
		if (u & 0x40000000)
			b += a;
		else
			a += b;
		if (u & 0x80000000)
			b += a;
		else
			a += b;
#ifdef LONG_IS_64BIT
		if (u & 0x100000000)
			b += a;
		else
			a += b;
		if (u & 0x200000000)
			b += a;
		else
			a += b;
		if (u & 0x400000000)
			b += a;
		else
			a += b;
		if (u & 0x800000000)
			b += a;
		else
			a += b;
		if (u & 0x1000000000)
			b += a;
		else
			a += b;
		if (u & 0x2000000000)
			b += a;
		else
			a += b;
		if (u & 0x4000000000)
			b += a;
		else
			a += b;
		if (u & 0x8000000000)
			b += a;
		else
			a += b;
		if (u & 0x10000000000)
			b += a;
		else
			a += b;
		if (u & 0x20000000000)
			b += a;
		else
			a += b;
		if (u & 0x40000000000)
			b += a;
		else
			a += b;
		if (u & 0x80000000000)
			b += a;
		else
			a += b;
		if (u & 0x100000000000)
			b += a;
		else
			a += b;
		if (u & 0x200000000000)
			b += a;
		else
			a += b;
		if (u & 0x400000000000)
			b += a;
		else
			a += b;
		if (u & 0x800000000000)
			b += a;
		else
			a += b;
		if (u & 0x1000000000000)
			b += a;
		else
			a += b;
		if (u & 0x2000000000000)
			b += a;
		else
			a += b;
		if (u & 0x4000000000000)
			b += a;
		else
			a += b;
		if (u & 0x8000000000000)
			b += a;
		else
			a += b;
		if (u & 0x10000000000000)
			b += a;
		else
			a += b;
		if (u & 0x20000000000000)
			b += a;
		else
			a += b;
		if (u & 0x40000000000000)
			b += a;
		else
			a += b;
		if (u & 0x80000000000000)
			b += a;
		else
			a += b;
		if (u & 0x100000000000000)
			b += a;
		else
			a += b;
		if (u & 0x200000000000000)
			b += a;
		else
			a += b;
		if (u & 0x400000000000000)
			b += a;
		else
			a += b;
		if (u & 0x800000000000000)
			b += a;
		else
			a += b;
		if (u & 0x1000000000000000)
			b += a;
		else
			a += b;
		if (u & 0x2000000000000000)
			b += a;
		else
			a += b;
		if (u & 0x4000000000000000)
			b += a;
		else
			a += b;
		if (u & 0x8000000000000000)
			b += a;
		else
			a += b;
#endif
		xp=int_nextW(xp);
		u = *xp;
	}
	
	ulong tmp = u;
	while (tmp)
	{
		if (tmp&1)
			b += a;
		else
			a += b;
		tmp >>= 1;
	}
	
	avma = ltop;
	return b;
}


GEN
fusc_large(GEN n)
{
	pari_sp ltop = avma;
	GEN a = gen_1, b = gen_0;
	while (signe(n))
	{
		if (mod2(n))
			b = addii(b, a);
		else
			a = addii(a, b);
		n = shifti(n, -1);	// TODO: Step through number as in fusc_small rather than shifting
	}
	b = gerepileupto(ltop, b);
	return b;
}

/******************************************************************************/
/**												 Other number theory															**/
/******************************************************************************/

GEN
Faulhaber(long e, GEN a)
{
	pari_sp ltop = avma;
	GEN ret = gen_0;
	GEN x = pol_x(0);
	if (!a)
		a = x;
	else if (!gcmpX(a))
		pari_err(typeer, "Faulhaber; must be a variable");
	pari_sp btop = avma;
	long i = 0;
	for (; i <= e; ++i)
	{
		ret = gadd(ret, gmul(gmul(binomialuu(e + 1, i), bernfrac(i)), gpowgs(x, e + 1 - i)));
		ret = gerepileupto(btop, ret);
	}
	ret = gsubstpol(gadd(gdivgs(ret, e + 1), gpowgs(x, e)), x, a);
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
rp(long b)
{
	pari_sp ltop = avma;
	GEN ret;
	GEN B = int2n(b - 1);
	ret = gnextprime(addii(B, randomi(B)));
	ret = gerepileupto(ltop, ret);
	return ret;
}


// FIXME: Take pity and do some numerical analysis here!
ulong
cuberoot(ulong n)
{
	ulong ret = pow(n + 0.5, 1.0/3);
#ifdef LONG_IS_64BIT
	if (n < 100000000000001ULL)
#endif
		return ret;
	if (n >= 18446724184312856125ULL)
		return 2642245ULL;
	if (ret * ret * ret > n) {
		ret--;
		while (ret * ret * ret > n)
			ret--;
		return ret;
	}
	while ((ret + 1) * (ret + 1) * (ret + 1) <= n)
		ret++;
	return ret;
}

GEN
cuberootint(GEN x)
{
	pari_sp ltop = avma;
	long t = typ(x);
	GEN ret = NEVER_USED;
	if (t == t_INT) {
		ret = gfloor(powrfrac(itor(x, lgefint(x)), 1, 3));
	} else if (t == t_REAL) {
		ret = gfloor(powrfrac(x, 1, 3));
		x = gfloor(x);
	} else {
		pari_err(typeer, "cuberootint");
	}

	if (cmpii(powis(ret, 3), x) == 1) {
		ret = subis(ret, 1);
		while (cmpii(powis(ret, 3), x) == 1)
			ret = subis(ret, 1);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	while (cmpii(powis(addis(ret, 1), 3), x) < 1)
		ret = addis(ret, 1);
	ret = gerepileupto(ltop, ret);
	return ret;
}


long
issquarefree_small(ulong n)
{
#define CUTOFF 1627ULL
	long tmp = n&3;
	if (!tmp)
		return 0;
	if (tmp == 2)
		n >>= 1;
	
	long p = 0;
	byteptr primepointer = diffptr;
	NEXT_PRIME_VIADIFF(p, primepointer);	// Skip 2

	// First loop: remove tiny primes, don't calculate cube roots
	// 99.7% of non-squarefree numbers are detected by this loop (or the above)
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > 97)
			break;
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
		}
	}
	
	// Beyond this point, 99.89% of numbers are squarefree.
	ulong last = cuberoot(n);
	ulong last1 = minuu(last, CUTOFF);
	for (;;)
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = cuberoot(n);
			last1 = minuu(last, CUTOFF);
		}
		
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > last1)
			break;
	}

#ifdef LONG_IS_64BIT
	if (n < CUTOFF * CUTOFF * CUTOFF)
#elif CUTOFF <= 1621
	pari_warn(warner, "issquarefree: cutoff marginal, performace suffers");
	if (n < CUTOFF * CUTOFF * CUTOFF)
#endif
	return n == 1 || !uissquare(n);

	// n is at least CUTOFF^3 and is not divisible by any prime under CUTOFF
	if (last < 65536) {	// maxprime() > 65536
		for (;;)
		{
			if (n%p == 0) {
				n /= p;
				if (n%p == 0)
					return 0;
				last = cuberoot(n);
			}
			
			NEXT_PRIME_VIADIFF(p, primepointer);
			if (p > last)
				break;
		}
		return !uissquare(n);
	}
	
	// n is at least 49 bits
	pari_sp ltop = avma;
	GEN f = pollardbrent(utoipos(n));
	if (f == NULL) {
		// Do nothing
	} else if (typ(f) == t_INT) {
		ulong ff = itos(f);
		n /= ff;
		if (!issquarefree_small(ff) || ugcd(ff, n % ff) > 1) {
			avma = ltop;
			return 0;
		}
		last = cuberoot(n);
	} else {	// typ(f) == t_VEC
		long l = lg(f);
		int i = 1;
		while (i < l) {
			/*
			 // Exponent from pollardbrent is guaranteed to equal 1
			if (gel(f[i], i+1) != gen_1) {
				avma = ltop;
				return 0;
			}*/
			ulong ff = itos(gel(f[i], i));
			n /= ff;
			if (!issquarefree_small(ff) || ugcd(ff, n % ff) > 1) {
				avma = ltop;
				return 0;
			}
			i += 3;	// [factor, exponent, class, factor, exponent, class, ...]
		}
		last = cuberoot(n);
	}
	avma = ltop;
	
	// cP(2^80)
	// 750000	72.81
	// 700000	72.26
	// 600000	72.76
	// 500000	72.50 (was: 71.86)
	// 400000	72.41
	// 300000	72.06
	// 250000	72.33
	// 100000	72.51
	if (last > 300000 || last > maxprime())	// TODO: Find good breakover point here
	{
		long ret = Z_issquarefree(stoi(n));
		avma = ltop;
		return ret;
	}
	for (;;)
	{
		if (n%p == 0) {
			n /= p;
			if (n%p == 0)
				return 0;
			last = cuberoot(n);
		}
		
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > last)
			break;
	}
	
	return !uissquare(n);
#undef CUTOFF
}


// TODO: All of the countPowerful functions could be improved greatly
// by using countSquarefree to check the # of squarefree numbers at the point
// that dividing will give exactly 1.	Cost: two invocations of countSquarefree
// at size ~ cuberoot(n).	Savings: ~0.63 cuberoot(n) invocations of
// issquarefree and ~0.38 cuberoot(n) powerings, square roots, and divisions.
// Big win!
ulong
ucountPowerfulu(ulong n)
{
#if 1
	// About 33% faster
	ulong k, breakpoint = cuberoot(n >> 2);
	ulong res = ucountSquarefree(cuberoot(n)) - ucountSquarefree(breakpoint);
#else
	ulong k, breakpoint = cuberoot(n), res = 0;
#endif
	for (k = 1; k <= breakpoint; k++)
		if (issquarefree_small(k))
			res += usqrtsafe(n / k) / k;
	return res;
}


ulong
ucountPowerfuli(GEN n)
{
	pari_sp ltop = avma;
	ulong cube_root = itou(cuberootint(n));
	ulong res = 0, k;
	for (k = 1; k <= cube_root; k++)
		if (issquarefree_small(k)) {
			res += itos(divis(sqrti(divis(n, k)), k));
			//res += itos(sqrti(divii(n, powuu(k, 3))));    // About 35% slower
			avma = ltop;
		}
	return res;
}


/*
// FIXME: Write this, I apparently needed it sometime before but never
// finished coding it.
GEN
listPowerful(GEN lim)
{
	// Should check for memer before attempting further...
	pari_sp ltop = avma;
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "listPowerful");

#ifdef LONG_IS_64BIT
	if (lim > mkintn(4, 227191940, 3022159750, 2605788960, 0))
		pari_err(overflower, "listPowerful");	// around 2^64 or more entries;
												// surely not enough memory
#else
	// 977971126747547528
	if (lim > uu32toi(227701646, 3932178312))
		pari_err(overflower, "listPowerful");	// 2^32 or more entries
#endif
	
}
*/


// FIXME: seems broken, at least for small values, unlike ucountPowerfuli;
// c.f. isPowerful
// Varies as zeta(3/2)/zeta(3) n^1/2 + zeta(2/3)/zeta(2) n^1/3 + o(n^1/6)
// estPowerful(n)=zeta(3/2)/zeta(3)*sqrt(n) + zeta(2/3)/zeta(2)*n^(1/3)
GEN
countPowerful(GEN n)
{
	pari_sp ltop = avma;
	GEN p1, ret;
pari_warn(warner, "Possibly broken, FIXME");
	
	if (lgefint(p1 = gfloor(n)) <= 3) {
		ret = utoi(ucountPowerfulu(itou(p1)));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	// ui, 32-bit: 3909962179575504899
	// ui, 64-bit: ~72047453657149422936422171552392422209
	if (typ(n) == t_INT)
		n = itor(n, lgefint(n));	// Is this the right amount of precision?
	else if (typ(n) != t_REAL)
		pari_err(typeer, "countPowerful");

	//p1 = gfloor(gpow(gadd(n, ghalf), ginv(stoi(3)), FAKE_PREC));
	p1 = gfloor(powrfrac(gadd(n, ghalf), 1, 3));
	n = gfloor(n);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k = gen_0;
	ret = gen_0;
	for (k = gen_1; cmpii(k, p1) <= 0; k = addis(k, 1))
	{
		if (Z_issquarefree(k))
			ret = addii(ret, sqrti(gdivent(n, powis(k, 3))));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &ret, &k);
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


INLINE long
moebiusu(ulong n)
{
	pari_sp ltop = avma;
	long ret = moebius(utoi(n));
	avma = ltop;
	return ret;
}


// TODO: Doesn't really save much time vs. the original.	To improve, a sieve
// would be needed (to calculate the values of moebius faster).
ulong
ucountSquarefree(ulong lim)
{
	ulong b = usqrtsafe(lim >> 1);
	ulong k;
	ulong ret = 0;
	for (k = 1; k <= b; k++)
		ret += moebiusu(k) * (lim / (k * k));
	ulong p3 = usqrtsafe(lim);
	for (k = b + 1; k <= p3; k++)
		ret += moebiusu(k);
	return ret;
}


GEN
countSquarefree(GEN lim)
{
	pari_sp ltop = avma;

	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "countSquarefree");

	GEN b, ret = gen_0, p3 = gen_0, p4 = gen_0;
	b = sqrti(shifti(lim, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN k = gen_0;
	for (k = gen_1; cmpii(k, b) <= 0; k = addis(k, 1))
	{
		// Clean up operations
		ret = addii(ret, mulis(gdivent(lim, gsqr(k)), moebius(k)));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &ret, &k);
	}
	p3 = sqrti(lim);
	btop = avma;
	st_lim = stack_lim(btop, 1);
	p4 = gen_0;
	for (k = addis(b, 1); cmpii(k, p3) <= 0; k = addis(k, 1))
	{
		p4 = addis(p4, moebius(k));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &p4, &k);
	}
	ret = addii(ret, p4);
	ret = gerepileupto(ltop, ret);
	return ret;
}

/******************************************************************************/
/**															 Factoring																	**/
/******************************************************************************/

GEN
Mfactor(GEN p, GEN lim, GEN start)
{
	pari_sp ltop = avma;

	// Check types
	if (typ(p) != t_INT)
		pari_err(arither1, "Mfactor");
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "Mfactor");
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err(arither1, "Mfactor");

	GEN v, k, p1, p2;
	v = cgetg(1, t_VEC);
	if (signe(p) < 1)
		pari_err(talker, "p must be positive");
	if (!(isprime(p)))
		pari_err(talker, "p must be prime");
	if (mod4(p) != 3)
		pari_warn(warner, "p must be a Mersenne exponent equal to 3 mod 4... I think");

	/* Really, only k in [0, 5, 8, 9] mod 12 need to be checked (at last for Mersenne prime exponents?). */
	/* So this should use a quick check, then a full loop with step [10p, 6p, 2p, 6p]. */
	/* Check for a factor of 3 also before starting loop. */
	p2 = shifti(p, 1);
	k = gceil(gdiv(subis(start, 1), p2));
	p1 = addis(mulii(p2, k), 1);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	GEN q = gen_0;
	for (q = p1; cmpii(q, lim) <= 0; q = addii(q, p2))
	{
		if (mod8(q) != 1 && mod8(q) != 7)
			continue;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &q, &v);
		if (!gequalgs(gpow(gmodulsg(2, q), modii(p, subis(q, 1)), FAKE_PREC), 1))
			continue;
		v = concat(v, q);
		long i = 2;
		while (gequal1(gpow(gmodulsg(2, powis(q, i)), modii(p, mulii(subis(q, 1), powis(q, i - 1))), FAKE_PREC)))
		{
			v = concat(v, q);
			i++;
		}
	}
	v = gerepileupto(ltop, v);
	return v;
}


// FIXME: Gives spurious factors of 2 sometimes
GEN
bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start)
{
	pari_sp ltop = avma;

	GEN v = cgetg(1, t_VEC);
	GEN p1 = gen_0, p2 = gen_0;
	if (typ(a) != t_INT || typ(b) != t_INT || typ(c) != t_INT)
		pari_err(arither1, "bigfactor");
	if (!start)
		start = gen_2;
	else if (typ(start) != t_INT)
		pari_err(arither1, "bigfactor");
	if (typ(lim) == t_REAL)
		lim = gfloor(lim);
	else if (typ(lim) != t_INT)
		pari_err(typeer, "bigfactor");
	long lm = itos(lim);
	if (lm > maxprime())
		pari_err(primer1, lim);

	if (signe(b) < 0)
	{
		// TODO: These two should have their formats changed to match that given below.
		if (equali1(a))
		{
			p1 = Z_factor(subsi(1, c));
			p1 = gerepileupto(ltop, p1);
			return p1;
		}
		if (equalim1(a))
		{
			p2 = Z_factor(subii(stoi(1 - 2 * mpodd(b)), c));
			p2 = gerepileupto(ltop, p2);
			return p2;
		}
		/* a^b not in Z */
		pari_err(talker, "not an integer power in bigfactor");
	}
	long p3 = minss(itos(a), lm);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long p = 0;
	byteptr primepointer = diffptr;
	GEN p5 = gen_0;		/* int */
	if (p3 > maxprime())
		pari_err(primer1, stoi(p3));

	// First loop -- deal with small primes
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		if (p > p3)
			break;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &v);
		if (cmpis(gcdii(a, stoi(p)), 1) > 0)
			p5 = b;	// What's the right way to handle this case?	It seems that something more efficient could be done.
		else
			p5 = stoi(smodis(b, p - 1));
		if (!gequal(powgi(gmodulo(a, stoi(p)), p5), c))
			continue;
		v = concat(v, stoi(p));
		long i = 2;
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN p6 = gen_0;
		for(;;)
		{
			if (cmpis(gcdii(a, stoi(p)), 1) > 0)
				p6 = b;
			else
				p6 = modii(b, mulis(powuu(p, i - 1), p - 1));
			if (!gequal(gpow(gmodulo(a, powuu(p, i)), p6, FAKE_PREC), c))
				break;
			v = concat(v, stoi(p));
			i++;
			if (low_stack(st_lim, stack_lim(btop, 1)))
				v = gerepileupto(btop, v);
		}
	}
	//long p4 = p;//itos(gprecprime(gmaxgs(start, p3))) + 1;
	v = gerepileupto(btop, v);
	btop = avma;	// Should this be reset (as currently) or removed?
	st_lim = stack_lim(btop, 1);	// Ditto
	p = 0;
	primepointer = diffptr;		/* bptr */

	//while (p < p4)
	//	NEXT_PRIME_VIADIFF(p, primepointer);

	// Second loop -- most of the work is done here
	for (;;)
	{
		NEXT_PRIME_VIADIFF(p, primepointer);
		//if (p < p4)
		//	continue;
		if (cmpsi(p, lim) > 0)
			break;
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &v);
		if (!gequal(gpowgs(gmodulo(a, stoi(p)), smodis(b, p - 1)), c))
			continue;
		v = concat(v, stoi(p));
		long i = 2;
		while (gequal(gpow(gmodulo(a, powis(stoi(p), i)), modii(b, mulis(powis(stoi(p), i - 1), p - 1)), FAKE_PREC), c))
		{
			v = concat(v, stoi(p));
			i++;
		}
	}
	v = gerepileupto(ltop, v);
	return v;
}


// Does d divide a^b - c?
long
bigdiv(GEN a, GEN b, GEN c, GEN d)
{
	pari_sp ltop = avma;
	long ret;
	if (typ(a) != t_INT || typ(b) != t_INT || typ(c) != t_INT || typ(d) != t_INT)
		pari_err(arither1, "bigdiv");
	
	if (signe(b) < 0)
	{
		if (equali1(a))
		{
			// Does d divide 1 - c?
			ret = !signe(modii(subsi(1, c), d));
			avma = ltop;
			return ret;
		}
		if (equalim1(a))
		{
			// Does d divide (-1)^b - c?
			ret = !signe(modii(subsi(mpodd(b) ? -1 : 1, c), d));
			avma = ltop;
			return ret;
		}
		/* a^b not in Z */
		pari_err(arither1, "bigdiv");
	} else if (!signe(b)) {
		// Does d divide a^0 - c?
		ret = !signe(modii(signe(a) ? subii(a, c) : subsi(1, c), d));
		avma = ltop;
		return ret;
	}
	
	if (signe(d) < 0)
		d = negi(d);
	if (cmpis(d, 2) <= 0)
	{
		if (cmpis(d, 0) == 0)
			pari_err(gdiver, "bigdiv");
		if (equali1(d))
		{
			// Does 1 divide a^b - c?
			avma = ltop;
			return 1;
		}
		
		// Does 2 divide a^b - c?
		ret = !mpodd(subii(a, c));
		avma = ltop;
		return ret;
	}
	
	if (cmpis(gcdii(a, d), 1) > 0)
		ret = gequal(powgi(gmodulo(a, d), b), c);	// Not as slow as it looks
	else
		ret = gequal(powgi(gmodulo(a, d), modii(b, eulerphi(d))), c);
	avma = ltop;
	return ret;
}

/******************************************************************************/
/**													Real and complex functions											**/
/******************************************************************************/

// Convenience function: binary logarithm of x
GEN
log_2(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN ret = NEVER_USED;	// to silence compiler, which doesn't know that pari_err never returns
	switch(typ(x)) {
		case t_INT:
			ret = mplog(itor(x, prec));
			break;
		case t_REAL:
			ret = mplog(x);
			break;
		case t_FRAC:
		case t_COMPLEX:
			ret = mplog(cxcompotor(x, prec));
			break;
		default:
			pari_err(typeer, "lg");
	}
	ret = divrr(ret, mplog2(prec));
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
contfracback(GEN v, GEN terms)
{
	pari_sp ltop = avma;
	GEN x = gen_0;
	long tterms = NEVER_USED;
	if (!terms)
		tterms = glength(v) - 1;
	else if (typ(terms) == t_INT)
		tterms = itos(terms);
	else
		pari_err(typeer, "contfracback");
	x = gcopy(gel(v, tterms + 1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i = 0;
	for (i = tterms; i >= 1; i--)
	{
		x = gadd(gel(v, i), ginv(x));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &x);
	}
	x = gerepileupto(ltop, x);
	return x;
}


GEN
W(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN e, t, w, ep, tmp;
	
	if (typ(x) == t_INT)
		x = itor(x, prec);
	else if (typ(x) != t_REAL)
		pari_err(typeer, "W");
	prec = precision(x);
	
	if (signe(x) <= 0)
	{
		if (!signe(x))
		{
			avma = ltop;
			return real_0(prec);
		}
		pari_sp btop = avma;
		GEN oneOverE = mpexp(stor(-1, prec));
		long c = absr_cmp(x, oneOverE);
		avma = btop;
		if (!c)
			return real_m1(prec);	// otherwise, sometimes sqrt becomes complex
		if (c > 0)	// x < -1/e
			pari_err(talker, "out of range");
	}
	t = real_1(prec);
	
	// Initial approximation for iteration
	if (cmprs(x, 1) < 0)
	{
		// t = 1 already, might as well use it for this calculation
		tmp = sqrtr(addrs(mulrr(shiftr(mpexp(t), 1), x), 2));
		w = subrs(mulrr(tmp, subir(gen_1, mulrr(tmp, addrr(invr(stor(3, prec)), divrs(mulrs(tmp, 11), 72))))), 1);
		// tmp = sqrt(2e * x + 2)
		// w = tmp * (1 - tmp * (1/3 + 11/72 * tmp)) - 1
		if (precision(w) < prec) {
			w = rtor(w, prec);
		}
	} else {
		w = mplog(x);	// Faster than the better approximation log + log log
	}
	if (cmprs(x, 3) > 0)
		w = subrr(w, mplog(w));

	ep = mulrr(eps(prec), addsr(1, absr(w)));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);

	while (cmprr(absr(t), ep) > 0)
	{
		// Halley loop
		e = mpexp(w);
		t = subrr(mulrr(w, e), x);
//pari_printf("	%Ps (off by %Ps)\n", w, t);
		if (cmprr(absr(t), ep) <= 0)
			break;	// Stops calculation when answer is already very close, to avoid division by 0
		tmp = addrs(w, 1);
		t = divrr(t, subrr(mulrr(e, tmp), divrr(mulrr(shiftr(addrs(w, 2), -1), t), tmp)));
		//if (!cmprs(tmp, 0))
		//	pari_printf("tmp is 0");
		//GEN dividebythis = subrr(mulrr(e, tmp), divrr(mulrr(shiftr(addrs(w, 2), -1), t), tmp));
		//if (!cmprs(dividebythis, 0))
		//	pari_printf("dividebythis is 0");
		//t = divrr(t, dividebythis);
		w = subrr(w, t);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &w);
	}

	w = gerepileupto(ltop, w);
if (precision(w) < prec)
	pari_warn(warner, "precision loss");
	return w;
}


/******************************************************************************/
/**														 Convenience																	**/
/******************************************************************************/

GEN
vecsum(GEN v)
{
	pari_sp ltop = avma;
	GEN p2 = gen_0;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "vecsum");
	long l1 = lg(v);
	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		p2 = gadd(p2, gel(v, i));
		p2 = gerepileupto(btop, p2);
	}
	p2 = gerepileupto(ltop, p2);
	return p2;
}


// TODO: Binary splitting for smaller subproducts.
GEN
vecprod(GEN v)
{
	pari_sp ltop = avma;
	GEN p2 = gen_1;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "vecprod");
	long l1 = lg(v);

	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		p2 = gmul(p2, gel(v, i));
		p2 = gerepileupto(btop, p2);
	}
	p2 = gerepileupto(ltop, p2);
	return p2;
}


GEN
vecgcd(GEN v)
{
	pari_sp ltop = avma;
	GEN l = gen_0;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "veclcm");
	long l1 = lg(v);
	pari_sp btop = avma;
	long i;
	for (i = 1; i < l1; ++i)
	{
		l = ggcd(l, gel(v, i));
		l = gerepileupto(btop, l);
	}
	l = gerepileupto(ltop, l);
	return l;
}


GEN
veclcm(GEN v)
{
	pari_sp ltop = avma;
	GEN l = gen_1;
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "veclcm");
	long l1 = lg(v);
	long i;
	for (i = 1; i < l1; ++i)
	{
		l = glcm(l, gel(v, i));
		l = gerepileupto(ltop, l);
	}
	return l;
}


GEN
oddres(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(typeer, "oddres");
	long v = vali(n);
	return v ? shifti(n, -v) : icopy(n);
}


// Faster than
//	 !cmpii(n, int2n(vali(n)))
//	 !cmpis(shifti(n, -vali(n)), 1)
//	 expi(n) == vali(n)
//	 hamming(n) == 1
// even for single-word values, and much faster for multiword values.
// If all you have is a word, you can just use n & !(n & (n-1)).
long
ispow2(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "ispow2");
	if (signe(n) < 1)
		return 0;
	pari_sp ltop = avma;
	
	GEN xp = int_LSW(n);
	long lx = lgefint(n);
	ulong u = *xp;
	long i = 3;
	for (; i < lx; ++i)
	{
		if (u != 0)
		{
			avma = ltop;
			return 0;
		}
		xp = int_nextW(xp);
		u = *xp;
	}
	avma = ltop;
	//return hamming_word(u) == 1;	// 14% slower
	return !(u & (u-1));
}


// TODO (possibly): give 32- and 64-bit specific code if it would differ, since
// sometimes code needs to be split anyway (or only one is relevant, etc.).
// TODO: Handle negatives, and possibly non-integer types?
void
toC(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(typeer, "toC");
	if (cmpis(n, 3) < 0)
	{
		if (cmpis(n, 2) == 0)
			pari_printf("gen_2\n");
		else if (equali1(n))
			pari_printf("gen_1\n");
		else if (signe(n) == 0)
			pari_printf("gen_0\n");
		else if (equalim1(n))
			pari_printf("gen_m1\n");
		else
			pari_err(alarmer, "can't handle negatives yet"); // white lie
		return;
	}
	if (ispow2(n)) {
		pari_printf("int2n(%ld)\n", expi(n));
		return;
	}
	
	pari_sp ltop = avma;
	GEN t;
	// words: number of 32-bit words in n.
#ifdef LONG_IS_64BIT
	long words = (lgefint(n) - 2) << 1;
	if ((ulong)*int_MSW(n) <= 0xFFFFFFFF)
		words--;
#else
	long words = lgefint(n) - 2;
#endif

	if (words == 1)
	{
		pari_printf("utoipos(%Ps)\n", n);
		avma = ltop;
		return;
	}
	if (words == 2)
	{
		pari_printf("uu32toi(%Ps, %Ps)\n", shifti(n, -32), remi2n(n, 32));
		avma = ltop;
		return;
	}

	// Large numbers
	// If efficiency mattered, walking through the binary representation
	// would be far more efficient.
	pari_printf("mkintn(%Ps", stoi(words));
	long i = words - 1;
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	for (; i >= 0; i--)
	{
		t = shifti(n, -(i * 32));
		pari_printf(", %Ps", t);
		n = subii(n, shifti(t, i * 32));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &n);
	}
	pari_printf(")\n");
	avma = ltop;
	return;
}


GEN
digits(GEN x)
{
	pari_sp ltop = avma;
	long s = sizedigit(x) - 1;
	if (gcmp(x, powis(stoi(10), s)) >= 0)
		s++;
	avma = ltop;
	return stoi(s);
}


GEN
eps(long prec)
{
	GEN ret = real_1(DEFAULTPREC);
	setexpo(ret, 1 - bit_accuracy(prec));
	return ret;
}

/******************************************************************************/
/**																	 I/O																		**/
/******************************************************************************/

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
		pari_err(openfiler, "input", name);
#else
		if (strlen(name) == 11 && name[0] == 'b' && name[1] >= '0' && name[1] <= '9' && name[2] >= '0' && name[2] <= '9' && name[3] >= '0' && name[3] <= '9' && name[4] >= '0' && name[4] <= '9' && name[5] >= '0' && name[5] <= '9' && name[6] >= '0' && name[6] <= '9' && name[7] == '.' && name[8] == 't' && name[9] == 'x' && name[10] == 't') {
			char command[65];
			sprintf(command, "wget http://oeis.org/A%c%c%c%c%c%c/%s", name[1],name[2],name[3],name[4],name[5],name[6],name);
			int result = system(command);
			if (result == -1)
				pari_warn(talker, "Download failed.");
		}
		f = fopen(name, "r");
		if (!f)
			pari_err(openfiler, "input", name);
#endif
	}
	
	GEN v = vectrunc_init(MAX_VECLEN + 1);
	char line[MAX_LINELEN];
	int i = 0;
	while(fgets(line, MAX_LINELEN, f) != NULL) {
		if (strlen(line) > MAX_LINELEN - 5)
			pari_err(talker, "Maximum line length exceeded; b-file probably not valid");
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
	// If no v is given, fine; but if it is it must be a vector.
	// Should this use is_vec_t(typ(v)) to allow t_COL as well?
	if (v && typ(v) != t_VEC)
		pari_err(typeer, "bfile");
	if (!offset)
		cur = 0;
	else if (typ(offset) == t_INT)
		cur = itos(offset) - 1;
	else
		pari_err(typeer, "bfile");
	
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
		pari_err(typeer, "bfile");
	}
	
	char* filename = GSTR(name);
	if (!v)
		return bfilein(filename);
	FILE *f = fopen(filename, "r");
	if (f) {
		pari_warn(warner, "File `%Ps' already exists. Moving to bfile.old...", name);
		fclose(f);
		rename(filename, "bfile.old");
	}
	
	f = fopen(filename, "w+");
	long l1 = lg(v);
	long i;
	for (i = 1; i < l1; ++i)
	{
		GEN e = gel(v, i);
		if (typ(e) != t_INT)
			pari_err(typeer, "bfile");
		if (cmpis(digits(e), 1000) > 0)
		{
			pari_warn(warner, "Next term has %Ps digits; exiting.\n", digits(e));
			break;
		}

		char *num = GENtostr(e);
		fprintf(f, "%ld %s\n", ++cur, num);
		pari_free(num);
	}
	fclose(f);
	if(offset)
		pari_printf("A%Ps: Terms %Ps..%ld written to %s\n", Anum, offset, cur, filename);
	else
		pari_printf("A%Ps: Terms 1..%ld written to %s\n", Anum, cur, filename);
	avma = ltop;
	return gnil;
}


GEN
fnice(GEN n)
{
	pari_sp ltop = avma;
	GEN f, s, s1, p1 = gen_0;
	long l2;
	if (typ(n) != t_INT)
		pari_err(arither1, "fnice");
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


GEN
tonice(GEN o, long prec)
{
	pari_sp ltop = avma;
	GEN s = gen_0, t = gen_0, t1 = gen_0, v = gen_0, p1 = gen_0, p2 = gen_0, p3 = gen_0;
	GEN p4 = gen_0;		/* genstr */
	if (typ(o) == t_POL)
	{
		/* Can handle only single-variable polynomials. */
		t = content(o);
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
				pari_err(user, "tonice: cannot display multivariate polynomials!");
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
		pari_err(typeer, "initial");
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
		pari_err(typeer, "medial");
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
		pari_err(typeer, "monomialnice");
	if (typ(degree) != t_INT)
		pari_err(typeer, "monomialnice");
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

/******************************************************************************/
/**															Set stuff																	 **/
/******************************************************************************/

// TODO: check if a = b and handle more efficiently (this is a common case)
GEN
sumset(GEN a, GEN b)		/* vecsmall */
{
	pari_sp ltop = avma;
	GEN c, p2, p4;
	long l1;
	long l3;
	l1 = glength(a)*glength(b);
	long l5;
	p2 = cgetg(l1+1, t_VEC);
	for (l5 = 1; l5 <= l1; ++l5)
		gel(p2, l5) = gen_0;
	c = p2;
	l3 = glength(a);
{
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i, l6;
	for (i = 1; i <= l3; ++i)
	{
		l6 = glength(b);
{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		long j;
		for (j = 1; j <= l6; ++j)
		{
			gel(c, ((i - 1)*glength(b)) + j) = gadd(gel(a, i), gel(b, j));
			if (low_stack(st_lim, stack_lim(btop, 1)))
				c = gerepilecopy(btop, c);
		}
}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			c = gerepilecopy(btop, c);
	}
}
	p4 = vecsort0(c, NULL, 8);
	p4 = gerepileuptoleaf(ltop, p4);
	return p4;
}


GEN
diffset(GEN a, GEN b)		/* vecsmall */
{
	pari_sp ltop = avma;
	GEN c = gen_0;
	long l1;
	GEN p2 = gen_0;		/* vec */
	long l3;
	GEN p4 = gen_0;		/* vecsmall */
	l1 = glength(a)*glength(b);
	{
		long l5;
		p2 = cgetg(l1+1, t_VEC);
		for (l5 = 1; l5 <= l1; ++l5)
			gel(p2, l5) = gen_0;
	}
	c = p2;
	l3 = glength(a);
	{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		long i, l6;
		for (i = 1; i <= l3; ++i)
		{
			l6 = glength(b);
			{
				pari_sp btop = avma, st_lim = stack_lim(btop, 1);
				long j;
				for (j = 1; j <= l6; ++j)
				{
					gel(c, ((i - 1)*glength(b)) + j) = gsub(gel(a, i), gel(b, j));
					if (low_stack(st_lim, stack_lim(btop, 1)))
						c = gerepilecopy(btop, c);
				}
			}
			if (low_stack(st_lim, stack_lim(btop, 1)))
				c = gerepilecopy(btop, c);
		}
	}
	p4 = vecsort0(geval(gtoset(c)), NULL, 0);
	p4 = gerepileuptoleaf(ltop, p4);
	return p4;
}

/******************************************************************************/
/**															 Statistics																 **/
/******************************************************************************/

long
infinite(GEN x)
{
	if (typ(x) != t_VEC || glength(x) != 1)
		return 0;
	GEN e = gel(x, 1);	// Nothing is created, so no garbage... right?
	
	// If e is gen_0, then is_pm1 is unpredicatable, but that doesn't matter
	// because then both branches return 0.
	return (typ(e) == t_INT && is_pm1(e)) ? signe(e) : 0;
}


long
isExtendedReal(GEN x)
{
	long t = typ(x);
	if (t == t_INT || t == t_FRAC || t == t_REAL)
		return 1;
	return infinite(x);
}


// FIXME: Infinities are broken?
GEN
normd(GEN a, GEN b, long prec)
{
	if (!isExtendedReal(a) || !isExtendedReal(b))
		pari_err(talker, "incorrect endpoint in normd");
	pari_sp ltop = avma;
	long tmp;
	GEN ret = NEVER_USED;
	
	/* Infinities */
	if ((tmp = infinite(a)))	// Assignment and test-if-0
	{
pari_warn(warner, "Doesn't work properly with infinities");
		if (tmp < 0)	// (-oo, b)
		{
			tmp = infinite(b);
			if (tmp > 0)
				ret = gen_1;
			else if (tmp < 0)
				ret = gen_0;
			else
				ret = gdivgs(mpneg(gerfc(mpdiv(b, gsqrt(gen_2, prec)), prec)), 2);
		} else {		// (oo, b)
			if (infinite(b) == 1)
				ret = gen_1;
			else
				pari_err(talker, "incorrect endpoint in normd");
		}
	} else if ((tmp = infinite(b))) {	// Assignment and test-if-0
pari_warn(warner, "Doesn't work properly with infinities");
		if (tmp < 0)	// (a, -oo)
			pari_err(talker, "incorrect endpoint in normd");
		ret = gdivgs(gerfc(mpdiv(a, gsqrt(gen_2, prec)), prec), 2);
	} else {
		GEN root2 = gsqrt(gen_2, prec);
		ret = gdivgs(gsub(gerfc(gdiv(a, root2), prec), gerfc(gdiv(b, root2), prec)), 2);
	}
	ret = gerepileupto(ltop, ret);
	return ret;
}


// Use the Box-Muller transform to generate random normal variables. Caches
// values, so multiple calls at the same precision are fast.
GEN
rnormal(long prec)
{
	if (rnormal_cached) {
		if (precision(rnormal_cached) != prec) {
			gunclone(rnormal_cached);
		} else {
			GEN ret = gcopy(rnormal_cached);
			gunclone(rnormal_cached);
			rnormal_cached = 0;
			return ret;
		}
	}
	pari_sp ltop = avma;
	GEN u1, u2, ret, outside, inside, cos_inside;
	u1 = randomr(prec);
	u2 = randomr(prec);
	outside = sqrtr_abs(shiftr(mplog(u1), 1));
	inside = mulrr(shiftr(mppi(prec), 1), u2);
	cos_inside = mpcos(inside);
	
	ret = mulrr(outside, cos_inside);
	rnormal_cached = gclone(ret);	// Cache for later use
	ret = mulrr(outside, cos_inside);
	ret = gerepileupto(ltop, ret);
		return ret;
}

/******************************************************************************/
/**												Verbose monstrosities														 **/
/******************************************************************************/

void
pBounds(GEN n, GEN verbose, long prec)
{
	pari_sp ltop = avma;
	GEN lower = gen_0, upper = gen_0, appx = gen_0, l = gen_0, ll = gen_0;
	if (!verbose)
		verbose = gen_0;
	if (gcmpgs(n, 6548) < 0)
	{
		if (gcmpgs(n, 1) < 0)
			pari_printf("There are no negative primes.\n");
		else
		{
			n = gfloor(n);
			pari_printf("p_%Ps = %Ps (exactly)\n", n, prime(gtos(n)));
		}
		avma = ltop;
		return;
	}
	n = gfloor(n);
	l = glog(n, prec);
	ll = glog(l, prec);
	lower = gmul(n, gsubgs(gadd(l, ll), 1));
	/* Dusart, n >= 2 */
	if (gcmpgs(n, 13196) > 0)
		lower = gmul(n, gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)), gdiv(strtor("2.25", prec), l)));
	/* Dusart, n >= 2 */
	appx = gmul(n, gadd(gadd(gsub(gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)), gdivsg(2, l)), gdiv(gdivgs(gsqr(ll), 2), gsqr(l))), gdiv(gmulsg(3, ll), gsqr(l))), gdiv(gdivgs(stoi(11), 2), gsqr(l))));
	/* + O(ll^3/l^3) */
	upper = gmul(n, gadd(l, ll));
	/* ?, n >= 6 */
	if (gcmpgs(n, 27076) >= 0)
		upper = gmul(n, gsub(gadd(gsubgs(gadd(l, ll), 1), gdiv(ll, l)), gdiv(strtor("1.8", prec), l)));
	/* Dusart, n >= 27076 */
	if (gcmpgs(n, 39017) >= 0)
		upper = gmin(upper, gmul(n, gsub(gadd(l, ll), strtor(".9484", prec))));
	/* Dusart, n >= 39017 */
	lower = gceil(lower);
	upper = gfloor(upper);
	pari_printf("%Ps (lower bound)\n", lower);
	if ((gcmp(lower, appx) < 0) && (gcmp(appx, upper) < 0))
		pari_printf("%Ps (approximate)\n", appx);
	pari_printf("%Ps (upper bound)\n", upper);
	if (!gequal0(verbose))
	{
		pari_printf("\nPierre Dusart, 'Autour de la fonction qui compte le nombre de nombres\n");
		pari_printf("premiers', doctoral thesis for l'Université de Limoges (1998).\n");
		if ((gcmp(lower, appx) < 0) && (gcmp(appx, upper) < 0))
		{
			pari_printf("Ernest Cesàro (1894). \"Sur une formule empirique de M. Pervouchine\". Comptes\n");
			pari_printf("rendus hebdomadaires des séances de l'Académie des sciences 119, pp. 848-849.\n");
		}
	}
	avma = ltop;
	return;
}

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
		pari_err(typeer, "forodd");
	if (typ(b) == t_REAL)
		b = gfloor(b);
	else if (typ(b) != t_INT)
		pari_err(typeer, "forodd");
	
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
		pari_err(notpoler, "sumformal; can only handle polynomials, not arbitrary functions");

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
// TODO: Use gtofp for conversions?
void
forbigprime(GEN ga, GEN gb, GEN code)
{
	pari_sp av = avma;
	long t = typ(ga);
	if (t == t_REAL || t == t_FRAC)
		ga = gceil(ga);
	else if (t != t_INT)
		pari_err(typeer, "forbigprime");

	t = typ(gb);
	if (t == t_REAL || t == t_FRAC)
		gb = gfloor(gb);
	else if (t != t_INT)
		pari_err(typeer, "forbigprime");

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
			pari_err(talker, "Only works for single-word integers.");
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
	if (maxprime() < 4294967296ULL && b > maxprime() * maxprime())
		pari_err(primer1, stoi(usqrtsafe(b)));
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
	ulong chunk = maxuu(0x100000, usqrtsafe(b) * __builtin_ffsll(b));
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
forthinprime(ulong a, ulong b, GEN code)
{
	pari_sp ltop = avma;
	long p[] = {evaltyp(t_INT)|_evallg(3), evalsigne(1)|evallgefint(3), 0};	// A GEN representing a prime to be passed to the code; its value is in p[2].
	push_lex((GEN)p,code);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);

	a |= 1;	// Make a odd (round up)
	b = (b-1)|1;	// Make b odd (round down)
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

/******************************************************************************/
/**									 Works-in-progress / limited utility										**/
/******************************************************************************/

GEN
checkVDW(GEN vv, GEN verbose)
{
	pari_sp ltop = avma;
	GEN r = gen_0, k = gen_0, s = gen_0, c = gen_0;
	if (!verbose)
		verbose = gen_1;
	r = stoi(glength(vv));
	c = cgetg(1, t_VEC);
	{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN i = gen_0;
		for (i = gen_1; gcmp(i, r) <= 0; i = gaddgs(i, 1))
		{
			if (!gequal(gel(vv, gtos(i)), vecsort0(gel(vv, gtos(i)), NULL, 8)))
			{
				if (!gequal0(verbose))
					pari_printf("Not a partition: numbers repeated in color %Ps.\n", i);
				avma = ltop;
				return gen_0;
			}
			s = gaddgs(s, glength(gel(vv, gtos(i))));
			c = concat(c, gel(vv, gtos(i)));
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 3, &i, &s, &c);
		}
	}
	c = vecsort0(c, NULL, 8);
	if (gcmpgs(gel(c, 1), 1) < 0)
	{
		if (!gequal0(verbose))
			pari_printf("Not a natural number partition: negative numbers in a member array.\n");
		avma = ltop;
		return gen_0;
	}
	if ((gcmpgs(gel(c, 1), 1) > 0) || (gcmpgs(gel(c, glength(c)), glength(c)) > 0))
	{
		if (!gequal0(verbose))
			pari_printf("Not a partition of an initial segment: not all numbers {1, 2, ..., %Ps} appear.\n", gel(c, glength(c)));
		avma = ltop;
		return gen_0;
	}
	k = gaddgs(longestProgression(gel(vv, 1)), 1);
	{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN i = gen_0;
		for (i = gen_2; gcmp(i, r) <= 0; i = gaddgs(i, 1))
		{
			k = gmax(k, gaddgs(longestProgression(gel(vv, gtos(i))), 1));
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 2, &i, &k);
		}
	}
	if (!gequal0(verbose))
		pari_printf("W(%Ps, %Ps) > %ld\n", r, k, glength(c));
	k = gerepileupto(ltop, k);
	return k;
}

GEN
longestProgression(GEN v)
{
	pari_sp ltop = avma;
	GEN r = gen_0, s = gen_0, t = gen_0, d = gen_0;
	long l1, l2;
	if (glength(v) < 3)
	{
		l1 = glength(v);
		avma = ltop;
		return stoi(l1);
	}
	s = gtoset(v);
	l2 = glength(v) - 1;
	{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN i = gen_0, p3 = gen_0;
		long l4;
		for (i = gen_1; gcmpgs(i, l2) <= 0; i = gaddgs(i, 1))
		{
			p3 = gaddgs(i, 1);
			l4 = glength(v);
			{
				pari_sp btop = avma, st_lim = stack_lim(btop, 1);
				GEN j = gen_0;
				for (j = p3; gcmpgs(j, l4) <= 0; j = gaddgs(j, 1))
				{
					t = gen_2;
					d = gsub(gel(v, gtos(j)), gel(v, gtos(i)));
					{
						pari_sp btop = avma;
						while (setsearch(s, gadd(gel(v, gtos(i)), gmul(d, t)), 0))
						{
							t = gaddgs(t, 1);
							t = gerepileupto(btop, t);
						}
					}
					r = gmax(r, t);
					if (low_stack(st_lim, stack_lim(btop, 1)))
						gerepileall(btop, 4, &j, &t, &d, &r);
				}
			}
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 5, &i, &p3, &t, &d, &r);
		}
	}
	r = gerepileupto(ltop, r);
	return r;
}

GEN
longestProgression1(GEN v)
{
	pari_sp ltop = avma;
	GEN Lstar = gen_2, L = gen_0, i = gen_0, k = gen_0, tmp = gen_0;
	long l1, l2;
	GEN p3 = gen_0;		/* vec */
	long l4, l5;
	l1 = glength(v);
	l2 = glength(v);
	{
		long l6, l7;
		p3 = cgetg(l1+1, t_MAT);
		for (l7 = 1; l7 <= l1; ++l7)
		{
			gel(p3, l7) = cgetg(l2+1, t_COL);
			for (l6 = 1; l6 <= l2; ++l6)
				gcoeff(p3, l6, l7) = gen_0;
		}
	}
	L = p3;
	if (glength(v) < 3)
	{
		l4 = glength(v);
		avma = ltop;
		return stoi(l4);
	}
	l5 = glength(v) - 1;
	{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN j = gen_0;
		long l8 = -1 > 0;		/* bool */
		for (j = stoi(l5); l8?gcmpgs(j, 1) <= 0:gcmpgs(j, 1) >= 0; j = gaddgs(j, -1))
		{
			i = gsubgs(j, 1);
			k = gaddgs(j, 1);
			{
				pari_sp btop = avma, st_lim = stack_lim(btop, 1);
				while ((gcmpgs(i, 0) > 0) && (gcmpgs(k, glength(v)) <= 0))
				{
					tmp = gsub(gadd(gel(v, gtos(i)), gel(v, gtos(k))), gmulsg(2, gel(v, gtos(j))));
					if (gcmpgs(tmp, 0) < 0)
						k = gaddgs(k, 1);
					else
					{
						if (gcmpgs(tmp, 0) > 0)
						{
							gcoeff(L, gtos(i), gtos(j)) = gen_2;
							i = gsubgs(i, 1);
						}
						else
						{
							gcoeff(L, gtos(i), gtos(j)) = gaddgs(gcoeff(L, gtos(j), gtos(k)), 1);
							Lstar = gmax(Lstar, gcoeff(L, gtos(i), gtos(j)));
							i = gsubgs(i, 1);
							k = gaddgs(k, 1);
						}
					}
					if (low_stack(st_lim, stack_lim(btop, 1)))
						gerepileall(btop, 5, &tmp, &k, &L, &i, &Lstar);
				}
			}
			{
				pari_sp btop = avma, st_lim = stack_lim(btop, 1);
				while (gcmpgs(i, 0) > 0)
				{
					gcoeff(L, gtos(i), gtos(j)) = gen_2;
					i = gsubgs(i, 1);
					if (low_stack(st_lim, stack_lim(btop, 1)))
						gerepileall(btop, 2, &L, &i);
				}
			}
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 6, &j, &i, &k, &tmp, &L, &Lstar);
		}
	}
	Lstar = gerepileupto(ltop, Lstar);
	return Lstar;
}


////////////////////////////////////////////////////////////////////////////////////// Newcomers

long
checkmult(GEN v, long verbose)
{
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkmult");

	// I suppose [] is a multiplicative sequence... if nothing else guard
	// against checking v[1] below.
	long n, l1 = lg(v), l2;
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	pari_sp ltop = avma;
	GEN f, target;

	// Require all arguments to be integers
	for (n = 2; n < 6 && n < l1; ++n)
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
			
	for (n = 6; n < l1; ++n) {
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
		if (uisprimepower(n))
			continue;
		f = Z_factor(stoi(n));
		l2 = glength(gel(f, 1));
		long i;
		
		// Set target = prod v[p^e] for each p^e || n
		target = gen_1;
		for (i = 1; i <= l2; ++i)
			target = mulii(target, gel(v, itos(powii(gcoeff(f, i, 1), gcoeff(f, i, 2)))));
		
		// If v[n] is not equal to the target, the sequence is not multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not multiplicative at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	avma = ltop;
	return 1;
}


long
checkcmult(GEN v, long verbose)
{
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkmult");

	// I suppose [] is a (completely) multiplicative sequence... if nothing else
	// guard against checking v[1] below.
	long n, l1 = lg(v), l2;
	if (l1 == 1)
		return 1;
	
	// At the moment v[1] must be equal to 1; definitions vary but this seems
	// sensible.
	if (!gequalgs(gel(v, 1), 1))
		return 0;
		
	pari_sp ltop = avma;
	GEN f, target;

	// Require all arguments to be integers
	for (n = 2; n < 4 && n < l1; ++n)
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
			
	for (n = 4; n < l1; ++n) {
		if (typ(gel(v, n)) != t_INT)
			pari_err(arither1, "checkmult");
		if (uisprime(n))
			continue;
		f = Z_factor(stoi(n));
		l2 = glength(gel(f, 1));
		long i;
		
		// Set target = prod v[p^e] for each p^e || n
		target = gen_1;
		for (i = 1; i <= l2; ++i)
			target = mulii(target, gel(v, itos(powii(gcoeff(f, i, 1), gcoeff(f, i, 2)))));
		
		// If v[n] is not equal to the target, the sequence is not
		// completely multiplicative.
		if (!gequal(gel(v, n), target)) {
			if (verbose)
				pari_printf("Not completely multiplicative at n = %Ps = %ld.\n", fnice(stoi(n)), n);
			avma = ltop;
			return 0;
		}
		avma = ltop;
	}
	avma = ltop;
	return 1;
}


long
checkdiv(GEN v, long verbose/*=1*/)
{
	pari_sp ltop = avma;
	long l1;	  /* lg */
	if (!is_matvec_t(typ(v)))
		pari_err(typeer, "checkdiv");

	l1 = lg(v);
	pari_sp btop = avma;
	long n, l2;
	long l3;	  /* lg */
	for (n = 1; n < l1; ++n) {
		l2 = n + n;
		l3 = lg(v);
{
		pari_sp btop = avma;
		GEN i = gen_0;
		long l4 = n > 0;	  /* bool */
		for (i = stoi(l2); l4?gcmpgs(i, l3-1) <= 0:gcmpgs(i, l3-1) >= 0; i = gaddgs(i, n)) {
			if (!gequal0(gmod(gel(v, gtos(i)), gel(v, n)))) {
				if (verbose)
					pariprintf("Not a divisibility sequence: a(%ld) = %Ps does not divide a(%Ps) = %Ps.\n", n, gel(v, n), i, gel(v, gtos(i)));
				avma = ltop;
				return 0;
			}
			i = gerepileupto(btop, i);
		}
}
		avma = btop;
	}
	avma = ltop;
	return 1;
}


// TODO: Better algorithm than continued fractions?  Failing that, at least
// control number of digits?
GEN
solvePell(GEN n, long prec)
{
	pari_sp ltop = avma;
	GEN myprec = gen_0, C = gen_0, k = gen_1, t = gen_0, x = gen_0, y = gen_0;
	if (typ(n) != t_INT)
		pari_err(typeer, "solvePell");
	myprec = stoi(getrealprecision());
	setrealprecision(125, &prec);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while (1) {
		C = contfrac0(gsqrt(n, prec), NULL, 0);
{
		pari_sp btop = avma, st_lim = stack_lim(btop, 1);
		GEN p1 = gen_0, p2 = gen_0;	  /* vec */
		while (gcmpgs(k, glength(C)) <= 0) {
			long i;
			p1 = cgetg(gtos(k)+1, t_VEC);
			for (i = 1; gcmpsg(i, k) <= 0; ++i)
				gel(p1, i) = gcopy(gel(C, i));
			t = contfracback(p1, NULL);
			x = numer(t);
			y = denom(t);
			if (gequal1(gsub(gsqr(x), gmul(n, gsqr(y))))) {
				setrealprecision(gtos(myprec), &prec);
				p2 = cgetg(3, t_VEC);
				gel(p2, 1) = gcopy(x);
				gel(p2, 2) = gcopy(y);
				p2 = gerepileupto(ltop, p2);
				return p2;
			}
			k = gaddgs(k, 1);
			if (low_stack(st_lim, stack_lim(btop, 1)))
				gerepileall(btop, 6, &p1, &t, &x, &y, &p2, &k);
		}
}
		setrealprecision(2*getrealprecision(), &prec);
		if (low_stack(st_lim, stack_lim(btop, 1)))
		gerepileall(btop, 5, &C, &t, &x, &y, &k);
	}
	avma = ltop;
	return gen_0;
}


GEN
tetrMod(GEN a, GEN b, GEN M)
{
	pari_err(talker, "busted");
	// FIXME: Handle the case where gcd(a, M) > 1.
	if (typ(a) != t_INT || typ(b) != t_INT || typ(M) != t_INT)
		pari_err(typeer, "tetrMod");
	switch (signe(b)) {
		case -1:
			pari_err(talker, "negative argument");
		case 0:
			return gen_1;
	}
	pari_sp ltop = avma;
	GEN e = icopy(a), v;
	// __builtin_ffsll or expi
	long vlen = itos_or_0(b)-1;
	if (vlen < 0)
		vlen = expi(subis(M, 1)) + 1;	// Upper bound on A003434
	else
		vlen = minuu(vlen, expi(subis(M, 1)) + 1);
	v = cgetg(vlen+1, t_VEC);
	gel(v, 1) = M;
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 2; i <= vlen; ++i)
	{
		gel(v, i) = geulerphi(gel(v, i - 1));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			v = gerepilecopy(btop, v);
	}
	for (i = glength(v); i >= 1; i--)
	{
		e = lift(gpow(gmodulo(a, gel(v, i)), e, FAKE_PREC));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			e = gerepilecopy(btop, e);
	}
	e = gerepileuptoint(ltop, e);
	return e;
}


GEN
tetrMod_tiny(ulong a, ulong b, ulong M)
{
	/*
	if (typ(a) != t_INT || typ(b) != t_INT || typ(M) != t_INT)
		pari_err(typeer, "tetrMod");
	pari_sp ltop = avma;
	GEN e = icopy(a), v, p1;
	long l3, l4;
	p1 = subis(b, 1);
	v = cgetg(itos(p1)+1, t_VEC);
	gel(v, 1) = M;
	l3 = glength(v);
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 2; i <= l3; ++i)
	{
		gel(v, i) = geulerphi(gel(v, i - 1));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			v = gerepilecopy(btop, v);
	}
	l4 = glength(v);
	btop = avma, st_lim = stack_lim(btop, 1);
	GEN ii = gen_0;
	long l6 = -1 > 0;
	for (ii = stoi(l4); l6?gcmpgs(ii, 1) <= 0:gcmpgs(ii, 1) >= 0; ii = gaddgs(ii, -1))
	{
		e = lift(gpow(gmodulo(a, gel(v, gtos(ii))), e, prec));
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 2, &i, &e);
	}
	e = gerepileuptoint(ltop, e);
	return e;
	*/
	return NEVER_USED;
}


/*
Continued fraction:
        av = avma; lx = lg(x);
        e = bit_accuracy(lx)-1-expo(x);
        if (e < 0) pari_err(talker,"integral part not significant in gboundcf");
        c = trunc2nr_lg(x,lx,0);
        y = int2n(e);
        a = Qsfcont(c,y, NULL, 0);
        b = addsi(signe(x), c);
        return gerepilecopy(av, Qsfcont(b,y, a, 0));*/

// FIXME: Needs numerical analysis to determine stopping point.  Also needs to
// handle rational numbers and intgers.
GEN
Engel(GEN x, long prec)
{
	GEN v, t, ret;
	switch (typ(x)) {
		case t_INT:
			if (signe(x) < 0)
				pari_err(talker, "negative argument");
			long n = itos(x);
			long i = 1;
			v = cgetg(n + 1, t_VEC);
			for (; i <= n; ++i)
				gel(v, i) = gen_1;
			return v;
			
		case t_FRAC:
			pari_err(impl, "fractions");
		
		case t_REAL:
			break;
		
		default:
			pari_err(typeer, "Engel");
	}
	
	pari_sp ltop = avma;
	v = listcreate();
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	while (1)
	{
		if (cmprr(x, real_0(prec)) == 0)
		{
			ret = gtovec(v);
			ret = gerepileupto(ltop, ret);
			return ret;
		} else {
			t = ceilr(invr(x));
		}
		listput(v, t, 0);
		x = subrs(mulri(x, t), 1);
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &x);
	}
	return NEVER_USED;
}


GEN
Eng(GEN n)
{
  pari_sp ltop = avma;
  GEN tmp = gen_0, s = gen_0;
  GEN p1 = gen_0;	  /* vec */
  GEN p2 = gen_0;
  GEN p3 = gen_0;	  /* vec */
  GEN p4 = gen_0;
  if (typ(n) != t_INT)
    pari_err(typeer, "Eng");
  s = strtoGENstr("");
  if (cmpis(n, 1000000) >= 0)
  {
    pariprintf("tmp: %Ps\n", tmp);
    tmp = truedivis(n, 1000000);
    s = Str(mkvec2(Eng(tmp), strtoGENstr(" million")));
    n = subis(n, gtos(gmulgs(tmp, 1000000)));
    if (!signe(n))
    {
      s = gerepileupto(ltop, s);
      return s;
    }
    s = Str(mkvec2(s, strtoGENstr(" ")));
  }
  if (cmpis(n, 1000) >= 0)
  {
    tmp = truedivis(n, 1000);
    s = Str(mkvec3(s, Eng(tmp), strtoGENstr(" thousand")));
    n = subis(n, gtos(gmulgs(tmp, 1000)));
    if (!signe(n))
    {
      s = gerepileupto(ltop, s);
      return s;
    }
    s = Str(mkvec2(s, strtoGENstr(" ")));
  }
  if (cmpis(n, 100) >= 0)
  {
    tmp = truedivis(n, 100);
    s = Str(mkvec3(s, Edigit(tmp), strtoGENstr(" hundred")));
    n = subis(n, gtos(gmulgs(tmp, 100)));
    if (!signe(n))
    {
      s = gerepileupto(ltop, s);
      return s;
    }
    s = Str(mkvec2(s, strtoGENstr(" ")));
  }
  if (cmpis(n, 20) < 0)
  {
    p1 = cgetg(20, t_VEC);
    gel(p1, 1) = strtoGENstr("one");
    gel(p1, 2) = strtoGENstr("two");
    gel(p1, 3) = strtoGENstr("three");
    gel(p1, 4) = strtoGENstr("four");
    gel(p1, 5) = strtoGENstr("five");
    gel(p1, 6) = strtoGENstr("six");
    gel(p1, 7) = strtoGENstr("seven");
    gel(p1, 8) = strtoGENstr("eight");
    gel(p1, 9) = strtoGENstr("nine");
    gel(p1, 10) = strtoGENstr("ten");
    gel(p1, 11) = strtoGENstr("eleven");
    gel(p1, 12) = strtoGENstr("twelve");
    gel(p1, 13) = strtoGENstr("thirteen");
    gel(p1, 14) = strtoGENstr("fourteen");
    gel(p1, 15) = strtoGENstr("fifteen");
    gel(p1, 16) = strtoGENstr("sixteen");
    gel(p1, 17) = strtoGENstr("seventeen");
    gel(p1, 18) = strtoGENstr("eighteen");
    gel(p1, 19) = strtoGENstr("ninteen");
    p2 = Str(mkvec2(s, gel(p1, itos(n))));
    p2 = gerepileupto(ltop, p2);
    return p2;
  }
  tmp = truedivis(n, 10);
  p3 = cgetg(10, t_VEC);
  gel(p3, 1) = gen_0;
  gel(p3, 2) = strtoGENstr("twenty");
  gel(p3, 3) = strtoGENstr("thirty");
  gel(p3, 4) = strtoGENstr("forty");
  gel(p3, 5) = strtoGENstr("fifty");
  gel(p3, 6) = strtoGENstr("sixty");
  gel(p3, 7) = strtoGENstr("seventy");
  gel(p3, 8) = strtoGENstr("eighty");
  gel(p3, 9) = strtoGENstr("ninety");
  s = Str(mkvec2(s, gel(p3, gtos(tmp))));
  n = subis(n, gtos(gmulgs(tmp, 10)));
  if (signe(n))
    p4 = Str(mkvec3(s, strtoGENstr("-"), Edigit(n)));
  else
    p4 = s;
  p4 = gerepileupto(ltop, p4);
  return p4;
}

GEN
Edigit(GEN n)
{
  pari_sp ltop = avma;
  GEN p1 = gen_0;	  /* vec */
  GEN p2 = gen_0;
  p1 = cgetg(10, t_VEC);
  gel(p1, 1) = strtoGENstr("one");
  gel(p1, 2) = strtoGENstr("two");
  gel(p1, 3) = strtoGENstr("three");
  gel(p1, 4) = strtoGENstr("four");
  gel(p1, 5) = strtoGENstr("five");
  gel(p1, 6) = strtoGENstr("six");
  gel(p1, 7) = strtoGENstr("seven");
  gel(p1, 8) = strtoGENstr("eight");
  gel(p1, 9) = strtoGENstr("nine");
  p2 = gcopy(gel(p1, gtos(n)));
  p2 = gerepileupto(ltop, p2);
  return p2;
}


// TODO: Would be much more efficient to just walk through the primelist...
ulong
ucomposite(long n)
{
	if (n < 2) {
		if (n < 1)
			pari_err(talker, "n-th composite meaningless if n = %ld", n);
		return 4;
	}
	double l = 1.0/log(n);
	
	// Series apparently due to Bojarincev 1967; more terms are known if desired
	ulong c = (ulong)(n * (1 + l * (1 + 2*l * (1 + 2*l))));
	int i = 0;
	int limit = (int)(n - c + uprimepi(c));
	for (; i <= limit; i++)
		if (uisprime(++c))
			i--;
	return c;
}
GEN
composite(long n) { return utoipos(ucomposite(n)); }


GEN
deBruijnXi(GEN x)
{
	double xx = rtodbl(x), left, right;
	if (xx < 1)
		pari_err(talker, "deBruijnXi: Can't find a xi given x < 1.");
	if (xx > 1)
		left = log(xx);
	else
		left = DBL_EPSILON;
	right = 1.35 * log(xx) + 1;	// Heuristic

	// Bisection
	while (right - left > left * DBL_EPSILON) {
		double m = (left + right) / 2;
		if (expm1(m) > xx * m)
			right = m;
		else
			left = m;
	}
	return dbltor((left + right) / 2);
}

GEN
rhoest(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN xi, ret;
	x = gtor(x, "rhoest", prec);
	xi = deBruijnXi(x);
	ret = gexp(gsub(gneg(veceint1(negr(xi), NULL, prec)), gmul(x, xi)), prec);
	ret = gdiv(gdiv(ret, sqrtr(mulrr(mulsr(2, mppi(prec)), x))), xi);
	ret = gerepileupto(ltop, ret);
	return ret;
}


static double rhoTable[] = {
	NEVER_USED, 1, 3.068528194e-1, 4.860838829e-2, 4.910925648e-3,
	3.547247005e-4, 1.964969635e-5, 8.745669953e-7, 3.232069304e-8,
	1.016248283e-9,	2.770171838e-11, 6.644809070e-13, 1.419713165e-14,
	2.729189030e-16, 4.760639989e-18, 7.589908004e-20
};
static int rhoTableLen = 15;
static double rhoScale = 1.130709295873035782;	// Last table entry, divided by
// rhoest at that point


GEN
DickmanRho(GEN x, long prec)
{
	pari_sp ltop = avma;
	GEN ret, left, right, scale;
	x = gtor(x, "DickmanRho", prec);
	if (cmprs(x, 2) <= 0) {
		ret = gsubsg(1, glog(gmaxgs(x, 1), prec));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
	if (gcmpgs(x, 3) <= 0) {
		ret = gadd(gadd(gsubsg(1, mulrr(subsr(1, mplog(subrs(x, 1))), mplog(x))), greal(dilog(subsr(1, x), prec))), divrs(sqrr(mppi(prec)), 12));
		ret = gerepileupto(ltop, ret);
		return ret;
	}
  
	double xx = rtodbl(x);
  
	// Asymptotic estimate (scaled for continuity)
	if (xx > rhoTableLen) {
		double scale = rhoScale;
		scale = (scale - 1) * sqrt(sqrt(rhoTableLen / xx)) + 1;
		/* Let the scale factor dwindle away, since the estimate is (presumably) */
		/* better in the long run than any scaled version of it.  The exponent */
		/* of 0.25 has been chosen to give the best results for 10 < x < 100 */
		/* with a table size of 10. */

		ret = precision0(mulrr(rhoest(x, prec), dbltor(scale)), 9);
		ret = gerepileupto(ltop, ret);
		return ret;
	}
  
	// Scaling factors: the factor by which the true value of rho differs from
	// the estimates at the endpoints.
	left = divrr(dbltor(rhoTable[(int)floor(xx)]), rhoest(floorr(x), prec));
	right = divrr(dbltor(rhoTable[(int)ceil(xx)]), rhoest(ceilr(x), prec));
	
	// Linear interpolation on the scale factors.
	scale = gadd(left, gmul(gsub(right, left), mpsub(x, floorr(x))));
	
	// Return a result based on the scale factor and the asymptotic formula.
	ret = precision0(gmul(rhoest(x, prec), scale), 9);
	ret = gerepileupto(ltop, ret);
	return ret;
}


GEN
graeffe(GEN f)
{
  pari_sp ltop = avma;
  GEN d = gen_0, g = gen_0, h = gen_0, p1 = gen_0;
  GEN p2 = gen_0;	  /* vec */
  GEN p3 = gen_0;
  GEN p4 = gen_0;	  /* vec */
  GEN x = pol_x(fetch_user_var("x")), p5 = gen_0;
  if (typ(f) != t_POL)
    pari_err(typeer, "graeffe");
  d = stoi(degpol(f));
  p1 = gaddgs(gdiventgs(d, 2), 1);
  {
    long i;
    p2 = cgetg(gtos(p1)+1, t_VEC);
    for (i = 1; gcmpsg(i, p1) <= 0; ++i)
      gel(p2, i) = polcoeff0(f, (2*i) - 2, -1);
  }
  g = p2;
  p3 = gdiventgs(gaddgs(d, 1), 2);
  {
    long i;
    p4 = cgetg(gtos(p3)+1, t_VEC);
    for (i = 1; gcmpsg(i, p3) <= 0; ++i)
      gel(p4, i) = polcoeff0(f, (2*i) - 1, -1);
  }
  h = p4;
  p5 = gsub(gsqr(gtopolyrev(g, -1)), gmul(x, gsqr(gtopolyrev(h, -1))));
  p5 = gerepileupto(ltop, p5);
  return p5;
}


long
iscyclo(GEN f)
{
	pari_sp ltop = avma;
	GEN f1, f2, fn, x = pol_x(fetch_user_var("x"));
	long l1;
	if (typ(f) != t_POL)
		pari_err(typeer, "iscyclo");
	f1 = graeffe(f);
	if (gequal(f, f1)) {
		avma = ltop;
		return 1;
	}
	fn = gsubst(f, gvar(x), gneg(x));
	if (gequal(f1, fn) && iscyclo(fn)) {
		avma = ltop;
		return 1;
	}
	l1 = !gequal0(gissquareall(f1, &f2)) && iscyclo(f2);
	avma = ltop;
	return l1;
}


long
istotient(GEN n)
{
	if (typ(n) != t_INT)
		pari_err(arither1, "istotient");
	if (signe(n) < 1)
		return 0;
	if(mod2(n))
		return isint1(n);

	pari_sp ltop = avma;
	GEN k, p, d, p2;
	k = n;
	while (1) {
		if (totientHelper(k, gen_2)) {
			avma = ltop;
			return 1;
		}
		if (mod2(k))
			break;
		k = shifti(k, -1);
		k = gerepileuptoint(ltop, k);
	}
	p2 = divisors(shifti(n, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long i;
	for (i = 1; i < lg(p2); ++i) {
		d = shifti(gel(p2, i), 1);
		if (!(isprime(p = addis(d, 1))))
			continue;
		k = diviiexact(n, d);
		while (1) {
			if (totientHelper(k, p)) {
				avma = ltop;
				return 1;
			}
			if (!dvdii(k, p))
				break;
			k = diviiexact(k, p);
		}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &k);
	}
	avma = ltop;
	return 0;
}


/* This function should be called only internally.  n and m are positive t_INT
 * values.  Test with
#
sum(n=1,35214/2,istotient(2*n)==0)
forstep(n=2,1e6,2,istotient(2*n))
 */
long
totientHelper(GEN n, GEN m)
{
	if (mod2(n))
		return equali1(n);
	pari_sp ltop = avma;
	GEN k, p, d, p1;
	p1 = divisors(shifti(n, -1));
	pari_sp btop = avma, st_lim = stack_lim(btop, 1);
	long l2;
	for (l2 = 1; l2 < lg(p1); ++l2) {
		d = shifti(gel(p1, l2), 1);
		if ((cmpii(d, m) < 0) || !(isprime(p = addis(d, 1))))
			continue;
		k = diviiexact(n, d);
		while (1) {
			if (totientHelper(k, p)) {
				avma = ltop;
				return 1;
			}
			if (!dvdii(k, p))
				break;
			k = diviiexact(k, p);
			/* // Slower:
			GEN quotient, remainder;
			quotient = dvmdii(k, p, &remainder);
			if (remainder != gen_0)	// documentation guarantees that gen_0, not a copy, is returned if p|k.
				break;
			k = quotient;
			*/
		}
		if (low_stack(st_lim, stack_lim(btop, 1)))
			gerepileall(btop, 1, &k);
	}
	avma = ltop;
	return 0;
}



////////////////////////////////////////////////////////////////////////////////////// End newcomers

////////////// Problems with gp2c //////////////
// Uses				Should use
// pariprintf		pari_printf
////////////////////////////////////////////////

// Warning types:
// user		User-initiated warning
// warnmem	Garbage collection
// warner	Generic warning
// warnprec	Precision increase
// warnfile	File I/O

// Error types:
// 0			Generic error
// talker2		?
// bugparier	Bug, please report
// alarmer		Generic error
// openfiler	File I/O
// talker		Generic error
// flagerr		Invalid flag
// impl			Not implemented
// archer		Not available on this system
// notfuncer	Not a function in function call
// precer		Precision too low
// typeer		Incorrect type
// consister	Inconsistent data
// user			User-initiated error
// errpile		Stack overflow
// overflower	Overflow
// matinv1		Non-invertible matrix (in gauss)
// mattype1		Not a square matrix
// arither1		Not an integer argument in an arithmetic function
// primer1		Not enough precomputed primes
// invmoder		Impossible inverse
// constpoler	Constant polynomial
// notpoler		Not a polynomial
// redpoler		Reducible polynomial
// zeropoler	Zero polynomial
// operi		"Impossible"
// operf		"Forbidden"
// gdiver		Division by zero
// memer		Not enough memory
// negexper		Negative exponent
// sqrter5		Non quadratic residue (in gsqrt)
// noer			Not an error...

/*

Possible additions to the OEIS:

A138591
a(n) = n + A000523(n + A000523(n))

A063274
For all n > N, a(n) <= 3. Possibly N = 119.
Heath-Brown, D. R. "Ternary Quadratic Forms and Sums of Three Square-Full Numbers." In Séminaire de Theorie des Nombres, Paris 1986-87	(Ed. C. Goldstein). Boston, MA: Birkhäuser, pp. 137-163, 1988.
 
A045911
Crossref A******

A******
Neither a cube nor the sum of a positive cube and a prime.
Hardy & Littlewood's conjecture that this sequence is finite.
H-L, ... (Conjecture L, p. 51)
Crossref A045911

*/

