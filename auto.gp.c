/*-*- compile-command: "/usr/bin/gcc -c -o auto.gp.o -O3 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -I"/usr/local/include" auto.gp.c && /usr/bin/gcc -o auto.gp.so -shared -O3 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari"; -*-*/
#include <pari/pari.h>
#include <pari/paripriv.h>
#include <float.h>
#include <limits.h>
/*
GP;install("isfactorial","lG","isfactorial","./auto.gp.so");
GP;addhelp(isfactorial, "isfactorial(n): Is n a factorial? Sloane's A012245; characteristic function of Sloane's A000142.");
GP;install("facmod","LL","facmod","./auto.gp.so");
GP;install("glnBell","L","lnBell", "./auto.gp.so");
GP;install("issm3","lL","issm3","./auto.gp.so");
GP;install("Bell","L","Bell","./auto.gp.so");
GP;addhelp(Bell, "Bell(n): Returns the n-th Bell or exponential number, Sloane's A000110.");
GP;install("log_2","Gp","lg","./auto.gp.so");
GP;addhelp(lg, "lg(x): Binary logarithm of x.");
GP;install("rad","G","rad","./auto.gp.so");
GP;addhelp(rad, "rad(n): Radical of n, the largest squarefree number dividing n. Sloane's A007947.");
GP;install("prp","lGDG","prp","./auto.gp.so");
GP;addhelp(prp, "prp(n,b=2): Is n a b-probable prime?");
GP;install("sprp","lGDG","sprp","./auto.gp.so");
GP;addhelp(sprp, "sprp(n,b=2): Is n a b-strong probable prime?");
GP;install("sopf","G","sopf","./auto.gp.so");
GP;addhelp(sopf, "sopf(n): Sum of distinct prime factors of n. Sloane's A008472.");
GP;install("sopfr","G","sopfr","./auto.gp.so");
GP;addhelp(sopfr, "sopfr(n): Sum of prime factors of n (with multiplicity). Sloane's A001414.");
GP;install("primorial","G","primorial","./auto.gp.so");
GP;addhelp(primorial, "primorial(n): Returns the product of each prime less than or equal to n. Sloane's A034386.");
GP;install("gpf","G","gpf","./auto.gp.so");
GP;addhelp(gpf, "gpf(n): The greatest prime factor of a number. Sloane's A006530.");
GP;install("lpf","G","lpf","./auto.gp.so");
GP;addhelp(lpf, "lpf(n): The least prime factor of a number. Sloane's A020639.");
GP;install("isTriangular","lG","isTriangular","./auto.gp.so");
GP;addhelp(isTriangular, "isTriangular(n): Is n a triangular number? Sloane's A010054; characteristic function of Sloane's A000217.");
GP;install("isHexagonal","lG","isHexagonal","./auto.gp.so");
GP;addhelp(isHexagonal, "isHexagonal(n): Is n a hexagonal number? Characteristic function of Sloane's A000384.");
GP;install("isFibonacci","lG","isFibonacci","./auto.gp.so");
GP;addhelp(isFibonacci, "isFibonacci(n): Is n a Fibonacci number? Sloane's A010056; characteristic function of Sloane's A000045.");
GP;install("largestSquareFactor","G","largestSquareFactor","./auto.gp.so");
GP;addhelp(largestSquareFactor, "largestSquareFactor(n): Largest square dividing n. Sloane's A008833.");
GP;install("ispow2","lG","ispow2","./auto.gp.so");
GP;addhelp(ispow2, "ispow2(n): Is n a power of two? Characteristic function of Sloane's A000079.");
GP;install("ispow3","lG","ispow3","./auto.gp.so");
GP;addhelp(ispow3, "ispow3(n): Is n a power of three? Characteristic function of Sloane's A000244.");
GP;install("issemiprime","lG","issemi","./auto.gp.so");
GP;addhelp(issemi, "issemi(n): Is n a semiprime? Sloane's A064911; characteristic function of Sloane's A001358.");
GP;install("istwo","lG","istwo","./auto.gp.so");
GP;addhelp(istwo, "istwo(n): Is the number a sum of two squares? Characteristic function of Sloane's A001481.");
GP;install("ways2","G","ways2","./auto.gp.so");
GP;addhelp(ways2, "ways2(n): Number of ways that n can be represented as a sum of two squares. Sloane's A000161.");
GP;install("isthree","lG","isthree","./auto.gp.so");
GP;addhelp(isthree, "isthree(n): Is the number the sum of three squares? Sloane's A071374; characteristic function of Sloane's A000378.");
GP;install("ways3","G","ways3","./auto.gp.so");
GP;addhelp(ways3, "ways3(n): Number of ways that n can be represented as a sum of three squares. Sloane's A000164.");
GP;install("msb","G","msb","./auto.gp.so");
GP;addhelp(msb, "msb(n): Most significant bit of n: returns the greatest power of 2 <= the number. Sloane's A053644.");
GP;install("Faulhaber","LDG","Faulhaber","./auto.gp.so");
GP;addhelp(Faulhaber, "Faulhaber(e,{a='x}): Returns the polynomial for the sum 1^e + 2^e + ... + x^e, evaluated at a.");
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
GP;install("oddres","G","oddres","./auto.gp.so");
GP;addhelp(oddres, "oddres(n): Returns the greatest odd number dividing n.");
GP;install("toC","vG","toC","./auto.gp.so");
GP;addhelp(toC, "toC(n): Format n for use with the PARI library (e.g., with gp2c programs).");
GP;install("digits","lG","digits","./auto.gp.so");
GP;addhelp(digits, "digits(n): Number of decimal digits in n. Sloane's A055642.");
GP;install("eps","p","eps","./auto.gp.so");
GP;addhelp(eps, "eps(): Returns machine epsilon for the current precision.");
GP;install("fnice","G","fnice","./auto.gp.so");
GP;addhelp(fnice, "fnice(n): Returns a string with a 'nice' factorization of n.");
GP;install("tonice","D0,G,p","nice","./auto.gp.so");
GP;addhelp(nice, "nice(o): Reformats the object o 'nicely' when possible. Currently chokes on multivariable polynomials.");
GP;install("sumset","D0,G,D0,G,","sumset","./auto.gp.so");
GP;addhelp(sumset, "sumset(A, B): Set of all numbers of the form a+b, a in A, b in B.");
GP;install("diffset","D0,G,D0,G,","diffset","./auto.gp.so");
GP;addhelp(diffset, "diffset(A, B): Set of all numbers of the form a-b, a in A, b in B.");
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
GP;install("forthinprime","vV=GGI","forthinprime","./auto.gp.so");
GP;addhelp(forthinprime, "forthinprime(X=a,b,seq): the sequence is evaluated, X running over the primes between a and b, even if b > primelimit. EXPERIMENTAL!");
GP;install("forbigprime","vV=GGI","forbigprime","./auto.gp.so");
GP;addhelp(forbigprime, "forbigprime(X=a,b,seq): the sequence is evaluated, X running over the primes between a and b.");
GP;install("fibmod","GG","fibmod","./auto.gp.so");
GP;addhelp(fibmod, "fibmod(n,m): Returns the nth Fibonacci number mod m. Same as finonacci(n)%m, but faster for large n.");
GP;install("bfile","GDGDG","bfile","./auto.gp.so");
GP;addhelp(bfile, "bfile(name, v, offset=1): If v is given, creates a b-file with the values of v, using name as the A-number (given as a number or a filename). If v is not given, open the b-file name (again, as a filename or number) and return a vector of its values.");
GP;install("checkmult","lGD1,L,","checkmult","./auto.gp.so");
GP;addhelp(checkmult, "checkmult(v,{verbose=1}): Is the sequence v multiplicative?");
GP;install("checkcmult","lGD1,L,","checkcmult","./auto.gp.so");
GP;addhelp(checkcmult, "checkcmult(v,{verbose=1}): Is the sequence v completely multiplicative?");
GP;install("checkadd","lGD1,L,","checkadd","./auto.gp.so");
GP;addhelp(checkadd, "checkadd(v,{verbose=1}): Is the sequence v additive?");
GP;install("checkcadd","lGD1,L,","checkcadd","./auto.gp.so");
GP;addhelp(checkcadd, "checkcadd(v,{verbose=1}): Is the sequence v completely additive?");
GP;install("checkdiv","lGD1,L,","checkdiv","./auto.gp.so");
GP;addhelp(checkdiv, "checkdiv(v,{verbose=1}): Is v a divisibility sequence?");
GP;install("solvePell","G","solvePell","./auto.gp.so");
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
GP;addhelp(rhoest, "rhoest(x): de Bruijn's asymptotic approximation for rho(x), rewritten as in van de Lune and Wattel 1969.  Curiously, their paper shows values for this estimate that differ from those calculated by this function, often as soon as the second decimal place -- but as the difference is in the direction of the true value, I have not looked further into this.");
GP;install("DickmanRho","Gp","DickmanRho","./auto.gp.so");
GP;addhelp(DickmanRho, "DickmanRho(x): Estimates the value of the Dickman rho function. For x <= 3 the exact values are used, up to rounding; up to 15 the value is interpolated using known values and rhoest; after 15 rhoest is used, along with a correction factor based on the last value in rhoTable.");
GP;alias(isPowerful, ispowerful);
GP;alias(hamming, hammingweight);
GP;alias(graeffe,polgraeffe);
GP;alias(dsum, sumdigits);
// New \/
GP;install("tetrMod","GGG","tetrMod","./auto.gp.so");
GP;addhelp(tetrMod, "tetrMod(a,b,M): Returns a^^b mod M.");
* // New /\
// No associated help
GP;install("init_auto","v","init_auto","./auto.gp.so");
GP;install("consistency","l","consistency","./auto.gp.so");
GP;install("ucountPowerfuli","lD0,G,","cP","./auto.gp.so");
GP;install("ucountSquarefree","lL","cS","./auto.gp.so");
*/

// Removed: //GP;install("pBounds","vD0,G,DGp","pBounds","./auto.gp.so");
// Removed: //GP;addhelp(pBounds, "pBounds(n, verbose=0): Estimates the nth prime. Set verbose=1 to get a list of sources for the results.");


//////////////////////////////////////////////////////////// New
GEN Bell(long n);
//////////////////////////////////////////////////////////// New


GEN listtovec_shallow(GEN v);
long checkmult(GEN v, long verbose);
long checkcmult(GEN v, long verbose);
long checkdiv(GEN v, long verbose);
GEN solvePell(GEN n);
GEN tetrMod(GEN a, GEN b, GEN M);
GEN Engel(GEN x, long prec);
GEN Eng(GEN n);
GEN Eng_small(long n);
GEN Eng_tiny(long n);
GEN Edigit(long n);
GEN composite(long n);
GEN deBruijnXi(GEN x);
GEN rhoest(GEN x, long prec);
GEN DickmanRho(GEN x, long prec);
long issemiprime(GEN n);
long uissemiprime(ulong n);
GEN phiset(GEN v);
GEN rad(GEN n);
INLINE long valu(ulong n);
long isPowerful_small(ulong n);
long prp(GEN n, GEN b);
long sprp(GEN n, GEN b);
GEN sopf(GEN n);
GEN sopfr(GEN n);
GEN primorial(GEN n);
GEN prodtree(GEN A, long start, long stop);
GEN gpf(GEN n);
GEN lpf(GEN n);
long isFibonacci(GEN n);
long isSmallFib(long n);
GEN fibmod(GEN n, GEN m);
long Pisano(long p, long e);
GEN largestSquareFactor(GEN n);
INLINE long hamming_word(ulong w);
long ispow2(GEN n);
long ispow3(GEN n);
long ispow3_tiny(ulong n);
long istwo(GEN n);
GEN ways2(GEN n);
long isthree(GEN n);
long sways3s(ulong n);
GEN ways3(GEN n);
GEN msb(GEN n);
GEN Faulhaber(long e, GEN a);
GEN rp(long b);
GEN countPowerful(GEN lim);
GEN countSquarefree(GEN lim);
ulong ucountSquarefree(ulong lim);
GEN Mfactor(GEN p, GEN lim, GEN start);
GEN bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start);
long bigdiv(GEN a, GEN b, GEN c, GEN d);
GEN contfracback(GEN v, GEN terms);
double W_small(double x);
GEN W(GEN x, long prec);
GEN vecsum(GEN v);
GEN vecprod(GEN v);
GEN vecgcd(GEN v);
GEN veclcm(GEN v);
GEN oddres(GEN n);
void toC(GEN n);
long digits(GEN x);
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
void bfileout(char* filename, GEN name, GEN v, GEN Anum, long offset);
void forbigprime(GEN ga, GEN gb, GEN code);
void forbigprime_sieve(ulong a, ulong b, GEN code);
void forthinprime(GEN ga, GEN gb, GEN code);
void forthinprimeuu(ulong a, ulong b, GEN code);
INLINE GEN gtor(GEN x, const char* funcName, long prec);
void init_auto(void);
/*End of prototype*/

#define FAKE_PREC 0		// Used when a precision is required but will not be used
#define NEVER_USED 0	// Used to initialize values so the compiler doesn't complain
GEN rnormal_cached;

#include "othergpincludes.h"
#include "prime.gp.c"
#include "arith.gp.c"
#include "numth.gp.c"
#include "rc.gp.c"
#include "conv.gp.c"
#include "io.gp.c"
#include "loops.gp.c"
#include "other.gp.c"

void
init_auto(void)
{
	rnormal_cached = 0;
}



/*
default(debug,4)
consistency()
*/




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

