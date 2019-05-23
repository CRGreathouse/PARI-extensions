#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Weverything"
#endif
#include <pari/pari.h>
#include <pari/paripriv.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
#include "othergpincludes.h"
#include <float.h>
#include <limits.h>

GEN Bell(long n);
long checkmult(GEN v, long verbose);
long checkcmult(GEN v, long verbose);
long checkdiv(GEN v, long verbose);
GEN solvePell(GEN n);
GEN Engel(GEN x);
GEN Eng(GEN n);
GEN Eng_small(long n);
GEN Eng_tiny(long n);
GEN Edigit(long n);
long countdigits(GEN x);
GEN composite(long n);
GEN deBruijnXi(GEN x);
GEN rhoest(GEN x, long prec);
GEN DickmanRho(GEN x, long prec);
long issemiprime(GEN n);
long uissemiprime(ulong n);
GEN rad(GEN n);
long prp(GEN n, GEN b);
long sprp(GEN n, GEN b);
GEN sopf(GEN n);
GEN sopfr(GEN n);
GEN primorial(GEN n);
GEN prodtree(GEN A, long start, long stop);
GEN gpf(GEN n);
GEN lpf(GEN n);
long isFibonacci(GEN n);
GEN fibmod(GEN n, GEN m);
long Pisano(long p, long e);
long ispow2(GEN n);
long ispow3(GEN n);
long istwo(GEN n);
GEN ways2(GEN n);
long isthree(GEN n);
long sways3s(ulong n);
GEN ways3(GEN n);
GEN msb(GEN n);
GEN Faulhaber(long e, GEN a);
GEN countPowerful(GEN lim);
GEN countSquarefree(GEN lim);
ulong ucountSquarefree(ulong lim);
GEN Mfactor(GEN p, GEN lim, GEN start);
GEN bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start);
long bigdiv(GEN a, GEN b, GEN c, GEN d);
GEN contfracback(GEN v, GEN terms);
double W_small(double x);
GEN oddres(GEN n);
void gToC(GEN n);
const char* toC(GEN n);
GEN eps(long prec);
GEN fnice(GEN n);
GEN tonice(GEN o, long prec);
GEN initial(GEN n, char *s);
GEN medial(GEN n, char *s);
GEN monomialnice(GEN coeff, GEN degree, GEN v);
GEN sumset(GEN a, GEN b);
GEN sumset_lim(GEN a, GEN b, GEN lim);
GEN diffset(GEN a, GEN b);
GEN normd(GEN a, GEN b, long prec);
GEN rnormal(long prec);
void pBounds(GEN n, GEN verbose, long prec);
GEN checkVDW(GEN vv, GEN verbose);
GEN longestProgression(GEN v);
GEN longestProgression1(GEN v);
GEN fusc(GEN n);
GEN fusc_large(GEN n);
ulong ucountPowerfulu(ulong lim);
long issquarefree_small(ulong n);
void forodd(GEN a, GEN b, GEN code);
char* getBValue(char*);
GEN bfile(GEN name, GEN v, GEN offset);
GEN bfilein(char* name);
void bfileout(char* filename, GEN name, GEN v, GEN Anum, long offset);
GEN primezeta(GEN s, long prec);
GEN primezeta_complex(GEN s);
GEN primezeta_real(GEN s);
GEN gtor(GEN x, const char* funcName, long prec);
void init_auto(void);
GEN taup_small(ulong p);
GEN HurwitzClassNumber_small(ulong n);
GEN HurwitzClassNumber(GEN n);
GEN taup_big(GEN p);
GEN taup(GEN p, long e);
GEN tau(GEN n);
GEN poleval_denseint(GEN x, GEN y);
long checkadd(GEN v, long verbose);
long checkcadd(GEN v, long verbose);
GEN cuberootint(GEN x);
GEN fibomod(long n, GEN m);
ulong fibomod_tiny(long n, ulong m);
void forpal(GEN a, GEN b, GEN code);
GEN glnBell(long n);
long infinite(GEN x);
long isExtendedReal(GEN x);
long isfactorial(GEN n);
double lnBell(long n);
GEN log_2(GEN x, long prec);
ulong lpfu(ulong n);
int palhelper(long digits, GEN a, GEN b, GEN code);
GEN prodtree_small(GEN A, long start, long stop);
GEN rev(GEN n, long B);
ulong ucomposite(long n);
ulong ucountPowerfuli(GEN n);
GEN listtovec_shallow(GEN v);
void listput_shallow(GEN L, GEN x);
long Collatz(GEN n);
long Collatz_tiny(ulong n);
long valu(ulong n) __attribute__ ((const));
long issm3(long n) __attribute__ ((const));
long ispow3_tiny(ulong n) __attribute__ ((const));
long isSmallFib(long n) __attribute__ ((const));
ulong cuberoot(ulong n) __attribute__ ((const));
ulong fusc_small(GEN n) __attribute__ ((pure));
void assume (int expr);

#define NEVER_USED 0	// Used to initialize values so the compiler doesn't complain
GEN rnormal_cached;
