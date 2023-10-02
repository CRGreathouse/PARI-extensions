#ifndef EXT_H
#define EXT_H
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

GEN Bell(long n);
long checkmult(GEN v, long verbose);
long checkcmult(GEN v, long verbose);
long checkdiv(GEN v, long verbose);
GEN solvePell(GEN n);
GEN Engel(GEN x);
GEN Eng(GEN n);
GEN composite(long n);
GEN rhoest(GEN x, long prec);
GEN DickmanRho(GEN x, long prec);
long issemiprime(GEN n);
GEN rad(GEN n);
long prp(GEN n, GEN b);
long sprp(GEN n, GEN b);
GEN sopf(GEN n);
GEN sopfr(GEN n);
GEN gpf(GEN n);
GEN lpf(GEN n);
long isFibonacci(GEN n);
GEN fibmod(GEN n, GEN m);
long ispow2(GEN n);
long ispow3(GEN n);
long istwo(GEN n);
GEN ways2(GEN n);
long isthree(GEN n);
GEN ways3(GEN n);
GEN msb(GEN n);
GEN Faulhaber(long e, GEN a);
GEN countPowerful(GEN lim);
GEN countSquarefree(GEN lim);
GEN Mfactor(GEN p, GEN lim, GEN start);
GEN bigfactor(GEN a, GEN b, GEN c, GEN lim, GEN start);
long bigdiv(GEN a, GEN b, GEN c, GEN d);
GEN contfracback(GEN v, GEN terms);
GEN oddres(GEN n);
void gToC(GEN n);
GEN eps(long prec);
GEN fnice(GEN n);
GEN tonice(GEN o, long prec);
GEN sumset(GEN a, GEN b);
GEN sumset_lim(GEN a, GEN b, GEN lim);
GEN diffset(GEN a, GEN b);
GEN normd(GEN a, GEN b, long prec);
GEN rnormal(long prec);
void pBounds(GEN n, GEN verbose, long prec);
long checkVDW(GEN vv, GEN verbose);
long longestProgression(GEN v);
long longestProgression1(GEN v);
GEN fusc(GEN n);
void forodd(GEN a, GEN b, GEN code);
GEN bfile(GEN name, GEN v, GEN offset);
GEN primezeta(GEN s, long prec);
GEN HurwitzClassNumber(GEN n);
GEN tau(GEN n);
long checkadd(GEN v, long verbose);
long checkcadd(GEN v, long verbose);
GEN cuberootint(GEN x);
void forpal(GEN a, GEN b, GEN code);
GEN glnBell(long n);
long isfactorial(GEN n);
GEN log_2(GEN x, long prec);
GEN rev(GEN n, long B);
long Collatz(GEN n);
long issm3(long n) __attribute__ ((const));
long isSmallFib(long n) __attribute__ ((const));

#define NEVER_USED 0	// Used to initialize values so the compiler doesn't complain
GEN rnormal_cached;

#endif
