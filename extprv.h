#include "ext.h"

GEN poleval_denseint(GEN x, GEN y);
void listput_shallow(GEN L, GEN x);
GEN gtor(GEN x, const char* funcName, long prec);
long countdigits(GEN x);
GEN listtovec_shallow(GEN v);
ulong cuberoot(ulong n) __attribute__ ((const));
void init_auto(void);
void assume (int expr);
GEN vecsum_zv(GEN x);
long uissemiprime(ulong n);
ulong lpfu(ulong n);
GEN primezeta_complex(GEN s);
GEN primezeta_real(GEN s);
double lnBell(long n);
double W_small(double x);
long isExtendedReal(GEN x) __attribute__ ((pure));
long valu(ulong n) __attribute__ ((const));
long words_free(void) __attribute__ ((pure));
