#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Weverything"
#endif
#include <pari/pari.h>
#include <pari/paripriv.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
#include <float.h>
#include <limits.h>
//#include <inttypes.h>

void
init_auto(void)
{
    rnormal_cached = 0;
}
