#!/bin/bash
cd ~/mth/pari

# 'Basic' warnings.
W='-Wall -Wextra -Wcast-align -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wunused-macros -Wswitch-enum -Wold-style-definition -Wpointer-arith -Wnested-externs -Wlogical-op -Winline -Winit-self -Wformat-security -Wformat-nonliteral -Wformat-y2k'

# Extra warnings; these may give lots of output but are probably worthwhile. Disable if needed.
W="$W -Wfloat-equal -Wdisabled-optimization -Wwrite-strings -Wtrampolines -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wshadow"
W="$W -Wvector-operation-performance -Wabi -Wpadded -Wpacked -Waggregate-return -Woverlength-strings -Wc++-compat "

# Not used: some warnings I have found to be unhelpful.
ignore="-Wbad-function-cast -Wconversion -Wswitch-default -Wunsuffixed-float-constants -Wjump-misses-init -Wunsafe-loop-optimizations -Wdouble-promotion -Wsystem-headers -Wnoexcept -Wbad-function-cast -Wdeclaration-after-statement -Wpedantic"

# Warnings which do not apply to C (but to C++, Objective C++, D, etc.).
ignore="$ignore -Wuseless-cast -Weffc++ -Wc++0x-compat -Wzero-as-null-pointer-constant -Wreorder -Wcast-result -Wabi-tag -Wsign-promo -Wnon-virtual-dtor -Wctor-dtor-privacy -Wctor-dtor-privacy -Wdelete-non-virtual-dtor -Wassign-intercept -Wdelete-non-virtual-dtor -Wc++0x-compat -Wundeclared-selector -Wold-style-cast -Wselector -Woverloaded-virtual"

# Basic optimizations.
O='-O2 -fno-strict-aliasing -fomit-frame-pointer -fPIC'

# Powerful Graphite optimizations.
O="$O -ftree-loop-linear -floop-interchange -floop-strip-mine -floop-block -fgraphite-identity -floop-nest-optimize -floop-parallelize-all"

# All options together here
GCCOPT="-march=native -m64 $W $O"

# What warnings are missed?
gcc $GCCOPT -Q --help=warnings | fgrep '[disabled]' | egrep -v `echo "'$ignore'" | tr ' ' '|' | sed s/-/[-]/ | sed s/+/[+]/g`

# Compile auto file
gcc -o auto.gp.o $GCCOPT -c -I"/usr/local/include" auto.gp.c
gcc -o auto.gp.so $GCCOPT -shared -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari

# Create run file to install auto and add associated help entries
egrep '^GP;' auto.gp.c | sed 's/^GP;//' > auto.gp.run

## Create documentation -- takes a long time
#doxygen Doxyfile
#{ for f in html/*.png ; do pngout $f ; done } &

# Run gp itself
gp
