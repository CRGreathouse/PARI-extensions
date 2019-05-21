#!/bin/bash
cd ~/mth/PARI-extensions
export GPDOCDIR="/home/charles/mth/pari/doc"

GCC=1
clang=0

#CC='gcc'
#CC='/usr/local/bin/gcc5.1'
#CC='/usr/bin/gcc-6'
#CC='gcc-7'
CC='gcc-9'
#CC='/usr/bin/clang'
#CC='/usr/lib/llvm-3.8/libexec/ccc-analyzer'

# 'Basic' warnings.
W='-Wall -Wextra -Wshadow -Wcast-align=strict -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wunused-macros -Wswitch-enum -Wold-style-definition -Wpointer-arith -Wnested-externs -Winline -Winit-self -Wformat-security -Wformat-nonliteral -Wformat-y2k'
#W="$W -Werror"

# GCC-only warnings
if [ $GCC -eq 1 ]
then
	W="$W -Wlogical-op -Wtrampolines -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format -Wvector-operation-performance"

	# New warnings for gcc 6, disable if using something else maybe
	W="$W -Wformat-signedness"

	# New warnings for gcc 7, disable if using something else maybe
	W="$W -Wnonnull -Wformat-overflow=2 -Wformat-truncation=2 -Wstringop-overflow=3 -Wduplicated-branches -Wrestrict"
	W="$W -Walloc-zero"
fi

# Extra warnings; these may give lots of output but are probably worthwhile. Disable if needed.
W="$W -Wfloat-equal -Wdisabled-optimization -Wwrite-strings"
W="$W -Wpadded -Wpacked -Waggregate-return -Woverlength-strings -Wc++-compat -Wstack-protector -Wundef"
W="$W -Wabi=11"	# warn about changes from GCC 7

# Warnings which should be errors
W="$W -Werror=int-conversion -Werror=incompatible-pointer-types"

if [ $clang -eq 1 ]
then
	# Clang gives way too many spurious warnings about LOCAL_HIREMAINDER
	#W="-Weverything -Wno-shadow"
	W="$W -Wno-padded -Wno-vla -Wno-reserved-id-macro -Wno-documentation-unknown-command"
fi

# Not used: some warnings I have found to be unhelpful.
ignore="-Walloca -Wbad-function-cast -Wsign-conversion -Wswitch-default -Wunsuffixed-float-constants -Wjump-misses-init -Wunsafe-loop-optimizations -Wdouble-promotion -Wsystem-headers -Wnoexcept -Wdeclaration-after-statement -Wpedantic -Wtraditional-conversion -Wtraditional"

# Warnings which do not apply to C (but to C++, Objective C++, D, etc.).
ignore="$ignore -Wuseless-cast -Wzero-as-null-pointer-constant -Wreorder -Wcast-result -Wabi-tag -Wsign-promo -Wnon-virtual-dtor -Wctor-dtor-privacy -Wctor-dtor-privacy -Wdelete-non-virtual-dtor -Wassign-intercept -Wdelete-non-virtual-dtor -Wc++0x-compat -Wundeclared-selector -Wold-style-cast -Wselector -Woverloaded-virtual -Wstrict-null-sentinel -Wstrict-selector-match -Wsynth -Wc++11-compat -Wc++1z-compat -Wc++14-compat -Wregister"

# C++/Objective-C++ only
ignore="$ignore -Wconditionally-supported"

# C++ only.
ignore="$ignore -Wpessimizing-move -Weffc++ -Wc++0x-compat"

# Fortran only.
ignore="$ignore -Waliasing -Wampersand -Warray-temporaries -Wc-binding-type -Wcharacter-truncation -Wcompare-reals -Wconversion"
ignore="$ignore -Winteger-division -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation"
ignore="$ignore -Wreal-q-constant -Wsurprising -Wtabs -Wtarget-lifetime -Wundefined-do-loop -Wunused -Wunused-dummy-argument -Wuse-without-only -Wzerotrip"

# Basic optimizations.
O='-O2'

if [ $GCC -eq 1 ]
then
	O="$O -fsplit-loops"
fi

# Non-optimization I'm hiding with optimizations for some reason.
O="$O -fno-strict-aliasing -fomit-frame-pointer -fPIC"

# Powerful Graphite optimizations.
O="$O -ftree-loop-linear -floop-interchange -floop-strip-mine -floop-block -fgraphite-identity -floop-nest-optimize -floop-unroll-and-jam"

# All options together here
ALLOPT="-march=native -m64 $W $O -fdiagnostics-color=always"

# What warnings are missed?
#$CC $ALLOPT -Q --help=warnings | fgrep '[disabled]' | egrep -v `echo "\s($ignore)\s" | tr ' ' '|' | sed s/+/[+]/g`

# What optimizations are turned off?
#$CC $ALLOPT -Q --help=optimizers | fgrep '[disabled]'

# Compile and link auto file
echo Compiling
#$CC -ftime-report -o auto.gp.o $ALLOPT -c -I"/usr/local/include" auto.gp.c
$CC -o auto.gp.o $ALLOPT -c -I"/usr/local/include" auto.gp.c || exit # exit with status of last command on failure
echo Linking
$CC -o auto.gp.so $ALLOPT -shared -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari || exit # exit with status of last command on failure

# Create run file to install auto and add associated help entries
egrep '^GP;' auto.gp.c | sed 's/^GP;//' > auto.gp.run

## Create documentation -- takes a long time
#doxygen Doxyfile
#{ for f in html/*.png ; do pngout $f ; done } &

# Run gp itself

mv auto.gp.run ../pari
mv auto.gp.so ../pari
cd ../pari
gp
cd ../PARI-extensions
