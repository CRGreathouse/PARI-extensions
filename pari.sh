#!/bin/bash
cd ~/mth/pari

# Compile auto file
GCCOPT='-O2 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC'
gcc -o auto.gp.o $GCCOPT -c -I"/usr/local/include" auto.gp.c
gcc -o auto.gp.so $GCCOPT -shared -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari

# Create run file to install auto and add associated help entries
egrep '^GP;' auto.gp.c | sed 's/^GP;//' > auto.gp.run

## Create documentation -- takes a long time
#doxygen Doxyfile
#{ for f in html/*.png ; do pngout $f ; done } &

# Run gp itself
gp
