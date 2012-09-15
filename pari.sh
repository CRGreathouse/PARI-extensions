#!/bin/bash
cd ~/mth/pari
gcc -o auto.gp.o  -O2 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -c -I"/usr/local/include" auto.gp.c
gcc -o auto.gp.so -O2 -Wall -Werror -fno-strict-aliasing -fomit-frame-pointer -fPIC -shared -Wl,-shared auto.gp.o -lc -lm -L"/usr/local/lib" -lpari
gp
