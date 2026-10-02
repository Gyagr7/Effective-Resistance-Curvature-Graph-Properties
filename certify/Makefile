# certify/ -- the exact RN/RP/SRN certifier used by ../resistance.py
# Needs GLPK and GMP:  glpk gmp (Arch) / libglpk-dev libgmp-dev (Debian).

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c11
LDLIBS   = -lglpk -lgmp -lm

rn: rn.c
	$(CC) $(CFLAGS) -o $@ rn.c $(LDLIBS)

check: rn
	./rn --exact examples.json

clean:
	rm -f rn

.PHONY: check clean
