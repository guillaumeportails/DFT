

CC      = gcc
CFLAGS  = -std=c99 -Wall -Wextra -O2 -g2
LDLIBS  = -lm


all: test

clean:
	/bin/rm -f dft.o dft check.m

test: dft
	./dft  2  0 > check.m && octave --no-gui check.m
	./dft  8  0 > check.m && octave --no-gui check.m
	./dft 13  0 > check.m && octave --no-gui check.m
	./dft  0  1 > check.m && octave --no-gui check.m
	./dft  0  3 > check.m && octave --no-gui check.m
	./dft  0  4 > check.m && octave --no-gui check.m
	./dft  3  2 > check.m && octave --no-gui check.m


dft: dft.o
dft.o: dft.c
