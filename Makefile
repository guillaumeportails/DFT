

CC      = gcc
CFLAGS  = -std=c99 -Wall -Wextra -O2 -g2
LDLIBS  = -lm


all: test

clean:
	/bin/rm -f *.o dft fft check.m

test: fft dft
	@echo
	./fft  0    > check.m && octave --no-gui check.m
	./fft  1    > check.m && octave --no-gui check.m
	./fft  2    > check.m && octave --no-gui check.m
	./fft 11    > check.m && octave --no-gui check.m
	./fft 12    > check.m && octave --no-gui check.m
	@echo
	./dft  2  0 > check.m && octave --no-gui check.m
	./dft  8  0 > check.m && octave --no-gui check.m
	./dft 13  0 > check.m && octave --no-gui check.m
	./dft  0  1 > check.m && octave --no-gui check.m
	./dft  0  3 > check.m && octave --no-gui check.m
	./dft  0  4 > check.m && octave --no-gui check.m
	./dft  3  2 > check.m && octave --no-gui check.m


fft: fft.o
fft.o: fft.c

dft: dft.o
dft.o: dft.c
