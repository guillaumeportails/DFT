
# -ftree-vectorize est nuisible, comme -O3
# clang moins bien

CC      = gcc
CFLAGS  = -std=c99 -Wall -Wextra -g2 -O2
LDLIBS  = -lm


all: test

clean:
	/bin/rm -f *.o dft fft check.m tmp* a.out

test: fft dft dft12.py
	python3 dft12.py > tmp.m && octave --no-gui tmp.m
	@echo
	./fft  2    > check.m && octave --no-gui check.m
	./fft 11    > check.m && octave --no-gui check.m
	@echo
	./dft  2  0 > check.m && octave --no-gui check.m
	./dft  8  0 > check.m && octave --no-gui check.m
	./dft 15  0 > check.m && octave --no-gui check.m
	./dft  0  1 > check.m && octave --no-gui check.m
	./dft  0  3 > check.m && octave --no-gui check.m
	./dft  0  4 > check.m && octave --no-gui check.m
	@echo
	./dft  2 -1 > check.m && octave --no-gui check.m
	./dft  4 -1 > check.m && octave --no-gui check.m
	./dft 15 -1 > check.m && octave --no-gui check.m
	@echo
	./dft  3  2 > check.m && octave --no-gui check.m



fft: fft.o
fft.o: fft.c

dft: dft.o
dft.o: dft.c
