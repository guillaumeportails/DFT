
# -ftree-vectorize est nuisible, comme -O3
# clang moins bien

CC              = gcc -std=c99
CXX             = g++ -std=c++11
CPPFLAGS        = -D_XOPEN_SOURCE  # -DNDEBUG
CFLAGS          = -Wall -Wextra -g2 -O0
CXXFLAGS        = -Wall -Wextra -g2 -O0
LDLIBS          = -lm


all: test

clean:
	/bin/rm -f *.o dft fft check.m tmp* a.out dft.asm

test: fft dft dft.asm dft12.py
	python3 dft12.py > tmp.m && octave --no-gui tmp.m
	@echo
	./fft  2 0 0   > check.m && octave --no-gui check.m
	./fft 11 0 0   > check.m && octave --no-gui check.m
	@echo
	./dft -2  1     > check.m && octave --no-gui check.m
	./dft -2  2     > check.m && octave --no-gui check.m
	./dft -2  8     > check.m && octave --no-gui check.m
	./dft -2  2 -S  > check.m && octave --no-gui check.m
	./dft -2  8 -S  > check.m && octave --no-gui check.m
	./dft -2 15     > check.m && octave --no-gui check.m
	./dft -2 15 -S  > check.m && octave --no-gui check.m
	./dft -2 16     > check.m && octave --no-gui check.m
	./dft -2 16 -S  > check.m && octave --no-gui check.m
	./dft -3  1     > check.m && octave --no-gui check.m
	./dft -3  3     > check.m && octave --no-gui check.m
	./dft -3  4     > check.m && octave --no-gui check.m
	@echo
	./dft -8  1   > check.m && octave --no-gui check.m
	./dft -8  2   > check.m && octave --no-gui check.m
	./dft -8  3   > check.m && octave --no-gui check.m
	./dft -8  4   > check.m && octave --no-gui check.m

#	@echo
#	./dft  2 -1 > check.m && octave --no-gui check.m
#	./dft  4 -1 > check.m && octave --no-gui check.m
#	./dft 15 -1 > check.m && octave --no-gui check.m
#	@echo
#	./dft  3  2 > check.m && octave --no-gui check.m

asm: dft.asm


fft: fft.o
fft.o: fft.c

dft: dft.o
dft.o: dft.cpp

# Code genere : gcc -S est moche. Il vaut mieux un objdump -S :
dft.asm: dft.o
	objdump -S -C -M intel $< > $@

