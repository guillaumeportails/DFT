
# -ftree-vectorize est nuisible, comme -O3
# clang moins bien

CC              = gcc -std=c99
CXX             = g++ -std=c++11
CPPFLAGS        = -D_XOPEN_SOURCE  # -DNDEBUG
optim			= -O2 -ffast-math -funroll-loops
CFLAGS          = -Wall -Wextra -g2 $(optim)
CXXFLAGS        = -Wall -Wextra -g2 $(optim)
LDLIBS          = -lm

SHELL           = bash	# {a..b}

all: test asm

clean:
	/bin/rm -f *.o dft fft check.m tmp* a.out *.asm *.log

test: fft dft dft12.py
	python3 dft12.py > tmp.m && octave --no-gui tmp.m
	rm -f fft.log dft.log
	@echo
	./fft  0    > check.m && octave --no-gui check.m | tee -a fft.log
	./fft  1    > check.m && octave --no-gui check.m | tee -a fft.log
	./fft  2    > check.m && octave --no-gui check.m | tee -a fft.log
	./fft 11    > check.m && octave --no-gui check.m | tee -a fft.log
	./fft 12    > check.m && octave --no-gui check.m | tee -a fft.log
	@echo
	set -e; for p in {1..16}; do \
		./dft -2 $$p     > check.m && octave --no-gui check.m | tee -a dft.log; \
		./dft -2 $$p -S  > check.m && octave --no-gui check.m | tee -a dft.log; \
		done
	@echo
	set -e; for p in {1..5}; do \
		./dft -3 $$p     > check.m && octave --no-gui check.m | tee -a dft.log; \
		done
	@echo
	set -e; for p in {1..8}; do \
		./dft -4 $$p     > check.m && octave --no-gui check.m | tee -a dft.log; \
		done
	@echo
	set -e; for p in {1..5}; do \
		./dft -8 $$p     > check.m && octave --no-gui check.m | tee -a dft.log; \
		done
	@echo; cat fft.log
	@echo; sort -n -k3,13 dft.log

#	@echo
#	./dft  2 -1 > check.m && octave --no-gui check.m
#	./dft  4 -1 > check.m && octave --no-gui check.m
#	./dft 15 -1 > check.m && octave --no-gui check.m
#	@echo
#	./dft  3  2 > check.m && octave --no-gui check.m

asm: dft4.asm dft8.asm


fft: fft.o
fft.o: fft.c

dft: dft.o dft2.o dft3.o dft4.o dft8.o
dft.o: dft.cpp dft.h
dft2.o: dft2.cpp dft.h
dft3.o: dft3.cpp dft.h
dft4.o: dft4.cpp dft.h
dft8.o: dft8.cpp dft.h


# Code genere : gcc -S est moche. Il vaut mieux un objdump -S :

dft4.asm: dft4.o
	objdump -S -C -M intel $< > $@

dft8.asm: dft8.o
	objdump -S -C -M intel $< > $@

