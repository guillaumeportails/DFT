# Colonnes de dft.log
#        3            5                 7       9        11         13
#   N = 16  OK 3.51e-15,  fft_base2,  199 io,  43 zmul,  64 zadd,  386 flop


flop(zmul,zadd) = zmul*6 + zadd*2

set style data linespoints
#set logscale x 2

#set y2tics
#plot '<grep -w fft_base2 dft.log' u 3:7 t 'FFT2-rw', '' u 3:9 t 'FFT2-z*' axes x1y2, '' u 3:11 t 'FFT2-z+' axes x1y2, \
#     '<grep -w fft_base4 dft.log' u 3:7 t 'FFT4-rw', '' u 3:9 t 'FFT4-z*' axes x1y2, '' u 3:11 t 'FFT4-z+' axes x1y2, \
#     '<grep -w fft_base8 dft.log' u 3:7 t 'FFT8-rw', '' u 3:9 t 'FFT8-z*' axes x1y2, '' u 3:11 t 'FFT8-z+' axes x1y2

set multiplot layout 1,2

set title "Memory Read/Write"
plot '<grep -w fft_base2 dft.log' u 3:7 t 'FFT2-rw', \
     '<grep -w fft_base4 dft.log' u 3:7 t 'FFT4-rw',  \
     '<grep -w fft_base8 dft.log' u 3:7 t 'FFT8-rw'

set title "FLOPs"
plot '<grep -w fft_base2 dft.log' u 3:(flop($9,$11)) t 'FFT2', \
     '<grep -w fft_base4 dft.log' u 3:(flop($9,$11)) t 'FFT4',  \
     '<grep -w fft_base8 dft.log' u 3:(flop($9,$11)) t 'FFT8'

#unset multiplot
pause -1
