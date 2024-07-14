// Cooley-Tukey + Rader : DFT de taille 2^P2 * 3^P3
//
// Cf https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
//
// "hus, for more moderate n, FFTW uses depth-rst recursion with a bounded radix (...)
//  but with much larger radices (radix 32 is common) and base cases (size 32 or 64 is common)"
// => Cas de base generable avc un template recusrsif ?
//
// Bit reversal algorithms : Cf https://upcommons.upc.edu/bitstream/handle/2117/1705/Porrata.pdf
// + Gold-Rader
// + Rodriguez

#include "dft.h"


// FFT pour N = 3^P3
void fft_base3 (cplx_t *X, myuint N, myuint P3)
{
    C.name(__func__);
    constexpr myuint R = 3;
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P3);

    const cplx_t  tw_1_1 =             zexpiy(-2*PI/R * 1 * 1);
    const cplx_t  tw_2_1 =             zexpiy(-2*PI/R * 2 * 1);
    const cplx_t& tw_1_2 = tw_2_1;  // zexpiy(-2*PI/R * 1 * 2);
    const cplx_t& tw_2_2 = tw_1_1;  // zexpiy(-2*PI/R * 1 * 2);  // (2*2) % R = (1*1)

    for (myuint s = 1; s <= P3; s++) {
        myuint h = ipow(R, s-1);                               // Read from table ? mem IO is slow
        myuint m = h * R;
        cplx_t wm = zexpiy(-2.0 * PI / m);                  // read from table
        for (myuint k = 0; k < N; k += m) {
            cplx_t w = 1.0;
            for (myuint j = 0; j < h; j++) {
                cplx_t const t0 = temp[k + j + 0*h];
                cplx_t const t1 = temp[k + j + 1*h] * w;
                cplx_t const t2 = temp[k + j + 2*h] * w*w;  // read w*w from table 
                temp[k + j + 0*h] = t0 + t1 + t2;
                temp[k + j + 1*h] = t0 + t1*tw_1_1 + t2*tw_2_1;
                temp[k + j + 2*h] = t0 + t1*tw_1_2 + t2*tw_2_2;
                w *= wm;
                C.zmul += 7; C.zadd += 6; C.mio += 6;
            }
        }
    }

    for (myuint i = 0; i < N; i++) X[i] = temp[i];
    C.mio += 2*N;
}

