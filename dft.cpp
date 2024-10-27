// Cooley-Tukey + Rader : DFT de taille 2^P2 * 3^P3
//
// Cf https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
//
// "... for more moderate n, FFTW uses depth-first recursion with a bounded radix (...)
//  but with much larger radices (radix 32 is common) and base cases (size 32 or 64 is common)"
// => Cas de base generable avc un template recursif ?
//
// Bit reversal algorithms : Cf https://upcommons.upc.edu/bitstream/handle/2117/1705/Porrata.pdf
// + Gold-Rader
// + Rodriguez

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "dft.h"



struct Counters C;


// Fonction d'affichage de tableau, syntaxe Octave
void print_complex_array(const char *v, cplx_t *X, myuint N)
{
    printf("%s = zeros(%d,1);\n", v, N);
    for (myuint i = 0; i < N; i++) {
        printf("%s(%2d) = %18.15f + %18.15fi;\n", v, i+1, std::real(X[i]), std::imag(X[i]));
    }
}






// FFT hybride pour N = 2^P2 * 3^P3
void fft_hybrid(cplx_t *X, myuint N, myuint P2, myuint P3)
{
    C.name(__func__);
    assert((P2 >= 1) && (P3 >= 1));

    cplx_t X2[N];
    cplx_t X3[N];    
    for(myuint i = 0; i < N; i++) X2[i] = X[i];
    for(myuint i = 0; i < N; i++) X3[i] = X[i];
    
    myuint const size2 = ipow(2, P2);
    for (myuint i = 0; i < N; i += size2)
    {
        cplx_t subarray2[size2];
        for (myuint j = 0; j < size2; j++) subarray2[j] = X2[i + j];
        fft_base2 (subarray2, size2, P2);
        for (myuint j = 0; j < size2; j++) X2[i + j] = subarray2[j];
    }

    myuint const size3 = ipow(3, P3);
    for (myuint i = 0; i < N; i += size3)
    {
        cplx_t subarray3[size3];
        for (myuint j = 0; j < size3; j++) subarray3[j] = X3[i + j];
        fft_base3 (subarray3, size3, P3);
        for (myuint j = 0; j < size3; j++) X3[i + j] = subarray3[j];
    }

    // for (myuint k = 0; k < N; k++)
    // {
    //     X[k] = 0;
    //     for (myuint p = 0; p < 2; p++) for (myuint q = 0; q < 3; q++)
    //     {
    //         X[k] += twiddle_factors[(p * q * k) % N] * subfft_p[q * 2 + p] * subfft_q[p * 3 + q];
    //     }
    // }

    //    for (int k = 0; k < a; k++) {
    //         double complex w2 = cexp(-2 * I * M_PI * k / a);
    //         double complex X1k = X1[k] + w1 * X1[k + a / 2];
    //         double complex X2k = X2[k] + w2 * X2[k + b / 2];
    //         X[k    ] = X1k + w2 * X2k;
    //         X[k + a] = X1k - w2 * X2k;
    //     }


}




static double urand()
{
  double const r = rand();
  double const k = RAND_MAX *0.5;
  return r/k - 1.0;
  
}

int main(int argc, char *argv[])
{
    myuint P2 = 0, P3 = 0, P4 = 0, P8 = 0;
    int o;
    void (*opt_fft_base2)(cplx_t *x, myuint N, myuint P2) = fft_base2;
//  bool opt_false = false;     //  const en fait, pour battre -O2
    while ((o = getopt(argc,argv,"2:3:4:8:SRn")) > 0) switch (o)
    {
        case '2' : P2 = atoi(optarg); break;
        case '3' : P3 = atoi(optarg); break;
        case '4' : P4 = atoi(optarg); break;
        case '8' : P8 = atoi(optarg); break;
        case 'S' : opt_fft_base2 = fft_stockham_base2; break;
        case 'R' : opt_fft_base2 = fft_base2_r4; break;
//      case 'n' : opt_false = true; break; // Option jamais donnee. Battre l'optimiseur.
    }
    assert((P2 > 0) || (P3 > 0) || (P4 > 0)|| (P8 > 0));
    assert((P2 == 0) || ((P4 == 0) && (P8 == 0)));
    assert((P4 == 0) || ((P2 == 0) && (P8 == 0)));
    assert((P8 == 0) || ((P2 == 0) && (P4 == 0)));

    myuint const N  = ipow(2,P2) * ipow(3,P3) * ipow(4,P4)* ipow(8,P8);
    printf("%% fprintf('N = %d\\n');\n", N);
    printf("isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;\n");

    // sanity check
    assert(ipow(7,0) ==  1);
    assert(ipow(5,1) ==  5);
    assert(ipow(3,2) ==  9);
    assert(ipow(2,3) ==  8);
    assert(ipow(2,4) == 16);
//  assert(ipow(1,6) ==  1);

    // Exemple de signal de taille N
    cplx_t x[N];
    for (myuint i = 0; i < N; i++) x[i] = cplx_t(urand(), urand());

    printf("%% Entree:\n");
    print_complex_array("X", x, N);

    C.raz();
    clock_t tic = clock();
    if (P3 == 0)
    {
        if (P2 != 0)
            opt_fft_base2(x, N, P2);
        else if (P4 != 0)
            fft_base4(x, N, P4);
        else if (P8 != 0)
            fft_base8(x, N, P8);
    }
    else if ((P2 == 0) && (P4 == 0) && (P8 == 0))
      fft_base3(x, N, P3);
//  else
//    fft_hybrid(x, N, P2, P3);
    clock_t toc = clock() - tic;

    printf("\n%% Sortie (FFT):  %.6fs clock\n", (double) toc/CLOCKS_PER_SEC);
    print_complex_array("Y", x, N);

    printf("\n");
//  if (N > 10000) printf("\nclk_ms = %.3f %% ms\n", (double) 1e3 * toc/CLOCKS_PER_SEC);
    printf("F = fft(X);\n");
    printf("err = max(abs(Y-F));\n");
    printf("if (err > 1e-6); fprintf('KO %%.2e\\n',err); error('KO');"
           " else fprintf('N = %5d  OK %%.2e,  %20s, %7u io, %7u zmul, %7u zadd, %7u flop\\n',err); end\n",
                        N, C.fn, C.mio, C.zmul, C.zadd, 6*C.zmul + 2*C.zadd);

    return 0;
}

