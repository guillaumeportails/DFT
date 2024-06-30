// Cooley-Tukey + Rader : DFT de taille 2^P2 * 3^P3

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <assert.h>

#define PI 3.141592653589793238462

static int ipow (int b, int n)         // pow(b,n)
{
  int r = 1;
  while (n-- > 0) r *= b;
  return r;
}

// Fonction pour permuter les elements selon l'ordre de bit inverse
static void base2_reverse_copy(complex double *source, complex double *dest, int N, int P2) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < P2; j++) {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        assert((rev >= 0) && (rev < N));
        dest[rev] = source[i];
    }
}

// Fonction pour permuter les elements selon l'ordre de base-3 inversee
static void base3_reverse_copy(complex double *source, complex double *dest, int N, int P3) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        int tmp = i;
        for (int j = 0; j < P3; j++) {
            rev = rev * 3 + tmp % 3;
            tmp /= 3;
        }
        assert((rev >= 0) && (rev < N));
        dest[rev] = source[i];
    }
}

// Fonction pour permuter les elements selon l'ordre de base-8 inversee
static void base8_reverse_copy(complex double *source, complex double *dest, int N, int P8) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        int tmp = i;
        for (int j = 0; j < P8; j++) {
            rev = rev * 8 + tmp % 8;
            tmp /= 8;
        }
        assert((rev >= 0) && (rev < N));
        dest[rev] = source[i];
    }
}

// FFT pour N = 2^P2
static void fft_base2(complex double *X, int N, int P2) {
    complex double temp[N];
    base2_reverse_copy(X, temp, N, P2);

    for (int s = 1; s <= P2; s++) {
        int m = 1 << s;
        complex double wm = cexp(-2.0 * I * PI / m);
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < m / 2; j++) {
                complex double t = w * temp[k + j + m / 2];
                complex double u = temp[k + j];
                temp[k + j] = u + t;
                temp[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}

// FFT pour N = 3^P3
static void fft_base3(complex double *X, int N, int P3) {
    complex double temp[N];
    base3_reverse_copy(X, temp, N, P3);

    for (int s = 1; s <= P3; s++) {
        int m = ipow(3, s);
        int h = m / 3;
        complex double wm = cexp(-2.0 * I * PI / m);
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < h; j++) {
                complex double const t0 = temp[k + j + 0*h];
                complex double const t1 = temp[k + j + 1*h] * w;
                complex double const t2 = temp[k + j + 2*h] * w*w;
                temp[k + j + 0*h] = t0 + t1 + t2;
                temp[k + j + 1*h] = t0 + t1*cexp(-2*1*1*I*PI/3)
                                       + t2*cexp(-2*2*1*I*PI/3);
                temp[k + j + 2*h] = t0 + t1*cexp(-2*1*2*I*PI/3)
                                       + t2*cexp(-2*2*2*I*PI/3);
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}

// FFT pour N = 8^P8
static void fft_base8(complex double *X, int N, int P8) {
    complex double temp[N];
    base8_reverse_copy(X, temp, N, P8);

    for (int s = 1; s <= P8; s++) {
        int m = ipow(8, s);
        complex double wm = cexp(-2.0 * I * PI / m);
        assert((m % 8) == 0);
        int h = m / 8;
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < h; j++) {
                complex double const t0 = temp[k + j + 0 * h];
                complex double const t1 = temp[k + j + 1 * h] * w;
                complex double const t2 = temp[k + j + 2 * h] * w*w;
                complex double const t3 = temp[k + j + 3 * h] * w*w*w;
                complex double const t4 = temp[k + j + 4 * h] * w*w*w*w;
                complex double const t5 = temp[k + j + 5 * h] * w*w*w*w*w;
                complex double const t6 = temp[k + j + 6 * h] * w*w*w*w*w*w;
                complex double const t7 = temp[k + j + 7 * h] * w*w*w*w*w*w*w;
                temp[k + j + 0 * h] =   t0 * cexp(-2*0*0*I*PI/8)
                                      + t1 * cexp(-2*0*1*I*PI/8)
                                      + t2 * cexp(-2*0*2*I*PI/8)
                                      + t3 * cexp(-2*0*3*I*PI/8)
                                      + t4 * cexp(-2*0*4*I*PI/8)
                                      + t5 * cexp(-2*0*5*I*PI/8)
                                      + t6 * cexp(-2*0*6*I*PI/8)
                                      + t7 * cexp(-2*0*7*I*PI/8);
                temp[k + j + 1 * h] =   t0 * cexp(-2*1*0*I*PI/8)
                                      + t1 * cexp(-2*1*1*I*PI/8)
                                      + t2 * cexp(-2*1*2*I*PI/8)
                                      + t3 * cexp(-2*1*3*I*PI/8)
                                      + t4 * cexp(-2*1*4*I*PI/8)
                                      + t5 * cexp(-2*1*5*I*PI/8)
                                      + t6 * cexp(-2*1*6*I*PI/8)
                                      + t7 * cexp(-2*1*7*I*PI/8);
                temp[k + j + 2 * h] =   t0 * cexp(-2*2*0*I*PI/8)
                                      + t1 * cexp(-2*2*1*I*PI/8)
                                      + t2 * cexp(-2*2*2*I*PI/8)
                                      + t3 * cexp(-2*2*3*I*PI/8)
                                      + t4 * cexp(-2*2*4*I*PI/8)
                                      + t5 * cexp(-2*2*5*I*PI/8)
                                      + t6 * cexp(-2*2*6*I*PI/8)
                                      + t7 * cexp(-2*2*7*I*PI/8);
                temp[k + j + 3 * h] =   t0 * cexp(-2*3*0*I*PI/8)
                                      + t1 * cexp(-2*3*1*I*PI/8)
                                      + t2 * cexp(-2*3*2*I*PI/8)
                                      + t3 * cexp(-2*3*3*I*PI/8)
                                      + t4 * cexp(-2*3*4*I*PI/8)
                                      + t5 * cexp(-2*3*5*I*PI/8)
                                      + t6 * cexp(-2*3*6*I*PI/8)
                                      + t7 * cexp(-2*3*7*I*PI/8);
                temp[k + j + 4 * h] =   t0 * cexp(-2*4*0*I*PI/8)
                                      + t1 * cexp(-2*4*1*I*PI/8)
                                      + t2 * cexp(-2*4*2*I*PI/8)
                                      + t3 * cexp(-2*4*3*I*PI/8)
                                      + t4 * cexp(-2*4*4*I*PI/8)
                                      + t5 * cexp(-2*4*5*I*PI/8)
                                      + t6 * cexp(-2*4*6*I*PI/8)
                                      + t7 * cexp(-2*4*7*I*PI/8);
                temp[k + j + 5 * h] =   t0 * cexp(-2*5*0*I*PI/8)
                                      + t1 * cexp(-2*5*1*I*PI/8)
                                      + t2 * cexp(-2*5*2*I*PI/8)
                                      + t3 * cexp(-2*5*3*I*PI/8)
                                      + t4 * cexp(-2*5*4*I*PI/8)
                                      + t5 * cexp(-2*5*5*I*PI/8)
                                      + t6 * cexp(-2*5*6*I*PI/8)
                                      + t7 * cexp(-2*5*7*I*PI/8);
                temp[k + j + 6 * h] =   t0 * cexp(-2*6*0*I*PI/8)
                                      + t1 * cexp(-2*6*1*I*PI/8)
                                      + t2 * cexp(-2*6*2*I*PI/8)
                                      + t3 * cexp(-2*6*3*I*PI/8)
                                      + t4 * cexp(-2*6*4*I*PI/8)
                                      + t5 * cexp(-2*6*5*I*PI/8)
                                      + t6 * cexp(-2*6*6*I*PI/8)
                                      + t7 * cexp(-2*6*7*I*PI/8);
                temp[k + j + 7 * h] =   t0 * cexp(-2*7*0*I*PI/8)
                                      + t1 * cexp(-2*7*1*I*PI/8)
                                      + t2 * cexp(-2*7*2*I*PI/8)
                                      + t3 * cexp(-2*7*3*I*PI/8)
                                      + t4 * cexp(-2*7*4*I*PI/8)
                                      + t5 * cexp(-2*7*5*I*PI/8)
                                      + t6 * cexp(-2*7*6*I*PI/8)
                                      + t7 * cexp(-2*7*7*I*PI/8);
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}



// FFT hybride pour N = 2^P2 * 3^P3
void fft_hybrid(complex double *X, int N, int P2, int P3) {
    if (P2 > 0) {
        int const size2 = ipow(2, P2);
        complex double subarray2[size2];

        for (int i = 0; i < N; i += size2) {
            for (int j = 0; j < size2; j++) subarray2[j] = X[i + j];
            fft_base2(subarray2, size2, P2);
            for (int j = 0; j < size2; j++) X[i + j] = subarray2[j];
        }
    }

    if (P3 > 0) {
        int const size3 = ipow(3, P3);
        complex double subarray3[size3];

        for (int i = 0; i < N; i += size3) {
            for (int j = 0; j < size3; j++) subarray3[j] = X[i + j];
            fft_base3(subarray3, size3, P3);
            for (int j = 0; j < size3; j++) X[i + j] = subarray3[j];
        }
    }
}

// Stockham FFT
// Cf http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham2.html
static void stockham (int n, int s, bool eo, complex double *x, complex double *y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    while (n >= 4)
    {
        const int m = n/2;
        const double theta0 = 2*PI/n;
        for (int p = 0; p < m; p++) {
            const complex double wp = cos(p*theta0) -I*sin(p*theta0);
            assert( (s == 1) || ((s%2) == 0) ); // Par q+=2 pour favoriser SIMD
            for (int q = 0; q < s; q++) {
                const complex double  a = x[q + s*(p + 0)]; // bi-sequentiel en read
                const complex double  b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;                // bi-sequentiel en write
                y[q + s*(2*p + 1)] = (a - b) * wp;
            }
        }
        n = m;
        s = 2*s;
        eo = !eo;
        complex double *z = y; y = x; x = z;
    }

    assert(n == 2);
    {
        complex double * z = eo ? y : x;
        for (int q = 0; q < s; q++) {
            const complex double  a = x[q + 0];
            const complex double  b = x[q + s];
            z[q + 0] = a + b;
            z[q + s] = a - b;
        }
    }
}

static void fft_stockham_base2 (complex double *X, int N)
{
    complex double Y[N];
    stockham(N, 1, false, X, Y);
}



// Fonction d'affichage du resultat
void print_complex_array(const char *v, complex double *X, int N) {
    for (int i = 0; i < N; i++) {
        printf("%s(%2d) = %18.15f + %18.15fi;\n", v, i+1, creal(X[i]), cimag(X[i]));
    }
}

static double urand()
{
  double const r = rand();
  double const k = RAND_MAX *0.5;
  return r/k - 1.0;
}

int main(int argc, char *argv[])
{
    int const P2 = (argc > 1) ? atoi(argv[1]) : 0;
    int const P3 = (argc > 2) ? atoi(argv[2]) : 0;
    int const P8 = (argc > 3) ? atoi(argv[3]) : 0;
    assert((P2 == 0) || (P8 == 0));

    int const N  = ipow(2,P2) * ipow(3,P3) * ipow(8,P8);

    printf("fprintf('N = %d\\n');\n", N);

    // Exemple de signal de taille N
    complex double x[N];
    for (int i = 0; i < N; i++) x[i] = urand() + I*urand();

    printf("%% Entree:\n");
    print_complex_array("X", x, N);

    clock_t tic = clock();
    if      (P2 > 0)
      fft_base2(x, N, P2);
//    fft_stockham_base2(x, N);
    else if (P3 > 0)
      fft_base3(x, N, P3);
    else if (P8 > 0)
      fft_base8(x, N, P8);
//  else
//    fft_hybrid(x, N, P2, P3);
    clock_t toc = clock() - tic;

    printf("\n%% Sortie (FFT):  %.6fs clock\n", (double) toc/CLOCKS_PER_SEC);
    print_complex_array("Y", x, N);

    printf("\n");
    if (N > 10000) printf("\nclocks = %.6f\n", (double) toc/CLOCKS_PER_SEC);
    printf("F = fft(X);\n");
    printf("err = max(abs(Y-F))\n");
    printf("if (err > 1e-6); fprintf('KO\\n'); else fprintf('OK\\n'); end\n");

    return 0;
}

