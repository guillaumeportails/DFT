// Cooley-Tukey + Rader : DFT de taille 2^N2 * 3^N3

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
static void base2_reverse_copy(complex double *source, complex double *dest, int N, int N2) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < N2; j++) {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        dest[rev] = source[i];
    }
}

// Fonction pour permuter les elements selon l'ordre de base-3 inversee
static void base3_reverse_copy(complex double *source, complex double *dest, int N, int N3) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        int tmp = i;
        for (int j = 0; j < N3; j++) {
            rev = rev * 3 + tmp % 3;
            tmp /= 3;
        }
        dest[rev] = source[i];
    }
}

// FFT pour N = 2^N2
static void fft_cooleytukey_base2(complex double *X, int N, int N2) {
    complex double *temp = malloc(N * sizeof(complex double));
    base2_reverse_copy(X, temp, N, N2);

    for (int s = 1; s <= N2; s++) {
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

    for (int i = 0; i < N; i++) {
        X[i] = temp[i];
    }

    free(temp);
}

// FFT pour N = 3^N3
static void fft_base3(complex double *X, int N, int N3) {
    complex double *temp = malloc(N * sizeof(complex double));
    base3_reverse_copy(X, temp, N, N3);

    for (int s = 1; s <= N3; s++) {
        int m = ipow(3, s);
        complex double wm = cexp(-2.0 * I * PI / m);
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < m / 3; j++) {
                complex double const t1 =     w * temp[k + j + 1 * m / 3];
                complex double const t2 = w * w * temp[k + j + 2 * m / 3];
                complex double const u  = temp[k + j];
                temp[k + j + 0 * m / 3] = u + t1 + t2;
                temp[k + j + 1 * m / 3] = u + t1 * cexp(-2.0 * I * PI / 3) + t2 * cexp(-4.0 * I * PI / 3);
                temp[k + j + 2 * m / 3] = u + t1 * cexp(-4.0 * I * PI / 3) + t2 * cexp(-8.0 * I * PI / 3);
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) {
        X[i] = temp[i];
    }

    free(temp);
}

// FFT hybride pour N = 2^N2 * 3^N3
void fft_hybrid(complex double *X, int N, int N2, int N3) {
    if (N2 > 0) {
        int const size2 = ipow(2, N2);
        complex double subarray2[size2];

        for (int i = 0; i < N; i += size2) {
            for (int j = 0; j < size2; j++) {
                subarray2[j] = X[i + j];
            }
            fft_cooleytukey_base2(subarray2, size2, N2);
            for (int j = 0; j < size2; j++) {
                X[i + j] = subarray2[j];
            }
        }
    }

    if (N3 > 0) {
        int const size3 = ipow(3, N3);
        complex double subarray3[size3];

        for (int i = 0; i < N; i += size3) {
            for (int j = 0; j < size3; j++) {
                subarray3[j] = X[i + j];
            }
            fft_base3(subarray3, size3, N3);
            for (int j = 0; j < size3; j++) {
                X[i + j] = subarray3[j];
            }
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

int main(int argc, char *argv[]) {
    int const N2 = (argc > 1) ? atoi(argv[1]) : 1;
    int const N3 = (argc > 2) ? atoi(argv[2]) : 0;
    int const N  = ipow(2,N2) * ((N3 > 0) ? ipow(3,N3) : 1);

    printf("fprintf('N = %d\\n');\n", N);
    
    // Exemple de signal de taille N
    complex double x[N];
    for (int i = 0; i < N; i++) {
        x[i] = sin(1+i*0.11) + I*sin(2-i*0.22);
    }

    printf("%% Entree:\n");
    print_complex_array("X", x, N);

    clock_t tic = clock();
    if (N3 < 0)
      fft_stockham_base2(x, N);
    else
      fft_hybrid(x, N, N2, N3);
    clock_t toc = clock() - tic;

    printf("\n%% Sortie (FFT):  %.6fs clock\n", (double) toc/CLOCKS_PER_SEC);
    print_complex_array("Y", x, N);

    printf("\n");
    if (N > 10000) printf("\nclocks = %.6f\n", (double) toc/CLOCKS_PER_SEC);
    printf("F = fft(X);\n");
    printf("e = max(abs(Y-F))\n");
    printf("if (e > 1e-6); fprintf('KO\\n');quit(1); else fprintf('OK\\n');quit(0); end\n");

    return 0;
}

