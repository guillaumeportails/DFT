// Cooley-Tukey : DFT de taille 2^N, out-of-place

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.141592653589793238462

static int ipow (int b, int n)         // pow(b,n)
{
  int r = 1;
  while (n-- > 0) r *= b;
  return r;
}

// Permuter les elements selon l'ordre de bit inverse des index
static void base2_reverse_copy(const complex double *source, complex double *dest, int N, int N2) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < N2; j++) {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        dest[rev] = source[i];
    }
}

// FFT pour N = 2^N2
static void fft_base2(const complex double *X, complex double *Y, int N, int N2) {
    base2_reverse_copy(X, Y, N, N2);

    for (int s = 1; s <= N2; s++) {
        int const m = 1 << s;
        complex double wm = cexp(-2.0 * I * PI / m);
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < m / 2; j++) {
                complex double t = w * Y[k + j + 1 * m / 2];
                complex double u =     Y[k + j + 0 * m / 2];
                Y[k + j + 0 * m / 2] = u + t;
                Y[k + j + 1 * m / 2] = u - t;
                w *= wm;
            }
        }
    }
}


//----------------------------------------------------------------

// Fonction d'affichage du resultat
void print_complex_array(const char *v, const complex double *X, int N) {
    for (int i = 0; i < N; i++) {
        printf("%s(%2d) = %18.15f + %18.15fi;\n", v, i+1, creal(X[i]), cimag(X[i]));
    }
}

int main(int argc, char *argv[]) {
    int const N2 = (argc > 1) ? atoi(argv[1]) : 1;
    int const N  = ipow(2,N2);

    printf("fprintf('N = %d = 2^%d\\n');\n", N, N2);
    
    // Exemple de signal de taille N
    complex double x[N];
    for (int i = 0; i < N; i++) {
        x[i] =    (1.0*sin(1+i*0.11) + 0.6*sin(1+i*0.21))
              + I*(1.0*sin(2-i*0.22) + 0.7*sin(2+i*0.22));
    }

    printf("%% Entree:\n");
    print_complex_array("X", x, N);

    complex double y[N];
    fft_base2(x, y, N, N2);

    printf("\n%% Sortie (FFT):\n");
    print_complex_array("Y", y, N);

    printf("\n");
    printf("F = fft(X);\n");
    printf("e = max(abs(Y-F))\n");
    printf("if (e > 1e-6); fprintf('KO\\n');quit(1); else fprintf('OK\\n');quit(0); end\n");

    return 0;
}

