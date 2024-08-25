// Cooley-Tukey : DFT de taille 2^N, Cooley-Tukey - Radix2 - DIT, out-of-place

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#define PI 3.141592653589793238462

// Table des racines 2^N iemes de l'unite (jusqu'a FFT de 65536 points)
static complex double w_p2[16];     // [0] indefini

static void init_onetime()
{
    for (int s = 1; s < 16; s++) {
        int const m = 1 << s;
        w_p2[s] = cexp(-2.0 * I * PI / m);
    }
}



static int ipow (int b, int n)         // pow(b,n)
{
    assert (b == 2);
    return 1 << n;
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

    // Shuffle X vers Y
    base2_reverse_copy(X, Y, N, N2);

    // La suite est in-place dans Y
    assert(N2 <= 16);
    for (int s = 1; s <= N2; s++) {
        int const m = 1 << s;
        complex double wm = w_p2[s];        // = cexp(-2.0 * I * PI / 2^s);
        for (int k = 0; k < N; k += m) {
            int const hm = m / 2;           // At least 1
            complex double w = 1.0; {
                complex double t = w * Y[k + 0 + 1 * hm];
                complex double u =     Y[k + 0 + 0 * hm];
                Y[k + 0 + 0 * hm] = u + t;
                Y[k + 0 + 1 * hm] = u - t;
                w *= wm;
            }
            for (int j = 1; j < hm; j++) {
                complex double t = w * Y[k + j + 1 * hm];
                complex double u =     Y[k + j + 0 * hm];
                Y[k + j + 0 * hm] = u + t;
                Y[k + j + 1 * hm] = u - t;
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

    init_onetime();

    int const N2 = (argc > 1) ? atoi(argv[1]) : 1;
    int const N  = ipow(2,N2);

    printf("fprintf('N = %d = 2^%d\\n');\n", N, N2);
    
    // Exemple de signal de taille N
    complex double X[N];
    for (int i = 0; i < N; i++) {
        X[i] =    (1.0*sin(1+i*0.11) + 0.6*sin(1+i*0.21))
              + I*(1.0*sin(2-i*0.22) + 0.7*sin(2+i*0.22));
    }

    printf("%% Entree:\n");
    print_complex_array("X", X, N);

    complex double Y[N];
    fft_base2(X, Y, N, N2);

    printf("\n%% Sortie (FFT):\n");
    print_complex_array("Y", Y, N);

    printf("\n");
    printf("F = fft(X);\n");
    printf("e = max(abs(Y-F))\n");
    printf("if (e > 1e-6); fprintf('KO\\n');quit(1); else fprintf('OK\\n');quit(0); end\n");

    return 0;
}

