// Cooley-Tukey + Rader : DFT de taille 2^N2 * 3^N3

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.141592653589793238462

int ipow (int b, int n)         // pow(b,n)
{
  int r = 1;
  while (n-- > 0) r *= b;
  return r;
}

// Fonction pour permuter les éléments selon l'ordre de bit inversé
void base2_reverse_copy(complex double *source, complex double *dest, int N, int N2) {
    for (int i = 0; i < N; i++) {
        int rev = 0;
        for (int j = 0; j < N2; j++) {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        dest[rev] = source[i];
    }
}

// Fonction pour permuter les éléments selon l'ordre de base-3 inversée
void base3_reverse_copy(complex double *source, complex double *dest, int N, int N3) {
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
void fft_base2(complex double *X, int N, int N2) {
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
void fft_base3(complex double *X, int N, int N3) {
    complex double *temp = malloc(N * sizeof(complex double));
    base3_reverse_copy(X, temp, N, N3);

    for (int s = 1; s <= N3; s++) {
        int m = ipow(3, s);
        complex double wm = cexp(-2.0 * I * PI / m);
        for (int k = 0; k < N; k += m) {
            complex double w = 1.0;
            for (int j = 0; j < m / 3; j++) {
                complex double t1 = w * temp[k + j + m / 3];
                complex double t2 = w * w * temp[k + j + 2 * m / 3];
                complex double u = temp[k + j];
                temp[k + j] = u + t1 + t2;
                temp[k + j + m / 3] = u + t1 * cexp(-2.0 * I * PI / 3) + t2 * cexp(-4.0 * I * PI / 3);
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
        int size2 = ipow(2, N2);
        complex double *subarray2 = malloc(size2 * sizeof(complex double));

        for (int i = 0; i < N; i += size2) {
            for (int j = 0; j < size2; j++) {
                subarray2[j] = X[i + j];
            }
            fft_base2(subarray2, size2, N2);
            for (int j = 0; j < size2; j++) {
                X[i + j] = subarray2[j];
            }
        }

        free(subarray2);
    }

    if (N3 > 0) {
        int size3 = ipow(3, N3);
        complex double *subarray3 = malloc(size3 * sizeof(complex double));

        for (int i = 0; i < N; i += size3) {
            for (int j = 0; j < size3; j++) {
                subarray3[j] = X[i + j];
            }
            fft_base3(subarray3, size3, N3);
            for (int j = 0; j < size3; j++) {
                X[i + j] = subarray3[j];
            }
        }

        free(subarray3);
    }
}

// Fonction d'affichage du résultat
void print_complex_array(const char *v, complex double *X, int N) {
    for (int i = 0; i < N; i++) {
        printf("%s(%2d) = %15.9f + %15.8fi;\n", v, i+1, creal(X[i]), cimag(X[i]));
    }
}

int main(int argc, char *argv[]) {
    int const N2 = (argc > 1) ? atoi(argv[1]) : 1;
    int const N3 = (argc > 2) ? atoi(argv[2]) : 0;
    int const N  = ipow(2,N2) * ipow(3,N3);

    printf("fprintf('N = %d\\n');\n", N);
    
    // Exemple de signal de taille N
    complex double x[N];
    for (int i = 0; i < N; i++) {
        x[i] = sin(1+i*0.11) + I*sin(2-i*0.22);
    }

    printf("%% Entrée:\n");
    print_complex_array("X", x, N);

    fft_hybrid(x, N, N2, N3);

    printf("\n%% Sortie (FFT):\n");
    print_complex_array("Y", x, N);

    printf("\n");
    printf("F = fft(X);\n");
    printf("e = max(abs(Y-F))\n");
    printf("if (e > 1e-3); fprintf('KO\\n');quit(1); else fprintf('OK\\n');quit(0); end\n");

    return 0;
}

