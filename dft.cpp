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

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <assert.h>

#define PI 3.141592653589793238462

// Sucre
typedef std::complex<double> cplx_t;
// Le compilo sait-il que cette fonction est pure ?
static inline cplx_t zexpiy (double y)     { return cplx_t(cos(y), sin(y)); }

// Fonction d'affichage de tableau, syntaxe Octave
void print_complex_array(const char *v, cplx_t *X, int N)
{
    printf("%s = zeros(%d,1);\n", v, N);
    for (int i = 0; i < N; i++) {
        printf("%s(%2d) = %18.15f + %18.15fi;\n", v, i+1, std::real(X[i]), std::imag(X[i]));
    }
}


static int ipow (int b, int n)         // pow(b,n)
{
    assert((b >= 2) && (n >= 0));
    int r = 1;
    while (n > 0) if (n%2 == 0) { b *= b; n /= 2; } else { r *= b; n--; }
    return r;
}

// Copier/permuter les elements selon l'ordre de bit inverse : algo explicite bit a bit
static void base2_reverse_copy (cplx_t *source, cplx_t *dest, int N, int P2)
{
    for (int i = 0; i < N; i++)
    {
        int rev = 0;
        for (int j = 0; j < P2; j++)
        {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        assert((rev >= 0) && (rev < N));
        dest[rev] = source[i];
    }
}

// Permuter les elements selon l'ordre de bit inverse : algo Gold-Rader
static void base2_reverse_shuffle (cplx_t *X, int N)
{
    if (N <= 2) return;
    int i = 0, j = 0;
    for (;;)
    {
        if (i < j) { assert(i<N && j<N); const cplx_t tmp = X[i]; X[i] = X[j]; X[j] = tmp; }
        if (++i >= N-1) break;
        int k = N/2;
        while (k <= j) { j -= k; k /= 2; }
        j += k;
    }
}

// Algo Rubio-Gomez-Drouiche
static void base2_reverse_rgd (cplx_t *X, int N, int P2)
{
    (void) X;
    int rev[N];
    rev[0] = 0;
    rev[1] = 1 << (P2-1);
    rev[2] = 1 << (P2-2);
    for (int np = (1 << 1) - 1, k = 2; k <= P2; k++)
    {
        int nk = (1 << k) - 1;
        rev[nk] = rev[np] + (1 << (P2-k));
        for (int j = 1; j <= np; j++) rev[nk-j] = rev[nk] - rev[j];
        np = nk;
    }

    // Check
    for (int i = 0; i < N; i++)
    {
        int r = 0;
        for (int j = 0; j < P2; j++)
        {
            r <<= 1;
            r |= (i >> j) & 1;
        }
        assert(rev[i] == r);
    }
}

#if 0
static void base2_reverse_shuffle_rius (cplx_t *X, int N)
{
   int k, 1,r, i;
  J[O]=O;
  J[l]=N/2;
  swap(l,N/2,X);
  for(i=2,1=2,r=N/4; r: 1*=2,r/=2)
  for (k=O; k<l; k++, i++) {
    J [il=J [kl+r;
    if (i<J[il) swap(i, J[i] ,X);
  }
}
static void base2_reverse_copy_rius_yong (N,J, X)
int N, J[]; complex X[ 1 ;
{ int kIl,r,i,j,N2,N4,N21;
N2=N/2; N4=N2/2; N21=N2+1; j = O ;
J[O]=O; swap(l,N2,X);
for (i=2,1=2, r=N4; r>l; 1*=2, r/=2)
for (k=O; k<l; k+=2, i+=2)
{ j=J[i]=J[k]+r;
if(i<j) { swap(i, j,x);
swap(it1, j+N2,X) ;
swap(i+N2l,]+NZl,X);
}
#endif


// Fonction pour permuter les elements selon l'ordre inverse des digits base-R
template <int R>     // N = R^P
void baseR_reverse_copy (cplx_t *source, cplx_t *dest, int N, int P)
{
    assert(ipow(R,P) == N);
    for (int i = 0; i < N; i++)
    {
        int rev = 0;
        int tmp = i;
        for (int j = 0; j < P; j++)
        {
            rev = rev * R + tmp % R;
            tmp /= R;
        }
        assert((rev >= 0) && (rev < N));
        dest[rev] = source[i];
    }
}


// FFT in-place pour N = 2^P2
static void fft_base2 (cplx_t *X, int N, int P2)
{
    constexpr int R = 2;
    base2_reverse_shuffle(X, N);
 
    for (int s = 1; s <= P2; s++)
    {
        int m = 1 << s;
        int h = m / R;
        cplx_t wm = zexpiy(-2.0 * PI / m);
        for (int k = 0; k < N; k += m)
        {
            cplx_t w = 1.0;
            for (int j = 0; j < h; j++)
            {
                cplx_t t0 = X[k + j + 0*h];
                cplx_t t1 = X[k + j + 1*h] * w;
                X[k + j + 0*h] = t0 + t1; // t1*zexpiy(-2*PI/2 * 0*1);
                X[k + j + 1*h] = t0 - t1; // t1*zexpiy(-2*PI/2 * 1*1);
                w *= wm;
            }
        }
    }
}

// FFT pour N = 3^P3
static void fft_base3 (cplx_t *X, int N, int P3)
{
    constexpr int R = 3;
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P3);

    const cplx_t  tw_1_1 =             zexpiy(-2*PI/R * 1 * 1);
    const cplx_t  tw_2_1 =             zexpiy(-2*PI/R * 2 * 1);
    const cplx_t& tw_1_2 = tw_2_1;  // zexpiy(-2*PI/R * 1 * 2);
    const cplx_t& tw_2_2 = tw_1_1;  // zexpiy(-2*PI/R * 1 * 2);  // (2*2) % R = (1*1)

    for (int s = 1; s <= P3; s++) {
        int h = ipow(R, s-1);
        int m = h * R;
        cplx_t wm = zexpiy(-2.0 * PI / m);
        for (int k = 0; k < N; k += m) {
            cplx_t w = 1.0;
            for (int j = 0; j < h; j++) {
                cplx_t const t0 = temp[k + j + 0*h];
                cplx_t const t1 = temp[k + j + 1*h] * w;
                cplx_t const t2 = temp[k + j + 2*h] * w*w;
                temp[k + j + 0*h] = t0 + t1 + t2;
                temp[k + j + 1*h] = t0 + t1*tw_1_1 + t2*tw_2_1;
                temp[k + j + 2*h] = t0 + t1*tw_1_2 + t2*tw_2_2;
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}


// FFT pour N = 4^P4
static void fft_base4 (cplx_t *X, int N, int P4)
{
    constexpr int R = 4;
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P4);

    // Le compilo saura t-il optimiser l'usage de ces valeurs canoniques ?
    const cplx_t  tw_1_1 =                     zexpiy(-2*PI/R * 1 * 1);
    const cplx_t  tw_2_1 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 2 * 1);    // == -1
    const cplx_t  tw_3_1 =                     zexpiy(-2*PI/R * 3 * 1);
    const cplx_t  tw_1_2 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 1 * 2);    // == -1
    const cplx_t  tw_2_2 = cplx_t( 1, 0);   // zexpiy(-2*PI/R * 2 * 2);    // == +1
    const cplx_t  tw_3_2 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 3 * 2);    // (3*2)%R = (2*1)

    for (int s = 1; s <= P4; s++) {
        int h = ipow(R, s-1);
        int m = h * R;
        cplx_t wm = zexpiy(-2.0 * PI / m);
        for (int k = 0; k < N; k += m) {
            cplx_t w = 1.0;
            for (int j = 0; j < h; j++) {
                cplx_t const t0 = temp[k + j + 0*h];
                cplx_t const t1 = temp[k + j + 1*h] * w;
                cplx_t const t2 = temp[k + j + 2*h] * w*w;
                cplx_t const t3 = temp[k + j + 3*h] * w*w*w;
                temp[k + j + 0*h] = t0 + t1 + t2 + t3;
                temp[k + j + 1*h] = t0 + t1*tw_1_1 + t2*tw_2_1 + t3*tw_3_1;
                temp[k + j + 2*h] = t0 + t1*tw_1_2 + t2*tw_2_2 + t3*tw_3_2;
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}


// FFT pour N = 8^P8
static void fft_base8 (cplx_t *X, int N, int P8)
{
    constexpr int R = 8;
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P8);

    for (int s = 1; s <= P8; s++)
    {
        int h = ipow(R, s-1);
        int m = h * R;
        cplx_t wm = zexpiy(-2.0 * PI / m);
        for (int k = 0; k < N; k += m)
        {
            cplx_t w = 1.0;
            for (int j = 0; j < h; j++)
            {
                cplx_t const t0 = temp[k + j + 0 * h];
                cplx_t const t1 = temp[k + j + 1 * h] * w;
                cplx_t const t2 = temp[k + j + 2 * h] * w*w;
                cplx_t const t3 = temp[k + j + 3 * h] * w*w*w;
                cplx_t const t4 = temp[k + j + 4 * h] * w*w*w*w;
                cplx_t const t5 = temp[k + j + 5 * h] * w*w*w*w*w;
                cplx_t const t6 = temp[k + j + 6 * h] * w*w*w*w*w*w;
                cplx_t const t7 = temp[k + j + 7 * h] * w*w*w*w*w*w*w;
                // Les zexpiy sont largement simplifiables ...
                temp[k + j + 0 * h] =   t0 * zexpiy(-2*PI/R * 0 * 0)
                                      + t1 * zexpiy(-2*PI/R * 0 * 1)
                                      + t2 * zexpiy(-2*PI/R * 0 * 2)
                                      + t3 * zexpiy(-2*PI/R * 0 * 3)
                                      + t4 * zexpiy(-2*PI/R * 0 * 4)
                                      + t5 * zexpiy(-2*PI/R * 0 * 5)
                                      + t6 * zexpiy(-2*PI/R * 0 * 6)
                                      + t7 * zexpiy(-2*PI/R * 0 * 7);
                temp[k + j + 1 * h] =   t0 * zexpiy(-2*PI/R * 1 * 0)
                                      + t1 * zexpiy(-2*PI/R * 1 * 1)
                                      + t2 * zexpiy(-2*PI/R * 1 * 2)
                                      + t3 * zexpiy(-2*PI/R * 1 * 3)
                                      + t4 * zexpiy(-2*PI/R * 1 * 4)
                                      + t5 * zexpiy(-2*PI/R * 1 * 5)
                                      + t6 * zexpiy(-2*PI/R * 1 * 6)
                                      + t7 * zexpiy(-2*PI/R * 1 * 7);
                temp[k + j + 2 * h] =   t0 * zexpiy(-2*PI/R * 2 * 0)
                                      + t1 * zexpiy(-2*PI/R * 2 * 1)
                                      + t2 * zexpiy(-2*PI/R * 2 * 2)
                                      + t3 * zexpiy(-2*PI/R * 2 * 3)
                                      + t4 * zexpiy(-2*PI/R * 2 * 4)
                                      + t5 * zexpiy(-2*PI/R * 2 * 5)
                                      + t6 * zexpiy(-2*PI/R * 2 * 6)
                                      + t7 * zexpiy(-2*PI/R * 2 * 7);
                temp[k + j + 3 * h] =   t0 * zexpiy(-2*PI/R * 3 * 0)
                                      + t1 * zexpiy(-2*PI/R * 3 * 1)
                                      + t2 * zexpiy(-2*PI/R * 3 * 2)
                                      + t3 * zexpiy(-2*PI/R * 3 * 3)
                                      + t4 * zexpiy(-2*PI/R * 3 * 4)
                                      + t5 * zexpiy(-2*PI/R * 3 * 5)
                                      + t6 * zexpiy(-2*PI/R * 3 * 6)
                                      + t7 * zexpiy(-2*PI/R * 3 * 7);
                temp[k + j + 4 * h] =   t0 * zexpiy(-2*PI/R * 4 * 0)
                                      + t1 * zexpiy(-2*PI/R * 4 * 1)
                                      + t2 * zexpiy(-2*PI/R * 4 * 2)
                                      + t3 * zexpiy(-2*PI/R * 4 * 3)
                                      + t4 * zexpiy(-2*PI/R * 4 * 4)
                                      + t5 * zexpiy(-2*PI/R * 4 * 5)
                                      + t6 * zexpiy(-2*PI/R * 4 * 6)
                                      + t7 * zexpiy(-2*PI/R * 4 * 7);
                temp[k + j + 5 * h] =   t0 * zexpiy(-2*PI/R * 5 * 0)
                                      + t1 * zexpiy(-2*PI/R * 5 * 1)
                                      + t2 * zexpiy(-2*PI/R * 5 * 2)
                                      + t3 * zexpiy(-2*PI/R * 5 * 3)
                                      + t4 * zexpiy(-2*PI/R * 5 * 4)
                                      + t5 * zexpiy(-2*PI/R * 5 * 5)
                                      + t6 * zexpiy(-2*PI/R * 5 * 6)
                                      + t7 * zexpiy(-2*PI/R * 5 * 7);
                temp[k + j + 6 * h] =   t0 * zexpiy(-2*PI/R * 6 * 0)
                                      + t1 * zexpiy(-2*PI/R * 6 * 1)
                                      + t2 * zexpiy(-2*PI/R * 6 * 2)
                                      + t3 * zexpiy(-2*PI/R * 6 * 3)
                                      + t4 * zexpiy(-2*PI/R * 6 * 4)
                                      + t5 * zexpiy(-2*PI/R * 6 * 5)
                                      + t6 * zexpiy(-2*PI/R * 6 * 6)
                                      + t7 * zexpiy(-2*PI/R * 6 * 7);
                temp[k + j + 7 * h] =   t0 * zexpiy(-2*PI/R * 7 * 0)
                                      + t1 * zexpiy(-2*PI/R * 7 * 1)
                                      + t2 * zexpiy(-2*PI/R * 7 * 2)
                                      + t3 * zexpiy(-2*PI/R * 7 * 3)
                                      + t4 * zexpiy(-2*PI/R * 7 * 4)
                                      + t5 * zexpiy(-2*PI/R * 7 * 5)
                                      + t6 * zexpiy(-2*PI/R * 7 * 6)
                                      + t7 * zexpiy(-2*PI/R * 7 * 7);
                w *= wm;
            }
        }
    }

    for (int i = 0; i < N; i++) X[i] = temp[i];
}



// FFT hybride pour N = 2^P2 * 3^P3
void fft_hybrid(cplx_t *X, int N, int P2, int P3)
{
    assert((P2 >= 1) && (P3 >= 1));

    cplx_t X2[N];
    cplx_t X3[N];    
    for(int i = 0; i < N; i++) X2[i] = X[i];
    for(int i = 0; i < N; i++) X3[i] = X[i];
    
    int const size2 = ipow(2, P2);
    for (int i = 0; i < N; i += size2)
    {
        cplx_t subarray2[size2];
        for (int j = 0; j < size2; j++) subarray2[j] = X2[i + j];
        fft_base2 (subarray2, size2, P2);
        for (int j = 0; j < size2; j++) X2[i + j] = subarray2[j];
    }

    int const size3 = ipow(3, P3);
    for (int i = 0; i < N; i += size3)
    {
        cplx_t subarray3[size3];
        for (int j = 0; j < size3; j++) subarray3[j] = X3[i + j];
        fft_base3 (subarray3, size3, P3);
        for (int j = 0; j < size3; j++) X3[i + j] = subarray3[j];
    }

    // for (int k = 0; k < N; k++)
    // {
    //     X[k] = 0;
    //     for (int p = 0; p < 2; p++) for (int q = 0; q < 3; q++)
    //     {
    //         X[k] += twiddle_factors[(p * q * k) % N] * subfft_p[q * 2 + p] * subfft_q[p * 3 + q];
    //     }
    // }


}

// Stockham FFT
// Cf http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham2.html
static void stockham (int n, int s, bool eo, cplx_t *x, cplx_t *y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    while (n >= 4)
    {
        const int m = n/2;
        const double theta0 = -2*PI/n;
        for (int p = 0; p < m; p++) {
            const cplx_t wp = zexpiy(p*theta0);
            assert( (s == 1) || ((s%2) == 0) ); // Par q+=2 pour favoriser SIMD
            for (int q = 0; q < s; q++) {
                const cplx_t  a = x[q + s*(p + 0)]; // bi-sequentiel en read
                const cplx_t  b = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  a + b;                // bi-sequentiel en write
                y[q + s*(2*p + 1)] = (a - b) * wp;
            }
        }
        n = m;
        s = 2*s;
        eo = !eo;
        cplx_t *z = y; y = x; x = z;
    }

    assert(n == 2);
    {
        cplx_t * z = eo ? y : x;
        for (int q = 0; q < s; q++)
        {
            const cplx_t  a = x[q + 0];
            const cplx_t  b = x[q + s];
            z[q + 0] = a + b;
            z[q + s] = a - b;
        }
    }
}

static void fft_stockham_base2 (cplx_t *X, int N, int P2)
{
    (void) P2;
    cplx_t Y[N];
    stockham(N, 1, false, X, Y);
}




static double urand()
{
  double const r = rand();
  double const k = RAND_MAX *0.5;
  return r/k - 1.0;
}

int main(int argc, char *argv[])
{
    int P2 = 0;
    int P3 = 0;
    int P8 = 0;

    int o;
    void (*opt_fft_base2)(cplx_t *x, int N, int P2) = fft_base2;
    bool opt_false = false;
    while ((o = getopt(argc,argv,"2:3:8:Sn")) > 0) switch (o)
    {
        case '2' : P2 = atoi(optarg); break;
        case '3' : P3 = atoi(optarg); break;
        case '8' : P8 = atoi(optarg); break;
        case 'S' : opt_fft_base2 = fft_stockham_base2; break;
        case 'n' : opt_false = true; break; // Option jamais donnee. Battre l'optimiseur.
    }
    assert((P2 > 0) || (P3 > 0) || (P8 > 0));
    assert((P2 == 0) || (P8 == 0));

    int const N  = ipow(2,P2) * ipow(3,P3) * ipow(8,P8);
    printf("fprintf('N = %d\\n');\n", N);
    printf("isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;\n");

    // sanity check
    assert(ipow(7,0) ==  1);
    assert(ipow(5,1) ==  5);
    assert(ipow(3,2) ==  9);
    assert(ipow(2,3) ==  8);
    assert(ipow(2,4) == 16);
//  assert(ipow(1,6) ==  1);
    for (int i = 2; i <= 16; i++) base2_reverse_rgd (nullptr, 1<<i, i);

    // Exemple de signal de taille N
    cplx_t x[N];
    for (int i = 0; i < N; i++) x[i] = cplx_t(urand(), urand());

    printf("%% Entree:\n");
    print_complex_array("X", x, N);

    clock_t tic = clock();
    if ((P3 == 0) && (P8 == 0))
    {
        if ((P2 % 2 == 0) && opt_false)
            fft_base4(x, N, P2/2);
        else
            opt_fft_base2(x, N, P2);
    }
    else if ((P2 == 0) && (P8 == 0))
      fft_base3(x, N, P3);
    else if ((P3 == 0) && (P2 == 0))
      fft_base8(x, N, P8);
//  else
//    fft_hybrid(x, N, P2, P3);
    clock_t toc = clock() - tic;

    printf("\n%% Sortie (FFT):  %.6fs clock\n", (double) toc/CLOCKS_PER_SEC);
    print_complex_array("Y", x, N);

    printf("\n");
    if (N > 10000) printf("\nclk_ms = %.3f %% ms\n", (double) 1e3 * toc/CLOCKS_PER_SEC);
    printf("F = fft(X);\n");
    printf("err = max(abs(Y-F));\n");
    printf("if (err > 1e-6); fprintf('KO %%f\\n',err); error('KO'); else fprintf('OK %%f\\n',err); end\n");

    return 0;
}

