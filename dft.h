// Explorations/pesages de DFT
//

#include <math.h>
#include <complex.h>
#include <assert.h>

// Ceci change avec chaque version du std C++ : M_PI is POSIX not std, C++20 <numbers>,
// ... what else next ?
static constexpr double PI = 3.141592653589793238462;

// size_t is possibly too wide  (64 bits)
typedef unsigned int myuint;

// Sucre
typedef std::complex<double> cplx_t;

// Pour les Twiddles
static inline constexpr cplx_t zexpiy (double y)
{
    return cplx_t(std::cos(y), std::sin(y));
}


// Static integer pow(base,exp)
template <myuint base, myuint exp>
struct kipow
{
    static_assert((base >= 1) /*&& (exp >= 0)*/, "positive");
    static constexpr myuint result = base * kipow<base, exp-1>::result;
};
// Bottom of recursion
template <myuint base>
struct kipow<base, 0>
{
    static_assert((base >= 1), "positive");
    static constexpr myuint result = 1;
};

static inline myuint ipow (myuint b, myuint n)         // pow(b,n)
{
    assert((b >= 2) /*&& (n >= 0)*/);
    myuint r = 1;
    while (n > 0) if (n%2 == 0) { b *= b; n /= 2; } else { r *= b; n--; }
    return r;
}



// Pesages
struct Counters 
{                           // Comptes :
    unsigned mio;           // + des I/O mem (1 = un scalaire, 32 ou 64 bits)
    unsigned zmul;          // + des multiplications complexes
    unsigned zadd;          // + des additions complexes
    const char *fn;         // Nom de la FFT pesee
    void raz ()                 { mio = 0; zmul = 0; zadd = 0; fn = nullptr; }
    void name (const char *n)   { if (fn == nullptr) fn = n; }
};
extern Counters C;



// Fonction pour permuter les elements selon l'ordre inverse des digits base-R
template <myuint R>                                        // N = R^P
void baseR_reverse_copy (cplx_t *source, cplx_t *dest, myuint N, myuint P)
{
    assert(ipow(R,P) == N);
    for (myuint i = 0; i < N; i++)
    {
        myuint rev = 0;
        myuint tmp = i;
        for (myuint j = 0; j < P; j++)
        {
            rev = rev * R + tmp % R;
            tmp /= R;
        }
        assert(/*(rev >= 0) &&*/ (rev < N));
        dest[rev] = source[i];                  C.mio += 2;
    }
}


// FFT de taille N = 2^P
extern void fft_base2 (cplx_t *X, myuint N, myuint P2);
extern void fft_stockham_base2 (cplx_t *X, myuint N, myuint P2);


// FFT de taille N = 3^P
extern void fft_base3 (cplx_t *X, myuint N, myuint P8);

// FFT de taille N = 4^P
extern void fft_base4 (cplx_t *X, myuint N, myuint P8);

// FFT de taille N = 8^P
extern void fft_base8 (cplx_t *X, myuint N, myuint P8);
