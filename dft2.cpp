// Cooley-Tukey + Rader : DFT de taille 2^P2
//
// 2iPI/2*k => Tous canoniques, aucune multiplication dans ce papillon

#include "dft.h"



// Copier/permuter les elements selon l'ordre de bit inverse : algo explicite bit a bit
static void base2_reverse_copy (cplx_t *source, cplx_t *dest, myuint N, myuint P2)
{
    for (myuint i = 0; i < N; i++)
    {
        myuint rev = 0;
        for (myuint j = 0; j < P2; j++)
        {
            rev <<= 1;
            rev |= (i >> j) & 1;
        }
        assert(/*(rev >= 0) &&*/ (rev < N));
        dest[rev] = source[i];                  C.mio += 2;
    }
}

// Permuter les elements selon l'ordre de bit inverse : algo Gold-Rader
static void base2_reverse_shuffle (cplx_t *X, myuint N)
{
    if (N <= 2) return;
    myuint i = 0, j = 0;
    for (;;)
    {
        if (i < j) { assert(i<N && j<N); const cplx_t tmp = X[i]; X[i] = X[j]; X[j] = tmp;      C.mio += 4; }
        if (++i >= N-1) break;
        myuint k = N/2;
        while (k <= j) { j -= k; k /= 2; }
        j += k;
    }
}

// Algo Rubio-Gomez-Drouiche : ne produit que les echanges necessaires (environ N/2 echanges)
static void base2_reverse_rgd (cplx_t *X, myuint N, myuint P2)
{
    (void) X;
    myuint rev[N];
    rev[0] = 0;
    rev[1] = 1 << (P2-1);
    rev[2] = 1 << (P2-2);
    for (myuint np = (1 << 1) - 1, k = 2; k <= P2; k++)
    {
        myuint nk = (1 << k) - 1;
        rev[nk] = rev[np] + (1 << (P2-k));
        for (myuint j = 1; j <= np; j++) { rev[nk-j] = rev[nk] - rev[j]; }
        np = nk;
    }

    // Check
    for (myuint i = 0; i < N; i++)
    {
        myuint r = 0;
        for (myuint j = 0; j < P2; j++)
        {
            r <<= 1;
            r |= (i >> j) & 1;
        }
        assert(rev[i] == r);
    }
}

// Faire un CHECK au boot du process
static struct ototest { ototest(); } oto;
ototest::ototest()
{
    for (myuint i = 2; i <= 16; i++) base2_reverse_rgd (nullptr, 1<<i, i);
}

#if 0
static void base2_reverse_shuffle_rius (cplx_t *X, myuint N)
{
   myuint k, 1,r, i;
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
myuint N, J[]; complex X[ 1 ;
{ myuint kIl,r,i,j,N2,N4,N21;
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


// FFT in-place pour N = 2^P2
void fft_base2 (cplx_t *X, myuint N, myuint P2)
{
    C.name(__func__);
    constexpr myuint R = 2;
    base2_reverse_shuffle(X, N);
 
    for (myuint s = 1; s <= P2; s++)
    {
        myuint m = 1 << s;
        myuint h = m / R;
        cplx_t wm = zexpiy(-2.0 * PI / m);      // exp(-2iPI / 2^s)

        cplx_t wmj[h];                          // exp(-2iPI * j / 2^s)   j=0..2^(s-1)  == 2^P values
        wmj[0] = cplx_t(1.0f, 0.0f);
        for (myuint j = 1; j < h; j++) wmj[j] = wmj[j-1] * wm;
        C.zmul += h-1; C.mio += h;

        for (myuint k = 0; k < N; k += m)
        {
            for (myuint j = 0; j < h; j++)
            {
                cplx_t t0 = X[j + k + 0*h];
                cplx_t t1 = X[j + k + 1*h] * wmj[j];
                X[j + k + 0*h] = t0 + t1;
                X[j + k + 1*h] = t0 - t1;
                C.zmul += 1; C.zadd += 2; C.mio += 5;
            }
        }
    }
}

// FFT-centered in-place pour N = 2^P2
void fft_base2_centered (cplx_t *X, myuint N, myuint P2)
{
    C.name(__func__);
    constexpr myuint R = 2;
    base2_reverse_shuffle(X, N);
 
    {
        myuint const s = 1;
        myuint const m = 1 << s;    // == 2
        myuint const h = m / R;     // == 1
//      cplx_t wm(-1,0);            // == exp(-2iPI / 2^s)
//      cplx_t wmj[1](1,0);         // exp(-2iPI * j / 2^s)   j=0..2^(s-1)

        for (myuint k = 0; k < N; k += m)
        {
            myuint const j = 0;     // loop j=[0,h[
            {
                cplx_t t0 =  X[j + k + 0*h];
                cplx_t t1 = -X[j + k + 1*h];        // minus on odd-index entries
                X[j + k + 0*h] = t0 + t1;
                X[j + k + 1*h] = t0 - t1;
                C.zmul += 1; C.zadd += 2; C.mio += 5;
            }
        }
    }

    for (myuint s = 2; s <= P2; s++)
    {
        myuint m = 1 << s;
        myuint h = m / R;
        cplx_t wm = zexpiy(-2.0 * PI / m);      // exp(-2iPI / 2^s)

        cplx_t wmj[h];                          // exp(-2iPI * j / 2^s)   j=0..2^(s-1)  == 2^P values
        wmj[0] = cplx_t(1.0f, 0.0f);
        for (myuint j = 1; j < h; j++) wmj[j] = wmj[j-1] * wm;
        C.zmul += h-1; C.mio += h;

        for (myuint k = 0; k < N; k += m)
        {
            for (myuint j = 0; j < h; j++)
            {
                cplx_t t0 = X[j + k + 0*h];
                cplx_t t1 = X[j + k + 1*h] * wmj[j];
                X[j + k + 0*h] = t0 + t1;
                X[j + k + 1*h] = t0 - t1;
                C.zmul += 1; C.zadd += 2; C.mio += 5;
            }
        }
    }
}


// Stockham FFT
// Cf http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham2.html
//
// Reputation de "moins amical" avec le cache sur Cooley-Tukey
// + les R/W de t0 et T1 a chaque boucle interne soint lointaines l'une de l'autre
//   ... mais c'est idem pour Cooley-Tukey ? Qui par contre ecrit la ou il a lu,
//       au contraire de Stockham

static void stockham (myuint n, myuint s, bool eo, cplx_t *x, cplx_t *y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    myuint const N = n;     // save for assertion
    while (n >= 4)
    {
        myuint const m = n / 2;
        double const theta0 = -2*PI/n;
        for (myuint p = 0; p < m; p++)
        {
            cplx_t const wp = zexpiy(p*theta0); // Tabule
            assert( (s == 1) || ((s%2) == 0) ); // Par q+=2 pour favoriser SIMD
            for (myuint q = 0; q < s; q++)
            {
                cplx_t const t0 = x[q + s*(p + 0)];     // bi-sequentiel en read
                cplx_t const t1 = x[q + s*(p + m)];
                y[q + s*(2*p + 0)] =  t0 + t1;          // bi-sequentiel en write
                y[q + s*(2*p + 1)] = (t0 - t1) * wp;
                C.zmul += 1; C.zadd += 2; C.mio += 4;
            }
        }
        n = m;
        s = 2*s;
        eo = !eo;
        cplx_t *z = y; y = x; x = z;
    }

    assert((n == 2) && (s == N/2)); // Dernier arrangement, sur les deux demi-tableaux initiaux
    {
        cplx_t * const z = eo ? y : x;
        for (myuint q = 0; q < s; q++)
        {
            const cplx_t t0 = x[q + 0];
            const cplx_t t1 = x[q + s];
            z[q + 0] = t0 + t1;
            z[q + s] = t0 - t1;
            C.zmul += 0; C.zadd += 2; C.mio += 4;
        }
    }
}

void fft_stockham_base2 (cplx_t *X, myuint N, myuint P2)
{
    C.name(__func__);
    (void) P2;
    cplx_t Y[N];
    stockham(N, 1, false, X, Y);
}
