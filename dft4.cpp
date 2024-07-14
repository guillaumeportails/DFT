// Cooley-Tukey + Rader : DFT de taille 4^P
//
// 2iPI/4*k => Tous canoniques, aucune multiplication dans ce papillon

#include "dft.h"


static constexpr myuint R = 4;

// Le compilo saurait-il optimiser l'usage de ces valeurs canoniques ?
// static constexpr cplx_t tw_1_1 = cplx_t(0, -1);   // zexpiy(-2*PI/R * 1 * 1);
// static constexpr cplx_t tw_2_1 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 2 * 1);    // == -1
// static constexpr cplx_t tw_3_1 = cplx_t(0,  1);   // zexpiy(-2*PI/R * 3 * 1);
// static constexpr cplx_t tw_1_2 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 1 * 2);    // == -1
// static constexpr cplx_t tw_2_2 = cplx_t( 1, 0);   // zexpiy(-2*PI/R * 2 * 2);    // == +1
// static constexpr cplx_t tw_3_2 = tw_2_1;          // zexpiy(-2*PI/R * 3 * 2);    // (3*2)%R = (2*1)
// static constexpr cplx_t tw_1_3 = tw_3_1;          // zexpiy(-2*PI/R * 1 * 3);
// static constexpr cplx_t tw_2_3 = tw_3_2;          // zexpiy(-2*PI/R * 2 * 3);
// static constexpr cplx_t tw_3_3 = tw_1_1;          // zexpiy(-2*PI/R * 3 * 3);    // (3*3)%R = (1*1)

// t *  exp(2*i*PI/R   *a*b)
static inline cplx_t mw_1_1 (const cplx_t & t)  { return cplx_t( std::imag(t), -std::real(t)); }
static inline cplx_t mw_2_1 (const cplx_t & t)  { return -t; }
static inline cplx_t mw_3_1 (const cplx_t & t)  { return cplx_t(-std::imag(t),  std::real(t)); }
static inline cplx_t mw_1_2 (const cplx_t & t)  { return mw_2_1(t); }
static inline cplx_t mw_2_2 (const cplx_t & t)  { return t; }
static inline cplx_t mw_3_2 (const cplx_t & t)  { return mw_2_1(t); }   // 3*2 == 2*1  % R
static inline cplx_t mw_1_3 (const cplx_t & t)  { return mw_3_1(t); }
static inline cplx_t mw_2_3 (const cplx_t & t)  { return mw_3_2(t); }
static inline cplx_t mw_3_3 (const cplx_t & t)  { return mw_1_1(t); }   // 3*3 == 1*1  % R



// FFT pour N = 4^P4
void fft_base4 (cplx_t *X, myuint N, myuint P4)
{
    C.name(__func__);
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P4);

    for (myuint s = 1; s <= P4; s++)
    {
        myuint h = ipow(R, s-1);
        myuint m = h * R;
        cplx_t wm = zexpiy(-2.0 * PI / m);         // exp(-2iPI/2^(2s))     Read from table

        // Exchange zmul for memIO 
        cplx_t wmj[h];                            // exp(-2iPI * j / 2^2s)   j=0..2^(2s-2)  == 4^P values
        wmj[0] = cplx_t(1.0f, 0.0f);
        for (myuint j = 1; j < h; j++) wmj[j] = wmj[j-1] * wm;
        C.zmul += h-1; C.mio += h;  // Si la relecture de wmj[j-1] est epargnee

        for (myuint k = 0; k < N; k += m)
        {
            for (myuint j = 0; j < h; j++)
            {                                      // w^i tabule
                cplx_t const w = wmj[j];
                cplx_t const t0 = temp[k + j + 0*h];
                cplx_t const t1 = temp[k + j + 1*h] * w;
                cplx_t const t2 = temp[k + j + 2*h] * w*w;      // wmj[ (2*j)%m + symetrie PI ]
                cplx_t const t3 = temp[k + j + 3*h] * w*w*w;    // wmj[ (3*j)%m ]
                temp[k + j + 0*h] = t0 + t1         + t2         + t3;
                temp[k + j + 1*h] = t0 + mw_1_1(t1) + mw_2_1(t2) + mw_3_1(t3);
                temp[k + j + 2*h] = t0 + mw_1_2(t1) + mw_2_2(t2) + mw_3_2(t3);
                temp[k + j + 3*h] = t0 + mw_1_3(t1) + mw_2_3(t2) + mw_3_3(t3);
                C.zmul += 6; C.zadd += 4*3; C.mio += 4*2+1;
            }
        }
    }

    for (myuint i = 0; i < N; i++) X[i] = temp[i];
    C.mio += 2*N;
}

