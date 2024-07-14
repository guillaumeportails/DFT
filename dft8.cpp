// DFT de taille 8^P
//
// Compilation separee pour voir l'ASM sans inliner ce code dans le main
//
// 2iPI/8*k => 4 cas(k impair) non canoniques en sqrt(2) imposant des multiplications

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "dft.h"


static constexpr myuint R = 8;


// Le compilo saura t-il optimiser l'usage de ces valeurs canoniques ?
// Pq pas une lambda ?
template <myuint a, myuint b>
static inline cplx_t multw (const cplx_t & t)
{ 
#if 0
    return t * zexpiy(-2*PI * (a * b) / R);
#else
    // Ces valeurs seront-telle en RAM data ou text, i.e. dans le cache L1 avec les instructions ?
//  static const cplx_t  tw_0 = cplx_t( 1, 0);      zexpiy(-2*PI/R * 0);    // == +1
    static const cplx_t  tw_1 =                     zexpiy(-2*PI/R * 1);
//  static const cplx_t  tw_2 = cplx_t( 0,-1);   // zexpiy(-2*PI/R * 2);    // == -i
    static const cplx_t  tw_3 =                     zexpiy(-2*PI/R * 3);
//  static const cplx_t  tw_4 = cplx_t(-1, 0);   // zexpiy(-2*PI/R * 4);    // == -1
    static const cplx_t  tw_5 =                     zexpiy(-2*PI/R * 5);
//  static const cplx_t  tw_6 = cplx_t(0, 1);   //  zexpiy(-2*PI/R * 6);    // == +i
    static const cplx_t  tw_7 =                     zexpiy(-2*PI/R * 7);
    
    switch ((a*b) % R) {
        case 0:           return t;
        case 1: C.zmul++; return t * tw_1;
        case 2:           return cplx_t(std::imag(t), -std::real(t));
        case 3: C.zmul++; return t * tw_3;
        case 4:           return -t;
        case 5: C.zmul++; return t * tw_5;
        case 6:           return cplx_t(-std::imag(t), std::real(t));;
        case 7: C.zmul++; return t * tw_7;
    };
#endif
}

void fft_base8 (cplx_t *X, myuint N, myuint P8)
{
    C.name(__func__);
    cplx_t temp[N];
    baseR_reverse_copy<R>(X, temp, N, P8);

    for (myuint s = 1; s <= P8; s++)
    {
        myuint h = ipow(R, s-1);                // R^(s-1)
        myuint m = h * R;                       // R^s
        cplx_t wm = zexpiy(-2.0 * PI / m);      // Read table

        // Exchange zmul for memIO 
        cplx_t wmj[m];                          // exp(-2iPI * j / R^s)   j=0..R^s  == R^P values
        wmj[0] = cplx_t(1.0f, 0.0f);
        for (myuint j = 1; j < m; j++) wmj[j] = wmj[j-1] * wm;
        C.zmul += m-1; C.mio += m;  // Si la relecture de wmj[j-1] est epargnee

        for (myuint k = 0; k < N; k += m)
        {
            for (myuint j = 0; j < h; j++)
            {                                       // Tous les w^i sont tabules
                cplx_t const t0 = temp[k + j + 0 * h];
                cplx_t const t1 = temp[k + j + 1 * h] * wmj[ j*1];
                cplx_t const t2 = temp[k + j + 2 * h] * wmj[(j*2)%m];
                cplx_t const t3 = temp[k + j + 3 * h] * wmj[(j*3)%m];
                cplx_t const t4 = temp[k + j + 4 * h] * wmj[(j*4)%m];
                cplx_t const t5 = temp[k + j + 5 * h] * wmj[(j*5)%m];
                cplx_t const t6 = temp[k + j + 6 * h] * wmj[(j*6)%m];
                cplx_t const t7 = temp[k + j + 7 * h] * wmj[(j*7)%m];
                C.mio += 7; C.zmul += 7;
                temp[k + j + 0 * h] =    multw<0,0>(t0)
                                      +  multw<0,1>(t1)
                                      +  multw<0,2>(t2)
                                      +  multw<0,3>(t3)
                                      +  multw<0,4>(t4)
                                      +  multw<0,5>(t5)
                                      +  multw<0,6>(t6)
                                      +  multw<0,7>(t7);
                temp[k + j + 1 * h] =    multw<1,0>(t0)
                                      +  multw<1,1>(t1)
                                      +  multw<1,2>(t2)
                                      +  multw<1,3>(t3)
                                      +  multw<1,4>(t4)
                                      +  multw<1,5>(t5)
                                      +  multw<1,6>(t6)
                                      +  multw<1,7>(t7);
                temp[k + j + 2 * h] =    multw<2,0>(t0)
                                      +  multw<2,1>(t1)
                                      +  multw<2,2>(t2)
                                      +  multw<2,3>(t3)
                                      +  multw<2,4>(t4)
                                      +  multw<2,5>(t5)
                                      +  multw<2,6>(t6)
                                      +  multw<2,7>(t7);
                temp[k + j + 3 * h] =    multw<3,0>(t0)
                                      +  multw<3,1>(t1)
                                      +  multw<3,2>(t2)
                                      +  multw<3,3>(t3)
                                      +  multw<3,4>(t4)
                                      +  multw<3,5>(t5)
                                      +  multw<3,6>(t6)
                                      +  multw<3,7>(t7);
                temp[k + j + 4 * h] =    multw<4,0>(t0)
                                      +  multw<4,1>(t1)
                                      +  multw<4,2>(t2)
                                      +  multw<4,3>(t3)
                                      +  multw<4,4>(t4)
                                      +  multw<4,5>(t5)
                                      +  multw<4,6>(t6)
                                      +  multw<4,7>(t7);
                temp[k + j + 5 * h] =    multw<5,0>(t0)
                                      +  multw<5,1>(t1)
                                      +  multw<5,2>(t2)
                                      +  multw<5,3>(t3)
                                      +  multw<5,4>(t4)
                                      +  multw<5,5>(t5)
                                      +  multw<5,6>(t6)
                                      +  multw<5,7>(t7);
                temp[k + j + 6 * h] =    multw<6,0>(t0)
                                      +  multw<6,1>(t1)
                                      +  multw<6,2>(t2)
                                      +  multw<6,3>(t3)
                                      +  multw<6,4>(t4)
                                      +  multw<6,5>(t5)
                                      +  multw<6,6>(t6)
                                      +  multw<6,7>(t7);
                temp[k + j + 7 * h] =    multw<7,0>(t0)
                                      +  multw<7,1>(t1)
                                      +  multw<7,2>(t2)
                                      +  multw<7,3>(t3)
                                      +  multw<7,4>(t4)
                                      +  multw<7,5>(t5)
                                      +  multw<7,6>(t6)
                                      +  multw<7,7>(t7);
                C.mio += 8*2; C.zadd += 8*7;
            }
        }
    }

    for (myuint i = 0; i < N; i++) X[i] = temp[i];
    C.mio += 2*N;
}
