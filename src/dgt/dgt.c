#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "gauss.h"
#include "const.h"

int64_t get_invk_modp(const int dim) {

    int64_t inv;

    switch(dim) {
        case 256: inv = -72057594021150720; break;
        case 512: inv = -36028797010575360; break;
        case 1024: inv = -18014398505287680; break;
        case 2048: inv = -9007199252643840; break;
        default: inv = 0; break;
    }

    return inv;

}

/* dgt() receives as input a folded signal of n/2 Gaussian integers */
void dgt(gauss_t* _x, const gauss_t* _input_signal) {
    
    uint64_t power;
    int i, j, k, l, m, n, stride, upperbound;
    gauss_t xi, xim, aux_mul, aux_power;

    k = n/2;

    for(i = 0; i < k; i++) {
        _x[i].re = _input_signal[i].re;
        _x[i].img = _input_signal[i].img;
    }

    upperbound = (int)(log(k)/log(2));
    
    for(stride = 0; stride < upperbound; stride++) {
        
        m = k/(2<<stride);
        
        for(l = 0; l < k/2; l++) {
            
            j = l/(k/(2*m));
            power = mod_uint64_t((uint64_t)(pow(gj[j], k >> (upperbound-stride))));
            i = j+(l%(k/(2*m)))*2*m;

            xi.re = _x[i].re;
            xi.img = _x[i].img;

            xim.re = _x[i+m].re;
            xim.img = _x[i+m].img;

            Gauss(&aux_power, (uint64_t)power, (uint64_t) 0);

            add(&_x[i], xi, xim);
            mul(&aux_mul, xi, xim);
            mul(&_x[i+m], aux_power, aux_mul);

        }

    }
}

void idgt(gauss_t *_output_signal, const gauss_t *_x) {

    uint64_t power;
    int i, j, k, l, m, n, stride, upperbound;
    int64_t inv;
    gauss_t xi, xim, aux_inv, aux_mul, aux_power;

    k = (int)n/2;
    upperbound = (int)(log(k)/log(2));

    gauss_t _copy[n/2];

    for(i = 0; i < k; i++) {
        _copy[i].re = _x[i].re;
        _copy[i].img = _x[i].img;
    }

    m = 1;

    for(stride = 0; stride < upperbound; stride++) {

        for(l = 0; l < (int)(k/2); l++) {

            j = (int)(l/(k/(2*m)));
            power = mod_uint64_t((uint64_t)(pow(invgj[j], k >> (stride+1))));
            i = j + (l % ((int)(k/(2*m))))*2*m;

            xi = _copy[i];
            xim = _copy[i+m];

            Gauss(&aux_power, (uint64_t)power, (uint64_t) 0);
            mul(&aux_mul, aux_power, xim);
            add(&_copy[i], xi, aux_mul);
            sub(&_copy[i+m], xi, aux_mul);

        }

        m = 2*m;

    }

    inv = get_invk_modp(k);

    aux_inv.re = inv;
    aux_inv.img = 0;

    for(i = 0; i < k; i++) {
        mul(&_output_signal[i], _copy[i], aux_inv);
    }

}

int main(int argc, char** argv[]) {
    srand(time(NULL)); /* For generating random 128-bit vectors */
    return 0;
}