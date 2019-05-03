#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "dgt.h"
#include "params.h"
#include "const.h"

unsigned less_than(const uint64_t x, const uint64_t y) {
    uint64_t less = (x-y);
    return (less >> sizeof(uint64_t)*8-1);
}

uint64_t select_uint64_t(uint64_t x, uint64_t y, uint64_t bit) {
    uint64_t mask = -bit;
    uint64_t output = mask & (x ^ y);
    return (output ^ x);
}

uint64_t mod_uint64_t(const uint64_t a) {   
    uint64_t output;
    output = select_uint64_t(output, output-p, !less_than(output,p));
    return output;
}

/* dgt() receives as input a folded signal of n/2 Gaussian integers */
void dgt(gauss_t* _x, gauss_t* _input_signal) {
    
    uint64_t power;
    int i, j, k, l, m, n, stride, upperbound;
    gauss_t xi, xim, aux;

    k = n/2;

    for(i = 0; i < k; i++) {
        _x[i].re = _input_signal[i].re;
        _x[i].img = _input_signal[i].img;
    }

    upperbound = (int)(log(k)/log(2));
    
    for(i = 0; i < upperbound; i++) {
        
        m = k/(2<<stride);
        
        for(l = 0; l < k/2; l++) {
            
            j = l/(k/(2*m));
            power = mod_uint64_t((uint64_t)(pow(gj[j], k >> (upperbound-stride))));
            i = j+(l%(k/(2*m)))*2*m;

            xi.re = _x[i].re;
            xi.img = _x[i].img;

            xim.re = _x[i+m].re;
            xim.img = _x[i+m].img;

            add(&_x[i], xi, xim);
            mul(&aux, xi, xim);
            mul(&_x[i+m], power, aux);

        }

    }
}

int main(int argc, char** argv[]) {
    srand(time(NULL)); /* For generating random 128-bit vectors */
    return 0;
}