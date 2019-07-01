#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "gaussian_integer.h"
#include "params.h"

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
    struct timespec now;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
    return now.tv_nsec + (timestamp_t)now.tv_sec * 1000000000.0;
}

#define MAX_RUN 5000

void print_signal(const gauss_t *_signal_a, const int length) {    
    int i;

    printf("{(%" PRId32 ", i%" PRId32 "), ", _signal_a[0].re, _signal_a[0].img); 
    for(i = 1; i < length-1; i++) {
        printf("(%" PRId32 ", i%" PRId32 "), ", _signal_a[i].re, _signal_a[i].img); 
    }
    printf("(%" PRId32 ", i%" PRId32 ")}\n", _signal_a[length-1].re, _signal_a[length-1].img); 
}

/* dgt() receives as input a folded signal of n/2 Gaussian integers */
void dgt(gauss_t *_x, const gauss_t *_input_signal) {    
    int i, j, l, m, stride;
    gauss_t xi, xim, aux_sub, aux_power;

    for(i = 0; i < PARAM_K; i++) {
        set_gauss(&_x[i], _input_signal[i].re, _input_signal[i].img);
    }
   
    for(stride = 0; stride < PARAM_K_LOG; stride++) {        
        m = PARAM_K/(2 << stride);
        
        for(l = 0; l < PARAM_K/2; l++) {            
            j = (2*m*l)/PARAM_K;
            i = j + (l%(PARAM_K/(2*m)))*2*m;

            set_gauss(&xi, _x[i].re, _x[i].img);
            set_gauss(&xim, _x[i+m].re, _x[i+m].img);
            set_gauss(&aux_power, __gj[j][stride], (int32_t) 0);
            add(&_x[i], xi, xim);
            sub(&aux_sub, xi, xim);
            mul(&_x[i+m], aux_power, aux_sub);
        }
    }
}

void idgt(gauss_t *_output_signal, const gauss_t *_x) {
    int i, j, l, m, stride;
    gauss_t xi, xim, aux_inv, aux_mul, aux_power;

    for(i = 0; i < PARAM_K; i++) {
        set_gauss(&_output_signal[i], _x[i].re, _x[i].img);
    }

    m = 1;
    for(stride = 0; stride < PARAM_K_LOG; stride++) {
        for(l = 0; l < PARAM_K/2; l++) {
            j = (2*m*l)/PARAM_K;
            i = j + (l % (PARAM_K/(2*m)))*2*m;

            set_gauss(&xi, _output_signal[i].re, _output_signal[i].img);
            set_gauss(&xim, _output_signal[i+m].re, _output_signal[i+m].img);
            set_gauss(&aux_power, __invgj[j][stride], (int32_t) 0);
            mul(&aux_mul, aux_power, xim);

            add(&_output_signal[i], xi, aux_mul);
            sub(&_output_signal[i+m], xi, aux_mul);
        }
        m = 2*m;
    }

    set_gauss(&aux_inv, invofkmodp, (int32_t) 0);

    for(i = 0; i < PARAM_K; i++) {
        mul(&_output_signal[i], _output_signal[i], aux_inv);
    }
}

void mul_poly(int32_t *output, const int32_t * _poly_a, const int32_t *_poly_b) {
    gauss_t _folded_a[PARAM_K], _folded_b[PARAM_K];
    gauss_t _dgt_a[PARAM_K], _dgt_b[PARAM_K];
    gauss_t _mul[PARAM_K], _output_gaussian[PARAM_K];
    gauss_t root;
    int i;

    for(i = 0; i < PARAM_K; i++) {
        set_gauss(&root, __nthroots[i][0], __nthroots[i][1]);
        set_gauss(&_folded_a[i], _poly_a[i], _poly_a[PARAM_K+i]);
        mul(&_folded_a[i], _folded_a[i], root);
        set_gauss(&_folded_b[i], _poly_b[i], _poly_b[PARAM_K+i]);
        mul(&_folded_b[i], _folded_b[i], root);
    }

    dgt(_dgt_a, _folded_a);
    dgt(_dgt_b, _folded_b);

    for(i = 0; i < PARAM_K; i++) {
        mul(&_mul[i], _dgt_a[i], _dgt_b[i]);
    }

    idgt(_output_gaussian, _mul);

    for(i = 0; i < PARAM_K; i++) {
        set_gauss(&root, __invnthroots[i][0], __invnthroots[i][1]);
        mul(&_output_gaussian[i], _output_gaussian[i], root);
        output[i] = _output_gaussian[i].re;
        output[i+PARAM_K] = _output_gaussian[i].img;
    }
}

void textbook_mul(int32_t *_output, const int32_t *_a, const int32_t *_b) {
    int i, j, power;
    int32_t aux, v;
    aux = 0; v = 0;

    for(i = 0; i < PARAM_N; i++) {
        _output[i] = (int32_t)0;
    }

    for(i = 0; i < PARAM_N; i++) {
        for(j = 0; j < PARAM_N; j++) {
            power = pow(-1, (int)((i+j)/PARAM_N));
            v = power;
            if(power < 0) {
                v = (int32_t)((int64_t)power + PARAM_Q);
            }
            aux = ((int64_t)((int64_t)_a[i]*_b[j] % PARAM_Q)*v) % PARAM_Q;
            _output[(i+j)%PARAM_N] = add_int32_t(_output[(i+j)%PARAM_N], aux);
        }
    }
}

void gen_poly(int32_t *_poly) {
    int i;
    
    for(i = 0; i < PARAM_N; i++) {
        _poly[i] = (int32_t) rand() % PARAM_Q;
    }    
}

void print_poly(int32_t *_poly) {
    int i;
    printf("[");
    for(i = 0; i < PARAM_N-1; i++) {
        printf("%" PRId32 ", ", _poly[i]);     
    }
    printf("%" PRId32 "]\n", _poly[PARAM_N-1]);     
}

int is_equal_poly(int32_t *_a, int32_t *_b) {
    int i, is_equal;

    is_equal = 1;

    for(i = 0; i < PARAM_N && is_equal == 1; i++) {
        if(_a[i] != _b[i]) {
            is_equal = 0;
        }
    }

    return is_equal;
}
int is_equal(const gauss_t *_signal_a, const gauss_t *_signal_b) {
    int i, is_equal;

    is_equal = 1;
    for(i = 0; i < PARAM_K && is_equal == 1; i++) {
        if(_signal_a[i].re != _signal_b[i].re || _signal_a[i].img != _signal_b[i].img) {
            is_equal = 0;
        }
    }

    return is_equal;
}

int test_dgt() {
    printf("\nTesting DGT/IDGT transform\n");
    gauss_t _signal_a[PARAM_K], _signal_b[PARAM_K], _output[PARAM_K];
    int i, n_run, counter;
    timestamp_t ts_start, ts_end, ts_min, ts_max, ts_avg, diff;

    counter = 0;
    ts_start = 0.0; ts_end = 0.0;
    ts_min = 10.0; ts_max = 0.0;
    ts_avg = 0.0; diff = 0.0;
    
    srand(time(NULL));

    for(n_run = 0; n_run < MAX_RUN; n_run++) {        
        for(i = 0; i < PARAM_K; i++) {
            set_gauss(&_signal_a[i], rand() % PARAM_Q, 0);
        }

        ts_start = get_timestamp();    
        dgt(_output, _signal_a);
        idgt(_signal_b, _output);
        ts_end = get_timestamp();

        if(is_equal(_signal_a, _signal_b)) counter++;
        
        diff = (ts_end-ts_start);

        ts_avg += diff;

        if(diff < ts_min) ts_min = diff;
        if(diff > ts_max) ts_max = diff;
    }
    
    printf("\tTotal running time: %.10lf s.\n", (double)(ts_avg)/(double)(1000000000.0));
    printf("\tAverage running time: %.10lf s.\n", (double)(ts_avg)/(double)(MAX_RUN*1000000000.0));
    printf("\tMinimum running time: %.10lf s.\n", (double)(ts_min)/(double)(1000000000.0));
    printf("\tMaximum running time: %.10lf s.\n", (double)(ts_max)/(double)(1000000000.0));

    if(counter == MAX_RUN) {
        printf("All %d tests have succeeded.\n.\n", MAX_RUN);
    } else
    printf("%d tests of %d were successful.\n..\n", counter, MAX_RUN);

    printf("--------------------------------------------------\n");
}

int test_poly_mul() {
    printf("\nTesting DGT polynomial multiplication\n");
    int32_t _a[PARAM_N], _b[PARAM_N], _c[PARAM_N], _c_textbook[PARAM_N];
    int n_run, counter;
    timestamp_t ts_start, ts_end, ts_min, ts_max, ts_avg, diff;

    counter = 0;
    ts_start = 0.0; ts_end = 0.0;
    ts_min = 10.0; ts_max = 0.0;
    ts_avg = 0.0; diff = 0.0;

    srand(time(NULL));

    for(n_run = 0; n_run < MAX_RUN; n_run++) {
        
        ts_start = get_timestamp();    
        gen_poly(_a);
        gen_poly(_b);        
        mul_poly(_c, _a, _b);   
        ts_end = get_timestamp();
        
        diff = (ts_end-ts_start);
        ts_avg += diff;

        if(diff < ts_min) ts_min = diff;
        if(diff > ts_max) ts_max = diff;        

        textbook_mul(_c_textbook, _a, _b);
 
        if(is_equal_poly(_c, _c_textbook)) {
            counter++;
        }
    }

    printf("\tTotal running time: %.10lf s.\n", (double)(ts_avg)/(double)(1000000000.0));
    printf("\tAverage running time: %.10lf s.\n", (double)(ts_avg)/(double)(MAX_RUN*1000000000.0));
    printf("\tMinimum running time: %.10lf s.\n", (double)(ts_min)/(double)(1000000000.0));
    printf("\tMaximum running time: %.10lf s.\n", (double)(ts_max)/(double)(1000000000.0));

    if(counter == MAX_RUN) {
        printf("All %d tests have succeeded.\n.\n", MAX_RUN);
    } else
    printf("%d tests of %d were successful.\n..\n", counter, MAX_RUN);

    printf("--------------------------------------------------\n");
}

int main(int argc, char** argv[]) {    
    
    test_dgt();
    test_poly_mul();

    return 0;
}