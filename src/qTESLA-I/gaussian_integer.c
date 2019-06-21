#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "gaussian_integer.h"
#include "params.h"
#include "poly.h"

#define RADIX32 32

void set_gauss(gauss_t *pointer, const int32_t re, const int32_t img) {
	pointer->re = re;
	pointer->img = img;
}

void print_gauss_t(const gauss_t x) {
	printf("%" PRId32, x.re);
	if(x.img >= 0) {
		printf(" + ");
	}
	printf("i%" PRId32 ", ", x.img);
}

int32_t sub_int32_t(const int32_t a, const int32_t b) {
	int32_t res;	
	res = (a-b);
    res += (res >> (RADIX32-1)) & PARAM_Q;    // If result[i] < 0 then add q
	return res;
}

int32_t add_int32_t(const int32_t a, const int32_t b) {
	int32_t res;
	res = a+b;
	res += (res >> (RADIX32-1)) & PARAM_Q;
	res -= PARAM_Q;
	res += (res >> (RADIX32-1)) & PARAM_Q;
	return res;
}

int32_t mul_int32_t(const int32_t a, const int32_t b) {
	int64_t x, y, res;
	x = (int64_t)a;
	y = (int64_t)b;
	res = x*y;
	return(reduce(res));
}

int overflow(const int32_t a, const int32_t b) {
	return((a+b) < a);
}

void add(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, add_int32_t(x.re,y.re), add_int32_t(x.img,y.img));
}

void sub(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, sub_int32_t(x.re,y.re), sub_int32_t(x.img,y.img));
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {	
	int32_t s1, s2, s3;

	s1 = mul_int32_t(x.re, y.re);
	s2 = mul_int32_t(x.img, y.img);	
	s3 = mul_int32_t(add_int32_t(x.re, x.img), add_int32_t(y.re, y.img));

	set_gauss(z, sub_int32_t(s1,s2), sub_int32_t(sub_int32_t(s3,s1),s2));
}