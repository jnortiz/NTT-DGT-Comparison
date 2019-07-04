#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "gaussian_integer.h"
#include "params.h"
#include "poly.h"

#define RADIX32 32

void set_gauss(gauss_t *pointer, const int64_t re, const int64_t img) {
	pointer->re = re;
	pointer->img = img;
}

void print_gauss_t(const gauss_t x) {
	printf("%" PRIu64, x.re);
	if(x.img >= 0) {
		printf(" + ");
	}
	printf("i%" PRIu64 ", ", x.img);
}

int64_t sub_int64_t(const int64_t a, const int64_t b) {
	return barr_reduce(a-b);
}

int64_t add_int64_t(const int64_t a, const int64_t b) {
	int64_t res;
	res = a+b;
	res -= PARAM_Q;
	res += (res >> (RADIX32-1)) & PARAM_Q;
	return res;
}

int64_t mul_int64_t(const int64_t a, const int64_t b) {
	// unsigned __int128 x, y;
	// int64_t res;	
	// x = (unsigned __int128) a;
	// y = (unsigned __int128) b;
	// res = reduce(a*b);	
	return reduce(a*b);
}

int overflow(const int64_t a, const int64_t b) {
	return((a+b) < a);
}

void add(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, add_int64_t(x.re,y.re), add_int64_t(x.img,y.img));
}

void sub(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, sub_int64_t(x.re,y.re), sub_int64_t(x.img,y.img));
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {	
	int64_t s1, s2, s3;

	s1 = mul_int64_t(x.re, y.re);
	s2 = mul_int64_t(x.img, y.img);	
	s3 = mul_int64_t(add_int64_t(x.re, x.img), add_int64_t(y.re, y.img));

	set_gauss(z, sub_int64_t(s1,s2), sub_int64_t(sub_int64_t(s3,s1),s2));
}