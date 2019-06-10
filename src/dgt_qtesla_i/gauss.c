#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "gauss.h"

uint64_t p = 0xFFFFFFFF00000001; 
int n = 512;
uint64_t gap = (uint64_t)(4294967295L); // UINT64_MAX-P+1 = 2^{64}-p

void set_gauss(gauss_t *pointer, const uint64_t re, const uint64_t img) {
	pointer->re = re;
	pointer->img = img;
}

void print_gauss_t(const gauss_t x) {
	printf("%" PRIu64, x.re);
	if(x.img >= 0) {
		printf("+");
	}
	printf("%" PRIu64 "i\n", x.img);
}

uint64_t sub_uint64_t(const uint64_t a, const uint64_t b) {
	uint64_t res;
	res = (a-b)+(b>a)*p; // Negative overflow
	return res;
}

uint64_t add_uint64_t(const uint64_t a, const uint64_t b) {
	uint64_t res;
	res = a+b;
	res += (res < a)*4294967295L; // Negative overflow
	res -= (res > p)*p; // (a+b) > p
	return res;
}

uint64_t mul_uint64_t(const uint64_t a, const uint64_t b) {
	unsigned __int128 x, y;
	uint64_t res;
	
	x = (unsigned __int128) a;
	y = (unsigned __int128) b;
	res = mod_uint128(x*y);
	
	return res;
}

int overflow(const uint64_t a, const uint64_t b) {
	return((a+b) < a);
}

uint64_t mod_uint128(const unsigned __int128 a) {
	uint64_t high_half;
	uint64_t low_half;
	
	high_half = (uint64_t) (a >> 64);
	low_half = (uint64_t) (a);
	
	/* It covers the case when (a*b) still fits in 64-bit and it is greater than p. */
	if(high_half == (uint64_t)(0) && low_half > p) {
		return (low_half-p);
	}
	
	uint64_t w0, w1, w2, w3;
	uint64_t W1, W2;
	uint64_t output, gap;
	int negative_overflow, greater_than_p, double_overflow;

	gap = (uint64_t)(4294967295L); // UINT64_MAX-P+1 = 2^{64}-p

	w0 = (uint32_t) a;
	w1 = (uint32_t)(a >> 32);
	w2 = (uint32_t)(a >> 64);
	w3 = (uint32_t)(a >> 96);

	W1 = (w1 << 32);
	W2 = (w2 << 32);

	output = W1+W2-w3-w2+w0;
	
	negative_overflow = (w2+w3 > W1+W2+w0) && !(overflow(W1,W2) || overflow(W1+W2,w0));
	greater_than_p = (output >= p);
	double_overflow = (overflow(W1,W2) || overflow(W1+W2,w0)) && (W1+W2+w0 > w2+w3);

	output = (p-output)*(negative_overflow)
			+ (output-p)*(!negative_overflow && greater_than_p)
			+ (output+gap)*(!negative_overflow && !greater_than_p && double_overflow)
		   	+ output*(!negative_overflow && !greater_than_p && !double_overflow);

	return output;
}

void add(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, add_uint64_t(x.re,y.re), add_uint64_t(x.img,y.img));
}

void sub(gauss_t *z, const gauss_t x, const gauss_t y) {
	set_gauss(z, sub_uint64_t(x.re,y.re), sub_uint64_t(x.img,y.img));
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {	
	uint64_t s1, s2, s3;

	s1 = mul_uint64_t(x.re, y.re);
	s2 = mul_uint64_t(x.img, y.img);	
	s3 = mul_uint64_t(add_uint64_t(x.re, x.img),  add_uint64_t(y.re, y.img));

	set_gauss(z, sub_uint64_t(s1,s2), sub_uint64_t(sub_uint64_t(s3,s1),s2));
}
