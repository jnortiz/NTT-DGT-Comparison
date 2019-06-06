#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

#include "gauss.h"

uint64_t p = 0xFFFFFFFF00000001; 
int64_t n = 512;

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

uint64_t mod_uint64_t(const uint64_t a) {
	uint64_t output = a;	
	if(a >= p) {
		output -= p;
	}
	return output;	
}

uint64_t sub_uint64_t(const uint64_t a, const uint64_t b) {
	uint64_t res;
	res = (a-b)+(b>a)*p;	
	return res;
}

uint64_t add_uint64_t(const uint64_t a, const uint64_t b) {
	uint64_t res = a+b;
	res += (res < a)*4294967295L;
	return res;
}

int overflow(const uint64_t a, const uint64_t b) {
	uint64_t sum = a+b;
	return(sum<a);
}

uint64_t mod_uint128(const unsigned __int128 a) {
	uint64_t half;
	half = (uint64_t)(a >> 64);

	if(half == (uint64_t)0) {
		return mod_uint64_t((uint64_t)a);
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
	unsigned __int128 a, b, xre, ximg, yre, yimg;
	uint64_t s1, s2, s3;

	xre = (unsigned __int128) x.re;
	yre = (unsigned __int128) y.re;
	s1 = mod_uint128(xre*yre);

	ximg = (unsigned __int128) x.img;
	yimg = (unsigned __int128) y.img;
	s2 = mod_uint128(ximg*yimg);

	a = (unsigned __int128) mod_uint128(xre + ximg);
	b = (unsigned __int128) mod_uint128(yre + yimg);
	s3 = mod_uint128(a*b);

	set_gauss(z, sub_uint64_t(s1,s2), sub_uint64_t(sub_uint64_t(s3,s1),s2));
}

void test_add() {
	gauss_t a, b, c;
	set_gauss(&a, 1, 0);
	set_gauss(&b, 2, 1);
	add(&c, a, b);
    printf("(%" PRIu64 ", %" PRIu64 ")\n", c.re, c.img);
}

void test_mod() {
	gauss_t a;
	a.re = mod_uint64_t((uint64_t)-1);
	a.img = mod_uint64_t((uint64_t)1);
    printf("(%" PRIu64 ", %" PRIu64 ")\n", a.re, a.img);
}

void test_large_mod() {
	unsigned __int128 c, a, b;
	uint64_t c_mod, a_mod, b_mod;

	a = (unsigned __int128)(p) << 1;
	a += 1;
	b = (unsigned __int128) p;
	b += 2;
	c = (unsigned __int128)(p) << 2;
	c += 3;

	a_mod = mod_uint128(a);
	b_mod = mod_uint128(b);
	c_mod = mod_uint128(c);

    printf("a: %" PRIu64 "\n", a_mod);
    printf("b: %" PRIu64 "\n", b_mod);
    printf("c: %" PRIu64 "\n", c_mod);
}

void test_sub() {
	gauss_t a, b, c;
	set_gauss(&a, 1, 0);
	set_gauss(&b, 2, 1);	
	sub(&c, a, b);
    printf("(%" PRIu64 ", %" PRIu64 ")\n", c.re, c.img);
}

void test_mul() {
	gauss_t a, b, c;
	set_gauss(&a, 3, 2);
	set_gauss(&b, 5, 1);	
	mul(&c, b, a);
    printf("(%" PRIu64 ", %" PRIu64 ")\n", c.re, c.img);
}

void test_large_mul() {
	gauss_t a, b, c;
	set_gauss(&a, p, 1);
	set_gauss(&b, 656189816218, 25679846546);	
	mul(&c, b, a);
    printf("(%" PRIu64 ", %" PRIu64 ")\n", c.re, c.img);
}

void test_add_uint64_t() {
	uint64_t a, b;

	a = 1024;
	//b = 18446744073709551615u; //2^64-1
	b = p;

    printf("a+b mod p: %" PRIu64 "\n", add_uint64_t(a,b));
}

void test_mod_128() {

	unsigned __int128 c;
	uint64_t a, b;
	a = 0xFFFFFFFF00000001u;
	b = rand() % 18446744073709551615u;
    printf("a: %" PRIu64 "\n", a);	
    printf("b: %" PRIu64 "\n", b);	

	c = ((unsigned __int128)(a) << 64) + b;
    printf("a mod p: %" PRIu64 "\n", mod_uint128(c));	

}


// int main(int argc, char** argv[]) {    

// 	test_mod_128();

// }