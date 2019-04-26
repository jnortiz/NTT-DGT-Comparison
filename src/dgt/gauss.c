#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

uint64_t p = 0xFFFFFFFF00000001; 
// 18 446 744 069 414 584 321
// 18 446 744 065 119 617 026
struct gauss {
	uint64_t re;
	uint64_t img;
};

typedef struct gauss gauss_t;

/* "Constructor" for type gauss_t */
void Gauss(gauss_t *pointer, uint64_t re, uint64_t img) {
	pointer->re = re;
	pointer->img = img;
}

/* Print a Gaussian integer */
void print_gauss_t(const gauss_t x) {
	printf("%" PRId64, x.re);
	if(x.img >= 0)
		printf("+");
	printf("%" PRId64 "i\n", x.img);
}

/* Compute the conjugate of a Gaussian integer */
gauss_t conjugate(const gauss_t x){
	gauss_t conj_x;
	Gauss(&conj_x, x.re, -x.img);
	return conj_x;
}

/* Return true if x is equal to y */
/*int is_equal(const gauss_t x, const gauss_t y) {
	return (x.re == y.re && x.img == y.img);
}*/

unsigned less_than(const uint64_t x, const uint64_t y) {
	uint64_t less = (x-y);
	return (less >> sizeof(uint64_t)*8-1);
}

int64_t select_int64_t(uint64_t x, uint64_t y, uint64_t bit) {
	uint64_t mask = -bit;
	int64_t output = mask & (x ^ y);
	return (output ^ x);
}

int64_t is_equal(uint64_t x, uint64_t y) {
	uint64_t zero = (x-y);
	return (1 ^ ((zero | -zero) >> 63) & 1);
}

/* Compute a mod b in Z_p */
uint64_t mod(unsigned __int128 a) {	

	uint64_t w0, w1, w2, w3, c;

	w0 = a & 0xFFFFFFFF;
	//printf("%" PRIu64 "\n", w0);
	w1 = ((uint64_t)(a >> 32)) & 0xFFFFFFFF;
	//printf("%" PRIu64 "\n", w1);
	w2 = ((uint64_t)(a >> 64)) & 0xFFFFFFFF;
	//printf("%" PRIu64 "\n", w2);
	w3 = (uint64_t)(a >> 96);
	//printf("%" PRIu64 "\n", w3);
	c = (w2+w1) << 32;
	//printf("%" PRIu64 "\n", c);

	c = c-w3-w2+w0;
	printf("Before select %" PRIu64 " | ", c);
	c = select_int64_t(c,c-p,less_than(p, c)||is_equal(p,c));
	printf("After select %" PRIu64 " ", c);

	return c;
}

/*
void add(gauss_t *z, const gauss_t x, const gauss_t y) {
	z->re = mod(x.re + y.re);
	z->img = mod(x.img + y.img);
}

void sub(gauss_t *z, const gauss_t x, const gauss_t y) {
	z->re = mod(x.re - y.re);
	z->img = mod(x.img - y.img);
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {
	__int128 s1, s2, s3;

	s1 = x.re * y.re;
	s2 = x.img * y.img;
	s3 = (x.re + x.img)*(y.re + y.img);

	z->re = mod(s1 - s2);
	z->img = mod(s3 - s1 - s2);
}*/

int main(int argc, char** argv[]) {
	//gauss_t x, y, z;
	//Gauss(&x, 1, -2);
	//Gauss(&y, 1, 2);

	//sub(&z, x, y);
	//print_gauss_t(z);
	unsigned __int128 a;
	int i;
	for(i = 0; i < 32; i++) {
		a = (p << i);
		printf("Reducing: p << %d | ", i);
		mod(a);
		printf("\n");
	}

	return 0;
}