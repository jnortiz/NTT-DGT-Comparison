#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>

uint64_t p = 0xFFFFFFFF00000001; 
unsigned __int128 max_128 = ((unsigned __int128)0xFFFFFFFFFFFFFFFFu << 64) + 0xFFFFFFFFFFFFFFFFu;

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

void print_gauss_t(const gauss_t x) {
	printf("%" PRId64, x.re);
	if(x.img >= 0)
		printf("+");
	printf("%" PRId64 "i\n", x.img);
}

gauss_t conjugate(const gauss_t x){
	gauss_t conj_x;
	Gauss(&conj_x, x.re, -x.img);
	return conj_x;
}

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

uint64_t mod_uint128(const unsigned __int128 a) {	
	uint64_t gap, w0, w1, w2, w3, W1, W2, sum, output;
	int8_t negative_overflow, greater_than_p, double_overflow;

	gap = (uint64_t)(18446744073709551615u)-p;

	w0 = (uint32_t) a;
	w1 = (uint32_t)(a >> 32);
	w2 = (uint32_t)(a >> 64);
	w3 = (uint32_t)(a >> 96);

	W1 = w1 << 32;
	W2 = w2 << 32;
	sum = W1+W2;

	output = sum-w3-w2+w0;

	negative_overflow = (w2+w3 > sum+w0) && !(less_than(W1,W2)||less_than(sum,w0));
	greater_than_p = (output >= p);
	double_overflow = (less_than(W1,W2) || less_than(sum,w0)) && (sum+w0 > w2+w3);

	output = (p-output)*negative_overflow
			+ (output-p)*(!negative_overflow && greater_than_p)
			+ (output+gap)*(!negative_overflow && !greater_than_p && double_overflow)
		   	+ output*(!(negative_overflow || greater_than_p || double_overflow));

	return output;
}

void add(gauss_t *z, const gauss_t x, const gauss_t y) {
	z->re = mod_uint64_t(x.re + y.re);
	z->img = mod_uint64_t(x.img + y.img);
}

void sub(gauss_t *z, const gauss_t x, const gauss_t y) {
	z->re = mod_uint64_t(x.re - y.re);
	z->img = mod_uint64_t(x.img - y.img);
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {
	__int128 s1, s2, s3;

	s1 = x.re * y.re;
	s2 = x.img * y.img;
	s3 = (x.re + x.img)*(y.re + y.img);

	z->re = mod_uint128(s1 - s2);
	z->img = mod_uint128(s3 - s1 - s2);
}

unsigned __int128 get_uint128_word() {
	__int128 a;
	a = rand() % max_128;
	return a;
}

int main(int argc, char** argv[]) {
	srand(time(NULL));
	int i;
	for(i = 0; i < 128; i++) {
		printf("%" PRIu64 "\n", mod_uint128(get_uint128_word()));
	}
	return 0;
}