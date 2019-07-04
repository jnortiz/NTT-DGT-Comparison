#include <inttypes.h>
#include "gaussian_integer.h"
#include "poly.h"

void set_gauss(gauss_t *pointer, const int32_t re, const int32_t img) {
	pointer->re = re;
	pointer->img = img;
}

void mul(gauss_t *z, const gauss_t x, const gauss_t y) {	
	int32_t s1, s2, s3;

	s1 = reduce((int64_t)x.re*y.re);
	s2 = reduce((int64_t)x.img*y.img);	
	s3 = reduce((int64_t)(x.re+x.img)*(y.re+y.img));
	set_gauss(z, barr_reduce(s1-s2), barr_reduce(s3-s1-s2));
}