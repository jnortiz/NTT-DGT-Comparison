CC = /usr/bin/gcc
LDFLAGS = -lcrypto

SOURCES = kem.c cpapke.c fips202.c dgt.c cpucycles.c poly.c reduce.c rng.c verify.c PQCgenKAT_kem.c
HEADERS = api.h cpapke.h fips202.h dgt.h cpucycles.h params.h poly.h reduce.h rng.h verify.h

SOURCES_MUL = test_mul.c poly.c reduce.c cpucycles.c rng.c dgt.c fips202.c
HEADERS_MUL = reduce.h cpucycles.h poly.h rng.h dgt.h fips202.h

PQCgenKAT_kem: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) -o $@ $(LDFLAGS)

test_mul: $(SOURCES_MUL) $(HEADERS_MUL)
	$(CC) $(SOURCES_MUL) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_kem PQCkemKAT_1792.req PQCkemKAT_1792.rsp test_mul
