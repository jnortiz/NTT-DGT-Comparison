CC = /usr/bin/gcc
LDFLAGS = -lcrypto

SOURCES = kem.c cpapke.c fips202.c cpucycles.c dgt.c poly.c reduce.c rng.c verify.c PQCgenKAT_kem.c
HEADERS = api.h cpapke.h fips202.h cpucycles.h dgt.h params.h poly.h reduce.h rng.h verify.h

PQCgenKAT_kem: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_kem PQCkemKAT_1888.req PQCkemKAT_1888.rsp
