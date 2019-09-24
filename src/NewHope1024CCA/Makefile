CC = /usr/bin/gcc
LDFLAGS = -lcrypto

SOURCES = kem.c cpapke.c fips202.c ntt.c dgt.c cpucycles.c poly.c precomp.c reduce.c rng.c verify.c PQCgenKAT_kem.c
HEADERS = api.h cpapke.h fips202.h ntt.h dgt.h cpucycles.c params.h poly.h reduce.h rng.h verify.h

PQCgenKAT_kem: $(SOURCES) $(HEADERS)
	$(CC) $(SOURCES) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	-rm PQCgenKAT_kem
