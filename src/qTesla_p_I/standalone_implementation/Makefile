CC = /usr/bin/gcc
LDFLAGS = -lcrypto

SOURCES = dgt.c poly.c reduce.c main.c
HEADERS = dgt.h poly.h reduce.h params.h

test_dgt: $(SOURCES) $(HEADERS)
	$(CC) -lm $(SOURCES) -o $@ $(LDFLAGS)

.PHONY: clean

clean:
	-rm test_dgt
