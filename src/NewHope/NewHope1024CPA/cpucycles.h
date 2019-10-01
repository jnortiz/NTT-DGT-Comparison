#ifndef CPUCYCLES_H
#define CPUCYCLES_H
    
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#define TARGET TARGET_AMD64
#define RADIX           64
#define RADIX32         32
typedef uint64_t        digit_t;        // Unsigned 64-bit digit
typedef int64_t         sdigit_t;       // Signed 64-bit digit

#define print_unit printf("cycles");
    
// Access system counter for benchmarking
int64_t cpucycles(void);

#endif
