#ifndef DGT_H
#define DGT_H

#include <stdint.h>
#include "params.h"

#define dgt DILITHIUM_NAMESPACE(dgt)
void dgt(uint32_t p[N]);

#define invdgt_tomont DILITHIUM_NAMESPACE(invdgt_tomont)
void invdgt_tomont(uint32_t p[N]);

#endif
