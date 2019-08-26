#ifndef DGT_H
#define DGT_H

#include "inttypes.h"

static uint16_t gj[];
static uint16_t invgj[];

void dgt(uint16_t *poly);
void idgt(uint16_t *poly);

#endif