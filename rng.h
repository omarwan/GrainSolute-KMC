#ifndef RNG_H
#define RNG_H

#include <stdint.h>

typedef struct {
    uint64_t state;
} RngState;

void rng_seed(RngState *rng, uint64_t seed);
uint32_t rng_next_u32(RngState *rng);
double rng_next_uniform(RngState *rng);       // [0,1)
double rng_next_exponential(RngState *rng, double rate); // Exponential with given rate

#endif
