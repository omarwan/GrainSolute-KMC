#include "rng.h"
#include <math.h>

// SplitMix64 for seeding xorshift64*
static uint64_t splitmix64_next(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

void rng_seed(RngState *rng, uint64_t seed) {
    rng->state = splitmix64_next(&seed);
    if (rng->state == 0) {
        rng->state = 0xdeadbeefULL; // avoid zero state
    }
}

uint32_t rng_next_u32(RngState *rng) {
    uint64_t x = rng->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng->state = x;
    return (uint32_t)((x * 0x2545F4914F6CDD1DULL) >> 32);
}

double rng_next_uniform(RngState *rng) {
    // 24 bits of mantissa precision is enough for uniform; scale to [0,1)
    const double inv = 1.0 / 16777216.0; // 2^24
    return (rng_next_u32(rng) >> 8) * inv;
}

double rng_next_exponential(RngState *rng, double rate) {
    double u = rng_next_uniform(rng);
    if (u <= 1e-12) {
        u = 1e-12; // avoid log(0)
    }
    return -log(u) / rate;
}
