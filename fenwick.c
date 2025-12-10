#include <stdlib.h>
#include <string.h>
#include "fenwick.h"

int fenwick_init(FenwickTree *ft, int n) {
    ft->tree = calloc((size_t)n + 1, sizeof(double));
    if (!ft->tree) {
        ft->size = 0;
        return -1;
    }
    ft->size = n;
    return 0;
}

void fenwick_free(FenwickTree *ft) {
    free(ft->tree);
    ft->tree = NULL;
    ft->size = 0;
}

void fenwick_add(FenwickTree *ft, int idx, double delta) {
    int i = idx + 1; // fenwick is 1-based
    while (i <= ft->size) {
        ft->tree[i] += delta;
        i += i & -i;
    }
}

void fenwick_build(FenwickTree *ft, double *values, int n) {
    // assumes ft already initialized with size n
    memset(ft->tree, 0, (size_t)(n + 1) * sizeof(double));
    for (int i = 0; i < n; i++) {
        fenwick_add(ft, i, values[i]);
    }
}

double fenwick_prefix_sum(const FenwickTree *ft, int idx) {
    double sum = 0.0;
    int i = idx + 1;
    while (i > 0) {
        sum += ft->tree[i];
        i -= i & -i;
    }
    return sum;
}

int fenwick_find_prefix(const FenwickTree *ft, double target) {
    // Binary lifting over Fenwick tree.
    int idx = 0;
    int bit;
    // Find largest power of two <= size.
    for (bit = 1; bit <= ft->size; bit <<= 1) {
        // loop to highest bit
    }
    bit >>= 1;

    double accumulated = 0.0;
    while (bit) {
        int next = idx + bit;
        if (next <= ft->size && accumulated + ft->tree[next] < target) {
            accumulated += ft->tree[next];
            idx = next;
        }
        bit >>= 1;
    }
    // idx is last position with prefix < target, so answer is idx
    // since we used idx starting at 0 (fenwick uses 1-based), return idx
    return idx; // caller typically uses idx as 0-based index where prefix >= target
}
