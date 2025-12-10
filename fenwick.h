#ifndef FENWICK_H
#define FENWICK_H

typedef struct {
    double *tree;
    int size; // number of elements (not tree length)
} FenwickTree;

// Allocate a Fenwick tree for n elements. Returns 0 on success.
int fenwick_init(FenwickTree *ft, int n);

// Free allocated memory.
void fenwick_free(FenwickTree *ft);

// Build the tree from an array of length n.
void fenwick_build(FenwickTree *ft, double *values, int n);

// Add delta to position idx (0-based).
void fenwick_add(FenwickTree *ft, int idx, double delta);

// Return prefix sum in [0, idx] (0-based).
double fenwick_prefix_sum(const FenwickTree *ft, int idx);

// Find the smallest index such that prefix sum >= target. Assumes target < total sum.
int fenwick_find_prefix(const FenwickTree *ft, double target);

// Return total sum.
static inline double fenwick_total(const FenwickTree *ft) {
    return fenwick_prefix_sum(ft, ft->size - 1);
}

#endif
