#ifndef COMMON_STRUCTS_H
#define COMMON_STRUCTS_H

typedef struct {
    int i;
    int j;
} AffectedSite;

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    Point left_neighbor;
    Point right_neighbor;
    Point top_neighbor;
    Point bottom_neighbor;
} Neighbors;

typedef struct {
    double uij;
    double Eo;
} EnergyResult;

typedef struct {
    int i;
    int j;
    int x;
    int y;
    double rate;
} TransitionRate;

typedef struct {
    int site_index;     // The 1D site index in the grid
    int flip_rate_index; // The index in the flip_rates array corresponding to this site
    double old_rate;    // The old rate before the update
    double new_rate;    // The new rate after the update
    double rate_diff;   // The difference (new_rate - old_rate)
} AffectedMaps;

typedef enum {
    PROC_FLIP = 0,
    PROC_SWAP = 1
} ProcessType;

#endif // COMMON_STRUCTS_H
