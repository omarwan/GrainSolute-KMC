// periodic_neighbor.h
#ifndef PERIODIC_NEIGHBOR_H
#define PERIODIC_NEIGHBOR_H
#include "common_structs.h"

Neighbors periodic_neighbor(int i, int j, int grid_size_x, int grid_size_y);
void init_periodic_neighbor_cache(int grid_size_x, int grid_size_y);
void free_periodic_neighbor_cache(void);

#endif 
