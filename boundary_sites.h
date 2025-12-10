// boundary_sites.h
#ifndef BOUNDARY_SITES_H
#define BOUNDARY_SITES_H
#include "common_structs.h"
#include "periodic_neighbor.h"

int **find_boundary_sites_and_neighbors(int **boundary_info, int **grains_grid, int grid_size_x, int grid_size_y);

#endif // BOUNDARY_SITES_H
