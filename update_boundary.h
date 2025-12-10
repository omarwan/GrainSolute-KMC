// update_boundary.h
#ifndef UPDATE_BOUNDARY_H
#define UPDATE_BOUNDARY_H
#include "common_structs.h"
#include "periodic_neighbor.h"

void update_boundary_site(int **grains_grid, int **current_boundary_info, int grid_size_x, int grid_size_y, int i, int j);

#endif // UPDATE_BOUNDARY_H
