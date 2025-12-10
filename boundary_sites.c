// boundary_sites.c
#include <stdlib.h>
#include "boundary_sites.h"

int **find_boundary_sites_and_neighbors(int **boundary_info, int **grains_grid, int grid_size_x, int grid_size_y) {
    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            Neighbors n = periodic_neighbor(i, j, grid_size_x, grid_size_y);
            int count = 0;
            if (grains_grid[n.left_neighbor.x][n.left_neighbor.y] != grains_grid[i][j]) count++;
            if (grains_grid[n.right_neighbor.x][n.right_neighbor.y] != grains_grid[i][j]) count++;
            if (grains_grid[n.top_neighbor.x][n.top_neighbor.y] != grains_grid[i][j]) count++;
            if (grains_grid[n.bottom_neighbor.x][n.bottom_neighbor.y] != grains_grid[i][j]) count++;
            boundary_info[i][j] = count;
        }
    }

    return boundary_info;
}
