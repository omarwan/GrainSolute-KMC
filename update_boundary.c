#include <stdlib.h>
#include "update_boundary.h"

void update_boundary_site(int **grains_grid, int **boundary_info, int grid_size_x, int grid_size_y, int i, int j) {
    // Iterate over the center site and its four neighbors
    Neighbors center_neighbors = periodic_neighbor(i, j, grid_size_x, grid_size_y);
    int xs[5] = {i, center_neighbors.left_neighbor.x, center_neighbors.right_neighbor.x, center_neighbors.top_neighbor.x, center_neighbors.bottom_neighbor.x};
    int ys[5] = {j, center_neighbors.left_neighbor.y, center_neighbors.right_neighbor.y, center_neighbors.top_neighbor.y, center_neighbors.bottom_neighbor.y};

    for (int idx = 0; idx < 5; idx++) {
        int x = xs[idx];
        int y = ys[idx];
        Neighbors n = periodic_neighbor(x, y, grid_size_x, grid_size_y);
        int count = 0;
        if (grains_grid[n.left_neighbor.x][n.left_neighbor.y] != grains_grid[x][y]) count++;
        if (grains_grid[n.right_neighbor.x][n.right_neighbor.y] != grains_grid[x][y]) count++;
        if (grains_grid[n.top_neighbor.x][n.top_neighbor.y] != grains_grid[x][y]) count++;
        if (grains_grid[n.bottom_neighbor.x][n.bottom_neighbor.y] != grains_grid[x][y]) count++;
        boundary_info[x][y] = count;
    }
}
