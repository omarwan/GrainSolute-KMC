// periodic_neighbor.c
#include "periodic_neighbor.h"
#include <stdlib.h>
#include <string.h>

static Neighbors *neighbor_cache = NULL;
static int cache_x = -1;
static int cache_y = -1;

void init_periodic_neighbor_cache(int grid_size_x, int grid_size_y) {
    size_t count = (size_t)grid_size_x * (size_t)grid_size_y;
    free(neighbor_cache);
    neighbor_cache = malloc(count * sizeof(Neighbors));
    if (!neighbor_cache) {
        cache_x = cache_y = -1;
        return;
    }
    cache_x = grid_size_x;
    cache_y = grid_size_y;

    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            Neighbors neighbors;

            neighbors.left_neighbor.x = (i - 1 + grid_size_x) % grid_size_x;
            neighbors.left_neighbor.y = j;

            neighbors.right_neighbor.x = (i + 1) % grid_size_x;
            neighbors.right_neighbor.y = j;

            neighbors.top_neighbor.x = i;
            neighbors.top_neighbor.y = (j - 1 + grid_size_y) % grid_size_y;

            neighbors.bottom_neighbor.x = i;
            neighbors.bottom_neighbor.y = (j + 1) % grid_size_y;

            neighbor_cache[i * grid_size_y + j] = neighbors;
        }
    }
}

void free_periodic_neighbor_cache(void) {
    free(neighbor_cache);
    neighbor_cache = NULL;
    cache_x = cache_y = -1;
}

Neighbors periodic_neighbor(int i, int j, int grid_size_x, int grid_size_y) {
    if (neighbor_cache && cache_x == grid_size_x && cache_y == grid_size_y) {
        return neighbor_cache[i * grid_size_y + j];
    }

    Neighbors neighbors;

    neighbors.left_neighbor.x = (i - 1 + grid_size_x) % grid_size_x;
    neighbors.left_neighbor.y = j;

    neighbors.right_neighbor.x = (i + 1) % grid_size_x;
    neighbors.right_neighbor.y = j;

    neighbors.top_neighbor.x = i;
    neighbors.top_neighbor.y = (j - 1 + grid_size_y) % grid_size_y;

    neighbors.bottom_neighbor.x = i;
    neighbors.bottom_neighbor.y = (j + 1) % grid_size_y;

    return neighbors;
}
