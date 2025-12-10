#include "affected_sites.h"
#include <stdio.h>

// Populate the list of 13 affected sites
void get_affected_sites_for_energy(AffectedSite *affected_sites, int i, int j, int grid_size_x, int grid_size_y) {
    int index = 0;

    // Sites along the same column (varying i)
    for (int di = -2; di <= 2; di++) {
        int ni = (i + di + grid_size_x) % grid_size_x;
        affected_sites[index++] = (AffectedSite){ni, j};
    }

    // Sites along the same row (varying j)
    for (int dj = -2; dj <= 2; dj++) {
        if (dj != 0) { // Skip (i, j) since it's already added
            int nj = (j + dj + grid_size_y) % grid_size_y;
            affected_sites[index++] = (AffectedSite){i, nj};
        }
    }

    // Diagonal neighbors
    affected_sites[index++] = (AffectedSite){(i + 1 + grid_size_x) % grid_size_x, (j + 1 + grid_size_y) % grid_size_y};
    affected_sites[index++] = (AffectedSite){(i - 1 + grid_size_x) % grid_size_x, (j + 1 + grid_size_y) % grid_size_y};
    affected_sites[index++] = (AffectedSite){(i + 1 + grid_size_x) % grid_size_x, (j - 1 + grid_size_y) % grid_size_y};
    affected_sites[index++] = (AffectedSite){(i - 1 + grid_size_x) % grid_size_x, (j - 1 + grid_size_y) % grid_size_y};
}


// Populate the list of affected sites in a 7x7 grid around (i, j)
void get_affected_sites(AffectedSite *affected_sites, int i, int j, int grid_size_x, int grid_size_y) {
    int index = 0;

    // Loop through a 7x7 grid centered around (i, j)
    for (int di = -3; di <= 3; di++) {
        for (int dj = -3; dj <= 3; dj++) {
            int ni = (i + di + grid_size_x) % grid_size_x; // Apply periodic boundary conditions for x
            int nj = (j + dj + grid_size_y) % grid_size_y; // Apply periodic boundary conditions for y

            affected_sites[index++] = (AffectedSite){ni, nj};
        }
    }
}

int get_swap_energy_sites(AffectedSite *affected_sites, int i, int j, int x, int y, int grid_size_x, int grid_size_y) {
    // Collect sites that contribute to energy when swapping solutes at (i,j) and (x,y)
    Point check_sites[9] = {
        {i, j},
        {(i - 1 + grid_size_x) % grid_size_x, j},
        {(i + 1) % grid_size_x, j},
        {i, (j - 1 + grid_size_y) % grid_size_y},
        {i, (j + 1) % grid_size_y},
        {(x - 1 + grid_size_x) % grid_size_x, y},
        {(x + 1) % grid_size_x, y},
        {x, (y - 1 + grid_size_y) % grid_size_y},
        {x, (y + 1) % grid_size_y}
    };

    int count = 0;
    for (int idx = 0; idx < 9; idx++) {
        int cx = check_sites[idx].x;
        int cy = check_sites[idx].y;
        int already = 0;
        for (int k = 0; k < count; k++) {
            if (affected_sites[k].i == cx && affected_sites[k].j == cy) {
                already = 1;
                break;
            }
        }
        if (!already) {
            affected_sites[count++] = (AffectedSite){cx, cy};
        }
    }
    return count;
}
