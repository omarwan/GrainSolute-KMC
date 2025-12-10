#ifndef AFFECTED_SITES_H
#define AFFECTED_SITES_H

#include "common_structs.h"
void get_affected_sites_for_energy(AffectedSite *affected_sites, int i, int j, int grid_size_x, int grid_size_y);
void get_affected_sites(AffectedSite *affected_sites, int i, int j, int grid_size_x, int grid_size_y);
int get_swap_energy_sites(AffectedSite *affected_sites, int i, int j, int x, int y, int grid_size_x, int grid_size_y);
#endif // AFFECTED_SITES_H
