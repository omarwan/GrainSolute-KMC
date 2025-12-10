// energy_calculation.c
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "energy_calculation.h"
#include "update_boundary.h"
#include "local_energy.h"
#include "printer.h"
#include "affected_sites.h"
#include "periodic_neighbor.h"

// Function to check if a point is already in the unique list
bool is_already_in(Point *unique_sites, int unique_count, Point p) {
    for (int i = 0; i < unique_count; i++) {
        if (unique_sites[i].x == p.x && unique_sites[i].y == p.y) {
            return true; // Duplicate found
        }
    }
    return false;
}

// Filter unique sites from check_sites
int get_unique_sites(Point *check_sites, int site_count, Point *unique_sites) {
    int unique_count = 0;

    for (int i = 0; i < site_count; i++) {
        if (!is_already_in(unique_sites, unique_count, check_sites[i])) {
            unique_sites[unique_count++] = check_sites[i];
        }
    }

    return unique_count; // Return count of unique sites
}

EnergyResult get_energy_and_Eo(int** grains_grid, int** solutes_grid, int** boundary_info, double **energy_cache, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y, int i, int j, int x, int y, 
                         int q_selected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                         double Ed_oo, double Eg_oo, double Ecorr, ProcessType process)
{
    (void)energy_cache;
    (void)solute_neighbor_counts;

    EnergyResult result = {0.0, 0.0}; // Initialize result structure

    double energy_before = 0, energy_after = 0;
    if (process == PROC_FLIP){
        int max_affected = 13;
        AffectedSite affected_sites[max_affected];
        get_affected_sites_for_energy(affected_sites,  i,  j,  grid_size_x,  grid_size_y);


        // 1) Compute local energy BEFORE
        for (int s = 0; s < max_affected; s++) {
            energy_before += compute_local_energy(
                affected_sites[s].i, affected_sites[s].j,
                grains_grid, solutes_grid, boundary_info,
                grid_size_x, grid_size_y,
                q_selected, Jgg, Js, Jsg, Jss, Jsgs, F
            );
        }

        // 2) Save old grain type to revert after
        int old_grain = grains_grid[i][j];

        // 3) Apply the flip in place
        if (x == -1 && y == -1) {
            grains_grid[i][j] = -1; 
        }
        else if (x == -2 && y == -2) {
            grains_grid[i][j] = q_selected;
        }
        else {
            grains_grid[i][j] = grains_grid[x][y];
        }

        // 4) Update boundary only for the flipped site (match v1 behaviour)
        update_boundary_site(grains_grid, boundary_info, grid_size_x, grid_size_y, i, j);

        // 5) Compute local energy AFTER
        for (int s = 0; s < max_affected; s++) {
            energy_after += compute_local_energy(
                affected_sites[s].i, affected_sites[s].j,
                grains_grid, solutes_grid, boundary_info,
                grid_size_x, grid_size_y,
                q_selected, Jgg, Js, Jsg, Jss, Jsgs, F
            );
        }

        // 6) Revert | this is just a check
        grains_grid[i][j] = old_grain;
        update_boundary_site(grains_grid, boundary_info, grid_size_x, grid_size_y, i, j);

        result.uij = energy_after - energy_before;
        result.Eo = Eg_oo;
    }
    else if (process == PROC_SWAP) {

        // Define all sites to check (neighbors of both (i, j) and (x, y))
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

        Point sites[8];
        int unique_count = get_unique_sites(check_sites, 9, sites);

        // for (int i = 0; i < unique_count; i++) {
        //         printf("Unique Site: (%d, %d) , %d \n", sites[i].x, sites[i].y, unique_count);
        //     }

        for (int s = 0; s < unique_count; s++) {
            int cx = sites[s].x;
            int cy = sites[s].y;
            energy_before += compute_local_energy(
                cx, cy,
                grains_grid, solutes_grid, boundary_info,
                grid_size_x, grid_size_y,
                q_selected, Jgg, Js, Jsg, Jss, Jsgs, F
            );
        }
        // printf("%d \n",solutes_grid[x][y]);

        // 2) Apply the swap in place
        int temp = solutes_grid[i][j];
        solutes_grid[i][j] = solutes_grid[x][y];
        solutes_grid[x][y] = temp;
        // print_matrix(solutes_grid, grid_size_x, grid_size_y);
        for (int s = 0; s < unique_count; s++) {
            int cx = sites[s].x;
            int cy = sites[s].y;
            energy_after += compute_local_energy(
                cx, cy,
                grains_grid, solutes_grid, boundary_info,
                grid_size_x, grid_size_y,
                q_selected, Jgg, Js, Jsg, Jss, Jsgs, F
            );
        }

        // 4) Reapply the swap in place
        temp = solutes_grid[i][j];
        solutes_grid[i][j] = solutes_grid[x][y];
        solutes_grid[x][y] = temp;
        // print_matrix(solutes_grid, grid_size_x, grid_size_y);

        result.uij = energy_after - energy_before;
        double n_avg = (boundary_info[i][j] + boundary_info[x][y]) / 2.0; 
        result.Eo = Ed_oo * (1 - Ecorr * phi_formula(n_avg));
    }
    return result;
}

void update_energy_cache_for_sites(AffectedSite *sites, int num_sites, double **energy_cache, int **grains_grid, int **solutes_grid, int **boundary_info, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y,
                         int q_selected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F) {
    for (int s = 0; s < num_sites; s++) {
        int cx = sites[s].i;
        int cy = sites[s].j;
        // update solute neighbor counts for this site
        Neighbors n = periodic_neighbor(cx, cy, grid_size_x, grid_size_y);
        Point neighbors[4] = {n.left_neighbor, n.right_neighbor, n.top_neighbor, n.bottom_neighbor};
        unsigned char cnt = 0;
        for (int k = 0; k < 4; k++) {
            if (solutes_grid[neighbors[k].x][neighbors[k].y] > 0) cnt++;
        }
        solute_neighbor_counts[cx][cy] = cnt;
        energy_cache[cx][cy] = compute_local_energy(
            cx, cy,
            grains_grid, solutes_grid, boundary_info,
            grid_size_x, grid_size_y,
            q_selected, Jgg, Js, Jsg, Jss, Jsgs, F
        );
    }
}


double get_transition_probability(double uij, double Eo, double T) {
    if (T <= 0) {
        fprintf(stderr, "Error: Temperature (T) must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (Eo <= 0) {
        fprintf(stderr, "Error: Energy offset (Eo) must be positive.\n");
        exit(EXIT_FAILURE);
    }

    double Eij;
    if (uij > 0) {
        Eij = uij + Eo * exp(-uij / (2.0 * Eo));
    } else {
        Eij = Eo * exp(uij / (2.0 * Eo));
    }

    // Avoid numerical overflow/underflow
    double Wij = exp(-Eij / T);
    if (Wij < 1e-300) {
        Wij = 1e-300; // Avoid underflow
    }
    return Wij;
}
