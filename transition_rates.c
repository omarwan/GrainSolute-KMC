#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <stdbool.h>
#include "transition_rates.h"
#include "energy_calculation.h"  
#include "periodic_neighbor.h"
#include "boundary_sites.h"
#include "update_boundary.h"

#define MAX_GRAIN_ID 128 // Max possible grain ID
#define MAX_UNIQUE_GRAINS 7
#define FORCE_TERM_X -2
#define FORCE_TERM_Y -2
#define N4_TERM_X -1
#define N4_TERM_Y -1

TransitionRate* get_rates(int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y,
                        double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                        double T, double Ed_oo, double Eg_oo, double Ecorr, ProcessType process, int q_no, int q_selected){

    TransitionRate *transition_rates = NULL; // Initialize to NULL
    int num_rates;

    if (process == PROC_FLIP) {

        num_rates = grid_size_x * grid_size_y * 6;
        transition_rates = calloc(num_rates, sizeof(TransitionRate));
        if (transition_rates == NULL) {
            fprintf(stderr, "Transition rates array is not properly allocated or initialized.\n");
            return NULL; // Safe exit for a pointer type
        }

        int index = 0;
        for (int i = 0; i < grid_size_x; i++) {
            for (int j = 0; j < grid_size_y; j++) {

                Neighbors neighbors_ij = periodic_neighbor(i, j, grid_size_x, grid_size_y);
                Point points_ij[4] = {
                    neighbors_ij.left_neighbor,
                    neighbors_ij.right_neighbor,
                    neighbors_ij.top_neighbor,
                    neighbors_ij.bottom_neighbor};

                int num_unique_grains = 0;
                bool grain_exists[MAX_GRAIN_ID] = {0};

                // Mark the grain at (i, j) as unique
                int grain_type_ij = grains_grid[i][j];
                if (!grain_exists[grain_type_ij]) {
                    grain_exists[grain_type_ij] = true;
                    num_unique_grains++;
                }

                for (int l = 0; l < 4; l++) {
                    int x = points_ij[l].x;
                    int y = points_ij[l].y;

                    transition_rates[index].i = i;
                    transition_rates[index].j = j;
                    transition_rates[index].x = x;
                    transition_rates[index].y = y;

                    int grain_type = grains_grid[x][y];

                    // Check if grain already exists in the list
                    if (grain_exists[grain_type]) {
                        transition_rates[index].rate = 0;
                    } else {
                        EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, 
                                        i, j, x, y, q_selected, Jgg, Js, Jsg, Jss, Jsgs, F, 
                                        Ed_oo, Eg_oo, Ecorr, PROC_FLIP);

                        // printf("uij: %f, Eo: %f\n", result.uij, result.Eo);
                        double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                        transition_rates[index].rate = transition_probability;
                        // printf("Transition probability: %f\n", transition_probability);

                        grain_exists[grain_type] = true;
                        num_unique_grains++;
                    }
                    index++;
                }

                ////////////////////////////////////////////////
                //////// Force term
                ////////////////////////////////////////////////
                transition_rates[index].i = i;
                transition_rates[index].j = j;
                transition_rates[index].x = FORCE_TERM_X;
                transition_rates[index].y = FORCE_TERM_Y;

                // Check if q_selected is unique
                if (!grain_exists[q_selected]) {
                    grain_exists[q_selected] = true;
                    num_unique_grains++;

                    EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, i, j, FORCE_TERM_X, FORCE_TERM_Y, q_selected, 
                                    Jgg, Js, Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_FLIP);
                    double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                    transition_rates[index].rate = transition_probability;

                } else {
                    transition_rates[index].rate = 0;
                }
                index++;
                ////////////////////////////////////////////////
                //////// n=4 term
                ////////////////////////////////////////////////
                transition_rates[index].i = i;
                transition_rates[index].j = j;
                transition_rates[index].x = N4_TERM_X;
                transition_rates[index].y = N4_TERM_Y;
                int scale_factor =q_no - num_unique_grains;

                EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, i, j, N4_TERM_X, N4_TERM_Y, q_selected, Jgg, Js, 
                                Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_FLIP);

                double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                transition_rates[index].rate = transition_probability * scale_factor;

                index++;
            }
        }
    }

    else if (process == PROC_SWAP) {
        num_rates = grid_size_x * grid_size_y * 4;
        transition_rates = malloc(num_rates * sizeof(TransitionRate));
        if (transition_rates == NULL) {
            fprintf(stderr, "Failed to allocate memory for swap transition rates\n");
            return NULL;
        }

        int index = 0;
        for (int i = 0; i < grid_size_x; i++) {
            for (int j = 0; j < grid_size_y; j++) {

                Neighbors neighbors_ij = periodic_neighbor(i, j, grid_size_x, grid_size_y);
                Point points_ij[4] = {neighbors_ij.left_neighbor, neighbors_ij.right_neighbor, neighbors_ij.top_neighbor,neighbors_ij.bottom_neighbor};

                for (int l = 0; l < 4; l++) {
                    int x = points_ij[l].x;
                    int y = points_ij[l].y;

                    transition_rates[index].i = i;
                    transition_rates[index].j = j;
                    transition_rates[index].x = x;
                    transition_rates[index].y = y;
                    if (solutes_grid[i][j] > 0 && solutes_grid[x][y] == 0){
                        EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, i, j, x, y, q_selected, Jgg, Js,
                                            Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_SWAP);
                        double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                        transition_rates[index].rate = transition_probability;
                    }else{
                        transition_rates[index].rate = 0;
                    }
                    index++;
                }
            }
        }
    }
    return transition_rates;
}


double update_flip_rates_for_affected_sites(TransitionRate *transition_rates, int num_rates, int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_counts, int grid_size_x, int grid_size_y,
                                     AffectedSite *affected_sites, int num_affected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                                     double T, double Ed_oo, double Eg_oo, double Ecorr, int q_no, int q_selected, FenwickTree *tree) {

    (void)num_rates; // unused, kept for signature compatibility
    double initial_sum = 0.0, final_sum = 0.0;

    // Loop through affected sites
    for (int index = 0; index < num_affected; index++) {
        int i = affected_sites[index].i;
        int j = affected_sites[index].j;

        Neighbors neighbors_ij = periodic_neighbor(i, j, grid_size_x, grid_size_y);
        Point points_ij[4] = {neighbors_ij.left_neighbor, neighbors_ij.right_neighbor, neighbors_ij.top_neighbor, neighbors_ij.bottom_neighbor};

        int num_unique_grains = 0;
        bool grain_exists[MAX_GRAIN_ID] = {0}; // Bitset to track which grains exist

        // Add (i, j) grain to unique grains
        int grain_type_ij = grains_grid[i][j];
        grain_exists[grain_type_ij] = true;
        num_unique_grains++;

        int rate_index;

        for (int l = 0; l < 4; l++) {
            int x = points_ij[l].x;
            int y = points_ij[l].y;

            rate_index = (i * grid_size_y + j) * 6 + l;

            transition_rates[rate_index].i = i;
            transition_rates[rate_index].j = j;
            transition_rates[rate_index].x = x;
            transition_rates[rate_index].y = y;

            double old_rate = transition_rates[rate_index].rate;
            initial_sum += old_rate;
            int grain_type = grains_grid[x][y];

            if (grain_exists[grain_type]) {
                transition_rates[rate_index].rate = 0;
                if (tree) {
                    fenwick_add(tree, rate_index, -old_rate);
                }
            } else {
                EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_counts, grid_size_x, grid_size_y, i, j, x, y, q_selected,
                                Jgg, Js, Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_FLIP);
                double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                transition_rates[rate_index].rate = transition_probability;
                if (tree) {
                    fenwick_add(tree, rate_index, transition_probability - old_rate);
                }

                grain_exists[grain_type] = true;
                num_unique_grains++;
                final_sum += transition_probability;
            }
        }

        // ////////////////////////////////////////////////
        // //////// force term
        // ////////////////////////////////////////////////
        rate_index = (i * grid_size_y + j) * 6 + 4;

        transition_rates[rate_index].i = i;
        transition_rates[rate_index].j = j;
        transition_rates[rate_index].x = FORCE_TERM_X;
        transition_rates[rate_index].y = FORCE_TERM_Y;
        
        double old_rate = transition_rates[rate_index].rate;
        initial_sum += old_rate;

        // Check if q_selected is unique
        if (!grain_exists[q_selected]) {
            grain_exists[q_selected] = true;
            num_unique_grains++;

                EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_counts, grid_size_x, grid_size_y, i, j, FORCE_TERM_X, FORCE_TERM_Y, q_selected, 
                                Jgg, Js, Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_FLIP);
            double transition_probability = get_transition_probability(result.uij, result.Eo, T);
            transition_rates[rate_index].rate = transition_probability;
            if (tree) {
                fenwick_add(tree, rate_index, transition_probability - old_rate);
            }

            // Add q_selected to unique_grains
            final_sum += transition_probability;
        } else {
            transition_rates[rate_index].rate = 0;
            if (tree) {
                fenwick_add(tree, rate_index, -old_rate);
            }
        }
        ////////////////////////////////////////////////
        //////// n=4 term
        ////////////////////////////////////////////////
        rate_index = (i * grid_size_y + j) * 6 + 5;

        transition_rates[rate_index].i = i;
        transition_rates[rate_index].j = j;
        transition_rates[rate_index].x = N4_TERM_X;
        transition_rates[rate_index].y = N4_TERM_Y;

        old_rate = transition_rates[rate_index].rate;
        initial_sum += old_rate;

        int scale_factor = q_no - num_unique_grains;
        if (scale_factor < 0) {
            printf("scale_factor is negative !!");
        }
        EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_counts, grid_size_x, grid_size_y, i, j, 
                        N4_TERM_X, N4_TERM_Y, q_selected, Jgg, Js, Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_FLIP);
        // printf("uij: %f, Eo: %f\n", result.uij, result.Eo);

        double transition_probability = get_transition_probability(result.uij, result.Eo, T);
        transition_rates[rate_index].rate = transition_probability * scale_factor;
        // printf("Transition probability: %.50f\n", transition_probability);

        if (tree) {
            fenwick_add(tree, rate_index, transition_rates[rate_index].rate - old_rate);
        }

        final_sum += transition_probability * scale_factor;
    }
    return final_sum - initial_sum;
}


double update_swap_rates_for_affected_sites(TransitionRate *transition_rates, int num_rates, int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_counts, int grid_size_x, int grid_size_y,
                                     AffectedSite *affected_sites, int num_affected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                                     double T, double Ed_oo, double Eg_oo, double Ecorr, int q_selected, FenwickTree *tree) {

    double initial_sum = 0.0, final_sum = 0.0;

    // // Loop through affected sites
    for (int index = 0; index < num_affected; index++) {
        int i = affected_sites[index].i;
        int j = affected_sites[index].j;
    // for (int i = 0; i < grid_size_x; i++) {
    //     for (int j = 0; j < grid_size_y; j++) {

            // Calculate swap transition rates for each affected site with its neighbors
            Neighbors neighbors = periodic_neighbor(i, j, grid_size_x, grid_size_y);
            Point neighbor_points[4] = {neighbors.left_neighbor, neighbors.right_neighbor, neighbors.top_neighbor, neighbors.bottom_neighbor};

            for (int k = 0; k < 4; k++) {
                int x = neighbor_points[k].x;
                int y = neighbor_points[k].y;


                // Index for the transition rate in the array should be based on both positions (i, j) and (x, y)
                int rate_index = (i * grid_size_y + j) * 4 + k;
                if (rate_index < num_rates) {
                    initial_sum += transition_rates[rate_index].rate;
                    transition_rates[rate_index].i = i;
                    transition_rates[rate_index].j = j;
                    transition_rates[rate_index].x = x;
                    transition_rates[rate_index].y = y;
                    double old_rate = transition_rates[rate_index].rate;
                    if (solutes_grid[i][j] > 0 && solutes_grid[x][y] == 0){
                        EnergyResult result = get_energy_and_Eo(grains_grid, solutes_grid, boundary_info, energy_cache, solute_counts, grid_size_x, grid_size_y, i, j, x, y, q_selected, Jgg, Js,
                                            Jsg, Jss, Jsgs, F, Ed_oo, Eg_oo, Ecorr, PROC_SWAP);
                        double transition_probability = get_transition_probability(result.uij, result.Eo, T);
                        transition_rates[rate_index].rate = transition_probability;
                        final_sum += transition_probability;
                    }else{
                        transition_rates[rate_index].rate = 0;
                    }
                    if (tree) {
                        fenwick_add(tree, rate_index, transition_rates[rate_index].rate - old_rate);
                    }
                }else{
                    printf("Index out of bounds in swap rates array\n");
                }
            }
        // }
    }
    return final_sum - initial_sum;
}
