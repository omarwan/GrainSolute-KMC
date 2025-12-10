#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdint.h>
#include "common_structs.h"
#include "read_params.h"
#include "grain_init.h"
#include "solute_init.h"
#include "printer.h"
#include "csv_saver.h"
#include "boundary_sites.h"
#include "transition_rates.h"
#include "affected_sites.h"
#include "update_boundary.h"
#include "allocation.h"
#include "energy_calculation.h"
#include "fenwick.h"
#include "rng.h"
#include "periodic_neighbor.h"
#include "local_energy.h"

static void assert_tree_matches_total(const char *label, double total, const FenwickTree *tree) {
    if (tree->size <= 0) return;
    double tree_total = fenwick_total(tree);
    double tol = fmax(1e-9, 1e-6 * fmax(tree_total, total));
    assert(fabs(tree_total - total) <= tol && label);
}

void print_flip_rates(TransitionRate *flip_rates, int num_flip_rates) {
    printf("Flip Transition Rates:\n");
    for (int k = 0; k < num_flip_rates; k++) {
        printf("index %d Site (%d, %d) (%d, %d): Rate = %.20f\n", k, flip_rates[k].i, flip_rates[k].j,flip_rates[k].y, flip_rates[k].x, flip_rates[k].rate);
    }
}

void save_simulated_time(const char *output_dir, double simulated_time, int flip_itr, double simulated_time_flip, int swap_itr, double simulated_time_swap, int step) {
    char filename[200];
    snprintf(filename, sizeof(filename), "%s/simulated_time.txt",output_dir);

    struct stat st;
    int file_exists = (stat(filename, &st) == 0) && (st.st_size > 0);

    FILE *file = fopen(filename, "a"); 
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    if (!file_exists) {
        fprintf(file, "Step, Flip Iteration, Flip Time, Swap Iteration, Swap Time, Simulated Time\n");
    }

    fprintf(file, "%d, %d, %f, %d, %f, %f\n", step, flip_itr, simulated_time_flip, swap_itr, simulated_time_swap, simulated_time);
    fclose(file);
}

static uint64_t get_seed_from_env(void) {
    const char *env = getenv("SEED");
    if (!env || !*env) return 0;
    char *endptr = NULL;
    uint64_t val = strtoull(env, &endptr, 10);
    if (endptr == env) return 0; // failed to parse
    return val;
}


////////////////////////////////////////
// MAIN
////////////////////////////////////////
int main(int argc, char *argv[]) {
    (void)argc; // unused
    clock_t start_time = clock(); // Start the timer
    char params_file[200];
    snprintf(params_file, sizeof(params_file), "%s/params.txt", argv[1]);

    int grid_size_x, grid_size_y, init_type, q_no, q_selected ;
    int long long step, n_step, e_step;
    double T, Jgg, Js, density, Jss, Jsgs, Jsg, F, Ed_oo, Eg_oo, Ecorr;
    char output_dir[200], grain_structure_file[200], solute_structure_file[200];
    
    read_params(params_file, output_dir, &grid_size_x, &grid_size_y, &init_type, &q_no, &q_selected, &Jgg, &Js, &Jsg, &Jss, &Jsgs, &F, &density, &n_step, &e_step, &T, &Ed_oo, &Eg_oo, &Ecorr, grain_structure_file, solute_structure_file);
    init_periodic_neighbor_cache(grid_size_x, grid_size_y);
    
    if (q_selected > q_no || q_selected < 1) {
        printf("Error: q_selected = %d is out of bounds. It should be between 1 and %d.\n", q_selected, q_no);
        exit(EXIT_FAILURE); // Terminates the program
    }

    uint64_t seed = get_seed_from_env();
    if (seed == 0) {
        seed = 12345ULL; // deterministic default; override with SEED env var
    }
    RngState rng;
    rng_seed(&rng, seed);


    // Dynamically allocate grids
    int** grains_grid = allocate_grid(grid_size_x, grid_size_y);
    int** solutes_grid = allocate_grid(grid_size_x, grid_size_y);
    int** boundary_info = allocate_grid(grid_size_x, grid_size_y);
    double** energy_cache = allocate_double_grid(grid_size_x, grid_size_y);
    unsigned char** solute_neighbor_counts = allocate_uchar_grid(grid_size_x, grid_size_y);
    if (!grains_grid || !solutes_grid || !boundary_info || !energy_cache || !solute_neighbor_counts) {
        fprintf(stderr, "Allocation failed\n");
        return 1;
    }

    initialize_grains_grid(grains_grid, grid_size_x, grid_size_y, init_type, q_no, grain_structure_file, &rng);
    initialize_solutes_grid(solutes_grid, grid_size_x, grid_size_y, init_type, density, solute_structure_file, &rng);
    find_boundary_sites_and_neighbors(boundary_info, grains_grid, grid_size_x, grid_size_y);
    // build initial energy cache and solute neighbor counts
    for (int ii = 0; ii < grid_size_x; ii++) {
        for (int jj = 0; jj < grid_size_y; jj++) {
            Neighbors n = periodic_neighbor(ii, jj, grid_size_x, grid_size_y);
            Point neighbors[4] = {n.left_neighbor, n.right_neighbor, n.top_neighbor, n.bottom_neighbor};
            unsigned char cnt = 0;
            for (int k = 0; k < 4; k++) {
                if (solutes_grid[neighbors[k].x][neighbors[k].y] > 0) cnt++;
            }
            solute_neighbor_counts[ii][jj] = cnt;
            energy_cache[ii][jj] = compute_local_energy(
                ii, jj, grains_grid, solutes_grid, boundary_info,
                grid_size_x, grid_size_y, q_selected, Jgg, Js, Jsg, Jss, Jsgs, F);
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Initiate the transition rates arrays
    //////////////////////////////////////////////////////////////////////
    int num_flip_rates = grid_size_x * grid_size_y * 6;
    TransitionRate *flip_rates = get_rates(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, Jgg, Js, Jsg, Jss, Jsgs, F, T, Ed_oo, Eg_oo, Ecorr, PROC_FLIP, q_no, q_selected);
    int num_swap_rates = grid_size_x * grid_size_y * 4;
    TransitionRate *swap_rates = get_rates(grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y, Jgg, Js, Jsg, Jss, Jsgs, F, T, Ed_oo, Eg_oo, Ecorr, PROC_SWAP, q_no, q_selected);
    FenwickTree flip_tree = {0};
    FenwickTree swap_tree = {0};
    if (fenwick_init(&flip_tree, num_flip_rates) != 0 || fenwick_init(&swap_tree, num_swap_rates) != 0) {
        fprintf(stderr, "Failed to allocate fenwick trees\n");
        return 1;
    }
    
    double total_flip_rate = 0.0;
    for (int r = 0; r < num_flip_rates; r++) {
        total_flip_rate += flip_rates[r].rate;
        fenwick_add(&flip_tree, r, flip_rates[r].rate);
    }
    double total_swap_rate = 0.0;
    for (int r = 0; r < num_swap_rates; r++) {
        total_swap_rate += swap_rates[r].rate;
        fenwick_add(&swap_tree, r, swap_rates[r].rate);
    }
    assert_tree_matches_total("init flip fenwick", total_flip_rate, &flip_tree);
    assert_tree_matches_total("init swap fenwick", total_swap_rate, &swap_tree);
    //////////////////////////////////////////////////////////////////////
    double simulated_time_flip = 0.0;
    double simulated_time_swap = 0.0;
    double simulated_time = 0.0;
    int flip_itr = 0;
    int swap_itr = 0;

    for (step = 0; step < n_step + 1; step++) {
        if (total_flip_rate <= 0.0) {
            fprintf(stderr, "Error: total_flip_rate must be equal or greater than zero.\n");
            exit(EXIT_FAILURE);
        }
        if (total_swap_rate <= 0.0) {
            fprintf(stderr, "Error: total_swap_rate must be equal or greater than zero.\n");
            exit(EXIT_FAILURE);
        }

        double tau_flip = rng_next_exponential(&rng, total_flip_rate);
        double tau_swap = rng_next_exponential(&rng, total_swap_rate);

        double r_flip = rng_next_uniform(&rng) * total_flip_rate;
        double r_swap = rng_next_uniform(&rng) * total_swap_rate;
        if (r_flip == 0.0) { r_flip = 1e-12; }
        if (r_swap == 0.0) { r_swap = 1e-12; }
        ///////////////// initiate an array for affected sites //////////////////////
        int max_affected = 49;
        AffectedSite affected_sites[max_affected];
        
        int i, j, x, y, flip_chosen_idx, swap_chosen_idx;
        ///////////////// flip //////////////////////
        int energy_sites_count = 0;
        AffectedSite energy_sites[50]; // big enough for 7x7 and small swap lists
        int is_flip = (tau_flip < tau_swap);
        if (is_flip){
            flip_itr += 1;
            simulated_time_flip += tau_flip;
            ////////////////// find the site to flip /////////////////////
            flip_chosen_idx = fenwick_find_prefix(&flip_tree, r_flip);
            ////////////////// flip the site & update boundary matrix /////////////////////
            i = flip_rates[flip_chosen_idx].i;
            j = flip_rates[flip_chosen_idx].j;
            x = flip_rates[flip_chosen_idx].x;
            y = flip_rates[flip_chosen_idx].y;

            if (x == -1 && y == -1) {
                int x1 = flip_rates[flip_chosen_idx-2].x;
                int y1 = flip_rates[flip_chosen_idx-2].y;
                int x2 = flip_rates[flip_chosen_idx-3].x;
                int y2 = flip_rates[flip_chosen_idx-3].y;
                int x3 = flip_rates[flip_chosen_idx-4].x;
                int y3 = flip_rates[flip_chosen_idx-4].y;
                int x4 = flip_rates[flip_chosen_idx-5].x;
                int y4 = flip_rates[flip_chosen_idx-5].y;

                int old_qs[] = {grains_grid[x1][y1], grains_grid[x2][y2], grains_grid[x3][y3], grains_grid[x4][y4], q_selected};
                int num_old_qs = 5; // Total number of unique old grain types
                int new_q;
                int is_valid;

                do {
                    new_q = (int)(rng_next_u32(&rng) % q_no) + 1; // Assuming q_no is the upper limit for grain types
                    is_valid = 1;

                    for (int k = 0; k < num_old_qs; k++) {
                        if (new_q == old_qs[k]) {
                            is_valid = 0;
                            break;
                        }
                    }
                } while (!is_valid);
                grains_grid[i][j] = new_q; // Update the grid with the new grain type
            }
            else if (x == -2 && y == -2){grains_grid[i][j] = q_selected;}
            else{grains_grid[i][j] = grains_grid[x][y];}
            update_boundary_site(grains_grid, boundary_info, grid_size_x, grid_size_y, i, j);
        }

        else{
            swap_itr += 1;
            simulated_time_swap += tau_swap;
            ////////////////// find the site to swap /////////////////////
            swap_chosen_idx = fenwick_find_prefix(&swap_tree, r_swap);
            ////////////////// swap solutes at the sites  /////////////////////
            i = swap_rates[swap_chosen_idx].i;
            j = swap_rates[swap_chosen_idx].j;
            x = swap_rates[swap_chosen_idx].x;
            y = swap_rates[swap_chosen_idx].y;
            
            int temp = solutes_grid[i][j];
            solutes_grid[i][j] = solutes_grid[x][y];
            solutes_grid[x][y] = temp;
            energy_sites_count = get_swap_energy_sites(energy_sites, i, j, x, y, grid_size_x, grid_size_y);
        }
        double rand_time = rng_next_uniform(&rng);
        if (rand_time <= 1e-12) { rand_time = 1e-12; }
        double tau_tot = -log(rand_time) / (total_swap_rate + total_flip_rate);
        simulated_time += tau_tot;
        // printf("%d, %d)\n", i, j);
        ////////////// Find the affected sites /////////////////////
        get_affected_sites(affected_sites, i, j, grid_size_x, grid_size_y);
        // update energy cache for affected sites before using it in rate updates
        if (is_flip) {
            get_affected_sites_for_energy(energy_sites, i, j, grid_size_x, grid_size_y);
            energy_sites_count = 13;
        } else if (energy_sites_count == 0) {
            energy_sites_count = get_swap_energy_sites(energy_sites, i, j, x, y, grid_size_x, grid_size_y);
        }
        update_energy_cache_for_sites(energy_sites, energy_sites_count, energy_cache, grains_grid, solutes_grid, boundary_info, solute_neighbor_counts, grid_size_x, grid_size_y,
                                      q_selected, Jgg, Js, Jsg, Jss, Jsgs, F);
        //////////////////// Update the rates /////////////////////
        double total_flip_rate_diff = update_flip_rates_for_affected_sites(flip_rates, num_flip_rates, grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y,
                                            affected_sites, max_affected, Jgg, Js, Jsg, Jss, Jsgs, F, T, Ed_oo, Eg_oo, Ecorr, q_no, q_selected, &flip_tree);
        total_flip_rate += total_flip_rate_diff;

        double total_swap_rate_diff = update_swap_rates_for_affected_sites(swap_rates, num_swap_rates, grains_grid, solutes_grid, boundary_info, energy_cache, solute_neighbor_counts, grid_size_x, grid_size_y,
                                            affected_sites, max_affected, Jgg, Js, Jsg, Jss, Jsgs, F, T, Ed_oo, Eg_oo, Ecorr, q_selected, &swap_tree);
        total_swap_rate += total_swap_rate_diff;
        assert_tree_matches_total("loop flip fenwick", total_flip_rate, &flip_tree);
        assert_tree_matches_total("loop swap fenwick", total_swap_rate, &swap_tree);
        
        if (step % e_step == 0){
            save_grid_to_csv(output_dir, grains_grid, grid_size_x, grid_size_y, step, "grains_grid");
            save_grid_to_csv(output_dir, solutes_grid, grid_size_x, grid_size_y, step, "solutes_grid");
            save_simulated_time(output_dir, simulated_time, flip_itr, simulated_time_flip, swap_itr, simulated_time_swap, step);
        }

    }
    free_grid(grains_grid, grid_size_x);
    free_grid(solutes_grid, grid_size_x);
    free_grid(boundary_info, grid_size_x);
    free_double_grid(energy_cache, grid_size_x);
    free_uchar_grid(solute_neighbor_counts, grid_size_x);
    free(flip_rates);
    free(swap_rates);
    fenwick_free(&flip_tree);
    fenwick_free(&swap_tree);
    free_periodic_neighbor_cache();
    clock_t end = clock();
    double elapsed_time = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("%f ", elapsed_time);
    return 0;
}
