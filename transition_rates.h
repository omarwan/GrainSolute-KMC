#ifndef TRANSITION_RATES_H
#define TRANSITION_RATES_H
#include "common_structs.h"
#include "fenwick.h"

TransitionRate* get_rates(int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y,
                        double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                        double T, double Ed_oo, double Eg_oo, double Ecorr, ProcessType process, int q_no, int q_selected);

double update_flip_rates_for_affected_sites(TransitionRate *transition_rates, int num_rates, int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_counts, int grid_size_x, int grid_size_y,
                                     AffectedSite *affected_sites, int num_affected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                                     double T, double Ed_oo, double Eg_oo, double Ecorr, int q_no, int q_selected, FenwickTree *tree);

double update_swap_rates_for_affected_sites(TransitionRate *transition_rates, int num_rates, int **grains_grid, int **solutes_grid, int **boundary_info, double **energy_cache, unsigned char **solute_counts, int grid_size_x, int grid_size_y,
                                     AffectedSite *affected_sites, int num_affected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                                     double T, double Ed_oo, double Eg_oo, double Ecorr, int q_selected, FenwickTree *tree);


#endif // TRANSITION_RATES_H
