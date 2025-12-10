// energy_calculation.h
#ifndef ENERGY_CALCULATION_H
#define ENERGY_CALCULATION_H
#include "common_structs.h"

// Ensure the function prototype matches the definition
EnergyResult get_energy_and_Eo(int** grains_grid, int** solutes_grid, int** boundary_info, double **energy_cache, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y, int i, int j, int x, int y, 
                         int q_selected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F,
                         double Ed_oo, double Eg_oo, double Ecorr, ProcessType process);

double get_transition_probability(double uij, double Eo, double T);
double phi_formula(double n);
void update_energy_cache_for_sites(AffectedSite *sites, int num_sites, double **energy_cache, int **grains_grid, int **solutes_grid, int **boundary_info, unsigned char **solute_neighbor_counts, int grid_size_x, int grid_size_y,
                         int q_selected, double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F);

#endif // ENERGY_CALCULATION_H
