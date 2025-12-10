#ifndef LOCAL_ENERGY_H
#define LOCAL_ENERGY_H

#include <stdio.h>
#include <math.h>

// Forward declarations
double phi_formula(double n);

double compute_local_energy(
    int i, int j,
    int** grains_grid,
    int** solutes_grid,
    int** boundary_info,
    int grid_size_x, int grid_size_y,
    int q_selected,
    double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F
);

#endif // LOCAL_ENERGY_H
