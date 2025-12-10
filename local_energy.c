#include "local_energy.h"
#include "periodic_neighbor.h"

double phi_formula(double n) {
    return 1.0 - (0.25 * (n - 2.0) * (n - 2.0));
}

double compute_local_energy(
    int i, int j,
    int** grains_grid,
    int** solutes_grid,
    int** boundary_info,
    int grid_size_x, int grid_size_y,
    int q_selected,
    double Jgg, double Js, double Jsg, double Jss, double Jsgs, double F
) {
    double Egg = 0, Ess = 0, Esg = 0, Force = 0;
    Neighbors n = periodic_neighbor(i, j, grid_size_x, grid_size_y);
    Point neighbors[4] = {n.left_neighbor, n.right_neighbor, n.top_neighbor, n.bottom_neighbor};

    // Check if the flipped cell has the special q_selected
    if (grains_grid[i][j] == q_selected) {
        Force = F;
    }

    // 1) Grain boundary energy
    Egg = Jgg / 2.0 * boundary_info[i][j];

    // 2) Solute–grain boundary energy
    if (solutes_grid[i][j] > 0) {
        double n_solgb = boundary_info[i][j];
        Esg = Js + Jsg * phi_formula(n_solgb);
    }

    // 3) Solute–solute energy (sum over all 4 neighbors)
    for (int k = 0; k < 4; k++) {
        int nx = neighbors[k].x;
        int ny = neighbors[k].y;
        if (solutes_grid[i][j] > 0 && solutes_grid[nx][ny] > 0) {
            double n_avg = 0.5 * (boundary_info[i][j] + boundary_info[nx][ny]);
            Ess += (Jss / 2.0) + (Jsgs / 2.0) * phi_formula(n_avg);
        }
    }
    // printf("Egg %f Esg %f Ess %f force %f at n = %d \n ", Egg, Esg, Ess, Force, boundary_info[i][j]);

    return Egg + Esg + Ess + Force;
}
