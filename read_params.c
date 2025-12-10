#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_params.h"

void read_params(const char *filename, char *output_dir, int *grid_size_x, int *grid_size_y, int *init_type, int *q_no, int *q_selected, double *Jgg, double *Js, double *Jsg, double *Jss, double *Jsgs, double *F, double *density, int long long *n_step, int long long *e_step, double *T, double *Ed_oo, double *Eg_oo, double *Ecorr, char *grain_structure_file, char *solute_structure_file) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening parameter file");
        exit(1);
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "grid_size_x=%d", grid_size_x) == 1) continue;
        if (sscanf(line, "grid_size_y=%d", grid_size_y) == 1) continue;
        if (sscanf(line, "init_type=%d", init_type) == 1) continue;
        if (sscanf(line, "q_no=%d", q_no) == 1) continue;
        if (sscanf(line, "q_selected=%d", q_selected) == 1) continue;
        if (sscanf(line, "Jgg=%lf", Jgg) == 1) continue;
        if (sscanf(line, "Js=%lf", Js) == 1) continue;
        if (sscanf(line, "Jsg=%lf", Jsg) == 1) continue;
        if (sscanf(line, "Jss=%lf", Jss) == 1) continue;
        if (sscanf(line, "Jsgs=%lf", Jsgs) == 1) continue;
        if (sscanf(line, "F=%lf", F) == 1) continue;
        if (sscanf(line, "density=%lf", density) == 1) continue;
        if (sscanf(line, "n_step=%lld", n_step) == 1) continue;
        if (sscanf(line, "e_step=%lld", e_step) == 1) continue;
        if (sscanf(line, "T=%lf", T) == 1) continue;
        if (sscanf(line, "Ed_oo=%lf", Ed_oo) == 1) continue;
        if (sscanf(line, "Eg_oo=%lf", Eg_oo) == 1) continue;
        if (sscanf(line, "Ecorr=%lf", Ecorr) == 1) continue;
        if (sscanf(line, "output_dir=%s", output_dir) == 1) continue;
        if (sscanf(line, "grain_structure_file=%s", grain_structure_file) == 1) continue;
        if (sscanf(line, "solute_structure_file=%s", solute_structure_file) == 1) continue;
    }
    printf("Params: GS_X=%d, GS_Y=%d, Init=%d, Q_No=%d, Q_Select=%d, Jgg=%.2f, Js=%.2f, Jsg=%.2f, Jss=%.2f, Jsgs=%.2f, F=%.2f, Density=%.2f, N_Step=%lld, E_Step=%lld, T=%.2f, Ed_oo=%.2f, Eg_oo=%.2f, Ecorr=%.2f, OutDir=%s, GrainFile=%s, SoluteFile=%s\n", 
    *grid_size_x, *grid_size_y, *init_type, *q_no, *q_selected, *Jgg, *Js, *Jsg, *Jss, *Jsgs, *F, *density, *n_step, *e_step, *T, *Ed_oo, *Eg_oo, *Ecorr, output_dir, grain_structure_file, solute_structure_file);

    fclose(file);
}
