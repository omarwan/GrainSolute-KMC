#ifndef READ_PARAMS_H
#define READ_PARAMS_H

// Declare the read_params function
void read_params(const char *filename, char *output_dir, int *grid_size_x, int *grid_size_y, int *init_type, int *q_no, int *q_selected, double *Jgg, double *Js, double *Jsg, double *Jss, double *Jsgs, double *F, double *density, int long long *n_step, int long long *e_step, double *T, double *Ed_oo, double *Eg_oo, double *Ecorr, char *grain_structure_file, char *solute_structure_file);

#endif