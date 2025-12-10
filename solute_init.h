#ifndef SOLUTE_INIT_H
#define SOLUTE_INIT_H

#include <stdio.h>
#include "rng.h"

// Definition of initialization types
#define INIT_FULLY_BULK 1
#define INIT_RANDOM 2
#define INIT_FROM_FILE 3

// Function declaration
void initialize_solutes_grid(int **solutes_grid, int grid_size_x, int grid_size_y, int init_type, double density, const char *file_path, RngState *rng);

#endif
