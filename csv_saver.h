#ifndef CSV_SAVER
#define CSV_SAVER

#include <stdio.h> 

void save_grid_to_csv(const char *output_dir, int **grid, int grid_size_x, int grid_size_y, int long long step, const char *grid_type);

#endif 
