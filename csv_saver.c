#include <stdio.h>
#include "csv_saver.h"

void save_grid_to_csv(const char *output_dir, int **grid, int grid_size_x, int grid_size_y, int long long step, const char *grid_type) {
    char filename[300]; // Buffer for the filename
    snprintf(filename, sizeof(filename), "%s/%s_step_%lld.csv", output_dir, grid_type, step);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        fprintf(stderr, "Failed to open file: %s\n", filename);
        return;
    }


    for (int i = 0; i < grid_size_x; i++) { // Iterate over rows
        for (int j = 0; j < grid_size_y; j++) { // Iterate over columns
            fprintf(file, "%d", grid[i][j]); // Write the value to the file
            if (j < grid_size_y - 1) {
                fprintf(file, ","); // Add a comma after each value except the last one in a row
            }
        }
        fprintf(file, "\n"); // New line at the end of each row
    }
    fclose(file); // Close the file
}
