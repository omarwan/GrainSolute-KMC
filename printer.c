#include <stdio.h>
#include "printer.h"
/**
 * Prints a given matrix of integers.
 * 
 * @param matrix The matrix to be printed.
 * @param rows The number of rows in the matrix.
 * @param cols The number of columns in the matrix.
 */
void print_matrix(int **matrix, int grid_size_x, int grid_size_y) {
    for (int i = 0; i < grid_size_x; i++) {
        for (int j = 0; j < grid_size_y; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");  // New line after each row
    }
}
