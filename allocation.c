// allocation.c
#include <stdlib.h>
#include <stdio.h>
#include "allocation.h"

/**
 * Allocates a 2D grid of integers with the specified dimensions.
 * Initializes all elements to 0.
 * 
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 * @return A pointer to the allocated grid, or NULL on failure.
 */
int** allocate_grid(int rows, int cols) {
    int** grid = malloc(rows * sizeof(int*));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for grid pointers.\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        grid[i] = calloc(cols, sizeof(int)); // Initialize all elements to 0
        if (grid[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for grid row %d.\n", i);
            // Free previously allocated rows
            while (i-- > 0) free(grid[i]);
            free(grid);
            return NULL;
        }
    }
    return grid;
}

/**
 * Frees a 2D grid of integers allocated with allocate_grid.
 * 
 * @param grid Pointer to the grid to free.
 * @param rows Number of rows in the grid.
 */
void free_grid(int** grid, int rows) {
    if (grid == NULL) return; // Handle null grid safely
    for (int i = 0; i < rows; i++) {
        free(grid[i]);
    }
    free(grid);
}

double** allocate_double_grid(int rows, int cols) {
    double** grid = malloc(rows * sizeof(double*));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for double grid pointers.\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        grid[i] = calloc(cols, sizeof(double));
        if (grid[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for double grid row %d.\n", i);
            while (i-- > 0) free(grid[i]);
            free(grid);
            return NULL;
        }
    }
    return grid;
}

void free_double_grid(double** grid, int rows) {
    if (grid == NULL) return;
    for (int i = 0; i < rows; i++) {
        free(grid[i]);
    }
    free(grid);
}

unsigned char** allocate_uchar_grid(int rows, int cols) {
    unsigned char** grid = malloc(rows * sizeof(unsigned char*));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for uchar grid pointers.\n");
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        grid[i] = calloc(cols, sizeof(unsigned char));
        if (grid[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for uchar grid row %d.\n", i);
            while (i-- > 0) free(grid[i]);
            free(grid);
            return NULL;
        }
    }
    return grid;
}

void free_uchar_grid(unsigned char** grid, int rows) {
    if (grid == NULL) return;
    for (int i = 0; i < rows; i++) {
        free(grid[i]);
    }
    free(grid);
}
