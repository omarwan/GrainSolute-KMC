// allocation.h
#ifndef ALLOCATION_H
#define ALLOCATION_H

// Function to allocate a 2D grid dynamically
int** allocate_grid(int rows, int cols);

// Function to free a 2D grid
void free_grid(int** grid, int rows);

double** allocate_double_grid(int rows, int cols);
void free_double_grid(double** grid, int rows);

unsigned char** allocate_uchar_grid(int rows, int cols);
void free_uchar_grid(unsigned char** grid, int rows);

#endif // ALLOCATION_H
