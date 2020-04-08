
// utiltkm.cc --
// some simple useful utilities 
// uran --- generates U(0,1) r.v.
// gran -- generators N(0,1) r.v.s, or pairs of them
// sort2lfd -- sort a double array, and a corresponding integer (index) array
// sort1lf -- sort a double array
// sort2fd -- sort a float array, and a corresponding integer (index) array
// sort1lf -- sort a float array

// All sorts are into increasing order

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef UTILTKM_H
#define UTILTKM_H
double gran();
void gran(double &r1, double &r2);
double uran(void);

void sort2lfd(int num_elements, double *array, int *array2);
void sort1lf(int num_elements, double *array);
void sort2fd(int num_elements, float *array, int *array2);
void sort1f(int num_elements, float *array);

#endif
