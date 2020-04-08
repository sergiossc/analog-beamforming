/**********************************************************************/
/* Generate two unit-variance Gaussian random numbers.
   This function calls rand(), so you can control the seed using
   void srand(unsigned int seed);

   Following Numerical Recipes
*/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <stdlib.h>
#include <math.h>

void gran2(double &r1, double &r2)
{
   double rsq, v1, v2, fac;

   do {
	  v1 = 2*(rand()/(double)RAND_MAX) - 1;
	  v2 = 2*(rand()/(double)RAND_MAX) - 1;
	  rsq = v1*v1 + v2*v2;
   } while(rsq > 1 || rsq == 0);
   fac = sqrt(-2*log(rsq)/rsq);
   r1 = v1*fac;
   r2 = v2*fac;
}
