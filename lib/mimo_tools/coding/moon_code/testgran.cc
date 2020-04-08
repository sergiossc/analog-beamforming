/**************************************
*
*  Program: testgran --- test the random number generators
*
*
*  Todd K. Moon
*
***************************************/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <iostream>
using namespace std;
#include "rands.h"

main()
{
   int i;
   double r,r1,r2;
   for(i = 0; i < 10; i++) {
	  r = gran();
	  gran2(r1,r2);
	  cout << "rand=" << uran() << " r=" << r << "r1,r2=" << r1 << " "
		   << r2 << endl;
   }
}	  

/*
Local Variables:
compile-command:"g++ -o testgran testgran.cc gran.c uran.c gran2.cc"
End:
*/
