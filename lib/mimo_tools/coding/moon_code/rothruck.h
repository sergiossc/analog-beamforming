// Define the linked list which holds all the y-roots
// found by the roth-ruck algorithm

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef ROTHRUCK_H
#define ROTHRUCK_H

class rpolynode
{
public:
   rpolynode *next;
   polynomialT<TYPE> f;
   rpolynode() {
	  next = NULL;
	  f = TYPE(0);
   }
   ~rpolynode() {};
};


rpolynode *
rothruck(polynomialT<polynomialT<TYPE> >&Q,int D);

#endif
