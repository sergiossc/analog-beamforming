// ******************************************************************
// instantiate.cc -- instantiate templatized classes using
// explicit instantiation declarations
// Created by Todd K. Moon, Electrical and Computer Engineering Dept.
// Utah State University.    Copyright January 1999
// *****************************************************************

// This file indicates explicitly to the compiler which templates are
// instantiated.  Unfortunately, this requires explicit information about the 
// class definitions; hence the .cc files are included here.

#include "polynomialT.cc"
#include "ModAr.h"

template class polynomialT<double>;
template class polytemp<double>;


template class polynomialT<ModAr>;
template class polytemp<ModAr>;

// template ostream& operator<< (ostream &os, const polynomialT<double> &p1);
