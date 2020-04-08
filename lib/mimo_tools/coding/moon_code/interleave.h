// interleave.h -- A random interleaver
// Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#ifndef INTERLEAVE_H
#define INTERLEAVE_H

class interleave {
   int size;
   int *pi;
   int *piinv;
public:
   interleave(int in_size, unsigned int seed=0);
   ~interleave() { delete[] pi;  delete[] piinv; };
   void Pi(double *in, double *out);
   void Pi(unsigned char *in, unsigned char *out);
   void Pi(double **in, double **out, int nrow);
   void Piinv(double *in, double *out);
   void Piinv(unsigned char *in, unsigned char *out);
   void Piinv(double **in, double **out, int nrow);
   void PiinvTimesoverlay(double **in, double **out, int nrow);
};


#endif

/*
Local Variables:
compile-command: "g++ -c interleave.cc"
End:
*/
