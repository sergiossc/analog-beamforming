int computeLm(int n,int k, int m)
// Lm = computeLm(n,k,m);
//
// Compute the maximum length of a list for an (n,k)
// code using GS(m) decoding

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

// Todd K. Moon, Feb. 12, 2004
{
   double v = k-1;
   double t;
   if(m==0) return 1;
   t = (v+2)/(2*v);
   int Lm = int(floor(sqrt( n*m*(m+1)/v + t*t) - t));
   return int(Lm);
}
