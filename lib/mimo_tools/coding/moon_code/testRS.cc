//  Program: testRS.cc --- test the RS decoder
//  Todd K. Moon
// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "RSdec.h"
#include "polynomialT.cc"

double uran(void);

// create instantiations of the polynomial class of type GFNUM2m
template class polynomialT<GFNUM2m>;
template class polytemp<GFNUM2m>;

main()
{
   int i,j,l;
   int j1,j2,j3;
   int j0=1;                    // starting index

   GFNUM2m::initgf(8,0x11D);  //   1 0001 1101   x^8 + x^4 + x^3 + x^2 + 1
   int n = 255;        // codeword length
   int k = 249;        // message length
   int t = 3;          // error correction capability
   RSdec decoderb(t,n);// build a decoder object

   GFNUM2m mb[255];    // codeword array
   GFNUM2m decb[255];  // decode array
   int loc[3];         // random locations to flip
   int val[3];         // random values to flip to
   int nerror;
   for(i = 0; i < 100000; i++) {
      cout << "i=" << i << endl;

      nerror = int(3*uran()) + 1;  // choose a random number of errors
      for(j = 0; j < nerror; j++) {
         loc[j] = int(n*uran());   // choose a random error location
         val[j] = int(256*uran()); // and a random error value
         mb[loc[j]] = val[j];      // set error
      }

      decoderb.decode(mb,decb);    // call the decoder
      for(j = 0; j < n; j++) {     // check the decoder
         if(decb[j] != 0) {
            cout << "Undecoded error: i=" << i << " " << loc[0] << " " <<
               loc[1] << " " << loc[2] << " " << val[0] << " " << val[1] <<
               " " << val[2] << endl;
            break;
         }
      }
      for(j = 0; j < nerror; j++) {
         mb[loc[j]] = 0;
      }
   }
}

/*
Local Variables:
compile-command: "g++ -o testRS -g testRS.cc GFNUM2m.cc RSdec.cc ChienSearch.cc uran.cc"
End:
*/
