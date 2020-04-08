// galtest2.cc -- test the Gallager (low-density parity-check) code
// decoder in the context of an AWGN
// Todd K. Moon

// Copyright 2004 by Todd K. Moon
// Permission is granted to use this program/data
// for educational/research only

#include "galdec.h"
#include <iostream>
using namespace std;
#include <math.h>

void randvec(double *y,int N, double sigma);
double gran(void);

main()
{
   int maxnumloop = 1000;  // maximum number of decoding iterations
   char *fname;
   // these SNR ranges should be set according to the rate of the code
   // Rate 1/2:  Approximately in the range 1-1.4 dB
   // Rate 1/3:  Approximately in the range 0.4-1 dB
   // maxcountarray tells how many blocks to check.  In general, as the
   // SNR increases, more blocks are necessary, since the probability of error
   // decreases and more blocks must 

   // Rate 1/2 data:
   fname = "A1-2.txt";
   double EbN0dBlist[] = {1,1.1,1.2,1.3,1.4};
                           // list of SNRs for plot values
   int maxcountarray[] = {20,20,20,1000,5000};
                           // number of blocks to process
                           // number of blocks to process
   // Rate 1/4 data:
   // fname = "A1-4.txt";
   // double EbN0dBlist[] = {.8,.9,1,1.1,1.2,1.3,1.4};
   // int maxcountarray[] = {20,20,20,20,100,1000,5000};

   // Rate 1/3 data:
   // fname = "A1-3.txt";
   // double EbN0dBlist[] = {0.4,0.5,0.6,0.7,0.8,0.9,1};
                           // list of SNRs for plot values
   // int maxcountarray[] = {20,20,20,20,100,2000,5000};

   int neb = sizeof(EbN0dBlist)/sizeof(double); // number of steps
   GALDEC galdec(fname,1);

   int i;
   int numloops;				// number of loops actually used
   int maxcount;				// maximum number of blocks to code
   double R = double(galdec.K)/double(galdec.N); // code rate
   cout << "Code rate=" << R << endl;
   double *y = new double[galdec.N]; // vector of channel outputs
   double *p = new double[galdec.N]; // vector of channel posterior probs

   double EbN0dB, EbN0, sigma2, sigma;
   int decval;					// return value from decoder
   int ctr;

   // variables to hold decoder statistics
   unsigned long int ncount, ndederrcount,nundederrcount,nblockcount;
   int decfailcount;			// count number of decoder failures
   int decsuccesscount;			// number of decoder successes
   unsigned long numdeciter;	// number of iterations of decoder
								// used to compute average number
   int maxnumdeciter;			// maximum number of iterations on decode

   double *ebnolist = new double[neb];   // list of EbN0 values
   unsigned int *nundedlist = new unsigned int[neb]; // list of
                                         // number of undedcoded samples
   unsigned int *ndedlist = new unsigned int[neb];	// list of
                                         // number of decoded samples
   double *errlist = new double[neb];    // error probability
   int *decfaillist = new int[neb];	     // decoder failures
   int *decsuccesslist = new int[neb];   // decoder successes
   double *numdeciteravglist = new double[neb]; // average number of decoder 
                                         // iterations
   int *maxnumdeciterlist = new int[neb];// maximum number of decoder iter's

   for(ctr=0; ctr < neb; ctr++) {
	  EbN0dB = EbN0dBlist[ctr];	// get SNR in DB
	  maxcount = maxcountarray[ctr];
	  EbN0 = pow(10.,EbN0dB/10.);  // get SNR
	  sigma2 = 1/(2*R*EbN0);       // compute noise variance
	  sigma = sqrt(sigma2);        // and standard deviation
	  cout << "EbN0dB=" << EbN0dB << "  EbN0=" << EbN0 << "  sigma2="
		   << sigma2 << endl;
	  ncount = 0;
	  ndederrcount = 0;
	  nundederrcount = 0;
	  nblockcount = 0;
	  decfailcount = 0;
	  decsuccesscount = 0;
	  numdeciter = 0;
	  maxnumdeciter = 0;
	  while(nblockcount < maxcount) { // do for specified number of blocks
         // The looping might be more efficient if looping continues until
         // a specified number of errors is counted.
		 nblockcount++;
		 randvec(y,galdec.N,sigma);  // Set up received vector, 
                                     // assume all 0 codeword transmitted
		 galdec.meastoprobAWGN(y,p,1,sigma2); // convert measurement to
		                                      // probability
		 decval = galdec.decode(p, 0, maxnumloop, numloops);  // decode
		 cout << "nblockcount=" << nblockcount << "  devcal=" << decval <<
			"  numloops=" << numloops << endl;
		 ncount += galdec.N;
		 if(!decval) {  // not decoded successfully - detected errors
			decfailcount++;
			for(i = 0; i < galdec.N; i++) {  // count the bits in error
			   if(galdec.c[i] != 0) {
				  ndederrcount++;
			   }
			}
		 }
		 else {  // check undected errors
			decsuccesscount++;
			numdeciter += numloops;
			if(numloops > maxnumdeciter)
			   maxnumdeciter = numloops;
			for(i = 0; i < galdec.N; i++) {  // count the bits in error
			   if(galdec.c[i] != 0) {
				  nundederrcount++;
			   }
			}
		 }
	  }
	  cout << "total="  << ncount << "  ndederr=" << ndederrcount <<
		 "  nundederr=" << nundederrcount << "  errorrate=" << 
			 double(ndederrcount+nundederrcount)/double(ncount) << endl;
	  cout << "avgdeciter=" << double(numdeciter)/double(decsuccesscount) <<
            "   maxdeciter=" << maxnumdeciter;
	  cout << "decsuccesscount=" << decsuccesscount << endl;
	  
      ebnolist[ctr] = EbN0;
	  nundedlist[ctr] = nundederrcount;
	  ndedlist[ctr] = ndederrcount;
	  errlist[ctr] = (double)(ndederrcount+nundederrcount)/
		 (double)ncount;
	  decfaillist[ctr] = decfailcount;
	  decsuccesslist[ctr] = decsuccesscount;
	  if(decsuccesscount)
		 numdeciteravglist[ctr] = (double)numdeciter/
			(double)decsuccesscount;
	  else
		 numdeciteravglist[ctr] = 0;
	  maxnumdeciterlist[ctr] = maxnumdeciter;
   }
   // print out summary information
   for(i = 0; i < ctr; i++) {
	  cout << EbN0dBlist[i] <<"\t" << ebnolist[i] << "\t" <<
		 nundedlist[i] << "\t" << ndedlist[i] << "\t" <<
		 errlist[i] << "\t" << decfaillist[i] << "\t" << 
		 decsuccesslist[i] << "\t" << numdeciteravglist[i] << 
		 maxnumdeciterlist[i] << endl;
			 
   }
}


void randvec(double *y,int N, double sigma)
// fill y with mean= -1 and variance sigma^2 Gaussian
// This corresponds to the all-zero vector
{
   int i;
   for(i = 0; i < N; i++) {
	  y[i] = -1 + sigma*gran();
   }
}

/*
Local Variables:
compile-command:"g++ -o galtest2 -g galtest2.cc  galdec.cc gran.cc"
End:
*/
