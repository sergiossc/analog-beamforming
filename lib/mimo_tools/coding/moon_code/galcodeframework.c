/**************************************
*
*  Program:
*
*
*  Todd K. Moon
*
***************************************/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <stdio.h>
#include <stdlib.h>
#include <pmatlib.h>
#include <math.h>
#define MAXCOUNT2 100  // number of blocks to check for each EbN0

typedef int BIT32U;
typedef int *BIT32Up;
#define EPS 1.e-10
#define TINYDIV EPS
#define CLIPONE 1-EPS
#define CLIPCHECK(q0,q1,one) if(q0>one){q0=one;q1=1-one;} \
                             else if(q1>one){q1=one;q0=1-one;}

/* decoder data -- this is kept here as static, rather than in a structure, 
   to slightly speed up the running code */

static int *na;			  /* [N] # of elements above in this column */
static double *rt;				/* [maxrowwt] temporary row info. */
static double **r1;				/* [maxcolwt][N] r1[m][n] */
static double **r0;				/* [maxcolwt][N] r0[m][n] */
static double **q0, **q1;		/* [maxcolwt][N] q[m][n] */
static double **deltar;			/* [maxcolwt][N] deltar[m][n] */
static double *q0p;			/* the pseudopriors */
/* data used in the computations is stored so that the 
   column index (the second index) goes from 0 to N-1,
   and the row index (the first index) goes from 0 to maxcolwt-1
*/

typedef struct {
  BIT32U N;						/* block length */
  BIT32U K;						/* message length */
  BIT32U M;						/* redudandcy */
  BIT32U *Nmlen;				/* [M] lengths of row weight vectors */
  BIT32U **Nm;				    /* [M][rowwt] 
								   set of bits n that participate in check m*/
  BIT32U *Mnlen;				/* [N] lengths of column weight vectors */
  BIT32U **Mn;				    /* [N][colwt] set of checks in which bit n 
								   participates */ 
  int maxcolwt;					/* maximum weight of columns */
  int maxrowwt;					/* maximum weight of rows */
} Astruct;

/* sparse file format (using Mackay's format)
   N M  -- block length (int) and redundancy (int)
   maxcolweight maxrowweight
   colweights  (with N elements in it)
   rowweights  (with M elements in it)
   column data:  (repeated N times)
      row(1) row(2) ... row
   row data:  (repeated M times)
      col(1) col(2) ... col
*/

   
allocdecodedat(Astruct *A) 
{
   CALLOCVECTOR(na,int,A->N);
   CALLOCVECTOR(rt,double,A->maxrowwt);
   CALLOCMATRIX(r1, double,A->maxcolwt,A->N);
   CALLOCMATRIX(r0,double,A->maxcolwt,A->N);
   CALLOCMATRIX(q0,double,A->maxcolwt,A->N);
   CALLOCMATRIX(q1,double,A->maxcolwt,A->N);
   CALLOCMATRIX(deltar,double,A->maxcolwt,A->N);
   CALLOCVECTOR(q0p,double,A->N);
}

freedecodedat()
{
   free(na);
   free(rt);
   pfree_matrix((void ***)r1);
   pfree_matrix((void ***)r0);
   pfree_matrix((void ***)q0);
   pfree_matrix((void ***)q1);
   pfree_matrix((void ***)deltar);
   free(q0p);
}

int decode(Astruct *A, double *fn, int maxnumloop, int *numloops,
		   double *q1p, char *x, int printstuff)
/* A = sparse matrix structure
   fn = prior probabilities (one for each N)
  
   maxnumloop = maximum number of decoding iterations
   numloops = number of decoding iterations actually used
   q1p = pseudoprior probability (allocate space before calling this function)
   x = decoded value (allocated space before calling this function)
   returns: 1 if decoding works; 0 for a decoding failure
*/
{
   int i,l,k,m,n,row;
   double prod, prod0, prod1,alpha;
   int paritycheck = 0;
   char z;
   int loopcount = 0;
   double sum;
   

}


int fread_ivector(int *w, int lo, int hi, FILE *fp, int offset)
{
  int i, status = 0 ;
  
  for (i=lo;i<=hi;i++) {
    if ( fscanf(fp,"%d ",&w[i]) == EOF ) {
	   w[i] -= offset;
      status = -1 ; 
      break ;
    }
  }

/*  if ( status < 0 ) 
    fprintf( stderr, 
	    "Warning: fread_ivector failed at component %d\n",i);
*/
  return status ; 
}

int fread_imatrix(int **b,int l1,int h1,int l2,int h2,FILE *fp, int offset)
{
  int	i, j , status = 0;

  for (i=l1; i<=h1; i++){
    /*    fprintf(stderr,"'");  fflush(stderr); */
    for (j=l2; j<=h2; j++){
      if ( fscanf(fp,"%d ",&b[i][j]) == EOF ) {
		 status -- ;
		 break;
      }
	  b[i][j] -= offset;
    }
    if ( status < 0 ) break ;
  }
  if ( status)
		  fprintf(stderr,"Warning: readinimatrix failed at component %d\n",i);
  return status ;
}


 
readgal(char *fname,Astruct **A, int offset)
/* Since Mackay's stuff has base index=1, to make it work with my code,
   use offset=1 when reading Mackay's files.
   To read something with baseindex=0, use offset=0
*/
{
  FILE *fp;
  BIT32U M;
  int status = 0;

  *A = (Astruct *)calloc(sizeof(Astruct),1);
  fp = fopen( fname, "r" );
  if( !fp )   fprintf( stderr, "No such file: %s\n", fname ), exit(0);

  do {
    if ( fscanf(fp,"%d %d " , &((*A)->N) , &((*A)->M) ) == EOF ) {
      status = -1 ;       break ;
    }
    if ( fscanf(fp,"%d %d " , &((*A)->maxcolwt) , &((*A)->maxrowwt))
		 == EOF ) {
      status = -1 ;       break ;
    }
	CALLOCVECTOR((*A)->Mnlen,int,(*A)->N);
	CALLOCVECTOR((*A)->Nmlen,int,(*A)->M);
    status += fread_ivector ( (*A)->Mnlen , 0 , (*A)->N-1 , fp, offset ) ; 
    status += fread_ivector ( (*A)->Nmlen , 0 , (*A)->M-1 , fp, offset ) ;
	CALLOCMATRIX((*A)->Mn,int,(*A)->N,(*A)->maxcolwt);
	CALLOCMATRIX((*A)->Nm,int,(*A)->M,(*A)->maxrowwt);
    status += fread_imatrix((*A)->Mn,0,(*A)->N-1,0,(*A)->maxcolwt-1,fp,offset);
    status += fread_imatrix((*A)->Nm,0,(*A)->M-1,0,(*A)->maxrowwt-1,fp,offset);
    /*    fprintf(stderr,";");  fflush(stderr); */
  } while ( 0 ) ; 

  if ( status < 0 ) { fprintf(stderr,"\nreadgal %d\n",status); }
  (*A)->K = (*A)->N - (*A)->M;
  fclose(fp);
  return status ; 
}
  

void freeA(Astruct *A)
/* Free the space allocated in A */
{
  int i;
  for(i = 0; i < A->M; i++) {
	free(A->Nm[i]);
  }
  free(A->Nm);
  free(A->Nmlen);
  for(i = 0; i < A->N; i++) {
	free(A->Mn[i]);
  }
  free(A->Mn);
  free(A->Mnlen);
}

printsparse(double **d, Astruct *A)
{
   int l,n, m,lastn,n1;
   int *na;

   CALLOCVECTOR(na,int,A->N);
   for(m = 0; m < A->M; m++) {
	  lastn = 0;
	  for(l = 0; l < A->Nmlen[m]; l++) {
		 n = A->Nm[m][l];
		 /* print 0s as necessary */
		 for(n1 = lastn; n1 < n; n1++) {
			printf("%.2f\t",0.0);
		 }
		 /* print the new number */
		 printf("%.2f\t",d[na[n]++][n]);
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < A->N; n1++) {
		 printf("%.2f\t",0.0);
	  }
	  printf("\n");
   }
   free(na);

}

printsparseA(Astruct *A)
{
   int l,n, m,lastn,n1;
   int lastm,m1;

   for(m = 0; m < A->M; m++) {
	  lastn = 0;
	  for(l = 0; l < A->Nmlen[m]; l++) {
		 n = A->Nm[m][l];
		 /* print 0s as necessary */
		 for(n1 = lastn; n1 < n; n1++) {
			printf("0 ");
		 }
		 /* print the new number */
		 printf("1 ");
		 lastn = n+1;
	  }
	  for(n1 = lastn; n1 < A->N; n1++) {
		 printf("0 ");
	  }
	  printf("\n");
   }
}

void randvec(double *y,int N, double sigma)
/* fill y with mean-1 and variance sigma^2 Gaussian */
{
   int i;
   for(i = 0; i < N; i++) {
	  y[i] = 1 + sigma*gran();
   }
}

void compute_like(double *y, double *f, int N, double sigma2)
{
   int i;
   double s2;
   s2 = -2/sigma2;
   for(i = 0; i < N; i++) {
	  f[i] = 1/(1+exp(s2*y[i]));
   }
}
	  
main()
{
   Astruct *A;  /* the sparse matrix */
   double *fn; /* the prior probabilities */
   double *fy;
   int maxnumloop = 1000;  // maximum number of decoding iterations
   int decval, numloops;
   double EbN0startdB = 0.4;
   double EbN0stepdB = 0.1;
   double EbN0enddB = 0.4;
   double EbN0dB;  // SNR in dB
   double EbN0;					   
   double R;
   double sigma2, sigma;
   double *y;
   double *q1p;
   char *x;
   int i;
   unsigned long int ncount, ndederrcount,nundederrcount,nblockcount;
   int decfailcount;			/* count number of decoder failures */
   int decsuccesscount;			/* number of decoder successes */
   unsigned long numdeciter;	/* number of iterations of decoder */
								/* used to compute average number */
   int maxnumdeciter;			/* maximum number of iterations on 
								   decode */
   /* Save information for the end */
   double *ebnolist;
   double *ebnodBlist;
   unsigned long int *nundedlist;
   unsigned long int *ndedlist;
   double *errlist;
   int *decfaillist;
   int *decsuccesslist;
   double *numdeciteravglist;
   int *maxnumdeciterlist;
   int ctr,neb;

   readgal("15000.10000.3.0.1523",&A,1);
   allocdecodedat(A);
   R = (double)A->K/(double)A->N;
   printf("R=%g\n",R);
   CALLOCVECTOR(q1p,double,A->N);
   CALLOCVECTOR(x,char,A->N);

   CALLOCVECTOR(y,double,A->N);
   CALLOCVECTOR(fy,double,A->N);
   ctr = 0;
   
   neb = floor((EbN0enddB - EbN0startdB)/EbN0stepdB+1);
   CALLOCVECTOR(ebnodBlist,double,neb);
   CALLOCVECTOR(ebnolist,double,neb);
   CALLOCVECTOR(nundedlist,unsigned long int,neb);
   CALLOCVECTOR(ndedlist,unsigned long int,neb);
   CALLOCVECTOR(errlist,double,neb);
   CALLOCVECTOR(decfaillist,int,neb);
   CALLOCVECTOR(decsuccesslist,int,neb);
   CALLOCVECTOR(numdeciteravglist,double,neb);
   CALLOCVECTOR(maxnumdeciterlist,int,neb);

   for(EbN0dB=EbN0startdB; EbN0dB <= EbN0enddB; EbN0dB += EbN0stepdB) {
	  EbN0 = pow(10.,EbN0dB/10.);
	  // sigma2 = R/(2*EbN0);
	  sigma2 = 1/(2*R*EbN0);
	  sigma = sqrt(sigma2);
	  printf("EbN0dB=%g  EbN0=%g  sigma2=%g\n",EbN0dB,EbN0,sigma2);
	  ncount = 0;
	  ndederrcount = 0;
	  nundederrcount = 0;
	  nblockcount = 0;
	  decfailcount = 0;
	  decsuccesscount = 0;
	  numdeciter = 0;
	  maxnumdeciter = 0;
	  while(nblockcount < MAXCOUNT2) {
		 nblockcount++;
		 randvec(y,A->N,sigma);
		 compute_like(y,fy,A->N,sigma2);
		 decval = decode(A, fy, maxnumloop, &numloops,q1p,x,0);
		 printf("%ld decval=%d   numloops=%d\n",nblockcount,decval,
				numloops);
		 ncount += A->N;
		 if(!decval) {  /* not decoded successfully - detected errors*/
			decfailcount++;
			for(i = 0; i < A->N; i++) {  /* count the bits in error */
			   if(x[i] != 1) {
				  // printf("Decoding problem: %d\n",x[i]);
				  ndederrcount++;
			   }
			}
		 }
		 else {  /* check undected errors */
			decsuccesscount++;
			numdeciter += numloops;
			if(numloops > maxnumdeciter)
			   maxnumdeciter = numloops;
			for(i = 0; i < A->N; i++) {  /* count the bits in error */
			   if(x[i] != 1) {
				  // printf("Decoding problem: %d\n",x[i]);
				  nundederrcount++;
			   }
			}
		 }
	  }
	  printf("total: %ld  ndederr=%ld  nundederr=%ld errorrate=%g\n"
			 ,ncount,ndederrcount,nundederrcount,
			 (double)(ndederrcount+nundederrcount)/(double)ncount);
	  printf("avgdeciter=%g maxdeciter=%d ",
			 (double)numdeciter/(double)decsuccesscount,maxnumdeciter);
	  printf("decsuccesscount=%d\n",decsuccesscount);
	  
      ebnodBlist[ctr] = EbN0dB;
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
	  ctr++;
   }
   /* print out summary information */
   for(i = 0; i < ctr; i++) {
	  printf("%g\t%g\t%ld\t%ld\t%g\t%d\t%d\t%g\t%d\n",
			 ebnodBlist[i],
			 ebnolist[i],nundedlist[i],ndedlist[i],
			 errlist[i],decfaillist[i],decsuccesslist[i],
			 numdeciteravglist[i],maxnumdeciterlist[i]);
   }

}


/*
Local Variables:
compile-command:"gcc -o galcodeframework -g galcodeframework.c -lpmatlib -lm"
End:
*/
