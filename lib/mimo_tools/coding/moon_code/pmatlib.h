/******************************
*
* matrix library using pointers 
*
* This library implements matrix arithmetic
*
* Storage is based on the array of pointers idea for
* flexible (and fast) indexing.
* 
* Todd K. Moon
* Utah State University
* Started: 4 September 1991
*********************************/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#ifndef _PMATLIB_H_
#define  _PMATLIB_H_

#include <stdio.h>      /* to define NULL */

#define THRESH 10e-30   /* An error threshold.  How close can two numbers be */
                        /* and still not be equal, > THRESH.*/

#define ONESVEC(mat,type,num) {int ionesvec; \
          mat = (type *)calloc(num,sizeof(type)); \
          if(mat==NULL) { \
           fprintf(stderr,"Cannot allocate space for vector %s\n",#mat);\
             exit(-1); } \
               for(ionesvec = 0; ionesvec < num; ionesvec++) \
               mat[ionesvec] = 1; }

#define ONESMAT(mat,type,m,n) {int ionesmat,jonesmat; \
          mat = (type **)pcalloc_matrix2(m,n,sizeof(type)); \
          if(mat==NULL) { \
           fprintf(stderr,"Cannot allocate space for matrix %s\n",#mat);\
             exit(-1); } \
               for(ionesmat = 0; ionesmat < m; ionesmat++) \
                  for(jonesmat = 0; jonesmat < n; jonesmat++) \
                      mat[ionesmat][jonesmat] = 1;}

/* allocate a vector of arbitrary type (even typdef'ed types) */
#define VECTOR(name,type,num) type *name; CALLOCVECTOR(name,type,num);

#define MATRIX(name,type,n1,n2) type **name; CALLOCMATRIX(name,type,n1,n2);
#define TENSOR(nam,type,n1,n2,n3)type ***nam;CALLOCTENSOR(nam,type,n1,n2,n3);

#define CALLOCVECTOR(name,type,num) {name = \
           (type *)pcalloc_vector2(num,sizeof(type)); \
           if(name == NULL) {  \
             fprintf(stderr,"Cannot allocate space for vector %s\n",#name);\
             exit(-1); }}

/* allocate a matrix of arbitrary type (even typdef'ed types) */
#define CALLOCMATRIX(name,type,rows,cols) {name = \
         (type **)pcalloc_matrix2(rows,cols,sizeof(type)); \
         if(name == NULL) \
           { fprintf(stderr,"Cannot allocate space for matrix %s\n",#name);\
           exit(-1); }}

/* allocate a tensor of arbitrary type (even typdef'ed types) */
#define CALLOCTENSOR(name,type,nmat,rows,cols) {name = \
         (type ***)pcalloc_tensor2(nmat,rows,cols,sizeof(type)); \
         if(name == NULL) \
           { fprintf(stderr,"Cannot allocate space for tensor %s\n",#name);\
           exit(-1); }}


/* declare a matrix of size (nr x nc) */
void **pcalloc_matrix(int type,int nr,int nc,char *name);

/* allocate a matrix of size nr x nc, with elements of given size */
/* used by macro CALLOCMATRIX */
void **pcalloc_matrix2(int length,int width, int sizeoftype);

/* declare a vector of size (length) */
void *pcalloc_vector(int type,int length,char *name);

/* allocate a matrix of size length with elements of a given size */
/* used by macro CALLOCVECTOR */
void *pcalloc_vector2(int length, int sizeoftype);

/* allocate a tensor of size nmat x nr x nc, with elements of given size */
/* used by macro CALLOCMATRIX */
void ***pcalloc_tensor2(int nmat, int length,int width, int sizeoftype);

/* free a matrix of size (nr x nc) */
void pfree_matrix(void ***array);

/* free an array of size (length) */
void pfree_vector(void **vector);

/* alloc a nstack x nrow x ncol tensor */
void ***pcalloc_tensor(int type,int nstack,int nrow,int ncol,char *s);
void pfree_tensor(void ****m,int nstack);


/* set a vector to 0 */
void pveczerod(int *v,int dim);
void pveczerof(float *v,int dim);
void pveczerolf(double *v,int dim);

/* set a matrix to 0 */
void pmatzerof(float **m, int r, int c);
void pmatzerolf(double **m, int r, int c);

/* set vector to a constant */
void pvecconstf(float *v, float c, int dim);
void pvecconstlf(double *v, double c, int dim);

/* set matrix to a constant */
void pmatconstf(float **mat, float c, int m, int n);
void pmatconstlf(double **mat, double c, int m, int n);


/* Add two vectors: c= a + b */
void pvecaddf(float *a, float *b, float *c, int n);
void pvecaddlf(double *a, double *b, double *c, int n);

/* Subtract two vectors:  c = a-b;*/
void pvecsubf(float *a, float *b, float *c, int n);
void pvecsublf(double *a, double *b, double *c, int n);

/* Multiply a vector times a scalar: c = b*a */
void pvecscalef(float *a,float b, float *c, int n);
void pvecscalelf(double *a, double b, double *c, int n);

/* Inner product: c = a'b */
float pvecinnerf(float *a, float *b, int n);
double pvecinnerlf(double *a, double *b, int n);

/* Multiply two vectors element by element */
void pvecvecmultf(float *a, float *b, float *c, int n);
void pvecvecmultlf(double *a, double *b, double *c, int n);

/* Copy vectors: a --> b */
void pveccopylf(double *a,double *b,int m);
void pveccopyf(float *a,float *b,int m);
void pveccopyd(int *a,int *b,int m);
void pveccopyc(char *a,char *b,int m);


/* Vector norms: */
double pvecnorm2lf(double *, int);
float pvecnorm2f(float *, int);

/********************************************************************/
/* convolve the sequence in1 (of n1 points) with the sequence in2
   (of n2 points) into out
   This does not do "fast" convolution (e.g., using FFTs).
   Return value is the length of the new sequence, n1+n2-1
*/
int pconvlf(double *in1, int n1, double *in2, int n2, double *out);
/********************************************************************/
/* return index of min (resp. max) of vectors */
int pminveclf(double *v, int n);
int pmaxveclf(double *v, int n);
int pminvecf(float *v, int n);
int pminvecf(float *v, int n);

/********************************************************************/

/* Add matrices of size (m x n): c = a + b */
void pmataddf(float **a,float **b,float ** c,int m,int n);
void pmataddlf(double **a,double **b,double ** c,int m,int n);

/* subtract matrices of size (m x n): c = a - b */
void pmatsubf(float **a,float **b,float ** c,int m,int n);
void pmatsublf(double **a,double **b,double ** c,int m,int n);

/* Multiply matrices element by element */
void pmatmatmultf(float **a,float **b,float ** c,int m,int n);
void pmatmatmultlf(double **a,double **b,double ** c,int m,int n);

/* Add and scale: c = d*(a+b), where d is a scalar constant */
void pmataddandscalelf(double **a,double **b,double d,double ** c,int m,int n);
void pmataddandscalef(float **a,float **b,float d,float ** c,int m,int n);

/********************************************************************/
/* Copy matrix: a --> b */
void pmatcopylf(double **a, double **b, int m, int n);
void pmatcopyf(float **a, float **b, int m, int n);

/* Multiply a matrix times a scalar: c = b*a, a is nxm */
void pmatscalef(float **a, float b, float **c, int n, int m);
void pmatscalelf(double **a, double b, double **c, int n, int m);

/* transpose a matrix (m x n) into b matrix, assumed to be declared as */
/* (n x m) */
void pmattposef(float **a, float **b, int m, int n);
void pmattpoself(double **a, double **b, int m, int n);


/* Multiply two matrices: c= a*b   (no transposes are assumed, so 
   it avoids the overhead of the switch)
*/
void pmatmultlf(double **a, double **b, double **c, int n, int m, int q);
void pmatmultf(float **a, float **b, float **c, int n, int m, int q);

/* Matrix multiplication with transposes: 
   Compute one of the following forms,
   depending on the value of type:
   type 0: c = a*b
   a is n x m
   b is m x q
   c is n x q
        
   type 1: c = a*b'
   a is n x m
   b is q x m
   c is n x q
        
   type 2: c = a'*b
   a is m x n
   b is m x q
   c is n x q
        
   type 3: c = a'*b'
   a is m x n
   b is q x m
   c is n x q
*/
void pmatmultTf(float **a,float **b,float **c,int n,int m,int q,int type);
void pmatmultTlf(double **a,double **b,double **c,int n,int m,int q,int type);


/* Multiply matrix x vector: a is nxm, b is mx1, c is nx1 */
void pmatvecmultf(float **a, float *b,float *c, int n, int m);
void pmatvecmultlf(double **a, double *b,double *c, int n, int m);

/* Multiply a*b', where a is mx1 and b is nx1, 
   to obtain the mxn outer product c */
void pvecouterlf(double *a, double *b, double **c, int m, int n);
void pvecouterf(float *a, float *b, float **c, int m, int n);

/* Multiply a vector (transpose) times a matrix: c = a'*b */
/* a is (nx1), b is (nxm), c is (1xm) (a vector) */
void pvecmatmultlf(double *a,double **b, double *c,int n,int m);
void pvecmatmultf(float *a,float **b, float *c,int n,int m);

/* compute v'Mw, and return the result. Generalized vector inner product */
double pvecinnerMlf(double *v, double **M, double *w,int n);
float pvecinnerMf(float *v, float **M, float *w, int n);


/* Multiply a vector (transpose) times a matrix transpose: c = a'*b' */
/* a is (nx1), b is (mxn), c is (1xm) (a vector) */
void pvecmatTmultlf(double *a,double **b, double *c,int n,int m);
void pvecmatTmultf(float *a,float **b, float *c,int n,int m);

/* Multiply a matrix transpose times a vector: c = A'*b, 
   A is (nxm), b is(nx1), c is (mx1) */
void pmatTvecmultf(float **A, float *b,float *c, int n, int m);
void pmatTvecmultlf(double **A, double *b,double *c, int n, int m);

/* form the product c = a * diag(d), where a is mxn and d is nx1 */
void pmatdiagmultf(float **a, float *d, float **c, int m, int n);
void pmatdiagmultlf(double **a, double *d, double **c, int m, int n);


/* Find the trace of a matrix. */
float pmattracef(float **A, int M);
double pmattracelf(double **A, int M);


/* compute the fobenius norm of a matrix */
double pfrobeniuslf( double **a, int rows, int cols );
float pfrobeniusf( float **a, int rows, int cols );

/* Copy tensors (space already allocated) a ---> b*/
void ptencopylf(double ***a, double ***b, int nstack, int r, int c);
void ptencopyf(float ***a, float ***b, int nstack, int r, int c);


/* Printing routines: */
void pmatprint(char *name,void **mat, int m, int n,int type);
void pmatprints(char *name,char **mat,int n);
void pmatprintsc(char *name,signed char **mat,int n,int m);
void pmatprintuc(char *name,unsigned char **mat,int n,int m);
void pmatprinthd(char *name,short int **mat, int m, int n);
void pmatprintd(char *name,int **mat, int m, int n);
void pmatprintx(char *name,int **mat,int n,int m);
void pmatprintf(char *name,float **mat, int m, int n);
void pmatprintlf(char *name,double **mat, int m, int n);
void pmatprintld(char *name,long int **mat, int m, int n);

/* control if a matrix is printed in a row or in matrix form */
/* default is to print in matrix form */
void setprintinrow();
void clearprintinrow();
/* set the format of integer printing for printf.  E.g., %4.3d is obtained with
    setintprintwidths(4,3);
   (default is 0.0)
*/
void setintprintwidths(int w1,int w2);
/* set the format of float printing for printf.  E.g., %4.3d is obtained with
    setintprintwidths(4,3);
   (default is 7.4)
*/
void setfloatprintwidths(int w1,int w2);
void setmatdelim(char *leftdelim, char *rightdelim);
void setmatrowdelim(char *rowdelim);
void setprintnamewnl(int p);


/* set output file stream */
void pprinttofile(FILE *fout);
void pprinttostdout();


void pvecprint(char *name,void *v,int n,int type);
void pvecprintd(char *name, int *v, int n);
void pvecprintx(char *name,int *v,int n);
void pvecprinthd(char *name, short int *v, int n);
void pvecprintld(char *name, long int *v, int n);
void pvecprintsc(char *name, signed char *v, int n);
void pvecprintuc(char *name, unsigned char *v, int n);
void pvecprints(char *name,char *v);
void vecprints(char *name,char *v);
void pvecprintf(char *name, float *v, int n);
void pvecprintlf(char *name, double *v, int n);

void ptenprint(char *name, void ***ten,int nstack,int n, int m, int type);
void ptenprints(char *name,char ***ten,int nstack,int m);
void ptenprintsc(char *name,signed char ***ten,int nstack,int n,int m);
void ptenprintuc(char *name,unsigned char ***ten,int nstack,int n,int m);
void ptenprinthd(char *name,short ***ten,int nstack,int n,int m);
void ptenprintd(char *name,int ***ten,int nstack,int n,int m);
void ptenprintx(char *name,int ***ten,int nstack,int n,int m);
void ptenprintf(char *name,float ***ten,int nstack,int n,int m);
void ptenprintlf(char *name,double ***ten,int nstack,int n,int m);
void ptenprintld(char *name,long ***ten,int nstack,int n,int m);


/* Compute the determinant */
double pdetlf(double **a, int n);
float pdetf(float **a, int n);


/* invert a matrix and return the inverse and the determinant.
   The method used is the LU decomposition.  (see pludcmp below).  
   The numerical analysts hate to invert matrices!!  The fundamental
   rule of numerical analysis is NEVER INVERT A MATRIX.
   (The alternative is to solve directly the set of equations that you
   need to solve using a decomposition instead of an inverse.)
   Unfortunately, most DSP algorithms are most easily expressed in
   terms of matrix inverse, so by default we invert matrices all over the
   place. 
   For general equation solution, however, see the discussion preceeding 
   pludcmp and plubksb).

   The determinant of a is returned.  If it is nonzero, then ainv 
   contains the inverse of a.
 */
float pmatinvf(float **a, float **ainv, int n);
double pmatinvlf(double **a, double **ainv, int n);

/* pludcmp -- compute the LU decomposition of a matrix.  That is, we
   can write
      A = LU
   where L is lower triangular and U is upper triangular.  The
   factorization is accomplished using Gauss elimination with
   pivoting.  The matrix A is overwritten in this routine.

   Once the factorization is accomplished, equations such as
      Ax = b
   can be solved readily from
      (LU)x = b
   and solving by backsubstitution from the triangular structore of L
   and U.  Note that if b changes, we can solve the equations without
   refactorizing A!

   Also, once the LU decomposition is performed, the determinant is
   easily calculated.

   In pludcmp, a is the matrix to factorize, n is the dimension (a
   must be square), indx is an array of integers used to hold the
   pivoting information, d is used to indicate the sign change of the
   pivot, eps indicates the smallest allowable pivot for nonsingular
   matrices.

   To solve a system of equations Ax = b, the following is the best
   way to proceed:

   pludcmpf(a,LUa,indx,n,&d,TINY);
   plubksubf(LUa,n,indx,b);

*/
void pludcmpf(float **a, float **LUa, int n, int *indx,float *d, float eps);
void pludcmplf(double **a, double **LU, int n, int *idx,double *d, double eps);


/* plubksb -- see discussion preceding pludcmp.  Perform
   backsubstition to solve Ax = b after do LU decomposition on A.
   a is the factored matrix
   n is the size
   indx is the array of pivots obtained from pludcmp
   b is the right hand side of the equation
*/
void plubksubf(float **LUa, int n, int *indx,float *b);
void plubksublf(double **LUa, int n, int *indx,double *b);

/* pmatsolve -- solve Ax = b, using ludcmp and lubksub.  Retains
   the factorization of A in internal storage, so that other 
   equations can be solved by calling pmatresolve, unless a function is
   called that modifies tmat1
*/
void pmatsolvef(float **A, float *x, float *b, int n);
void pmatresolvef(float *x, float *b, int n);
void pmatsolvelf(double **A, double *x, double *b, int n);
void pmatresolvelf(double *x, double *b, int n);

/* 
   lueps is a number such that if a pivot is smaller than it, the matrix
   is determined to be singular.  A small value (such as 1.e-20) will
   probably suffice for most work. 

   This number is used by all functions which use LU factorization, 
   including pludcmp, pmatsolve, pmatinv, pdet, and pmineig.

   The default value is 1e-20.

*/
void psetlueps(double lueps);



/* pcholf -- do a cholesky factorization.  That is, given a positive-definite
   symmetric matrix A, factor A as 
   A = LL', where L is a lower triangular matrix
   A is not modified

*/
/* Compute the Cholesky factorization of a symmetric matrix A */
void pcholf(float **A, float **L, int n);
void pchollf(double **A, double **L, int n);

/* Compute the  LDL decomposition for a symmetric matrix, 
   factoring A so that A = L*diag(D)*L' */
void pLDLlf( double **A, double **L, double *D, int n );
void pLDLf( float **A, float **L, float *D, int n );

/********************************************************************/

/* SVD Decomposition:  A = UWV', where
   A is mxn
   U is mxn
   W is an nxn diagonal matrix
   V is nxn
W is returned as a nx1 vector.  The matrices must all be allocated
before calling the routine.

Properties of the SVD:  V is orthogonal
                        U is column orthogonal
*/
void psvdf(float **a, float **u, float *w, float **v, int m, int n);
void psvdlf(double **a, double **u, double *w, double **v, int m, int n);

/********************************************************************/

/*********************************************************
Finds eigenvalues and eigenvectors of a square matrix of dimension dim.  
       Q -- the matrix whose columns are the eigenvectors.
       Lambda -- a diagonal matrix containing the eigenvalues.
       epsilon -- Specifies the precision desired for the eigenvalues,
          by specifying how small the relative adjustment per iteration \
      must be before the eigenvalue is considered found.  Note that
          this is relative and not absolute precision.  A precision of 5e-20 
      for instance will require that every bit in the mantissa of the 
      double precision floating point number be the same, 
      which is unreasonable.  
       maxits -- specifies the maximum number of iterations before giving
          up on the desired accuracy.
This routine uses the QR method of finding eigenvalues without any of 
the neat tricks to speed it up.  It could be better.
Return value:
        the number of iterations required.              
    Q, Lambda are modified as return values also.
*********************************************************/
int peigenlf(double **A, double **Q, double *Lambda, int dim);
int peigenf(float **A, float **Q, float *Lambda, int dim);

/* find the eigenvalues and eigenvectors of a symmetric 2x2 matrix */
/* S is the matrix, v has the eigenvectors as columns, and e 
   has the normalized eigenvalues */
void eigen2lf(double **S, double **v, double *e);
void eigen2f(float **S, float **v, float *e);

/* find the eigenvalues and eigenvectors of a symmetric 3x3 matrix */
/* S is the matrix, v has the eigenvectors as columns, and e 
   has the normalized eigenvalues */
void eigen3lf(double **S, double **v, double *e);


/* houseqr -- get the QR factorization of A, storing the
   Q matrix implicitly (in terms of normalized householder vectors)
   in QR 
*/
void houseqrlf(double **A, double **QR, int m, int n);
void houseqrf(float **A, float **QR, int m, int n);

/* findQlf -- using the matrix QR returned by houseqr, extract Q */
void findQlf(double **QR, double **Q, int m, int n);
void findQf(float **QR, float **Q, int m, int n);

/* findQandRlf -- using the matrix QR returned by houseqr, extract Q and R */
void findQandRlf(double **QR, double **Q, double **R, int m, int n);
void findQandRf(float **QR, float **Q, float **R, int m, int n);



/********************************************************************/

void pmatpsolve1lf(double **A, double *x, double *b, int m, int n);
/* Least-squares Solve Ax = b, where A is overdetermined, using 
   a QR factorization 
*/
void pmatpresolve1lf(double *x, double *b, int m, int n);
/* Least-squares Solve Ax = b, where A is the same as the last call to 
  pmatpsolve1  (and tmat1 has not been modified) */


void applyQlf(double **QR, double *b, int m, int n);
/* Given the compact householder representation in QR, compute
   Q' b, replacing the result in b
*/

void ptrisolupperlf(double **R, double *x, double *b, int n) ;
/* solve the triangular system Rx = b, where R is nxn upper triangular,
   assuming that R[i][i] != 0.
*/

void pmatpsolve1f(float **A, float *x, float *b, int m, int n);
/* Least-squares Solve Ax = b, where A is overdetermined, using 
   a QR factorization 
*/
void pmatpresolve1f(float *x, float *b, int m, int n);
/* Least-squares Solve Ax = b, where A is the same as the last call to 
  pmatpsolve1  (and tmat1 has not been modified) */


void applyQf(float **QR, float *b, int m, int n);
/* Given the compact householder representation in QR, compute
   Q' b, replacing the result in b
*/

void ptrisolupperf(float **R, float *x, float *b, int n) ;
/* solve the triangular system Rx = b, where R is nxn upper triangular,
   assuming that R[i][i] != 0.
*/


/*********************************************************************
Toeplitz stuff:
  ptoepinvlf -- computes the inverse of a symmetric toeplitz matrix
  (from Golub and Van Loan, 2nd ed. p. 190

  pdurbinlf -- solves the yule walker equations T_n y = -(r1 r2 ... rn)
  (this is more efficient than  computing the inverse using ptoepinvlf 
  and then multiplying)

   plevinsonlf -- solve the Toepliz equation T_n x = b for arbitrary right-hand
   side.  (this is more efficient than computing the inverse and multipling)

**********************************************************************/
void ptoepinvlf(double **H,     /* matrix to invert */
               double **B,      /* the inverse */
               int n);          /* matrix size */
void ptoepinvf(float **H,     /* matrix to invert */
               float **B,      /* the inverse */
               int n);          /* matrix size */

/* solve the yule-walker equations  T_n y = -(r1,...,rn)' */
void pdurbinlf(double *r, double *y, int n);
void pdurbinf(float *r, float *y, int n);

/*********************************************************
  Toeplitz solution:  Solve Tx = b, where T is a symmetric
  toeplitz matrix, determined by the vector r
**********************************************************/
void plevinsonlf(double *r, double *x, double *b, int n);
void plevinsonf(float *r, float *x, float *b, int n);


/*********************************************************
  Inverse of a triangular matrix.
  Computes  L^-1
  Where:
  L is rows x rows lower triangular.
  Linv is rows x rows lower triangular.
  Return value:  Linv is modified as a return value.
  *********************************************************/
void ptriinvlf(double **L, double **Linv, int rows);
void ptriinvf(float **L, float **Linv, int rows);



/*********************************************************/
/* This subroutine will attempt to determine the most significant
   eigenvalue and eigenvector for the square symmetric input matrix A.
   The rate of
   convergence depends on the ratio of the most significant eigenvalue
   and the second most significant eigenvalue.  (the larger the
   better)  The Power Method is used.  For more information see:
   "MATRIX Computions" second edition By Gene H. Golub and Charles F.
   Van Loan,ISBN 0-8018-3772-3 or 0-8018-3739-1(pbk.)

   NOTE: this algorithm is may not converge!!!  The number return value is
   the last_lamda - final_lamda.  It will somewhat reflect if the
   algorithm has converged.

   A:  is a MxM matrix
   v:  is where the eigenvector corresponding to lamda that is returned
   lamda: is the largest eigenvalue for A
   M:  is the dimension of A and v.
   Max_iterations:  is the maximuim number of iterations that the
   method will attempt before stopping.  {This guarrenties that the
   algorithm will stop.}
   Error_threshold:  this is the normal stopping codition.   i.e. the
   algorithm will stop when the:

   (last_lamda - final_lamda) < Error_threshold 
 */
float pmaxeigenf(float **A, float *v, float *lamda, int M, int
		 Max_iterations, float Error_threshold);
double pmaxeigenlf(double **A, double *v, double *lamda, int M, int
		  Max_iterations, double Error_threshold);


/*********************************************************/
/* This subroutine will attempt to determine the least significant
   eigenvalue and eigenvector for the input matrix A.  The rate of
   convergence depends on the ratio of the most significant eigenvalue   and 
   the second most significant eigenvalue.  (the larger the
   better)  The Power Method is used.  For more information see:
   "MATRIX Computions" second edition By Gene H. Golub and Charles F.
   Van Loan,ISBN 0-8018-3772-3 or 0-8018-3739-1(pbk.)
   
   NOTE: this algorithm is may not converge!!!  The number return value is
   the last_lamda - final_lamda.  It will somewhat reflect if the
   algorithm has converged.

   A:  is a MxM matrix
   v:  is where the eigenvector corresponding to lamda that is returned
   lamda: is the largest eigenvalue for A
   M:  is the dimension of A and v.
   Max_iterations:  is the maximuim number of iterations that the
   method will attempt before stopping.  {This guarrenties that the
   algorithm will stop.}
   Error_threshold:  this is the normal stopping codition.   i.e. the
   algorithm will stop when the:
                 (last_lamda - final_lamda) < Error_threshold 
 */
float pmineigenf(float **A, float *v, float *lamda, int M, int
		 Max_iterations, float Error_threshold);
double pmineigenlf(double **A, double *v, double *lamda, int M, int
		   Max_iterations, double Error_threshold);


/******************************************************************/
void pfreealltemps(void);   /* free all the temporary vectors and
                               matrices used by the linear algebra routines */


/******************************************************************/
/* Stuff for reading ascii from files */
/******************************************************************/

#ifndef TKUCHAR
#define TKUCHAR    0     /* Unsigned char */
#endif

#ifndef TKSHINT
#define TKSHINT    1     /* Short int */
#endif

#ifndef TKINT
#define TKINT      2     /* Int */
#endif

#ifndef TKFLOAT
#define TKFLOAT    3     /* Float */
#endif

#ifndef TKDOUBLE
#define TKDOUBLE   4     /* Double */
#endif

#ifndef TKLINT
#define TKLINT     5     /* Long int */
#endif

#ifndef TKSCHAR
#define TKSCHAR    6     /* Signed char */
#endif

#ifndef TKCHAR
#define TKCHAR     9    /* Ascii char */
#endif

#ifndef TKHEX
#define TKHEX     10    /* Hexadecimal int */
#endif

/* asread -- Read a scalar from the file fname with name varname into 
   varptr, with
   type, where type is specified by the standard bread library types.  Return
   1 if successful and 0 if not successful */
int asread(char *fname, char *varname,void *varptr, int type);

/* avread -- Read a vector from the file fname with name varname into p with
type, where type is specified by the standard bread library types.
The number of elements read is set in n.

If the variable is not found, 0 is returned; otherwise 1 is returned */
int avread(char *fname, char *varname,void **varptr,int *n, int type);

/* amread -- Read a matrix from the file fname with name varname into p with
type, where type is specified by the standard bread library types.
The number of rows and columns are set in m and n.
Returns  1 if successful and 0 if not successful */
int amread(char *fname, char *varname, void ***varptr, 
           int *m, int *n,int type) ;

/* read a (3-dimensional) tensor from the data file
  the tensor is assumed to be stored in the form  p[m][n][q] 

p[0][0][0] p[0][0][1] ...
p[0][1][0] p[0][1][1] ...
p[0][2][0] p[0][2][1] ...
p[1][0][0] p[1][0][1] ...
p[1][1][0] p[1][1][1] ...
p[1][2][0] p[1][2][1] ...
 ...

Note that m is NOT determined -- it must be passed in to determine how
to block the information.  If the number of rows detected is not a
multiple of m, then an error is declared.  Otherwise, the value of n
is determined by (number of rows)/m.

*/
int atread(char *fname, char *varname,void ****varptr,
            int m, int *n, int *q, int type);



/*
The routines asread, avread, amread, and atread allow for convenient input of
ascii (human-readable) data into programs. 
by which this can be accomplished. 

Variable-by-variable reading

Here is a sample data file:

; Comments can appear on lines by themselves
% Comments are delimited by ; or % or # or C++ or 
# C comments, as long as the C comments begin and end on the 
; same line

n: 2      ; colon after a variable name is ok
m 3       ; only white space after a variable name is ok
n = 2     ; = sign after a variable is ok
  j = 4   ; white space preceeding the variable does not matter
nvect = 2 3 4  ; vectors can be specified
f -1.34e-2     ; floating point formats are ok
g = 34.23
fvect 3.45 23.4

; blank lines don't matter
imat 2 3       ; elements are counted in the first line of matrix to
     1 2       ; determine matrix width
     5  ; incomplete data will be set to zero.

fmat = 2.34, 2
	3. 6

Note the following syntax rules:
1) In general white space is ignored
2) ;, # or % all indicate the beginning of a comment; 
    all text following the character is ignored
3) C or C++ comments can be used, provided that the C comment 
    fits all within one line
3) Any separator (that does not denote a comment) can be used between
data fields.
4) For vectors, the number of fields read determines the size of the
vector.  It may be necessary to pad with zeros in order to get a
desired size.
5) Blank lines are not necessary to separate the fields, even when
reading matrices.
6) Do not put a complete comment line in the middle of a matrix.  That
is, do not do this:

fmat 2 3
; bad place for a comment!
    4 5

The data above can be read with one of the following statements:
   retval = asread(fname,varname,&var,type);                 // read a scalar
   retval = avread(fname,varname,&var,&nr,type);             // read a vector
   retval = amread(fname,varname,&var,&nr,&nc,type);         // read a matrix
   retval = atread(fname,varnme,&var,nstack,&nr,&nc,type);   // read a tensor
The retval==0 indicates an error.

The type variable is defined as follows:
	0 - unsigned char    UCHAR
	1 - short int        SHINT
	2 - int              INT
	3 - float            FLOAT
	4 - double           DOUBLE
	5 - long int         LINT
	6 - signed char      SCHAR
	9 - character ascii CHAR
	10 - hexadecimal     HEX
             (can be prefaced by 0x)
There is not strong type checking between the unsigned char and
signed char.  Caveat Emptor.

Note the following specifics for various types:
ASCII CHARACTERS:
1) The white space preceeding the characters is skipped.
2) For strings, the string is read up until the end of line, 
   and any \n character is removed.
3) The strings are returned null terminated.
4) For character matrcies, the ^ character is used to denote continuation
(so the program can distinguish from variable identifiers).
For example:

strmat first_input
^      second_input
^      third_input

The ^ and preceeding white space for each line is stripped.
5) The longest line of input is used to determine the width of the array.
6) To preserve comment characters, use \ before them: \; \# \%  or quote them
   "#"  '#'  ";"  ';'  "%"  '%'
7) For reading single characters, no quotes are necessary.
   However, quotes can be used:
   c1 = $     ; this is ok
   c2 = '@'   ; a quoted character --- will read @
   c3 = "#"   ; a quoted character --- will read #
   To get the quoted itself, put an escape before it:
   c4 = \'    ; read '
   c5 = \"    ; read "
   or quote it:
   c6 = '"'  ; read "
   c7 = "'"  ; read '
   To get \, just put it by itself
   c8 = \
   or quote it
   c9 = '\'
   c10 = "\"

   To get a space, it must be quoted or used with \:
   c11 = ' '
   c12 = " "
   c13 = \ 

   Other than preceding comment or quote characters or space, the \ is
not regarded as special.

8) For reading character strings, comment indicators $ % ; inside quotes are  
   ignored.  However, removal of quoting characters ', ", or \ is not done.

9) Hex can be just the number:  3a
   or the number preceded by 0:    03a
   or the number preceded by 0x:   0x3a
   However, when forming a hex array, you must make sure that the first
   number on each row does not look like a C identifier, 
   preceding it with 0 or 0x if necessary.

One such read statement is required for each variable.  For example,
the preceding data could be read with the following piece of program
code:

   int n,m,j;        // declare the variables used
   int *nvect;
   float f,g;
   float *fvect;
   int **imat;
   int **fmat;
   char ac;
   char *str;					// ascii string

   char **strmat;				// character matrix

   int nvect_n;       // the number read into nvect
   int fvect_n;       // the number read into fvect
   int nvect_comp;    // the number read into vect_comp
   int imat_r,imat_c; // imat rows and columns
   int fmat_r,fmag_c; // fmat rows and columns

   char fname[20];    // the file name

   // ...

   asread(fname,"m",&m,TKINT);  // read scalars
   asread(fname,"j",&j,TKINT);
   asread(fname,"n",&n,TKINT);
   asread(fname,"f",&f,TKFLOAT);
   asread(fname,"g",&g,TKFLOAT);
   asread(fname,"ac",&ac,TKCHAR);	// read a character as ascii
   asread(fname,"scal_comp",&scal_comp,12);

   avread(fname,"nvect",&nvect,&nvect_n,TKINT);  // read vectors
   avread(fname,"fvect",&fvect,&fvect_n,TKFLOAT);
   avread(fname,"str",&str,&strn,CHAR); // character string
   vect_comp = (complx *)avread(fname,"vect_comp",&vect_comp,&nvect_comp,12);

   amread(fname,"imat",&imat,&imat_r,&imat_c,TKINT);  // read matrices
   amread(fname,"fmat",&fmat,&fmat_r,&fmat_c,TKFLOAT); 
   amread(fname,"strmat",&strmat,&nrow,&nchars,TKCHAR);


Note that the ORDER IN WHICH THE VARIABLES ARE READ DOES NOT MATTER.
This is one of the useful and flexible strengths of this library.

When there are multiple appearences of a variable in the data file, the
FIRST one found it used.  (After finding a variable, no more searching
is performed.)

If the named variable is not found in the file, then the variable pointer is 
NOT modified.  Thie means that you can set up default values which may be
overwritten by calling a*read. 

   int newvar=5;   // set up variable with default value 
   asread(fname,"newvar",&newvar,TKINT);  // if "newvar" not found,
                                          // the original value is retained

However, for avread, amread, and atread, the size parameters are set to zero
if the named variable is not found in the file.

 Also, if the variable is not found, the return
value is 0:

   if(!asread(fname,"newvar",&newvar,TKINT)) {
       newvar = 2;   // set up default value this way 
   }


Todd K. Moon.  Utah State University

*/

/* sort array into INCREASING order, and shuffle array2 at the same time
   if array2 = 0,1,...,N and A represents the original _unsorted_ array,
   then after the sort, A[array2[0]], A[array2[1]], ...
   represents the sorted data in array.
*/
void sort2lfd(int num_elements, double *array, int *array2);
void sort2fd(int num_elements, float *array, int *array2);

/* sort array into increasing order */
void sort1lf(int num_elements, double *array);
void sort1f(int num_elements, float *array);


/**********************************************************************/
/* Generate a unit-variance random number.
   This function calls rand(), so you can control the see using
   void srand(unsigned int seed);
*/
double gran(void);
/* Generate a uniform random number between 0 and 1 */
double uran(void);

#endif  /* end of _PMATLIB_H_ */

/*
Local Variables:
compile-command: "gcc -o testpmatlib -g testpmatlib.c pmatlib.c -lm"
End:
*/
