/*********************
*
*  pmatlib.c -- linear algebra library 
*  using pointer-allocated matrices and vectors
*
*  Todd K. Moon, Utah State University  (and other sources)
*********************************************/
/* Copyright 2004 by Todd K. Moon
 Permission is granted to use this program/data
 for educational/research only
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "pmatlib.h"

#define MAXSTR 5000
#define isCidstart(c) (isalpha(c) || (c == '_'))
#define isCid(c) (isalnum(c) || (c == '_'))
#define isnum(c) (isdigit(c) || (c=='-') || (c=='.'))
#define isintnum(c) (isdigit(c) || (c=='-'))
#define isfnum(c) (isdigit(c) || (c=='-') || (c=='.') || (c=='e') || (c=='E'))

#define ishexnum(c) (isdigit(c) || (c=='a') || (c=='A') || (c=='b') || 	(c=='B') || (c=='c') || (c=='C') || (c=='d') || (c=='D') || (c=='e') ||  (c=='E') || (c=='f') || (c=='F') || (c=='x') || (c== 'X'))

#define iscommentstart(c) (c=='%' || c == '#' || c == ';')
#define isquote(c) (c=='"' || c=='\'')
#define LPR '('
#define RPR ')'

static int *indx = NULL;        /* (nx1) used to record index changes in LU*/
static int mi=0;                /* length of indx */

static float *tvec1f = NULL;    /* (nx1) tempoarary vect. */
static float *tvec2f = NULL;    /* (nx1) tempoarary vect. */
static float *tvec3f = NULL;    /* (nx1) tempoarary vect. */
static double *tvec1lf = NULL;  /* (nx1) tempoarary vect. */
static double *tvec2lf = NULL;  /* (nx1) tempoarary vect. */
static double *tvec3lf = NULL;  /* (nx1) tempoarary vect. */

static int mv1f=0,mv2f=0,mv3f=0;  /* size of tvec1f, tvec2f, tvec3f */
static int mv1lf=0,mv2lf=0,mv3lf=0;  /* size of tvec1lf, tvec2lf,tvev3lf */

static int mm1f=0, nm1f=0;   /* size of tmat1f */
static int mm2f=0, nm2f=0;   /* size of tmat2f */
static int mm3f=0, nm3f=0;   /* size of tmat3f */
static int mm4f=0, nm4f=0;   /* size of tmat4f */
static int mm5f=0, nm5f=0;   /* size of tmat5f */
static int mm1lf=0, nm1lf=0;   /* size of tmat1lf */
static int mm2lf=0, nm2lf=0;   /* size of tmat2lf */
static int mm3lf=0, nm3lf=0;   /* size of tmat3lf */
static int mm4lf=0, nm4lf=0;   /* size of tmat4lf */
static int mm5lf=0, nm5lf=0;   /* size of tmat5lf */


static float **tmat1f = NULL;   /* (nxn) */
static float **tmat2f = NULL;   /* (nxn) */
static float **tmat3f = NULL;   /* (nxn) */
static float **tmat4f = NULL;   /* (mxn) */
static float **tmat5f = NULL;   /* (mxn) */


static double **tmat1lf = NULL; /* (nxn) */
static double **tmat2lf = NULL; /* (nxn) */
static double **tmat3lf = NULL; /* (nxn) */
static double **tmat4lf = NULL; /* (mxn) */
static double **tmat5lf = NULL; /* (mxn) */

static float *pf,*pf1,*pf2,*pf3; /* For pointer indexing. */
static double *pd,*pd1,*pd2,*pd3; /* For pointer indexing. */

static int pbytes(int type);
static char *pdtn(int type);
static void pmaterr(char *error_text);
static int m_foundf_( float a, float b, float e);
static int m_foundlf_( double a, double b, double e);



/* stuff to set print options */
static int intwid1 = 0,intwid2 = 0;
static int flwid1=7,flwid2=4;
static int matprintnamewnl = 0;  /* print a matrix name with a newline */
static char matleftdelim[5] = {0}, matrightdelim[5] = {0};
static char matrowdelim[5] = {'\n'};  /* row delimiter: \n or ; are typical */
static int matprintnlafter = 1;

/* set output file pointer */
static FILE *fout = NULL;

#define TEMPVEC(name,newsize,size,type) if(newsize>size) { pfree_vector((void **)&name); name = (type *)pcalloc_vector2(newsize,sizeof(type)); size = newsize; }
#define TEMPMAT(name,newm,newn,m,n,type) if(newm > m || newn > n) { pfree_matrix((void ***)&name); name = (type **)pcalloc_matrix2(newm,newn,sizeof(type)); m = newm;  n = newn;}
static double lueps = 1e-20;

/***********************************************************************/
static double at,bt,ct;
#define SVDPYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? (ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static double maxarg1,maxarg2;
#define SVDMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define SVDSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*****************************************************************/

void **pcalloc_matrix(int type,int length,int width,char *name)
{
   void **array;
   register int i;
   int nbytes,bwidth;
   
   if(length == 0 || width == 0)
      return NULL;
   nbytes = pbytes(type);
   bwidth = width * nbytes;
   array = (void **)calloc(length,sizeof(void *));
   if(array == NULL) {
      fprintf(stderr,"Can't allocate matrix \"%s\" of type %s.\n",name,
              pdtn(type));
      exit(-1);
   }
   
   array[0] = (char *)calloc(length,bwidth);
   if(array[0] == NULL) {
      fprintf(stderr,"Can't allocate matrix \"%s\" of type %s.\n",name,
              pdtn(type));
      exit(-1);
   }
   for(i=1 ; i<length ; i++)
      array[i] = array[0] + i * bwidth;
   return (void **)array;
}

void **pcalloc_matrix2(int length,int width, int sizeoftype)
{
   void **array;
   register int i;
   int nbytes,bwidth;
   
   if(length == 0 || width == 0)
      return NULL;
   array = (void **)calloc(length,sizeof(void *));
   if(array == NULL) return array;
   
   bwidth = width * sizeoftype;

   array[0] = (void *)calloc(length,bwidth);
   if(array[0] == NULL) {
      free(array);
      return NULL;
   }
   for(i=1 ; i<length ; i++)
      array[i] = array[0] + i * bwidth;
   return (void **)array;
}

/*****************************************************************/

void *pcalloc_vector(int type,int length,char *name)
{
   void *array;
   
   if(length == 0)
      return NULL;
   if(NULL == (array = (void *)calloc(length,pbytes(type)))){
      fprintf(stderr,"Can't allocate array \"%s\" of type %s.\n",name,
              pdtn(type));
      exit(-1);
   }
   return array;
}

void *pcalloc_vector2(int length, int sizeoftype) 
{
   void *array;
   if(length == 0) return NULL;
   if((array = (void *)calloc(length,sizeoftype)) == NULL) {
      return NULL;
   }
   return array;
}


/*****************************************************************/

void pfree_matrix(void ***array)
{
   if(*array != NULL) {
      free(*(array)[0]);
      free(*array);
      *array = NULL;
   }
}


void pfree_vector(void **array)
{
   if(*array != NULL){
      free(*array);
      *array = NULL;
   }
}

/*****************************************************************/
void ***pcalloc_tensor(int type,int nstack,int nrow,int ncol,char *s)
{
   int i;
   void ***m;
   char errmsg[100];
   
   if(nstack == 0 || ncol == 0 || nrow == 0)
      return NULL;
   m=(void ***)malloc( (unsigned) (nstack*sizeof(void***)));
   if (m==NULL){
      fprintf(stderr,"Can't allocate tensor \"%s\" of type %s.\n",s,
              pdtn(type));
      exit(2);
   }
   errmsg[0] = 0;
   strcat(errmsg,s);
   strcat(errmsg,"  (calloc_tensor)");
   for(i = 0; i < nstack; i++) {
      m[i] = (void **)pcalloc_matrix(type,nrow,ncol,errmsg);
   }
   return m;
}


void ***pcalloc_tensor2(int nmat,int nrow,int ncol,int sizeoftype)
{
   int i;
   void ***m;

   if(nmat == 0 || ncol == 0 || nrow == 0)
      return NULL;
   m=(void ***)malloc( (unsigned) (nmat*sizeof(void***)));
   if (m==NULL) return(m);

   for(i = 0; i < nmat; i++) {
      m[i] = (void **)pcalloc_matrix2(nrow,ncol,sizeoftype);
   }
   return m;
}

void pfree_tensor(void ****m,int nstack)
{
   int i;
   
   if(*m != NULL) {
      for(i = 0; i < nstack; i++) {
         pfree_matrix(*m +i);
      }
      free(*m);
      *m = NULL;
   }
}

/*****************************************************************/

void pveczerod(int *v,int dim)
{
   int i;
   for(i = 0; i < dim; i++)  v[i] = 0;
}

void pveczerof(float *v,int dim)
{
   int i;
   for(i = 0; i < dim; i++)  v[i] = 0;
}

void pveczerolf(double *v,int dim)
{
   int i;
   for(i = 0; i < dim; i++)  v[i] = 0;
}

void pmatzerof(float **m, int r, int c)
{
   int i,j;
   for(i = 0; i < r; i++) {
      pf = m[i];
      for(j = 0; j < c; j++) pf[j] = 0;
   }
}


void pmatzerolf(double **m, int r, int c)
{
   int i,j;
   for(i = 0; i < r; i++) {
      pd = m[i];
      for(j = 0; j < c; j++) pd[j] = 0;
   }
}


void pvecconstf(float *v,float c,int dim)
{
   int i;
   for(i = 0; i < dim; i++)  v[i] = c;
}


void pvecconstlf(double *v,double c,int dim)
{
   int i;
   for(i = 0; i < dim; i++)  v[i] = c;
}

void pmatconstf(float **mat, float c,int m, int n)
{
   int i,j;
   for(i = 0; i < m; i++) {
      pf=mat[i];
      for(j = 0; j < n; j++) pf[j] = c;
   }
}

void pmatconstlf(double **mat, double c,int m, int n)
{
   int i,j;
   for(i = 0; i < m; i++) {
      pd=mat[i];
      for(j = 0; j < n; j++) pd[j] = c;
   }
}

/*****************************************************************/
void pvecaddf(float *a,float *b, float *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i] + b[i];
}

void pvecaddlf(double *a,double *b, double *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i] + b[i];
}

/*****************************************************************/
void pvecsubf(float *a,float *b, float *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i] - b[i];
}


void pvecsublf(double *a,double *b, double *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i] - b[i];
}


/*****************************************************************/

void pvecscalef(float *a,float b, float *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i]*b;
}

void pvecscalelf(double *a,double b, double *c,int n)
{
   int i;
   for(i =0; i < n; i++) 
      c[i] = a[i]*b;
}

/*****************************************************************/
/* Inner product: c = a'b */
float pvecinnerf(float *a, float *b, int n)
{
   int i;
   
   float c = 0;
   for(i=0; i < n; i++) c += a[i] * b[i];
   return(c);
}

double pvecinnerlf(double *a, double *b, int n)
{
   int i;
   double c = 0; 
   for(i=0; i < n; i++) c += a[i] * b[i];
   return(c);
}


/*****************************************************************/
void pvecvecmultf(float *a, float *b, float *c, int n)
{
   int i;
   for(i = 0; i < n; i++) c[i] = a[i]*b[i];
}

void pvecvecmultlf(double *a, double *b, double *c, int n)
{
   int i;
   for(i = 0; i < n; i++) c[i] = a[i]*b[i];
}

/*****************************************************************/
/* copy: a --> b */
void pveccopylf(double *a,double *b,int m)
{
   int i;
   for(i = 0; i < m; i++)
      b[i] = a[i];
}

void pveccopyf(float *a,float *b,int m)
{
   int i;
   for(i = 0; i < m; i++)
      b[i] = a[i];
}

void pveccopyd(int *a,int *b,int m)
{
   int i;
   for(i = 0; i < m; i++)
      b[i] = a[i];
}


void pveccopyc(char *a,char *b,int m)
{
   int i;
   for(i = 0; i < m; i++)
      b[i] = a[i];
}


/*****************************************************************/
/* Vector norm */
double pvecnorm2lf(double *v, int n)
{
   int i;
   double sum = 0;
   for(i=0,pd=v; i < n; i++) 
      sum += v[i] * v[i];
   return(sqrt(sum));
}

float pvecnorm2f(float *v, int n)
{
   int i;
   double sum = 0;
   for(i = 0; i < n; i++) 
      sum += v[i] * v[i];
   return(sqrt(sum));
}

/*****************************************************************/
/* convolve the sequence in1 (of n1 points) with the sequence in2
   (of n2 points) into out
   This does not do "fast" convolution (e.g., using FFTs).
   Return value is the length of the new sequence
*/
int pconvlf(double *in1, int l1, double *in2, int l2, double *out)
{
   int i,j,nout;
   nout = l1+l2-1;
   /* clear out the destination */
   memset(out,0,sizeof(double)*(nout));
   for(i = 0; i < l1; i++) {
      for(j = 0; j < l2; j++) {
         out[i+j] += in1[i]*in2[j];
      }
   }
   return nout;
}

int pconvf(float *in1, int l1, float *in2, int l2, float *out)
{
   int i,j,nout;
   nout = l1+l2-1;
   /* clear out the destination */
   memset(out,0,sizeof(float)*(nout));
   for(i = 0; i < l1; i++) {
      for(j = 0; j < l2; j++) {
         out[i+j] += in1[i]*in2[j];
      }
   }
   return nout;
}


/*****************************************************************/
int pminveclf(double *v, int n)
{
   int i,j;
   for(i=1, j=0; i < n; i++) 
      if(v[i] < v[j]) j = i;
   return j;
}

int pmaxveclf(double *v, int n)
{
   int i,j;
   for(i=1, j=0; i < n; i++) 
      if(v[i] > v[j]) j = i;
   return j;
}

int pminvecf(float *v, int n)
{
   int i,j;
   for(i=1, j=0; i < n; i++) 
      if(v[i] < v[j]) j = i;
   return j;
}

int pmaxvecf(float *v, int n)
{
   int i,j;
   for(i=1, j=0; i < n; i++) 
      if(v[i] > v[j]) j = i;
   return j;
}

/*****************************************************************/
void pmataddf(float **a,float **b,float ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pf=c[i],pf1=a[i],pf2=b[i]; j < n; j++) {
         pf[j] = pf1[j] + pf2[j];
      }
   }
}

void pmataddlf(double **a,double **b,double ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pd=c[i],pd1=a[i],pd2=b[i]; j < n; j++) {
         pd[j] = pd1[j] + pd2[j];
      }
   }
}

/*****************************************************************/
void pmatsublf(double **a,double **b,double ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pd=c[i],pd1=a[i],pd2=b[i]; j < n; j++) {
         pd[j] = pd1[j] - pd2[j];
      }
   }
}

void pmatsubf(float **a,float **b,float ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pf=c[i],pf1=a[i],pf2=b[i]; j < n; j++) {
         pf[j] = pf1[j] - pf2[j];
      }
   }
}


/*****************************************************************/
void pmatmatmultlf(double **a,double **b,double ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pd=c[i],pd1=a[i],pd2=b[i]; j < n; j++) {
         pd[j] = pd1[j] * pd2[j];
      }
   }
}

void pmatmatmultf(float **a,float **b, float ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pf=c[i],pf1=a[i],pf2=b[i]; j < n; j++) {
         pf[j] = pf1[j] * pf2[j];
      }
   }
}

/*****************************************************************/
void pmataddandscalelf(double **a,double **b,double d,double ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pd=c[i],pd1=a[i],pd2=b[i]; j < n; j++) {
         pd[j] = d*(pd1[j] + pd2[j]);
      }
   }
}

void pmataddandscalef(float **a,float **b,float d,float ** c,int m,int n)
{
   int i,j;
   
   for(i = 0; i < m; i++) {
      for(j=0,pf=c[i],pf1=a[i],pf2=b[i]; j < n; j++) {
         pf[j] = d*(pf1[j] + pf2[j]);
      }
   }
}

/*****************************************************************/
/* Copy matrices: a --> B */
void pmatcopylf(double **a,double **b,int m,int n)
{
   int i,j;
   for(i = 0; i < m; i++)
      for(j=0,pd=b[i],pd1=a[i]; j < n; j++)
         pd[j] = pd1[j];
}

void pmatcopyf(float **a,float **b,int m,int n)
{
   int i,j;
   for(i = 0; i < m; i++)
      for(j=0,pf=b[i],pf1=a[i]; j < n; j++)
         pf[j] = pf1[j];
}

/*****************************************************************/
void pmatscalef(float **a, float b,float **c, int n, int m)
{
   int i,j;
        
   for(i = 0; i < n; i++) {
      for(j=0,pf=c[i],pf1=a[i]; j < m; j++)
         pf[j] = pf1[j] * b;
   }
}


void pmatscalelf(double **a, double b,double **c, int n, int m)
{
   int i,j;
   
   for(i = 0; i < n; i++) {
      for(j=0,pd=c[i],pd1=a[i]; j < m; j++)
         pd[j] = pd1[j] * b;
   }
}


/*****************************************************************/
void pmattpoself(double **a,double **b,int m,int n)
   /* transpose a matrix (m x n) into b matrix, assumed to be declared as */
         /* (n x m) */
{
   int i,j;
   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
         b[j][i] = a[i][j];
      }
   }
}


void pmattposef(float **a,float **b,int m,int n)
{
   int i,j;
   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
         b[j][i] = a[i][j];
      }
   }
}


/*****************************************************************/
void pmatmultlf(double **a, double **b, double **c, int n, int m, int q)
/* a is nxm, b is mxq, c is nxq */
{
   int i,j,k;
   for(i=0 ; i<n ; i++){
      for(k=0,pd1=c[i] ; k<q ; k++)
         pd1[k] = 0.0;
      for(k=0,pd2=a[i] ; k<m ; k++)
         if(pd2[k] != 0.0)
            for(j=0,pd1=c[i],pd3=b[k] ; j<q ; j++)
               pd1[j] += pd2[k] * pd3[j];
   }
}


void pmatmultf(float **a, float **b, float **c, int n, int m, int q)
/* a is nxm, b is mxq, c is nxq */
{
   int i,j,k;
   for(i=0 ; i<n ; i++){
      for(k=0,pf1=c[i] ; k<q ; k++)
         pf1[k] = 0.0;
      for(k=0,pf2=a[i] ; k<m ; k++)
         if(pf2[k] != 0.0)
            for(j=0,pf1=c[i],pf3=b[k] ; j<q ; j++)
               pf1[j] += pf2[k] * pf3[j];
   }
}


/*****************************************************************/
void pmatmultTlf(double **a,double **b,double **c,int n,int m,int q,int type)
   /* type 0: c = a*b
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
{
   int i,j,k;
   double temp;
   
   switch(type) {
   case 0:              /* c = a*b  */
      for(i=0 ; i<n ; i++){
         for(k=0,pd1=c[i] ; k<q ; k++)
            pd1[k] = 0.0;
         for(k=0,pd2=a[i] ; k<m ; k++)
            if(pd2[k] != 0.0)
               for(j=0,pd1=c[i],pd3=b[k] ; j<q ; j++)
                  pd1[j] += pd2[k] * pd3[j];
      }
      break;
   case 1:                      /* c = a*b' */
      for(i=0; i < n; i++) {
         for(k=0,pd=c[i]; k < q; k++) {
            temp = 0.0;
            for(j=0,pd1=a[i],pd2=b[k]; j < m; j++) {
               temp += pd1[j] * pd2[j];
            }
            pd[k] = temp;
         }
      }
      break;
   case 2:                      /* c = a'*b  */
      for(i = 0; i < n; i++)
         for(k = 0,pd=c[i]; k < q; k++)
            pd[k] = 0.0;
      for(j = 0; j < m; j++)
         for(i = 0,pd2=a[j] ; i < n; i++)
            if(pd2[i] != 0.0)
               for(k = 0,pd=c[i],pd1=b[j] ; k < q; k++)
                  pd[k] += pd2[i] * pd1[k];
      break;
   case 3:                      /* c = a'*b' */

/* type 3: c = a'*b'
a is m x n
b is q x m
c is n x q
*/
      pd2 = a[0];
      for(i = 0; i < n; i++) {
         for(k=0,pd=c[i]; k < q; k++) {
            temp = 0;
            for(j=0,pd1=b[k]; j < m; j++,pd2+=n)
               temp += pd2[i] * pd1[j];
            pd2 = a[0];
            pd[k] = temp;
         }
      }
      break;
   default:                     /* invalid type */
      pmaterr("invalid multiplication type\n");
   }
}

void pmatmultTf(float **a,float **b,float **c,int n,int m,int q,int type)
   /* type 0: c = a*b
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
{
   int i,j,k;
   float temp;
   
   switch(type) {
   case 0:              /* c = a*b  */
      for(i=0 ; i<n ; i++){
         for(k=0,pf1=c[i] ; k<q ; k++)
            pf1[k] = 0.0;
         for(k=0,pf2=a[i] ; k<m ; k++)
            if(pf2[k] != 0.0)
               for(j=0,pf1=c[i],pf3=b[k] ; j<q ; j++)
                  pf1[j] += pf2[k] * pf3[j];
      }
      break;
   case 1:                      /* c = a*b' */
      for(i=0; i < n; i++) {
         for(k=0,pf=c[i]; k < q; k++) {
            temp = 0.0;
            for(j=0,pf1=a[i],pf2=b[k]; j < m; j++) {
               temp += pf1[j] * pf2[j];
            }
            pf[k] = temp;
         }
      }
      break;
   case 2:                      /* c = a'*b  */
      for(i = 0; i < n; i++)
         for(k = 0,pf=c[i]; k < q; k++)
            pf[k] = 0.0;
      for(j = 0; j < m; j++)
         for(i = 0,pf2=a[j] ; i < n; i++)
            if(pf2[i] != 0.0)
               for(k = 0,pf=c[i],pf1=b[j] ; k < q; k++)
                  pf[k] += pf2[i] * pf1[k];
      break;
   case 3:                      /* c = a'*b' */
      pf2 = a[0];
      for(i = 0; i < n; i++) {
         for(k=0,pf=c[i]; k < q; k++) {
            temp = 0;
            for(j=0,pf1=b[k]; j < m; j++,pf2+=n)
               temp += pf2[i] * pf1[j];
            pf2 = a[0];
            pf[k] = temp;
         }
      }
      break;
   default:                     /* invalid type */
      pmaterr("invalid multiplication type\n");
   }
}



/*****************************************************************/
void pmatvecmultf(float **a, float *b,float *c, int n, int m)
   /* c = A*b, A is (nxm), b is(mx1), c is (nx1) */
{
   int i,j;
   float temp;
   
   for(i=0; i < n; i++) {
      temp = 0.0;
      for(j=0,pf=a[i]; j < m; j++)
         temp += pf[j] * b[j];
      c[i] = temp;
   }
}

void pmatvecmultlf(double **a, double *b,double *c, int n, int m)
   /* c = A*b, A is (nxm), b is(mx1), c is (nx1) */
{
   int i,j;
   double temp;
   
   
   for(i=0; i < n; i++) {
      temp = 0;
      for(j=0,pd2=a[i]; j < m; j++)
         temp += pd2[j] * b[j];
      c[i] = temp;
   }
}


/*****************************************************************/
/* Multiply a*b', where a is mx1 and b is nx1, 
   to obtain the mxn outer product c */
void pvecouterlf(double *a, double *b, double **c, int m, int n)
{
   int i, j;
   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
         c[i][j] = a[i]*b[j];
      }
   }
}

void pvecouterf(float *a, float *b, float **c, int m, int n)
{
   int i, j;
   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
         c[i][j] = a[i]*b[j];
      }
   }
}


/*****************************************************************/
/* Compute c = v'Mv (a generalized vector inner product */
double pvecinnerMlf(double *v, double **M, double *w, int n)
{
   double sum = 0;
   int i,j;
   
   for(i=0; i < n; i++)
      for(j=0,pd2=M[i]; j < n; j++)
         sum += v[i] * pd2[j] * w[j];
   return(sum);
}

float pvecinnerMf(float *v, float **M, float *w, int n)
{
   float sum = 0;
   int i,j;
   
   for(i=0; i < n; i++)
      for(j=0,pf2=M[i]; j < n; j++)
         sum += v[i] * pf2[j] * w[j];
   return(sum);
}

/*****************************************************************/
/* Multiply a vector (transpose) times a matrix: c = a'*B */
/* a is (nx1), b is (nxm), c is (1xm) (a vector) */
void pvecmatmultf(float *a, float **b,float *c, int n, int m)
{
   int i,j;
   
   for(j=0; j < n; j++)
      c[j] = 0.0;
   for(j=0; j < n; j++){
      for(i=0,pf=b[j]; i < m; i++)
         c[i] += pf[i] * a[j];
   }
}

void pvecmatmultlf(double *a, double **b,double *c, int n, int m)
{
   int i,j;
        
   for(j=0; j < n; j++)
      c[j] = 0.0;
   for(j=0; j < n; j++){
      for(i=0,pd=b[j]; i < m; i++)
         c[i] += pd[i] * a[j];
   }
}


/*****************************************************************/
/* Multiply a vector (transpose) times a matrix transpose: c = a'*B' */
/* a is (nx1), b is (mxn), c is (1xm) (a vector) */
void pvecmatTmultlf(double *a, double **b,double *c, int n, int m)
{
   int i,j;
   double temp;
   
   for(i=0; i < m; i++) {
      temp = 0;
      for(j=0,pd=b[i]; j < n; j++)
         temp += a[j] * pd[j];
      c[i] = temp;
   }
}

void pvecmatTmultf(float *a, float **b,float *c, int n, int m)
{
   int i,j;
   float temp;
   
   for(i=0; i < m; i++) {
      temp = 0;
      for(j=0,pf=b[i]; j < n; j++)
         temp += a[j] * pf[j];
      c[i] = temp;
   }
}



/*****************************************************************/

/* Multiply a matrix transpose times a vector: c = A'*b, 
   A is (nxm), b is(nx1), c is (mx1) */
void pmatTvecmultf(float **a, float *b,float *c, int n, int m)
{
   int i,j;
   
   for(j=0; j < m; j++)
      c[j] = 0.0;
   for(j=0; j < n; j++){
      for(i=0,pf=a[j]; i < m; i++)
         c[i] += pf[i] * b[j];
   }
}

void pmatTvecmultlf(double **a, double *b,double *c, int n, int m)
{
   int i,j;
   
   for(j=0; j < m; j++)
      c[j] = 0;
   for(j=0; j < n; j++) {
      for(i = 0, pd=a[j]; i < m; i++)
         c[i] += pd[i] * b[j];
   }
}


/* form the product c = a * diag(d), where a is nxn */
void pmatdiagmultf(float **a, float *d, float **c, int m, int n)
{
   int i,j;
   for(i = 0; i < m; i++) {
      for(j = 0,pf1 = a[i],pf2 = c[i]; j < n; j++,pf1++,pf2++) {
         *pf2 = *pf1 * d[j];
      }
   }
}

void pmatdiagmultlf(double **a, double *d, double **c, int m, int n)
{
   int i,j;
   for(i = 0; i < m; i++) {
      for(j = 0,pd1 = a[i],pd2 = c[i]; j < n; j++,pd1++,pd2++) {
         *pd2 = *pd1 * d[j];
      }
   }
}


/* Find the trace of a matrix. */
float pmattracef(float **A, int M)
{
   int i;
   float trace = 0;
   
   for(i = 0; i < M; i++) trace += A[i][i];
   return(trace);
}

double pmattracelf(double **A, int M)
{
   int i;
   double trace = 0;
   
   for(i = 0; i < M; i++) trace += A[i][i];
   return(trace);
}

/*********************************************************
Finds the frobenius norm of a matrix.
Where:
  a is rows x cols.
  Return value:
  double ||a|| where ||a|| is the frobenius norm.
*********************************************************/
double pfrobeniuslf( double **a, int rows, int cols )
{
   int i, j;         /*Counters.*/
   double sum;       /*The running frobenius norm value.*/
   
   sum = 0.0;
   for ( i = 0; i < rows; i++ ) {
      for ( j=0,pd=a[i]; j < cols; j++,pd++ ) {
         sum += (*pd)*(*pd);
      }
   }
   return( sqrt( sum ) );
}


float pfrobeniusf( float **a, int rows, int cols )
{
   int i, j;         /*Counters.*/
   float sum;       /*The running frobenius norm value.*/
   
   sum = 0.0;
   for ( i = 0; i < rows; i++ ) {
      for ( j=0,pf=a[i]; j < cols; j++,pf++ ) {
         sum += (*pf)*(*pf);
      }
   }
   return( sqrt( sum ) );
}

/*****************************************************************/
void ptencopylf(double ***a, double ***b, int nstack, int r, int c)
{
   int i,j,k;

   for(i = 0; i < nstack; i++) {
      for(j = 0; j < r; j++) {
	 for(k = 0; k < c; k++) {
	    b[i][j][k] = a[i][j][k];
	 }
      }
   }
}

void ptencopyf(float ***a, float ***b, int nstack, int r, int c)
{
   int i,j,k;

   for(i = 0; i < nstack; i++) {
      for(j = 0; j < r; j++) {
	 for(k = 0; k < c; k++) {
	    b[i][j][k] = a[i][j][k];
	 }
      }
   }
}


/*****************************************************************/


void pmatprint(char *name, void **mat, int m,int n, int type)
{
   switch(type) {
   case TKSCHAR:
      pmatprintsc(name,(signed char **)mat,m,n);
      break;
   case TKUCHAR:
      pmatprintuc(name,(unsigned char **)mat,m,n);
      break;
   case TKCHAR:
      pmatprints(name,(char **)mat,m);
      break;
   case TKSHINT:
      pmatprinthd(name,(short **)mat,m,n);
      break;
   case TKINT:
      pmatprintd(name,(int **)mat,m,n);
      break;
   case TKFLOAT:
      pmatprintf(name,(float **)mat,m,n);
      break;
   case TKDOUBLE:
      pmatprintlf(name,(double **)mat,m,n);
      break;
   case TKLINT:
      pmatprintld(name,(long **)mat,m,n);
      break;
   case TKHEX:
      pmatprintx(name,(int **)mat,m,n);
      break;
   default:
      fprintf(stderr,"undefined type %d\n",type);
   }
}



void pmatprints(char *name,char **mat,int n)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;

   fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }

   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      fprintf(fout,"%s",mat[i]);
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

void pmatprintsc(char *name,signed char **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*d\t",intwid1,intwid2,(int)mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

void pmatprintuc(char *name,unsigned char **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*ud\t",intwid1,intwid2,(unsigned int)mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pmatprinthd(char *name,short int **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*hd\t",intwid1,intwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}


/*****************************************************************/
void pmatprintd(char *name,int **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*d\t",intwid1,intwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

void pmatprintx(char *name,int **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*x\t",intwid1,intwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pmatprintf(char *name,float **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*g\t",flwid1,flwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pmatprintlf(char *name,double **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");

   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*g\t",flwid1,flwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pmatprintld(char *name,long int **mat,int n,int m)
{
   int i,j;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(mat == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++) {
      for(j = 0; j < m; j++) {
         fprintf(fout,"%*.*ld\t",intwid1,intwid2,mat[i][j]);
      }
      if(i < n-1) fprintf(fout,"%s",matrowdelim);
   }
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}


/*****************************************************************/
/* set the format of integer printing for printf.  E.g., %4.3d is obtained with
    setintprintwidths(4,3);
   (default is 0.0)
*/
void setintprintwidths(int w1,int w2)
{  intwid1 = w1;
 intwid2 = w2;
}

/* set the format of float printing for printf.  E.g., %4.3d is obtained with
    setintprintwidths(4,3);
   (default is 7.4)
*/
void setfloatprintwidths(int w1,int w2)
{  flwid1 = w1;
 flwid2 = w2;
}

void setmatdelim(char *leftdelim, char *rightdelim)
{
   strcpy(matleftdelim,leftdelim);
   strcpy(matrightdelim,rightdelim);
}

void setmatrowdelim(char *rowdelim)
{  strcpy(matrowdelim,rowdelim);
}


void setprintnamewnl(int p)
{  matprintnamewnl = p;
}

void setprintnlafter(int p)
{  matprintnlafter = p;
}


/*****************************************************************/
void pprinttofile(FILE *f)
{   fout = f;
}

void pprinttostdout()
{   fout = stdout;
}


/*****************************************************************/
void pvecprint(char *name, void *vec, int n, int type)
{
   switch(type) {
   case TKSCHAR:
      pvecprintsc(name,(signed char *)vec,n);
      break;
   case TKUCHAR:
      pvecprintuc(name,(unsigned char *)vec,n);
      break;
   case TKCHAR:
      pvecprints(name,(char *)vec);
      break;
   case TKSHINT:
      pvecprinthd(name,(short *)vec,n);
      break;
   case TKINT:
      pvecprintd(name,(int *)vec,n);
      break;
   case TKFLOAT:
      pvecprintf(name,(float *)vec,n);
      break;
   case TKDOUBLE:
      pvecprintlf(name,(double *)vec,n);
      break;
   case TKLINT:
      pvecprintld(name,(long *)vec,n);
      break;
   case TKHEX:
      pvecprintx(name,(int *)vec,n);
      break;
   default:
      fprintf(stderr,"undefined type %d\n",type);
   }
}

   
/*****************************************************************/
void pvecprintd(char *name,int *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*d\t",intwid1,intwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pvecprintx(char *name,int *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*x\t",intwid1,intwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pvecprinthd(char *name,short int *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*hd\t",intwid1,intwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}


/*****************************************************************/
void pvecprintld(char *name,long int *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*ld\t",intwid1,intwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}


/*****************************************************************/
void pvecprintsc(char *name,signed char *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*d\t",intwid1,intwid2,(int)v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}


void pvecprintuc(char *name,unsigned char *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*ud\t",intwid1,intwid2,(unsigned int)v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}

void pvecprints(char *name,char *v) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",v);
   if(matprintnlafter) fprintf(fout,"\n");
}

/*****************************************************************/
void pvecprintf(char *name,float *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*g\t",flwid1,flwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");

}


/*****************************************************************/
void pvecprintlf(char *name,double *v,int n) 
{
   int i;
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s\t",name);
   if(v == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   fprintf(fout,"%s",matleftdelim);
   for(i = 0; i < n; i++)
      fprintf(fout,"%*.*g\t",flwid1,flwid2,v[i]);
   fprintf(fout,"%s",matrightdelim);
   if(matprintnlafter) fprintf(fout,"\n");
}




/*****************************************************************/

void ptenprint(char *name, void ***ten,int nstack,int n, int m, int type)
{
   switch(type) {
   case TKCHAR:
      ptenprints(name,(char ***)ten,nstack,n);
      break;
   case TKSCHAR:
      ptenprintsc(name,(signed char ***)ten,nstack,n,m);
      break;
   case TKUCHAR: 
      ptenprintuc(name,(unsigned char ***)ten,nstack,n,m);
      break;
   case TKSHINT:
      ptenprinthd(name,(short ***)ten,nstack,n,m);
      break;
   case TKINT:
      ptenprintd(name,(int ***)ten,nstack,n,m);
      break;
   case TKFLOAT:
      ptenprintf(name,(float ***)ten,nstack,n,m);
      break;
   case TKDOUBLE:
      ptenprintlf(name,(double ***)ten,nstack,n,m);
      break;
   case TKLINT:
      ptenprintld(name,(long ***)ten,nstack,n,m);
      break;
   case TKHEX:
      ptenprintx(name,(int ***)ten,nstack,n,m);
      break;
   default:
      fprintf(stderr,"Type %d not implemented for ptenprint\n", type);
   }
}


void ptenprints(char *name,char ***ten,int nstack,int m)
{
   int i,j;

   if(fout == NULL) fout = stdout;   
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < m; j++) {
	 fprintf(fout,"%s\n",ten[i][j]);
      }
      fprintf(fout,"\n");
   }
}

/*****************************************************************/
void ptenprintsc(char *name,signed char ***ten,int nstack,int n,int m)
{
   int i,j,k;
   

   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*d\t",intwid1,intwid2,(int)ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}

/*****************************************************************/
void ptenprintuc(char *name,unsigned char ***ten,int nstack,int n,int m)
{
   int i,j,k;
   

   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*d\t",intwid1,intwid2,(int)ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}


/*****************************************************************/
void ptenprinthd(char *name,short ***ten,int nstack,int n,int m)
{
   int i,j,k;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*hd\t",intwid1,intwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}


/*****************************************************************/
void ptenprintd(char *name,int ***ten,int nstack,int n,int m)
{
   int i,j,k;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*d\t",intwid1,intwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}

/*****************************************************************/
void ptenprintx(char *name,int ***ten,int nstack,int n,int m)
{
   int i,j,k;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*x\t",intwid1,intwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}


/*****************************************************************/
void ptenprintf(char *name,float ***ten,int nstack,int n,int m)
{
   int i,j,k;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*g\t",flwid1,flwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}


/*****************************************************************/
void ptenprintlf(char *name,double ***ten,int nstack,int n,int m)
{
   int i,j,k;
   
   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*g\t",flwid1,flwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}



/*****************************************************************/
void ptenprintld(char *name,long ***ten,int nstack,int n,int m)
{
   int i,j,k;

   if(fout == NULL) fout = stdout;
   if(strlen(name)) fprintf(fout,"%s",name);
   if(matprintnamewnl) fprintf(fout,"\n");
   if(ten == NULL) {
      fprintf(fout,"(null)\n");
      return;
   }
   for(i = 0; i < nstack; i++) {
      for(j = 0; j < n; j++) {
	 for(k = 0; k < m; k++) {
	    fprintf(fout,"%*.*ld\t",intwid1,intwid2,ten[i][j][k]);
	 }
	 fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
   }
}


/*****************************************************************/
/*****************************************************************/
/* compute the determinant */
float pdetf(float **a,int n)
{
   float d;
   int j;
   
   TEMPVEC(indx,n,mi,int);
   TEMPMAT(tmat1f,n,n, mm1f, nm1f, float);

   pludcmpf(a,tmat1f,n,indx,&d,lueps);	/* Lu decompose the matrix */
   
   /* for compatibility with the previous library, also compute the */
   /* determinant */
   for(j = 0; j < n; j++) {
      d *= tmat1f[j][j];
   }
   return(d);
}

double pdetlf(double **a,int n)
{
   double d;
   int j;
   
   TEMPVEC(indx,n,mi,int);
   TEMPMAT(tmat1lf,n,n, mm1lf, nm1lf, double);

   pludcmplf(a,tmat1lf,n,indx,&d,lueps);	/* Lu decompose the matrix */
   
   /* for compatibility with the previous library, also compute the */
   /* determinant */
   for(j = 0; j < n; j++) {
      d *= tmat1lf[j][j];
   }
   return(d);
}



/*****************************************************************/
void pmatsolvef(float **A, float *x, float *b, int n)
/* solve Ax = b */
{
   float d;
   int i;
   for(i = 0; i < n; i++) {
      x[i] = b[i];
   }
   TEMPMAT(tmat1f,n,n,mm1f,nm1f,float);
   TEMPVEC(indx,n,mi,int);
   pludcmpf(A,tmat1f,n,indx,&d,lueps);
   plubksubf(tmat1f,n,indx,x);
}


void pmatresolvef(float *x, float *b, int n)
/* Solve Ax = b, where A is the same as the last call to pmatsolve
   (and tmat1 has not been modified) */
{
   int i;
   for(i = 0; i < n; i++) {
      x[i] = b[i];
   }
   plubksubf(tmat1f,n,indx,x);
}


void pmatsolvelf(double **A, double *x, double *b, int n)
/* solve Ax = b */
{
   double d;
   int i;
   for(i = 0; i < n; i++) {
      x[i] = b[i];
   }
   TEMPMAT(tmat1lf,n,n,mm1lf,nm1lf,double);
   TEMPVEC(indx,n,mi,int);
   pludcmplf(A,tmat1lf,n,indx,&d,lueps);
   plubksublf(tmat1lf,n,indx,x);
}


void pmatresolvelf(double *x, double *b, int n)
/* Solve Ax = b, where A is the same as the last call to pmatsolve
   (and tmat1 has not been modified) */
{
   int i;
   for(i = 0; i < n; i++) {
      x[i] = b[i];
   }
   plubksublf(tmat1lf,n,indx,x);
}

/*****************************************************************/
float pmatinvf(float **a,float **ainv, int n)
{
   float d;
   int i,j;
   
   TEMPVEC(indx,n,mi,int);
   /* TEMPVEC(tvec1,n,mv1f,float);  <-- done by pludcmp */
   TEMPMAT(tmat1f,n,n,mm1f,nm1f,float);

   pludcmpf(a,ainv,n,indx,&d,lueps);	/* Lu decompose the matrix */
   
   /* for compatibility with the previous library, also compute the */
   /* determinant */
   for(j = 0; j < n; j++) {
      d *= ainv[j][j];
   }
   if(d) {
      for(j = 0; j < n; j++) {
	 for(i = 0; i < n; i++) tvec1f[i] = 0;
	 tvec1f[j] = 1.; 
	 plubksubf(ainv,n,indx,tvec1f);
	 for(i = 0; i < n; i++)
	    tmat1f[i][j] = tvec1f[i];
      }
      /* now copy the inverse into the called matrix */
      for(i = 0; i < n; i++) 
	 for(j = 0,pf=ainv[i],pf1=tmat1f[i]; j < n; j++)
	    pf[j] = pf1[j];
   }
   return(d);
} /* end of inverse */

void psetlueps(double luepsin)
{  lueps = luepsin;
}

double pmatinvlf(double **a,double **ainv, int n)
{
   double d;
   int i,j;
   register double *p,*p1;

   TEMPVEC(indx,n,mi,int);
   /* TEMPVEC(tvec1lf,n,mv1lf,double);  <-- done by pludcmp */
   TEMPMAT(tmat1lf,n,n,mm1lf,nm1lf,double);
   
   pludcmplf(a,ainv,n,indx,&d,lueps);	/* Lu decompose the matrix */
   
   /* for compatibility with the previous library, also compute the */
   /* determinant */
   for(j = 0; j < n; j++) {
      d *= ainv[j][j];
   }
   for(j=0,p1=tvec1lf; j < n; j++,p1++) {
      for(i=0,p=tvec1lf; i < n; i++,p++) *p = 0;
      *p1 = 1.; 
      plubksublf(ainv,n,indx,tvec1lf);
      for(i=0,p=tvec1lf; i < n; i++,p++)
	 tmat1lf[i][j] = *p;
   }
   /* now copy the inverse into the called matrix */
   for(i = 0; i < n; i++) 
      for(j=0,p=ainv[i],p1=tmat1lf[i]; j < n; j++,p++,p1++)
	 *p = *p1;
   return(d);
} /* end of inverse */


/*****************************************************************/



void pludcmpf(float **a,float **LUa,int n,int *indx,float *d,float TINY)
{
   int i,imax,j,k;
   float big,dum,sum,temp;
   
   TEMPVEC(tvec1f,n,mv1f,float);
   
   imax = 0;
   *d=1.0;
   for (i=0,pf2=tvec1f;i<n;i++,pf2++) {
      big=0.0;
      for (j=0,pf=LUa[i],pf1=a[i];j<n;j++,pf++,pf1++) {
		 *pf = *pf1;
		 if ((temp=fabs(*pf1)) > big) big=temp;
      }
      if (big  < TINY) pmaterr("Error: singular matrix in pludcmpf");
      *pf2 = 1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
	 pf = LUa[i];
	 sum=pf[j];
	 for (k=0;k<i;k++) sum -= pf[k]*LUa[k][j];
	 pf[j]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) {
	 pf = LUa[i];
	 sum=pf[j];
	 for (k=0;k<j;k++)
	    sum -= pf[k]*LUa[k][j];
	 pf[j]=sum;
	 if ( (dum=tvec1f[i]*fabs(sum)) >= big) {
	    big=dum;
	    imax=i;
	 }
      }
      if (j != imax) {
	 for (k=0,pf=LUa[j];k<n;k++,pf++) {
	    dum=LUa[imax][k];
	    LUa[imax][k] = *pf;
	    *pf = dum;
	 }
	 *d = -(*d);
	 tvec1f[imax]=tvec1f[j];
      }
      indx[j]=imax;
      if (LUa[j][j] == 0) {
	 LUa[j][j] = TINY;
      }
      if (j != n - 1) {
	 dum=1.0/(LUa[j][j]);
	 for (i=j+1;i<n;i++) LUa[i][j] *= dum;
      }
   }
}



void pludcmplf(double **a,double **LUa,int n,int *indx,double *d,double TINY)
{
   int i,imax,j,k;
   double big,dum,sum,temp;
   
   TEMPVEC(tvec1lf,n,mv1lf,double);
   
   imax = 0;
   *d=1.0;
   for (i=0,pd2=tvec1lf;i<n;i++,pd2++) {
      big=0.0;
      for (j=0,pd=LUa[i],pd1=a[i];j<n;j++,pd++,pd1++) {
	 *pd = *pd1;
	 if ((temp=fabs(*pd1)) > big) big=temp;
      }
      if (big  < TINY) pmaterr("Error: singular matrix in pludcmplf");
      *pd2 = 1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
	 pd = LUa[i];
	 sum=pd[j];
	 for (k=0;k<i;k++) sum -= pd[k]*LUa[k][j];
	 pd[j]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) {
	 pd = LUa[i];
	 sum=pd[j];
	 for (k=0;k<j;k++)
	    sum -= pd[k]*LUa[k][j];
	 pd[j]=sum;
	 if ( (dum=tvec1lf[i]*fabs(sum)) >= big) {
	    big=dum;
	    imax=i;
	 }
      }
      if (j != imax) {
	 for (k=0,pd=LUa[j];k<n;k++,pd++) {
	    dum=LUa[imax][k];
	    LUa[imax][k] = *pd;
	    *pd = dum;
	 }
	 *d = -(*d);
	 tvec1lf[imax]=tvec1lf[j];
      }
      indx[j]=imax;
      if (LUa[j][j] == 0) {
	 LUa[j][j] = TINY;
      }
      if (j != n - 1) {
	 dum=1.0/(LUa[j][j]);
	 for (i=j+1;i<n;i++) LUa[i][j] *= dum;
      }
   }
}


/*****************************************************************/
void plubksubf(float ** a,int n,int *indx,float *b)
{
   int i,ii=-1,ip,j;
   float sum;
   
   for (i=0;i<n;i++) {
      pf = a[i];
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii>= 0)
	 for (j=ii;j<=i-1;j++) sum -= pf[j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n-1;i>=0;i--) {
      pf = a[i];
      sum=b[i];
      for (j=i+1;j<n;j++) sum -= pf[j]*b[j];
      b[i]=sum/pf[i];
   }
}

void plubksublf(double ** a,int n,int *indx,double *b)
{
   int i,ii=-1,ip,j;
   double sum;
   
   for (i=0;i<n;i++) {
      pd = a[i];
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii>= 0)
	 for (j=ii;j<=i-1;j++) sum -= pd[j]*b[j];
      else if (sum) ii=i;
      b[i]=sum;
   }
   for (i=n-1;i>=0;i--) {
      sum=b[i];
      pd = a[i];
      for (j=i+1;j<n;j++) sum -= pd[j]*b[j];
      b[i]=sum/pd[i];
   }
}



/*****************************************************************/
/* perform a cholesky factorization of the positive definite symmetric */
/*  matrix A, A = LL' , L a */
/* lower triangular matrix */
/* See Burden and Faires, third ed. p. 351 */
void pcholf(float **A, float **L, int n)
{
   int i,j,k;
   float temp;
   float l;
   
   if(A[0][0] < 0) {
      FILE *fo;
      fprintf(stderr,
	      "Warning: sqrt domain error (0) in Cholesky: %f\n",A[0][0]);
      fo = fout;
      fout = stderr;
      pmatprintf("A",A,n,n);
      fout = fo;
   }
   L[0][0] = sqrt(A[0][0]);	/* step 1 */
   /* clear out the upper triangle */
   for(i = 0; i < n; i++) {
      pf = L[i];
      for(j = i+1; j < n; j++)
	 pf[j] = 0;
   }
   
   l = L[0][0];
   for(j = 1; j < n; j++) {	/* step 2 */
      L[j][0] = A[j][0]/l;
   }
   /* step 3 */
   for(i = 1; i < n-1; i++) {	/* step 4 */
      temp = 0;
      for(k = 0,  pf = L[i]; k <= i-1; k++,pf++) {
	 temp += *pf * *pf;
      }
      temp = A[i][i] - temp;
      if(temp < 0) {
	 FILE *fo;
	 fprintf(stderr,
		 "Warning: sqrt domain error (1) in Cholesky: %f\n",temp);
	 fo = fout;
	 fout = stderr;
	 pmatprintf("A",A,n,n);
	 fout = fo;
      }
      L[i][i] = sqrt(temp);
      
      for(j = i+1; j < n; j++) { /* step 5 */
	 temp = 0;
	 for(k = 0,pf=L[j],pf1=L[i]; k <= i-1; k++,pf++,pf1++) {
	    temp += *pf1 * *pf;
	 }
	 L[j][i] = 1./L[i][i]*(A[j][i] - temp);
      }
   } /* end step 3 */
   
   temp = 0;			/* step 6 */
   for(k = 0,pf = L[n-1]; k < n -1 ; k++,pf++) {
      temp += *pf * *pf;
   }
   temp = A[n-1][n-1] - temp;
   if(temp < 0) {
      fprintf(stderr,"Warning: sqrt domain error in (2) Cholesksy: %f\n",temp);
   }
   L[n-1][n-1] = sqrt(temp);
}



void pchollf(double **A, double **L, int n)
{
   int i,j,k;
   double temp;
   double l;
   
   if(A[0][0] < 0) {
      FILE *fo;
      fprintf(stderr,
	      "Warning: sqrt domain error (0) in Cholesky: %lf\n",A[0][0]);
      fo = fout;
      fout = stderr;
      pmatprintlf("A",A,n,n);
      fout = fo;
   }
   L[0][0] = sqrt(A[0][0]);	/* step 1 */
   for(i = 0; i < n; i++) {
      pd = L[i];
      for(j = i+1; j < n; j++)
	 pd[j] = 0;
   }
   
   l = L[0][0];
   for(j = 1; j < n; j++) {	/* step 2 */
      L[j][0] = A[j][0]/l;
   }
   /* step 3 */
   for(i = 1; i < n-1; i++) {	/* step 4 */
      temp = 0;
      for(k=0,pd=L[i] ; k <= i-1; k++,pd++) {
	 temp += *pd * *pd;
      }
      temp = A[i][i] - temp;
      if(temp < 0) {
	 fprintf(stderr,
		 "Warning: sqrt domain error (1) in Cholesky: %f\n",temp);

      }
      L[i][i] = sqrt(temp);
      for(j = i+1; j < n; j++) { /* step 5 */
	 temp = 0;
	 for(k=0,pd=L[j]; k <= i-1; k++,pd++) {
	    temp += L[j][k]*L[i][k];     /* *pd * *pd; */
	 }
	 L[j][i] = 1./L[i][i]*(A[j][i] - temp);
      }
   } /* end step 3 */
   
   temp = 0;			/* step 6 */
   for(k=0,pd=L[n-1]; k < n-1 ; k++,pd++) {
      temp += *pd * *pd;
   }
   temp = A[n-1][n-1] - temp;
   if(temp < 0) {
      FILE *fo;
      fprintf(stderr,"Warning: sqrt domain error in (2) Cholesksy: %f\n",temp);
      fo = fout;
      fout = stderr;
      pmatprintlf("A",A,n,n);
      fout = fo;
   }
   L[n-1][n-1] = sqrt(temp);
}


/*********************************************************
Computes the  LDL decomposition for a symmetric matrix
 (Algorithm lifted from Golub, 2nd ed. pp. 137-138.)

 Factor A so that A = L*diag(D)*L'

Where:
  A is the matrix.
  L is the lower-triangular component.
  D is the diagonal component.
  n is the dimension of the matrix.
  Return value:
  (none).
  *********************************************************/
void pLDLlf( double **A, double **L, double *D, int n )
{
   int i, k, p;
   
   TEMPVEC(tvec1lf,n,mv1lf,double);
   
   if ( n == 1 ) {
      L[0][0] = 1.0;
      D[0] = A[0][0];
      return;
   }
   tvec1lf[0] = 0.0;
   for ( i = 0; i < n; i++ ) {
      for( k = i,pd = (L[i]+k); k < n; k++,pd++)
	 *pd = (i == k) ? 1.0 : 0.0;
   }
   for ( k = 0; k < n; k++ ) {
      pd = L[k];
      for ( p = 0; p < k; p++ )
	 tvec1lf[p] = D[p] * pd[p];
      D[k] = A[k][k];
      for ( p = 0; p < k; p++ )
	 D[k] -= pd[p] * tvec1lf[p];
      if ( D[k] < THRESH ) {
	 pmaterr("LDL error");
      }
      for ( i = k + 1; i < n; i++ ) {
	 pd = L[i];
	 pd[k] = A[i][k];
	 for (p = 0; p < k; p++ )
	    pd[k] -= pd[p] * tvec1lf[p];
	 pd[k] /= D[k];
      }
   }
}

void pLDLf( float **A,float **L, float *D, int n )
{
   int i, k, p;
   
   TEMPVEC(tvec1f,n,mv1f,float);
   
   if ( n == 1 ) {
      L[0][0] = 1.0;
      D[0] = A[0][0];
      return;
   }
   tvec1f[0] = 0.0;
   for ( i = 0; i < n; i++ ) {
      for( k = i,pf = (L[i]+k); k < n; k++,pf++)
	 *pf = (i == k) ? 1.0 : 0.0;
   }
   for ( k = 0; k < n; k++ ) {
      pf = L[k];
      for ( p = 0; p < k; p++ )
	 tvec1f[p] = D[p] * pf[p];
      D[k] = A[k][k];
      for ( p = 0; p < k; p++ )
	 D[k] -= pf[p] * tvec1f[p];
      if ( D[k] < THRESH ) {
	 pmaterr("LDL error");
      }
      for ( i = k + 1; i < n; i++ ) {
	 pf = L[i];
	 pf[k] = A[i][k];
	 for (p = 0; p < k; p++ )
	    pf[k] -= pf[p] * tvec1f[p];
	 pf[k] /= D[k];
      }
   }
}


/***********************************************************************/
/* SVD factorization: a = u diag(w) v' */

void psvdlf(double **a,	/* matrix to factor */
	       double **u,	/* u */
	       double *w,	/* diagonal matrix of singular values */
	       double **v,	/* v */
	       int m,int n)	/* rows and columns of a */
{
   int flag,i,its,j,jj,k,l,nm;
   double c,f,h,s,x,y,z;
   double anorm=0.0,g=0.0,scale=0.0;
   
   if (m < n) pmaterr("SVD: You must augment A with extra zero rows");
   
   TEMPVEC(tvec1lf,m,mv1lf,double);

   l = 1;
   nm = 0;
   for(i = 0; i < m; i++) 
      for(j = 0,pd = u[i],pd1=a[i]; j < n; j++,pd++,pd1++)
	 *pd = *pd1;
   for (i=0;i<n;i++) {
      l=i+1;
      tvec1lf[i]=scale*g;
      g=s=scale=0.0;
      if (i < m) {
	 for (k=i;k<m;k++) scale += fabs(u[k][i]);
	 if (scale) {
	    for (k=i;k<m;k++) {
	       u[k][i] /= scale;
	       s += u[k][i]*u[k][i];
	    }
	    f=u[i][i];
	    g = -SVDSIGN(sqrt(s),f);
	    h=f*g-s;
	    u[i][i]=f-g;
	    if (i != n-1) {
	       for (j=l;j<n;j++) {
		  for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
		  f=s/h;
		  for (k=i;k<m;k++) u[k][j] += f*u[k][i];
	       }
	    }
	    for (k=i;k<m;k++) u[k][i] *= scale;
	 }
      }
      w[i]=scale*g;
      g=s=scale=0.0;
      if (i < m && i != n-1) {
	 for (k=l;k<n;k++) scale += fabs(u[i][k]);
	 if (scale) {
	    for (k=l;k<n;k++) {
	       u[i][k] /= scale;
	       s += u[i][k]*u[i][k];
	    }
	    f=u[i][l];
	    g = -SVDSIGN(sqrt(s),f);
	    h=f*g-s;
	    u[i][l]=f-g;
	    for (k=l;k<n;k++) tvec1lf[k]=u[i][k]/h;
	    if (i != m-1) {
	       for (j=l;j<m;j++) {
		  for (s=0.0,k=l;k<n;k++) s += u[j][k]*u[i][k];
		  for (k=l;k<n;k++) u[j][k] += s*tvec1lf[k];
	       }
	    }
	    for (k=l;k<n;k++) u[i][k] *= scale;
	 }
      }
      anorm=SVDMAX(anorm,(fabs(w[i])+fabs(tvec1lf[i])));
   }
   for (i=n-1;i>=0;i--) {
      if (i < n-1) {
	 if (g) {
	    for (j=l;j<n;j++)
	       v[j][i]=(u[i][j]/u[i][l])/g;
	    for (j=l;j<n;j++) {
	       for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
	       for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	    }
	 }
	 for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
      }
      v[i][i]=1.0;
      g=tvec1lf[i];
      l=i;
   }
   for (i=n-1;i>=0;i--) {
      l=i+1;
      g=w[i];
      if (i < n-1)
	 for (j=l;j<n;j++) u[i][j]=0.0;
      if (g) {
	 g=1.0/g;
	 if (i != n-1) {
	    for (j=l;j<n;j++) {
	       for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
	       f=(s/u[i][i])*g;
	       for (k=i;k<m;k++) u[k][j] += f*u[k][i];
	    }
	 }
	 for (j=i;j<m;j++) u[j][i] *= g;
      } else {
	 for (j=i;j<m;j++) u[j][i]=0.0;
      }
      ++u[i][i];
   }
   for (k=n-1;k>=0;k--) {
      for (its=1;its<=30;its++) {
	 flag=1;
	 for (l=k;l>=0;l--) {
	    nm=l-1;
	    if (fabs(tvec1lf[l])+anorm == anorm) {
	       flag=0;
	       break;
	    }
	    if (fabs(w[nm])+anorm == anorm) break;
	 }
	 if (flag) {
	    c=0.0;
	    s=1.0;
	    for (i=l;i<=k;i++) {
	       f=s*tvec1lf[i];
	       if (fabs(f)+anorm != anorm) {
		  g=w[i];
		  h=SVDPYTHAG(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s=(-f*h);
		  for (j=0;j<m;j++) {
		     y=u[j][nm];
		     z=u[j][i];
		     u[j][nm]=y*c+z*s;
		     u[j][i]=z*c-y*s;
		  }
	       }
	    }
	 }
	 z=w[k];
	 if (l == k) {
	    if (z < 0.0) {
	       w[k] = -z;
	       for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
	    }
	    break;
	 }
	 if (its == 30) pmaterr("No convergence in 30 SVDCMP iterations");
	 x=w[l];
	 nm=k-1;
	 y=w[nm];
	 g=tvec1lf[nm];
	 h=tvec1lf[k];
	 f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	 g=SVDPYTHAG(f,1.0);
	 f=((x-z)*(x+z)+h*((y/(f+SVDSIGN(g,f)))-h))/x;
	 c=s=1.0;
	 for (j=l;j<=nm;j++) {
	    i=j+1;
	    g=tvec1lf[i];
	    y=w[i];
	    h=s*g;
	    g=c*g;
	    z=SVDPYTHAG(f,h);
	    tvec1lf[j]=z;
	    c=f/z;
	    s=h/z;
	    f=x*c+g*s;
	    g=g*c-x*s;
	    h=y*s;
	    y=y*c;
	    for (jj=0;jj<n;jj++) {
	       x=v[jj][j];
	       z=v[jj][i];
	       v[jj][j]=x*c+z*s;
	       v[jj][i]=z*c-x*s;
	    }
	    z=SVDPYTHAG(f,h);
	    w[j]=z;
	    if (z) {
	       z=1.0/z;
	       c=f*z;
	       s=h*z;
	    }
	    f=(c*g)+(s*y);
	    x=(c*y)-(s*g);
	    for (jj=0;jj<m;jj++) {
	       y=u[jj][j];
	       z=u[jj][i];
	       u[jj][j]=y*c+z*s;
	       u[jj][i]=z*c-y*s;
	    }
	 }
	 tvec1lf[l]=0.0;
	 tvec1lf[k]=f;
	 w[k]=x;
      }
   }
}




/* factor a = u diag(w) v' */
void psvdf(float **a,	/* matrix to factor */
	      float **u,	/* u */
	      float *w,		/* diagonal matrix of singular values */
	      float **v,	/* v */
	      int m,int n)	/* rows and columns of a */
{
   int flag,i,its,j,jj,k,l,nm;
   float c,f,h,s,x,y,z;
   float anorm=0.0,g=0.0,scale=0.0;
   
   if (m < n) pmaterr("SVDCMP: You must augment A with extra zero rows");

   TEMPVEC(tvec1f,m,mv1f,float);

   l = 1;
   nm = 0;
   for(i = 0; i < m; i++) for(j = 0; j < n; j++) u[i][j] = a[i][j];
   for (i=0;i<n;i++) {
      l=i+1;
      tvec1f[i]=scale*g;
      g=s=scale=0.0;
      if (i < m) {
	 for (k=i;k<m;k++) scale += fabs(u[k][i]);
	 if (scale) {
	    for (k=i;k<m;k++) {
	       u[k][i] /= scale;
	       s += u[k][i]*u[k][i];
	    }
	    f=u[i][i];
	    g = -SVDSIGN(sqrt(s),f);
	    h=f*g-s;
	    u[i][i]=f-g;
	    if (i != n-1) {
	       for (j=l;j<n;j++) {
		  for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
		  f=s/h;
		  for (k=i;k<m;k++) u[k][j] += f*u[k][i];
	       }
	    }
	    for (k=i;k<m;k++) u[k][i] *= scale;
	 }
      }
      w[i]=scale*g;
      g=s=scale=0.0;
      if (i < m && i != n-1) {
	 for (k=l;k<n;k++) scale += fabs(u[i][k]);
	 if (scale) {
	    for (k=l;k<n;k++) {
	       u[i][k] /= scale;
	       s += u[i][k]*u[i][k];
	    }
	    f=u[i][l];
	    g = -SVDSIGN(sqrt(s),f);
	    h=f*g-s;
	    u[i][l]=f-g;
	    for (k=l;k<n;k++) tvec1f[k]=u[i][k]/h;
	    if (i != m-1) {
	       for (j=l;j<m;j++) {
		  for (s=0.0,k=l;k<n;k++) s += u[j][k]*u[i][k];
		  for (k=l;k<n;k++) u[j][k] += s*tvec1f[k];
	       }
	    }
	    for (k=l;k<n;k++) u[i][k] *= scale;
	 }
      }
      anorm=SVDMAX(anorm,(fabs(w[i])+fabs(tvec1f[i])));
   }
   for (i=n-1;i>=0;i--) {
      if (i < n-1) {
	 if (g) {
	    for (j=l;j<n;j++)
	       v[j][i]=(u[i][j]/u[i][l])/g;
	    for (j=l;j<n;j++) {
	       for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
	       for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	    }
	 }
	 for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
      }
      v[i][i]=1.0;
      g=tvec1f[i];
      l=i;
   }
   for (i=n-1;i>=0;i--) {
      l=i+1;
      g=w[i];
      if (i < n-1)
	 for (j=l;j<n;j++) u[i][j]=0.0;
      if (g) {
	 g=1.0/g;
	 if (i != n-1) {
	    for (j=l;j<n;j++) {
	       for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
	       f=(s/u[i][i])*g;
	       for (k=i;k<m;k++) u[k][j] += f*u[k][i];
	    }
	 }
	 for (j=i;j<m;j++) u[j][i] *= g;
      } else {
	 for (j=i;j<m;j++) u[j][i]=0.0;
      }
      ++u[i][i];
   }
   for (k=n-1;k>=0;k--) {
      for (its=1;its<=30;its++) {
	 flag=1;
	 for (l=k;l>=0;l--) {
	    nm=l-1;
	    if (fabs(tvec1f[l])+anorm == anorm) {
	       flag=0;
	       break;
	    }
	    if (fabs(w[nm])+anorm == anorm) break;
	 }
	 if (flag) {
	    c=0.0;
	    s=1.0;
	    for (i=l;i<=k;i++) {
	       f=s*tvec1f[i];
	       if (fabs(f)+anorm != anorm) {
		  g=w[i];
		  h=SVDPYTHAG(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s=(-f*h);
		  for (j=0;j<m;j++) {
		     y=u[j][nm];
		     z=u[j][i];
		     u[j][nm]=y*c+z*s;
		     u[j][i]=z*c-y*s;
		  }
	       }
	    }
	 }
	 z=w[k];
	 if (l == k) {
	    if (z < 0.0) {
	       w[k] = -z;
	       for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
	    }
	    break;
	 }
	 if (its == 30) pmaterr("No convergence in 30 SVDCMP iterations");
	 x=w[l];
	 nm=k-1;
	 y=w[nm];
	 g=tvec1f[nm];
	 h=tvec1f[k];
	 f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	 g=SVDPYTHAG(f,1.0);
	 f=((x-z)*(x+z)+h*((y/(f+SVDSIGN(g,f)))-h))/x;
	 c=s=1.0;
	 for (j=l;j<=nm;j++) {
	    i=j+1;
	    g=tvec1f[i];
	    y=w[i];
	    h=s*g;
	    g=c*g;
	    z=SVDPYTHAG(f,h);
	    tvec1f[j]=z;
	    c=f/z;
	    s=h/z;
	    f=x*c+g*s;
	    g=g*c-x*s;
	    h=y*s;
	    y=y*c;
	    for (jj=0;jj<n;jj++) {
	       x=v[jj][j];
	       z=v[jj][i];
	       v[jj][j]=x*c+z*s;
	       v[jj][i]=z*c-x*s;
	    }
	    z=SVDPYTHAG(f,h);
	    w[j]=z;
	    if (z) {
	       z=1.0/z;
	       c=f*z;
	       s=h*z;
	    }
	    f=(c*g)+(s*y);
	    x=(c*y)-(s*g);
	    for (jj=0;jj<m;jj++) {
	       y=u[jj][j];
	       z=u[jj][i];
	       u[jj][j]=y*c+z*s;
	       u[jj][i]=z*c-y*s;
	    }
	 }
	 tvec1f[l]=0.0;
	 tvec1f[k]=f;
	 w[k]=x;
      }
   }
}


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
  This routine uses the QR method of finding eigenvalues.

  Return value:
  The number of iterations required.              
  Q, Lambda are modified as return values also.
  *********************************************************/
int peigenlf( double **A, double **Q, double *Lambda, int dim)
{
   int i, i1,j, k, found;
   double sum;
   double **save;  /* switch pointers to avoid copying time */
   int maxits = 30;
   double epsilon = 1e-10;

   TEMPMAT(tmat1lf,dim,dim,mm1lf,nm1lf,double);
   TEMPMAT(tmat2lf,dim,dim,mm2lf,nm2lf,double);
   TEMPMAT(tmat3lf,dim,dim,mm3lf,nm3lf,double);
   TEMPMAT(tmat4lf,dim,dim,mm4lf,nm4lf,double);
   TEMPMAT(tmat5lf,dim,dim,mm5lf,nm5lf,double); 


   /* tmat3 acts like Lambda */
   /* tmat5 acts like Qold */
   for ( i = 0; i < dim; i++ ) {
      for ( j = 0; j < dim; j++ ) {
	 tmat3lf[i][j] = A[i][j];
	 tmat5lf[i][j] = (i == j) ? 1.0:0.0; 
      }
   }
   for ( i = 0; i < maxits; i++ ) {
      houseqrlf(tmat3lf,tmat1lf,dim,dim);  /* tmat1lf has Q and R data in it */
      findQlf(tmat1lf,tmat2lf,dim,dim);   /* tmat2 has Q */

      /* form tmat4 = R*Q */
      for(i1 = 0; i1 < dim; i1++) {
	 for(j = 0; j < dim; j++) {
	    sum = 0;
	    for(k = i1; k < dim; k++) {
	       sum += tmat1lf[i1][k]*tmat2lf[k][j];
	    }
	    tmat4lf[i1][j] = sum;
	 }
      }

      found = 1;
      for( j = 0; j < dim; j++ ) {
	 if ( !m_foundlf_( tmat4lf[j][j], tmat3lf[j][j], epsilon ) ) {
	    found = 0;
	    break;
	 }
      }
      save = tmat3lf;
      /* RQ -> tmat3 */
      tmat3lf = tmat4lf;
      tmat4lf = save;


      /* RQ -> tmat3 */
      /* tmat4 = Qold * Q */
      pmatmultlf(tmat5lf,tmat2lf,tmat4lf,dim,dim,dim);
      /* tmat4 -> Q */
      save = tmat5lf;
      tmat5lf = tmat4lf;
      tmat4lf = save;
      if ( found ) {
	 break;
      }
   }
   for(j = 0; j < dim; j++) {
      Lambda[j] = tmat3lf[j][j];
      for(i = 0; i < dim; i++) {
	 Q[i][j] = tmat5lf[i][j];
      }
   }
   return(i);
}


int peigenf( float **A, float **Q, float *Lambda, int dim)

{
   int i, i1,j, k, found;
   float sum;
   float **save;  /* switch pointers to avoid copying time */
   int maxits = 30;
   float epsilon = 1e-10;

   TEMPMAT(tmat1f,dim,dim,mm1f,nm1f,float);
   TEMPMAT(tmat2f,dim,dim,mm2f,nm2f,float);
   TEMPMAT(tmat3f,dim,dim,mm3f,nm3f,float);
   TEMPMAT(tmat4f,dim,dim,mm4f,nm4f,float);
   TEMPMAT(tmat5f,dim,dim,mm5f,nm5f,float); 


   /* tmat3 acts like Lambda */
   /* tmat5 acts like Qold */
   for ( i = 0; i < dim; i++ ) {
      for ( j = 0; j < dim; j++ ) {
	 tmat3f[i][j] = A[i][j];
	 tmat5f[i][j] = (i == j) ? 1.0:0.0; 
      }
   }
   for ( i = 0; i < maxits; i++ ) {
      houseqrf(tmat3f,tmat1f,dim,dim);  /* tmat1f has Q and R data in it */
      findQf(tmat1f,tmat2f,dim,dim);   /* tmat2 has Q */

      /* form tmat4 = R*Q */
      for(i1 = 0; i1 < dim; i1++) {
	 for(j = 0; j < dim; j++) {
	    sum = 0;
	    for(k = i1; k < dim; k++) {
	       sum += tmat1f[i1][k]*tmat2f[k][j];
	    }
	    tmat4f[i1][j] = sum;
	 }
      }

      found = 1;
      for( j = 0; j < dim; j++ ) {
	 if ( !m_foundf_( tmat4f[j][j], tmat3f[j][j], epsilon ) ) {
	    found = 0;
	    break;
	 }
      }
      save = tmat3f;
      /* RQ -> tmat3 */
      tmat3f = tmat4f;
      tmat4f = save;


      /* RQ -> tmat3 */
      /* tmat4 = Qold * Q */
      pmatmultf(tmat5f,tmat2f,tmat4f,dim,dim,dim);
      /* tmat4 -> Q */
      save = tmat5f;
      tmat5f = tmat4f;
      tmat4f = save;
      if ( found ) {
	 break;
      }
   }
   for(j = 0; j < dim; j++) {
      Lambda[j] = tmat3f[j][j];
      for(i = 0; i < dim; i++) {
	 Q[i][j] = tmat5f[i][j];
      }
   }
   return(i);
}


/* compute the compact form of the householder transformation 
   from GvL, 2nd ed., p. 212 
*/
#define HOUSESIGN(x) ((x) >= 0 ? 1 : -1)

void houseqrlf(double **A, double **QR, int m, int n)
{
   int i,j,k,i1,j1,n1;
   double mu,a,beta;
   double *v;
   double vpv,x1,sum;
   TEMPVEC(tvec1lf,m,mv1lf,double);
   TEMPVEC(tvec2lf,m,mv2lf,double);

   /* copy over the answer */
   for(i = 0; i < m; i++) {
      for(j = 0,pd = QR[i], pd1 = A[i]; j < n; j++) {
	 *pd++ = *pd1++;
      }
   }
   if(m < n) pmaterr("Error: rows must exceed columns in houseqr");
   for(j = 0; j < n; j++) { /* work over each column */

      /* HOUSE:  find the householder vector for A(j:(m-1),j) --
         that is, it will compact the energy in column j into the
         jth position.  GvL, p. 196 */
      n1 = m-j;
      vpv = 1;          /* use to accumulate v'v */
      v = (tvec1lf+j);
      mu = 0;
      for(i = j; i < m; i++) {
	 a = QR[i][j];
	 mu += a*a;
	 tvec1lf[i] = a;
      }
      mu = sqrt(mu);
      if(mu != 0) {
	 x1 = tvec1lf[j];
	 beta = x1 + HOUSESIGN(x1)*mu;

	 for(i = 1; i < n1; i++) {
	    v[i] /= beta;
	    vpv += v[i]*v[i];
	 }
      }
      v[0] = 1;

      /* ROW.HOUSE(A(j:m,j:n)) overwrite A(j:m) with PA(j:m,j:n),
	 where P = I - 2vv'/v'v.  Note that this only deals with the
         last rows and columns.  GvL, p. 197 */
      beta = -2/vpv;
      for(i1 = j; i1 < n; i1++) {
	 sum = 0;
	 for(j1 = j; j1 < m; j1++) {
	    sum += beta*QR[j1][i1]*tvec1lf[j1];
	 }
	 tvec2lf[i1] = sum;
      }
      QR[j][j] += tvec1lf[j]*tvec2lf[j];
      for(j1 = j; j1 < m; j1++) {   /* operate on submatrix A(j:m-1,j:n-1) */
	 for(i1 = j+1; i1 < n; i1++){  /* < make j+1 to save time */
	    QR[j1][i1] += tvec1lf[j1]*tvec2lf[i1];
	 }
      }
      /* Save the householder vectors */
      if(j < m) {
	 for(i1 = j+1; i1 < m; i1++) {
	    QR[i1][j] = tvec1lf[i1];
	 }
      }
   }
}      



void houseqrf(float **A, float **QR, int m, int n)
{
   int i,j,k,i1,j1,n1;
   float mu,a,beta;
   float *v;
   float vpv,x1,sum;
   TEMPVEC(tvec1f,m,mv1f,float);
   TEMPVEC(tvec2f,m,mv2f,float);

   /* copy over the answer */
   for(i = 0; i < m; i++) {
      for(j = 0,pf = QR[i], pf1 = A[i]; j < n; j++) {
	 *pf++ = *pf1++;
      }
   }
   if(m < n) pmaterr("Error: rows must exceed columns in houseqr");
   for(j = 0; j < n; j++) { /* work over each column */

      /* HOUSE:  find the householder vector for A(j:(m-1),j) --
         that is, it will compact the energy in column j into the
         jth position.  GvL, p. 196 */
      n1 = m-j;
      vpv = 1;          /* use to accumulate v'v */
      v = (tvec1f+j);
      mu = 0;
      for(i = j; i < m; i++) {
	 a = QR[i][j];
	 mu += a*a;
	 tvec1f[i] = a;
      }
      mu = sqrt(mu);
      if(mu != 0) {
	 x1 = tvec1f[j];
	 beta = x1 + HOUSESIGN(x1)*mu;

	 for(i = 1; i < n1; i++) {
	    v[i] /= beta;
	    vpv += v[i]*v[i];
	 }
      }
      v[0] = 1;

      /* ROW.HOUSE(A(j:m,j:n)) overwrite A(j:m) with PA(j:m,j:n),
	 where P = I - 2vv'/v'v.  Note that this only deals with the
         last rows and columns.  GvL, p. 197 */
      beta = -2/vpv;
      for(i1 = j; i1 < n; i1++) {
	 sum = 0;
	 for(j1 = j; j1 < m; j1++) {
	    sum += beta*QR[j1][i1]*tvec1f[j1];
	 }
	 tvec2f[i1] = sum;
      }
      QR[j][j] += tvec1f[j]*tvec2f[j];
      for(j1 = j; j1 < m; j1++) {   /* operate on submatrix A(j:m-1,j:n-1) */
	 for(i1 = j+1; i1 < n; i1++){  /* < make j+1 to save time */
	    QR[j1][i1] += tvec1f[j1]*tvec2f[i1];
	 }
      }
      /* Save the householder vectors */
      if(j < m) {
	 for(i1 = j+1; i1 < m; i1++) {
	    QR[i1][j] = tvec1f[i1];
	 }
      }
   }
}      





void findQandRlf(double **QR, double **Q, double **R, int m, int n)
/* given the compact householder representation in QR,
   return Q and R */
{
   int i,j,i1;
   double vi,vpv,beta,sum;
   TEMPVEC(tvec1lf,m,mv1lf,double);  /* w is tvec1lf */
   TEMPVEC(tvec2lf,m,mv2lf,double);  /* v is tvec2lf */

   for(i = 0;i < m; i++) {
      for(j = 0; j < m; j++) {
	 Q[i][j] = 0;
      }
      Q[i][i] = 1;
   }
   for(i = 0;i < m; i++) {
      for(j = 0; j < n; j++) {
	 R[i][j] = 0;
      }
      for(j = i; j < n; j++) {
	 R[i][j] = QR[i][j];
      }
   }
   /* get the first block directly from column [n-1] */
   /* find v'v.  In this first loop, tvect1lf is v */
   vpv = 1;
   for(i = n; i < m; i++) {
      vi = tvec2lf[i] = QR[i][n-1];
      vpv += vi*vi;
   }
   beta = -2/vpv;
   Q[n-1][n-1] = 1+beta;
   for(j = n; j < m; j++) {
      Q[n-1][j] = vi = beta*tvec2lf[j];
      /* Q[j][n-1] = vi; */
   }
   for(i = n; i < m; i++) {
      pd = Q[i];
      for(j = i; j < m; j++) {
	 pd[j] += beta*tvec2lf[j]*tvec2lf[i];
      }
   }
   for(i = n; i < m; i++) { /* fill in the rest by symmetry */
      pd = Q[i];
      for(j = n-1; j < i; j++) {
	 pd[j] = Q[j][i];
      }
   }
   for(i = n-2; i >= 0; i--) {
      /* v[0] = 1,  v[j-i] = QR[j][i], j=i+1:m */
      vpv = 1;
      tvec2lf[i] = 1;
      for(j = i+1; j < m; j++) {
	 tvec2lf[j] = QR[j][i];
	 vpv += tvec2lf[j]*tvec2lf[j];
      }
      beta = -2/vpv;
      /* compute w = beta * A' * v */
      sum = beta;
      for(j = i; j < m; j++) {
	 for(i1 = i+1; i1 < m; i1++) {
	    sum += beta * Q[i1][j]*tvec2lf[i1];
	 }
	 tvec1lf[j] = sum;
	 sum = 0;
      }
      for(j = i; j < m; j++) {
	 for(i1 = i; i1 < m; i1++) {
	    Q[j][i1] += tvec2lf[j]*tvec1lf[i1];
	 }
      }
   }
}

void findQandRf(float **QR, float **Q, float **R, int m, int n)
/* given the compact householder representation in QR,
   return Q and R */
{
   int i,j,i1;
   float vi,vpv,beta,sum;
   TEMPVEC(tvec1f,m,mv1f,float);  /* w is tvec1f */
   TEMPVEC(tvec2f,m,mv2f,float);  /* v is tvec2f */

   for(i = 0;i < m; i++) {
      for(j = 0; j < m; j++) {
	 Q[i][j] = 0;
      }
      Q[i][i] = 1;
   }
   for(i = 0;i < m; i++) {
      for(j = 0; j < n; j++) {
	 R[i][j] = 0;
      }
      for(j = i; j < n; j++) {
	 R[i][j] = QR[i][j];
      }
   }
   /* get the first block directly from column [n-1] */
   /* find v'v.  In this first loop, tvect1f is v */
   vpv = 1;
   for(i = n; i < m; i++) {
      vi = tvec2f[i] = QR[i][n-1];
      vpv += vi*vi;
   }
   beta = -2/vpv;
   Q[n-1][n-1] = 1+beta;
   for(j = n; j < m; j++) {
      Q[n-1][j] = vi = beta*tvec2f[j];
      /* Q[j][n-1] = vi; */
   }
   for(i = n; i < m; i++) {
      pf = Q[i];
      for(j = i; j < m; j++) {
	 pf[j] += beta*tvec2f[j]*tvec2f[i];
      }
   }
   for(i = n; i < m; i++) { /* fill in the rest by symmetry */
      pf = Q[i];
      for(j = n-1; j < i; j++) {
	 pf[j] = Q[j][i];
      }
   }
   for(i = n-2; i >= 0; i--) {
      /* v[0] = 1,  v[j-i] = QR[j][i], j=i+1:m */
      vpv = 1;
      tvec2f[i] = 1;
      for(j = i+1; j < m; j++) {
	 tvec2f[j] = QR[j][i];
	 vpv += tvec2f[j]*tvec2f[j];
      }
      beta = -2/vpv;
      /* compute w = beta * A' * v */
      sum = beta;
      for(j = i; j < m; j++) {
	 for(i1 = i+1; i1 < m; i1++) {
	    sum += beta * Q[i1][j]*tvec2f[i1];
	 }
	 tvec1f[j] = sum;
	 sum = 0;
      }
      for(j = i; j < m; j++) {
	 for(i1 = i; i1 < m; i1++) {
	    Q[j][i1] += tvec2f[j]*tvec1f[i1];
	 }
      }
   }
}


void findQlf(double **QR, double **Q, int m, int n)
/* given the compact householder representation in QR,
   return Q */
{
   int i,j,i1;
   double vi,vpv,beta,sum;
   TEMPVEC(tvec1lf,m,mv1lf,double);  /* w is tvec1lf */
   TEMPVEC(tvec2lf,m,mv2lf,double);  /* v is tvec2lf */

   for(i = 0;i < m; i++) {
      for(j = 0; j < m; j++) {
	 Q[i][j] = 0;
      }
      Q[i][i] = 1;
   }
   /* get the first block directly from column [n-1] */
   /* find v'v.  In this first loop, tvect1lf is v */
   vpv = 1;
   for(i = n; i < m; i++) {
      vi = tvec2lf[i] = QR[i][n-1];
      vpv += vi*vi;
   }
   beta = -2/vpv;
   Q[n-1][n-1] = 1+beta;
   for(j = n; j < m; j++) {
      Q[n-1][j] = vi = beta*tvec2lf[j];
      /* Q[j][n-1] = vi; */
   }
   for(i = n; i < m; i++) {
      pd = Q[i];
      for(j = i; j < m; j++) {
	 pd[j] += beta*tvec2lf[j]*tvec2lf[i];
      }
   }
   for(i = n; i < m; i++) { /* fill in the rest by symmetry */
      pd = Q[i];
      for(j = n-1; j < i; j++) {
	 pd[j] = Q[j][i];
      }
   }
   for(i = n-2; i >= 0; i--) {
      /* v[0] = 1,  v[j-i] = QR[j][i], j=i+1:m */
      vpv = 1;
      tvec2lf[i] = 1;
      for(j = i+1; j < m; j++) {
	 tvec2lf[j] = QR[j][i];
	 vpv += tvec2lf[j]*tvec2lf[j];
      }
      beta = -2/vpv;
      /* compute w = beta * A' * v */
      sum = beta;
      for(j = i; j < m; j++) {
	 for(i1 = i+1; i1 < m; i1++) {
	    sum += beta * Q[i1][j]*tvec2lf[i1];
	 }
	 tvec1lf[j] = sum;
	 sum = 0;
      }
      for(j = i; j < m; j++) {
	 for(i1 = i; i1 < m; i1++) {
	    Q[j][i1] += tvec2lf[j]*tvec1lf[i1];
	 }
      }
   }
}


void pmatpsolve1lf(double **A, double *x, double *b, int m, int n)
/* Least-squares Solve Ax = b, where A is overdetermined, using 
   a QR factorization 
*/
{
   TEMPMAT(tmat1lf,m,n,mm1lf,nm1lf,double);  /* the QR factorization */
   TEMPVEC(tvec1lf,m,mv1lf,double);   /* the b vector (to save it) */
   houseqrlf(A,tmat1lf,m,n);
   pveccopylf(b,tvec1lf,m);
   applyQlf(tmat1lf,tvec1lf,m,n);
   ptrisolupperlf(tmat1lf, x, tvec1lf, n);
}


void pmatpresolve1lf(double *x, double *b, int m, int n)
/*Least-squares Solve Ax = b, where A is the same as the last call to 
  pmatpsolve1  (and tmat1 has not been modified) */
{
   TEMPVEC(tvec1lf,m,mv1lf,double);   /* the b vector (to save it) */
   pveccopylf(b,tvec1lf,m);
   applyQlf(tmat1lf,tvec1lf,m,n);
   ptrisolupperlf(tmat1lf, x, tvec1lf, n);
}

   


void applyQlf(double **QR, double *b, int m, int n)
/* Given the compact householder representation in QR, compute
   Q' b, replacing the result in b
*/
{
   int j,j1;
   double w,beta;
   double vpv;

   for(j = 0; j < n; j++) {
      /* apply the jth Householder transformation */
      vpv = 1;
      for(j1 = j+1; j1 < m; j1++) {
	 vpv += QR[j1][j]*QR[j1][j];
      }
      /* ROW.HOUSE(A(j:m,j:n)) overwrite A(j:m) with PA(j:m,j:n),
	 where P = I - 2vv'/v'v.  Note that this only deals with the
         last rows and columns.  GvL, p. 197 */
      beta = -2/vpv;
      w = b[j];
      for(j1 = j+1; j1 < m; j1++) {
	 w += b[j1]*QR[j1][j];
      }
      w = w*beta;
      b[j] += w;
      for(j1 = j+1; j1 < m; j1++) {   /* operate on submatrix A(j:m-1,j:n-1) */
	 b[j1] += QR[j1][j]*w;
      }
   }
}


void ptrisolupperlf(double **R, double *x, double *b, int n) 
/* solve the triangular system Rx = b, where R is nxn upper triangular,
   assuming that R[i][i] != 0.
*/
{
   int i, j;
   double sum;
   for(i = n-1; i>= 0; i--) {
      sum = 0;
      for(j = i+1; j < n; j++) {
	 sum += R[i][j]*x[j];
      }
      x[i] = (b[i] - sum)/R[i][i];
   }
}      
   

void pmatpsolve1f(float **A, float *x, float *b, int m, int n)
/* Least-squares Solve Ax = b, where A is overdetermined, using 
   a QR factorization 
*/
{
   TEMPMAT(tmat1f,m,n,mm1f,nm1f,float);  /* the QR factorization */
   TEMPVEC(tvec1f,m,mv1f,float);   /* the b vector (to save it) */
   houseqrf(A,tmat1f,m,n);
   pveccopyf(b,tvec1f,m);
   applyQf(tmat1f,tvec1f,m,n);
   ptrisolupperf(tmat1f, x, tvec1f, n);
}


void pmatpresolve1f(float *x, float *b, int m, int n)
/*Least-squares Solve Ax = b, where A is the same as the last call to 
  pmatpsolve1  (and tmat1 has not been modified) */
{
   TEMPVEC(tvec1f,m,mv1f,float);   /* the b vector (to save it) */
   pveccopyf(b,tvec1f,m);
   applyQf(tmat1f,tvec1f,m,n);
   ptrisolupperf(tmat1f, x, tvec1f, n);
}

   


void applyQf(float **QR, float *b, int m, int n)
/* Given the compact householder representation in QR, compute
   Q' b, replacing the result in b
*/
{
   int j,j1;
   double w,beta;
   double vpv;

   for(j = 0; j < n; j++) {
      /* apply the jth Householder transformation */
      vpv = 1;
      for(j1 = j+1; j1 < m; j1++) {
	 vpv += QR[j1][j]*QR[j1][j];
      }
      /* ROW.HOUSE(A(j:m,j:n)) overwrite A(j:m) with PA(j:m,j:n),
	 where P = I - 2vv'/v'v.  Note that this only deals with the
         last rows and columns.  GvL, p. 197 */
      beta = -2/vpv;
      w = b[j];
      for(j1 = j+1; j1 < m; j1++) {
	 w += b[j1]*QR[j1][j];
      }
      w = w*beta;
      b[j] += w;
      for(j1 = j+1; j1 < m; j1++) {   /* operate on submatrix A(j:m-1,j:n-1) */
	 b[j1] += QR[j1][j]*w;
      }
   }
}


void ptrisolupperf(float **R, float *x, float *b, int n) 
/* solve the triangular system Rx = b, where R is nxn upper triangular,
   assuming that R[i][i] != 0.
*/
{
   int i, j;
   float sum;
   for(i = n-1; i>= 0; i--) {
      sum = 0;
      for(j = i+1; j < n; j++) {
	 sum += R[i][j]*x[j];
      }
      x[i] = (b[i] - sum)/R[i][i];
   }
}      
   

	 
void findQf(float **QR, float **Q, int m, int n)
/* given the compact householder representation in QR,
   return Q */
{
   int i,j,i1;
   float vi,vpv,beta,sum;
   TEMPVEC(tvec1f,m,mv1f,float);  /* w is tvec1f */
   TEMPVEC(tvec2f,m,mv2f,float);  /* v is tvec2f */

   for(i = 0;i < m; i++) {
      for(j = 0; j < m; j++) {
	 Q[i][j] = 0;
      }
      Q[i][i] = 1;
   }
   /* get the first block directly from column [n-1] */
   /* find v'v.  In this first loop, tvect1f is v */
   vpv = 1;
   for(i = n; i < m; i++) {
      vi = tvec2f[i] = QR[i][n-1];
      vpv += vi*vi;
   }
   beta = -2/vpv;
   Q[n-1][n-1] = 1+beta;
   for(j = n; j < m; j++) {
      Q[n-1][j] = vi = beta*tvec2f[j];
      /* Q[j][n-1] = vi; */
   }
   for(i = n; i < m; i++) {
      pf = Q[i];
      for(j = i; j < m; j++) {
	 pf[j] += beta*tvec2f[j]*tvec2f[i];
      }
   }
   for(i = n; i < m; i++) { /* fill in the rest by symmetry */
      pf = Q[i];
      for(j = n-1; j < i; j++) {
	 pf[j] = Q[j][i];
      }
   }
   for(i = n-2; i >= 0; i--) {
      /* v[0] = 1,  v[j-i] = QR[j][i], j=i+1:m */
      vpv = 1;
      tvec2f[i] = 1;
      for(j = i+1; j < m; j++) {
	 tvec2f[j] = QR[j][i];
	 vpv += tvec2f[j]*tvec2f[j];
      }
      beta = -2/vpv;
      /* compute w = beta * A' * v */
      sum = beta;
      for(j = i; j < m; j++) {
	 for(i1 = i+1; i1 < m; i1++) {
	    sum += beta * Q[i1][j]*tvec2f[i1];
	 }
	 tvec1f[j] = sum;
	 sum = 0;
      }
      for(j = i; j < m; j++) {
	 for(i1 = i; i1 < m; i1++) {
	    Q[j][i1] += tvec2f[j]*tvec1f[i1];
	 }
      }
   }
}


/*********************************************************
  A utility for eigen().  This function calculates the relative
  error of a and b given the tolerance or error, e.  
  Return value:
  if a == b (within tolerance), TRUE is returned, else
  FALSE is returned.
  *********************************************************/
static int m_foundlf_( double a, double b, double e)
{
   if (a == 0.0)       /* if a == 0, cannot divide by a for relative error.*/
      if(b == 0.0)
	 return(1);   /*But b == 0 also, so a == b.*/
      else
	 return(0);  /*Not equal. This means it has to be exactly zero. */
   if ( fabs( (a-b)/a ) < e )   /*Within tolerance.*/
      return(1);
   else
      return(0);
}
static int m_foundf_( float a, float b, float e)
{
   if (a == 0.0)       /* if a == 0, cannot divide by a for relative error.*/
      if(b == 0.0)
	 return(1);   /*But b == 0 also, so a == b.*/
      else
	 return(0);  /*Not equal.  This means it has to be exactly zero.*/
   if ( fabs( (a-b)/a ) < e )   /*Within tolerance.*/
      return(1);
   else
      return(0);
}

/*********************************************************/
/* find the eigenvalues and eigenvectors of a symmetric 2x2 matrix */
/* S is the matrix, v has the eigenvectors as columns, and e 
   has the normalized eigenvalues */
void eigen2lf(double **S, double **v, double *e)
{
   double a,b,c,n,lambda,t,tv[2];
   a = S[0][0];
   b = S[0][1];
   c = S[1][1];
   if(b == 0) {
      e[0] = a;
      v[0][0] = 1; v[1][0] = 0;
      e[1] = c;
      v[0][1] = 0;  v[1][1] = 1;
   }
   else {
    lambda = e[0] = 0.5*((a+c)+ sqrt(a*a + c*c + 4*b*b - 2*a*c));
    t = (lambda-a)/b;
    tv[0] = 1;  tv[1] = t;
    n = pvecnorm2lf(tv,2);
    v[0][0] = tv[0]/n;  v[1][0] = tv[1]/n;

    lambda = e[1] = 0.5*((a+c) - sqrt(a*a + c*c + 4*b*b - 2*a*c));
    t = (lambda-a)/b;
    tv[0] = 1;  tv[1] = t;
    n = pvecnorm2lf(tv,2);
    v[0][1] = tv[0]/n;  v[1][1] = tv[1]/n;
   }
}


void eigen2f(float **S, float **v, float *e)
{
   float a,b,c,n,lambda,t,tv[2];
   a = S[0][0];
   b = S[0][1];
   c = S[1][1];
   if(b == 0) {
      e[0] = a;
      v[0][0] = 1; v[1][0] = 0;
      e[1] = c;
      v[0][1] = 0;  v[1][1] = 1;
   }
   else {
    lambda = e[0] = 0.5*((a+c)+ sqrt(a*a + c*c + 4*b*b - 2*a*c));
    t = (lambda-a)/b;
    tv[0] = 1;  tv[1] = t;
    n = pvecnorm2f(tv,2);
    v[0][0] = tv[0]/n;  v[1][0] = tv[1]/n;

    lambda = e[1] = 0.5*((a+c) - sqrt(a*a + c*c + 4*b*b - 2*a*c));
    t = (lambda-a)/b;
    tv[0] = 1;  tv[1] = t;
    n = pvecnorm2f(tv,2);
    v[0][1] = tv[0]/n;  v[1][1] = tv[1]/n;
   }
}



/*********************************************************/
static void findns3lf(int j, double lambda, double a, double b, double c, 
		    double d, double e, double f,
		    double **v);
static void findns3f(int j, float lambda, float a, float b, float c, 
		    float d, float e, float f,
		    float **v);

/*********************************************************/
/* find the eigenvalues and eigenvectors of a symmetric 3x3 matrix */
/* S is the matrix, v has the eigenvectors as columns, and e 
   has the normalized eigenvalues */


void eigen3lf(double **S, double **v, double *ev)
{
   double a,b,c,d,e,f,h,a2,a1,a0,q,r,g,x,y,theta,lambda,det,n;
   double x11,x12,x22,b1,b2,v1,v2,h3,sqrt32,st3;
   double s1ps2, s1ms2j;
   int cas;
   int j;			/* which dimension we are looking at */
   double vt[3];

   a = S[0][0];  b = S[0][1];  c = S[0][2];  d = S[1][1];  e = S[1][2];  
   f = S[2][2];
   a2 = -(a+d+f);
   /* Find the coefficients of the characteristic polynomial */
   /* (Mathematica helped here) */
   a1 = -(b*b + c*c - a*d + e*e - a*f - d*f);
   a0 = c*c*d - 2*b*c*e + a*e*e + b*b*f - a*d*f;
   /* Now find the solution to the cubic; See Abramowitz&Stegun, p. 17 */
   q = a1/3 - a2*a2/9;
   r = (a1*a2-3*a0)/6 - a2*a2*a2/27;
   g = q*q*q+r*r;   /* g will be <= 0 */
   /* complex arithmetic in Z&S has been converted to real for maximum root */
   x = sqrt(-g);  
   h = sqrt(r*r - g);
   theta = atan2(x,r);
   h3 = pow(h,1./3.);
   sqrt32 = sqrt(3.)/2;
   s1ps2 = 2*h3*cos(theta/3);
   s1ms2j = 2*h3*sin(theta/3);
   j = 0;
   lambda = ev[j] = s1ps2 - a2/3;
   findns3lf(j,lambda, a,b,c,d,e,f,v);
   j = 1;
   lambda = ev[j] = -.5*s1ps2 - a2/3 + sqrt32*s1ms2j;
   findns3lf(j,lambda, a,b,c,d,e,f,v);
   j = 2;
   lambda = ev[j] = -.5*s1ps2 - a2/3 - sqrt32*s1ms2j;
   findns3lf(j,lambda, a,b,c,d,e,f,v);
}

void eigen3f(float **S, float **v, float *ev)
{
   float a,b,c,d,e,f,h,a2,a1,a0,q,r,g,x,y,theta,lambda,det,n;
   float x11,x12,x22,b1,b2,v1,v2,h3,sqrt32,st3;
   float s1ps2, s1ms2j;
   int cas;
   int j;			/* which dimension we are looking at */
   float vt[3];

   a = S[0][0];  b = S[0][1];  c = S[0][2];  d = S[1][1];  e = S[1][2];  
   f = S[2][2];
   a2 = -(a+d+f);
   /* Find the coefficients of the characteristic polynomial */
   /* (Mathematica helped here) */
   a1 = -(b*b + c*c - a*d + e*e - a*f - d*f);
   a0 = c*c*d - 2*b*c*e + a*e*e + b*b*f - a*d*f;
   /* Now find the solution to the cubic; See Abramowitz&Stegun, p. 17 */
   q = a1/3 - a2*a2/9;
   r = (a1*a2-3*a0)/6 - a2*a2*a2/27;
   g = q*q*q+r*r;   /* g will be <= 0 */
   /* complex arithmetic in Z&S has been converted to real for maximum root */
   x = sqrt(-g);  
   h = sqrt(r*r - g);
   theta = atan2(x,r);
   h3 = pow(h,1./3.);
   sqrt32 = sqrt(3.)/2;
   s1ps2 = 2*h3*cos(theta/3);
   s1ms2j = 2*h3*sin(theta/3);
   j = 0;
   lambda = ev[j] = s1ps2 - a2/3;
   findns3f(j,lambda, a,b,c,d,e,f,v);
   j = 1;
   lambda = ev[j] = -.5*s1ps2 - a2/3 + sqrt32*s1ms2j;
   findns3f(j,lambda, a,b,c,d,e,f,v);
   j = 2;
   lambda = ev[j] = -.5*s1ps2 - a2/3 - sqrt32*s1ms2j;
   findns3f(j,lambda, a,b,c,d,e,f,v);
}


/* findns3 --- find the nullspace vector for a 3x3 matrix with a nullspace */
static void findns3lf(int j, double lambda, double a, double b, double c, 
		    double d, double e, double f,
		       double **v)
{
   double x11, x12, x22, det, b1, b2,n, v1, v2;
   int cas;
   double vt[3];

   /* Now find the eigenvector */
   a = a-lambda;  d = d-lambda; f = f-lambda;
   /* Assume first that eigenvector has a component in [1 0 0] direction */
   /* Set up matrix corresponding to setting v1=1: */
   x11 = (b*b + d*d + e*e);
   x12 = (b*c + d*e + e*f);
   x22 = (c*c + e*e + f*f);
   /* check that this works */
   det = x11*x22 - x12*x12;
   b1 = -(a*b + b*d + c*e);
   b2 = -(a*c + b*e + c*f);
   cas = 1;
   if(fabs(det) < 1e-10) { /* ok -- no component in [1 0 0] direction */
      /* assume eigenvector has a component in [0 1 0] direction */
      /* Set up matrix corresponding to v2=1 */
      cas = 2;
      x11 = a*a + b*b + c*c;
      x12 = a*c + b*e + c*f;
      x22 = c*c + e*e + f*f;
      b1 = -(a*b + b*d + c*e);
      b2 = -(c*b + e*d + e*f);
      det = x11*x22 - x12*x12;
      if(fabs(det) < 1e-10) {  /* must have [0 0 1] as eigenvector */
	 v[0][j] = 0;
	 v[1][j] = 0;
	 v[2][j] = 1;
	 cas = 3;
      }
   }
   /* Solve Xv = b */
   if(cas != 3) {
      v1 = (x22*b1 - x12*b2)/det;
      v2 = (-x12*b1 + x11*b2)/det;
      if(cas == 1) {
	 vt[0] = 1;
	 vt[1] = v1;
	 vt[2] = v2;
      } else {
	 vt[0] = v1;
	 vt[1] = 1;
	 vt[2] = v2;
      }
      n = pvecnorm2lf(vt,3);
      v[0][j] = vt[0]/n;
      v[1][j] = vt[1]/n;
      v[2][j] = vt[2]/n;
   }
}

static void findns3f(int j, float lambda, float a, float b, float c, 
		    float d, float e, float f,
		       float **v)
{
   float x11, x12, x22, det, b1, b2,n, v1, v2;
   int cas;
   float vt[3];

   /* Now find the eigenvector */
   a = a-lambda;  d = d-lambda; f = f-lambda;
   /* Assume first that eigenvector has a component in [1 0 0] direction */
   /* Set up matrix corresponding to setting v1=1: */
   x11 = (b*b + d*d + e*e);
   x12 = (b*c + d*e + e*f);
   x22 = (c*c + e*e + f*f);
   /* check that this works */
   det = x11*x22 - x12*x12;
   b1 = -(a*b + b*d + c*e);
   b2 = -(a*c + b*e + c*f);
   cas = 1;
   if(fabs(det) < 1e-10) { /* ok -- no component in [1 0 0] direction */
      /* assume eigenvector has a component in [0 1 0] direction */
      /* Set up matrix corresponding to v2=1 */
      cas = 2;
      x11 = a*a + b*b + c*c;
      x12 = a*c + b*e + c*f;
      x22 = c*c + e*e + f*f;
      b1 = -(a*b + b*d + c*e);
      b2 = -(c*b + e*d + e*f);
      det = x11*x22 - x12*x12;
      if(fabs(det) < 1e-10) {  /* must have [0 0 1] as eigenvector */
	 v[0][j] = 0;
	 v[1][j] = 0;
	 v[2][j] = 1;
	 cas = 3;
      }
   }
   /* Solve Xv = b */
   if(cas != 3) {
      v1 = (x22*b1 - x12*b2)/det;
      v2 = (-x12*b1 + x11*b2)/det;
      if(cas == 1) {
	 vt[0] = 1;
	 vt[1] = v1;
	 vt[2] = v2;
      } else {
	 vt[0] = v1;
	 vt[1] = 1;
	 vt[2] = v2;
      }
      n = pvecnorm2f(vt,3);
      v[0][j] = vt[0]/n;
      v[1][j] = vt[1]/n;
      v[2][j] = vt[2]/n;
   }
}


/****************************************************
 *
 *
 * compute the inverse of a symmetric Toeplitz matrix
 *
 *  This is based on Golub and Van Loan, 2nd ed. section 4.7
 *
 *  Todd K. Moon
 *
 *****************************************************/

void ptoepinvlf(double **H,		/* matrix to invert */
		double **B,		/* the inverse */
		int n)			/* matrix size */
{
	
   double alpha,beta,sum,gamma_r0,norm;
   int i,j,k;

   TEMPVEC(tvec1lf,n,mv1lf,double);
   TEMPVEC(tvec2lf,n,mv2lf,double);
   /* work1 -> tvec1 */
   /*
     First use the Durbin algorithm (Alg 4.7.1) to solve
     T(n-1)y = -(r(1), ..., r(n-1))'
     
     Note: H[0][0]  = r[0],  H[0][1] = r[1], H[0][k] = r[k]
     use tvec1lf for y
   */
   
   norm = 1./H[0][0];			/* 1/r(0) */
   tvec1lf[0] = -H[0][1]*norm;
   beta = 1.0;
   alpha = -H[0][1]*norm;
   for(k=0; k < n-2; k++) {/* n-2 because solving T_{n-1} equation */
      beta = (1.0-alpha*alpha)*beta;
      sum = H[0][k+2]*norm;
      for(i = 0; i <= k; i++) {
	 sum += H[0][k-i+1]*norm*tvec1lf[i];
      }
      alpha = -sum/beta;
      /* Here use tvec2lf for z */
      for(i = 0; i <= k; i++) {
	 tvec2lf[i] = tvec1lf[i] + alpha*tvec1lf[k-i];
      }
      for(i = 0; i <= k; i++) {
	 tvec1lf[i] = tvec2lf[i];
      }
      tvec1lf[k+1] = alpha;
   }
   /* tvec1lf[0] = y[0], tvec1lf[1] = y[1], ..., tvec1lf[n-2] = y[n-2]; */
   /*
     Now do rest of calculations
   */
   sum = 1;
   for(i = 0; i < n-1; i++) {
      sum += H[0][i+1]*norm*tvec1lf[i];
   }
   /*   gamma = 1.0/sum; */
   gamma_r0 = norm/sum;
   
   /* Here use tvec2lf for nu */
   for(i = 0; i < n-1; i++) {
      tvec2lf[i] = gamma_r0*tvec1lf[n-i-2];
   }
   B[0][0] = gamma_r0;
   for(j = 1; j < n; j++) {
      B[0][j] = tvec2lf[n-j-1];
   }
   for(i = 1; i < (n-1)/2 +1; i++) {
      for(j = i; j <= n-i+1; j++) {
	 B[i][j] = B[i-1][j-1] + (tvec2lf[n-j-1]*tvec2lf[n-i-1] - 
				  tvec2lf[i-1]*tvec2lf[j-1])/gamma_r0;
      }
   }
   /*
     The algorithm computes the upper wedge.  Fill in the rest
   */
   for(j = 0; j < n - (n+1)/2; j++) {
      for(i = j+1; i < n; i++) {
	 B[i][j] = B[j][i];
      }
   }
   for(i = 1; i < n; i++) {
      for(j = n-i; j < n; j++) {
	 B[i][j] = B[n-j-1][n-i-1];
      }
   }
}

void ptoepinvf(float **H,		/* matrix to invert */
		float **B,		/* the inverse */
		int n)			/* matrix size */
{
	
   float alpha,beta,sum,gamma_r0,norm;
   int i,j,k;

   TEMPVEC(tvec1f,n,mv1f,float);
   TEMPVEC(tvec2f,n,mv2f,float);
   /* work1 -> tvec1 */
   /*
     First use the Durbin algorithm (Alg 4.7.1) to solve
     T(n-1)y = -(r(1), ..., r(n-1))'
     
     Note: H[0][0]  = r[0],  H[0][1] = r[1], H[0][k] = r[k]
     use tvec1f for y
   */
   
   norm = 1./H[0][0];			/* 1/r(0) */
   tvec1f[0] = -H[0][1]*norm;
   beta = 1.0;
   alpha = -H[0][1]*norm;
   for(k=0; k < n-2; k++) {/* n-2 because solving T_{n-1} equation */
      beta = (1.0-alpha*alpha)*beta;
      sum = H[0][k+2]*norm;
      for(i = 0; i <= k; i++) {
	 sum += H[0][k-i+1]*norm*tvec1f[i];
      }
      alpha = -sum/beta;
      /* Here use tvec2f for z */
      for(i = 0; i <= k; i++) {
	 tvec2f[i] = tvec1f[i] + alpha*tvec1f[k-i];
      }
      for(i = 0; i <= k; i++) {
	 tvec1f[i] = tvec2f[i];
      }
      tvec1f[k+1] = alpha;
   }
   /* tvec1f[0] = y[0], tvec1f[1] = y[1], ..., tvec1f[n-2] = y[n-2]; */
   /*
     Now do rest of calculations
   */
   sum = 1;
   for(i = 0; i < n-1; i++) {
      sum += H[0][i+1]*norm*tvec1f[i];
   }
   /*   gamma = 1.0/sum; */
   gamma_r0 = norm/sum;
   
   /* Here use tvec2f for nu */
   for(i = 0; i < n-1; i++) {
      tvec2f[i] = gamma_r0*tvec1f[n-i-2];
   }
   B[0][0] = gamma_r0;
   for(j = 1; j < n; j++) {
      B[0][j] = tvec2f[n-j-1];
   }
   for(i = 1; i < (n-1)/2 +1; i++) {
      for(j = i; j <= n-i+1; j++) {
	 B[i][j] = B[i-1][j-1] + (tvec2f[n-j-1]*tvec2f[n-i-1] - 
				  tvec2f[i-1]*tvec2f[j-1])/gamma_r0;
      }
   }
   /*
     The algorithm computes the upper wedge.  Fill in the rest
   */
   for(j = 0; j < n - (n+1)/2; j++) {
      for(i = j+1; i < n; i++) {
	 B[i][j] = B[j][i];
      }
   }
   for(i = 1; i < n; i++) {
      for(j = n-i; j < n; j++) {
	 B[i][j] = B[n-j-1][n-i-1];
      }
   }
}



/***************************************************************************/
/* solve the yule-walker equations  T_n y = -(r1,...,rn)' */
   /*
     Use the Durbin algorithm (Alg 4.7.1) to solve
     T(n)y = -(r(1), ..., r(n))'
     
     use tvec1f for z
   */
void pdurbinlf(double *r, double *y, int n)
{
   double beta,alpha,sum;
   int k,i;
   
   TEMPVEC(tvec1lf,n,mv1lf,double);
   y[0] = -r[1]/r[0];
   beta = r[0];
   alpha = y[0];
   for(k=0; k < n-1; k++) {	/* n-2 because solving T_{n-1} equation */
      beta = (1.0-alpha*alpha)*beta;
      sum = r[k+2];
      for(i = 0; i <= k; i++) {
	 sum += r[k-i+1]*y[i];
      }
      alpha = -sum/beta;
      /* Here use tvec1lf for z */
      for(i = 0; i <= k; i++) {
	 tvec1lf[i] = y[i] + alpha*y[k-i];
      }
      for(i = 0; i <= k; i++) {
	 y[i] = tvec1lf[i];
      }
      y[k+1] = alpha;
   }
}


void pdurbinf(float *r, float *y, int n)
{
   float beta,alpha,sum;
   int k,i;
   
   TEMPVEC(tvec1f,n,mv1f,float);
   
   y[0] = -r[1]/r[0];
   beta = r[0];
   alpha = y[0];
   for(k=0; k < n-1; k++) {	/* n-2 because solving T_{n-1} equation */
      beta = (1.0-alpha*alpha)*beta;
      sum = r[k+2];
      for(i = 0; i <= k; i++) {
	 sum += r[k-i+1]*y[i];
      }
      alpha = -sum/beta;
      /* Here use tvec1f for z */
      for(i = 0; i <= k; i++) {
	 tvec1f[i] = y[i] + alpha*y[k-i];
      }
      for(i = 0; i <= k; i++) {
	 y[i] = tvec1f[i];
      }
      y[k+1] = alpha;
   }
}


/*********************************************************
  Toeplitz solution:  Solve Tx = b, where T is a symmetric
  toeplitz matrix, determined by the vector r
**********************************************************/
void plevinsonlf(double *r, double *x, double *b, int n)
{
   int i,k;
   double beta, alpha, sum, mu;
   TEMPVEC(tvec1lf,n,mv1lf,double);  /* tvec1lf is y */
   TEMPVEC(tvec2lf,n,mv2lf,double);  /* tvec12f is z */

   alpha = tvec1lf[0] = -r[1]/r[0];
   x[0] = b[0]/r[0];
   beta = r[0];
   for(k = 0; k < n-1; k++) {
      beta = (1-alpha*alpha)*beta;
      sum = b[k+1];
      for(i = 0; i <= k; i++)
	 sum -= r[k-i+1]*x[i];
      mu = sum/beta;
      for(i = 0; i <= k; i++) {
	 x[i] += mu * tvec1lf[k-i];
      }
      x[k+1] = mu;
      if(k < n-2) {
	 sum = r[k+2];
	 for(i = 0; i <= k; i++)
	    sum += r[k-i+1]*tvec1lf[i];
	 alpha = -sum/beta;
	 for(i = 0; i <= k; i++)
	    tvec2lf[i] = tvec1lf[i] + alpha*tvec1lf[k-i];
	 for(i = 0; i <= k; i++) 
	    tvec1lf[i] = tvec2lf[i];
	 tvec1lf[k+1] = alpha;
      }
   }
}


void plevinsonf(float *r, float *x, float *b, int n)
{
   int i,k;
   float beta, alpha, sum, mu;
   TEMPVEC(tvec1f,n,mv1f,float);  /* tvec1lf is y */
   TEMPVEC(tvec2f,n,mv2f,float);  /* tvec12f is z */

   alpha = tvec1f[0] = -r[1]/r[0];
   x[0] = b[0]/r[0];
   beta = r[0];
   for(k = 0; k < n-1; k++) {
      beta = (1-alpha*alpha)*beta;
      sum = b[k+1];
      for(i = 0; i <= k; i++)
	 sum -= r[k-i+1]*x[i];
      mu = sum/beta;
      for(i = 0; i <= k; i++) {
	 x[i] += mu * tvec1f[k-i];
      }
      x[k+1] = mu;
      if(k < n-2) {
	 sum = r[k+2];
	 for(i = 0; i <= k; i++)
	    sum += r[k-i+1]*tvec1f[i];
	 alpha = -sum/beta;
	 for(i = 0; i <= k; i++)
	    tvec2f[i] = tvec1f[i] + alpha*tvec1f[k-i];
	 for(i = 0; i <= k; i++) 
	    tvec1f[i] = tvec2f[i];
	 tvec1f[k+1] = alpha;
      }
   }
}
   


/*********************************************************
  Householder triangularization.  Very similar to the QR
  decomposition, but does not save Q.  Can be performed "in place."
  Computes AT = L
  Where:
  A is rows x cols
  L is rows x cols lower triangular.
  Return value:  L is modified as a return value.
  *********************************************************/

/* void ptrilf(double **A, double **L, int rows, int cols)
{
   int i, j, k;
   double snorm;
   
   if(cols < rows) {
      pmaterr("Error in ptrilf: cols < rows");
   }
   
   TEMPVEC(tvec1lf,rows,mv1lf,double);
   TEMPVEC(tvec2lf,rows,mv2lf,double);
   
   snorm = 0.0;
   for(i=0,pd=A[0] ; i<cols; i++) 
      snorm += pd[i] * pd[i];
   snorm = sqrt(snorm);
   snorm = (A[0][0] > 0.0) ? snorm : -snorm;
   if(snorm == 0.0){
      pmaterr("Error in ptrilf: A not full rank");
   }
   
   tvec1lf[0] = snorm + A[0][0];
   snorm = 1.0 / (snorm * tvec1lf[0]);
   for(i=1,pd=A[0] ; i<cols ; i++)
      tvec1lf[i] = pd[i];
   
   tvec2lf[0] = 1.0;
   for(i=1 ; i<rows ; i++){
      tvec2lf[i] = 0.0;
      for(j=0,pd=A[i] ; j<cols ; j++)
	 tvec2lf[i] += tvec1lf[j] * pd[j];
      tvec2lf[i] *= snorm;
   }
   
   for(i=0 ; i<rows ; i++)
      for(j=0,pd1=L[i],pd=A[i] ; j<cols ; j++){
	 pd1[j] = pd[j] - tvec1lf[j] * tvec2lf[i];
      }
   
   for ( k=1; k<rows; k++ ) {
      snorm = 0.0;
      for(i=k,pd=L[k] ; i<cols; i++) 
	 snorm += pd[i] * pd[i];
      snorm = sqrt(snorm);
      snorm = (L[k][k] > 0.0) ? snorm : -snorm;
      if(snorm == 0.0){
	 pmaterr("Error in ptrilf: A not full rank");
      }
      
      tvec1lf[k-1] = 0.0;
      tvec1lf[k] = snorm + L[k][k];
      snorm = 1.0 / (snorm * tvec1lf[k]);
      for(i=k+1,pd=L[k] ; i<cols ; i++)
	 tvec1lf[i] = pd[i];
      
      tvec2lf[k-1] = 0.0;
      tvec2lf[k] = 1.0;
      for(i=k+1 ; i<rows ; i++){
	 tvec2lf[i] = 0.0;
	 for(j=0,pd=L[i] ; j<cols ; j++)
	    tvec2lf[i] += tvec1lf[j] * pd[j];
	 tvec2lf[i] *= snorm;
      }
      
      for(i=0 ; i<rows ; i++)
	 for(j=0,pd=L[i] ; j<cols ; j++){
	    pd[j] -= tvec1lf[j] * tvec2lf[i];
	 }
   }
}


void ptrif(float **A, float **L, int rows, int cols)
{
   int i, j, k;
   double snorm;
   
   if(cols < rows) {
      pmaterr("Error in ptrilf: cols < rows");
   }
   
   TEMPVEC(tvec1f,rows,mv1f,float);
   TEMPVEC(tvec2f,rows,mv2f,float);
   
   snorm = 0.0;
   for(i=0,pf=A[0] ; i<cols; i++) 
      snorm += pf[i] * pf[i];
   snorm = sqrt(snorm);
   snorm = (A[0][0] > 0.0) ? snorm : -snorm;
   if(snorm == 0.0){
      pmaterr("Error in ptrilf: A not full rank");
   }
   
   tvec1f[0] = snorm + A[0][0];
   snorm = 1.0 / (snorm * tvec1f[0]);
   for(i=1,pf=A[0] ; i<cols ; i++)
      tvec1f[i] = pf[i];
   
   tvec2f[0] = 1.0;
   for(i=1 ; i<rows ; i++){
      tvec2f[i] = 0.0;
      for(j=0,pf=A[i] ; j<cols ; j++)
	 tvec2f[i] += tvec1f[j] * pf[j];
      tvec2f[i] *= snorm;
   }
   
   for(i=0 ; i<rows ; i++)
      for(j=0,pf1=L[i],pf=A[i] ; j<cols ; j++){
	 pf1[j] = pf[j] - tvec1f[j] * tvec2f[i];
      }
   
   for ( k=1; k<rows; k++ ) {
      snorm = 0.0;
      for(i=k,pf=L[k] ; i<cols; i++) 
	 snorm += pf[i] * pf[i];
      snorm = sqrt(snorm);
      snorm = (L[k][k] > 0.0) ? snorm : -snorm;
      if(snorm == 0.0){
	 pmaterr("Error in ptrilf: A not full rank");
      }
      
      tvec1f[k-1] = 0.0;
      tvec1f[k] = snorm + L[k][k];
      snorm = 1.0 / (snorm * tvec1f[k]);
      for(i=k+1,pf=L[k] ; i<cols ; i++)
	 tvec1f[i] = pf[i];
      
      tvec2f[k-1] = 0.0;
      tvec2f[k] = 1.0;
      for(i=k+1 ; i<rows ; i++){
	 tvec2f[i] = 0.0;
	 for(j=0,pf=L[i] ; j<cols ; j++)
	    tvec2f[i] += tvec1f[j] * pf[j];
	 tvec2f[i] *= snorm;
      }
      
      for(i=0 ; i<rows ; i++)
	 for(j=0,pf=L[i] ; j<cols ; j++){
	    pf[j] -= tvec1f[j] * tvec2f[i];
	 }
   }
}
*/

/*********************************************************
  Inverse of a triangular matrix.  Can be performed "in place."
  Computes  L^-1
  Where:
  L is rows x rows lower triangular.
  Linv is rows x rows lower triangular.
  Return value:  Linv is modified as a return value.
  *********************************************************/
#define TRIABS(a) ((a>0.0) ? (a) : -(a))

void ptriinvlf(double **L, double **Linv, int rows)
{
   int i, j, k;
   
   for(i=0 ; i<rows ; i++){
      pd1=L[i];
      if(TRIABS(pd1[i]) < lueps){
	 pmaterr("Error in ptriinvlf: Singular matrix");
      }
      for(j=0,pd=Linv[i] ; j<i ; j++){
	 pd[j] = pd1[j] * -Linv[j][j];
	 for(k=j+1 ; k<i ; k++)
	    pd[j] += pd1[k] * -Linv[k][j];
	 pd[j] /= pd1[i];
      }
      pd[j] = 1.0 / pd1[i];
      for(j++ ; j<rows ; j++)
	 pd[j] = 0.0;
   }
}

void ptriinvf(float **L, float **Linv, int rows)
{
   register int i, j, k;
   
   for(i=0 ; i<rows ; i++){
      pf1=L[i];
      if(TRIABS(pf1[i]) < lueps){
	 pmaterr("Error in ptriinvlf: Singular matrix");
      }
      for(j=0,pf=Linv[i] ; j<i ; j++){
	 pf[j] = pf1[j] * -Linv[j][j];
	 for(k=j+1 ; k<i ; k++)
	    pf[j] += pf1[k] * -Linv[k][j];
	 pf[j] /= pf1[i];
      }
      pf[j] = 1.0 / pf1[i];
      for(j++ ; j<rows ; j++)
	 pf[j] = 0.0;
   }
}

#undef TRIABS

/*********************************************************/
float pmaxeigenf(float **A, float *v, float *lamda, int M, int
             Max_iterations, float Error_threshold)
/* This subroutine will attempt to determine the most significant
   eigenvalue and eigenvector for the square symmetric input matrix A.
   The rate of convergence depends on the ratio of the most significant 
   eigenvalue
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
{
   int i;
   /* Power method interation variables */
   float last_lamda; 
   float error;      /* change in lamda between iterations */

   TEMPVEC(tvec1f,M,mv1f,float);
   /* randomly select a starting point for the eigenvector iteration */
   for(i = 0; i < M; i++) v[i] = drand48();
   
   i = 0;
   last_lamda = 0.0;
   do {
      /* update the estimation of the eigenvector */
      pmatvecmultf(A,v,tvec1f,M,M);
      /* normalize the eigenvector */
      error = 1/pvecnorm2f(tvec1f,M);
      pvecscalef(tvec1f,error,v,M);
      /* Estimate the eigenvalue */
      *lamda= pvecinnerMf(v,A,v,M);
      /* have we converged? */
      error = fabs(last_lamda - *lamda);
      last_lamda = *lamda;
      i++;
   } while( (i < Max_iterations) && ( error > Error_threshold) );
   return(error);
}   


double pmaxeigenlf(double **A, double *v, double *lamda, int M, int
             Max_iterations, double Error_threshold)
{
   int i;
   /* Power method interation variables */
   double last_lamda; 
   double error;      /* change in lamda between iterations */

   TEMPVEC(tvec1lf,M,mv1lf,double);

   /* randomly select a starting point for the eigenvector iteration */
   for(i = 0; i < M; i++) v[i] = drand48();
   
   i = 0;
   last_lamda = 0.0;
   do {
      /* update the estimation of the eigenvector */
      pmatvecmultlf(A,v,tvec1lf,M,M);
      
      /* normalize the eigenvector */
      error = 1/pvecnorm2lf(tvec1lf,M);
      pvecscalelf(tvec1lf,error,v,M);
      
      /* Estimate the eigenvalue */
      *lamda = pvecinnerMlf(v,A,v,M);
      
      /* have we converged? */
      error = fabs(last_lamda - *lamda);
      last_lamda = *lamda;
      i++;
   } while( (i < Max_iterations) && ( error > Error_threshold) );
   return(error);
}   


/*********************************************************/
float pmineigenf(float **A, float *v, float *lamda, int M, int
             Max_iterations, float Error_threshold)
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
{
   int i;
   float last_lamda; 
   float error;      /* change in lamda between iterations */

   /* q->tvec1   z->tvec2, p->tvec3   pi->indx   a->tmat1 */
   TEMPVEC(indx,M,mi,int);
   TEMPMAT(tmat1f,M,M,mm1f,nm1f,float);

   /* decompose A so that the problem of Ax(i+1)=x(i) can be solved */
   /* quickly */
   pludcmpf(A,tmat1f,M,indx,v,lueps);  /* uses tvec1 */
   
   /* randomly select a starting point for the eigenvector iteration */
   for(i = 0; i < M; i++) v[i] = drand48();
   /* The Power Method iteration */
   i = 0;
   last_lamda = 0.0;
   do {
      /* update the estimation of the eigenvector */
      plubksubf(tmat1f,M,indx,v);
      /* normalize the eigenvector */
      error = 1/pvecnorm2f(v,M);
      pvecscalef(v,error,v,M);
      /* Estimate the eigenvalue */
      *lamda = pvecinnerMf(v,A,v,M);
      /* have we converged? */
      error = fabs(last_lamda - *lamda);
      last_lamda = *lamda;
      i++;
   } while( (i < Max_iterations) && ( error > Error_threshold) );
   return(error);
}   


double pmineigenlf(double **A, double *v, double *lamda, int M, int
             Max_iterations, double Error_threshold)
{
   int i;
   double last_lamda; 
   double error;      /* change in lamda between iterations */

   /* q->tvec1   z->tvec2, p->tvec3   pi->indx   a->tmat1 */
   TEMPVEC(indx,M,mi,int);
   TEMPMAT(tmat1lf,M,M,mm1lf,nm1lf,double);
   /* decompose A so that the problem of Ax(i+1)=x(i) can be solved */
   /* quickly */
   pludcmplf(A,tmat1lf,M,indx,v,lueps);  /* uses tvec1 */
   /* randomly select a starting point for the eigenvector iteration */
   for(i = 0; i < M; i++) v[i] = drand48();
   /* The Power Method iteration */
   i = 0;
   last_lamda = 0.0;
   do {
      /* update the estimation of the eigenvector */
      plubksublf(tmat1lf,M,indx,v);
      /* normalize the eigenvector */
      error = 1/pvecnorm2lf(v,M);
      pvecscalelf(v,error,v,M);
      /* Estimate the eigenvalue */
      *lamda = pvecinnerMlf(v,A,v,M);
      /* have we converged? */
      error = fabs(last_lamda - *lamda);
      last_lamda = *lamda;
      i++;
   } while( (i < Max_iterations) && ( error > Error_threshold) );
   /* make copy of the estimated eigenvector to return to the caller */
   return(error);   

}   


/*******************************************************************
Utility routines:
*/

static char *pdtn(int dt)
{
   char	*cp;
   
   switch(dt) {
   case TKUCHAR: cp = "unsigned char";	break;
   case TKSHINT: cp = "short int";	break;
   case TKINT: cp = "int";		break;
   case TKFLOAT: cp = "float";		break;
   case TKDOUBLE:cp = "double";		break;
   case TKLINT:	cp = "long int";	break;
   case TKSCHAR:cp = "signed char";	break;
      /* case 12:cp = "Complx (8)";	break; */
      /* case 13:cp = "Dcomplx (16)";	break; */
   case TKCHAR:cp = "char (ascii)";	break;
   case TKHEX:cp = "hex (int)";	        break;
   default:cp = "Unknown";				break;
   }
   return(cp);
}

static int pbytes(int type)
{
   int size;
   
   switch(type){
   case TKUCHAR: size = sizeof(char); break; 	/* UNSIGNED CHAR */
   case TKSHINT: size = sizeof(short int); break;	/* SHORT INT */
   case TKINT: size = sizeof(int); break;	/* INT */
   case TKHEX: size = sizeof(int); break;       /* HEX INT */
   case TKFLOAT: size = sizeof(float); break;	/* FLOAT */
   case TKDOUBLE: size = sizeof(double); break;	/* DOUBLE */
   case TKLINT: size = sizeof(long int); break;	/* LONG INT */
   case TKSCHAR: size = sizeof(char); break;	/* SIGNED BYTE */
   case TKCHAR:	size = sizeof(char); break;	/* ascii char */
      /* case 12:size = 2*sizeof(float); break;	*/ /* COMPLX */
      /* case 13:size = 2*sizeof(double); break;	*/ /* DCOMPLX */
   default:size = 0; break;
   }
   return(size);
}

/*****************************************************************/
static void pmaterr(char *error_text)
{
   
   
   fprintf(stderr,"Error in pmatlib: ");
   fprintf(stderr,"%s\n",error_text);
   exit(0);
}



#undef TEMPVEC
#undef TEMPMAT


void pfreealltemps(void)
{
   pfree_vector((void **)&indx);
   pfree_vector((void **)&tvec1f);
   pfree_vector((void **)&tvec2f);
   pfree_vector((void **)&tvec3f);
   pfree_vector((void **)&tvec1lf);
   pfree_vector((void **)&tvec2lf);
   pfree_vector((void **)&tvec3lf);
   mi=0;                /* length of indx */
   mv1f=mv2f=mv3f=0;  /* size of tvec1f, tvec2f, tvec3f */
   mv1lf=mv2lf=mv3lf=0;  /* size of tvec1lf, tvec2lf,tvev3lf */

   pfree_matrix((void ***) &tmat1f);
   pfree_matrix((void ***) &tmat2f);
   pfree_matrix((void ***) &tmat3f);
   pfree_matrix((void ***) &tmat4f);
   pfree_matrix((void ***) &tmat5f);
   pfree_matrix((void ***) &tmat1lf);
   pfree_matrix((void ***) &tmat2lf);
   pfree_matrix((void ***) &tmat3lf);
   pfree_matrix((void ***) &tmat4lf);
   pfree_matrix((void ***) &tmat5lf);

   mm1f=nm1f=mm2f=nm2f=mm3f=nm3f=mm4f=nm4f=mm5f=nm5f=0;
   mm1lf=nm1lf=mm2lf=nm2lf=mm3lf=nm3lf=mm4lf=nm4lf=mm5lf=nm5lf=0;
}



/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/* Stuff for ascii reading routines */


static int countmargs(char *s,int type,int *m,int *n,FILE *fp);
static int countargs(char *s,int type);
static int parsestr(char *str,char *name, int *startvar);
void plclearcomments(char *str);
static int plgetname(char **p,char *name);
static void plskwh(char **p);
static int plskwhcolcom(char **p);
static int plskiptoint(char **p);
static int plskiptohex(char **p);
static int plskiptofloat(char **p);
static int plskipwhite(char **p);
static int plskiptocomp(char **p);
static int plchecknumstr(char *str,int type);
static int readnumsc(char *str,char *a,int n);        /* characters (as numbers) */
static int readnumshd(char *str,short int *a,int n);  /* short ints */
static int readnumsd(char *str,int *a,int n);         /* ints */
static int readnumsx(char *str,int *a,int n);         /* hex to ints */
static int readnumsld(char *str,long int *a,int n);   /* long ints */
static int readnumsf(char *str,float *a,int n);       /* floats */
static int readnumslf(char *str,double *a,int n);     /* doubles */
static int readnumss(char *str,char *a,int n);  /* reads character strings */
static int readnums(int type,char *str,char *a,int n); /* read a number of specified type */
static int parsereadmatc(char ***mat,int n,int m,FILE *fp,char *str,
			 int startvar);
static int parsereadmathd(short int ***mat,int n,int m,FILE *fp,char *str,
			  int startvar);
static int parsereadmatd(int ***mat,int n,int m,FILE *fp,char *str,
			 int startvar);
static int parsereadmatld(long int ***mat,int n,int m,FILE *fp,char *str,
			  int startvar);
static int parsereadmatf(float ***mat,int n,int m,FILE *fp,char *str,
			 int startvar);
static int parsereadmatlf(double ***mat,int n,int m,FILE *fp,char *str,
			  int startvar);
static int parsereadtensor(void ****ten,int m,int n,int q,FILE *fp,char *str,
			   int type);
/* static int readnumscompf(char *str,complx *a,int n); */
/* static int readnumscomplf(char *str,dcomplx *a,int n); */
/*static int parsereadmatClf(dcomplx ***mat,int n,int m,FILE *fp,char *str,int
  startvar); */
/* static int parsereadmatCf(complx ***mat,int n,int m,FILE *fp,char *str,int
   startvar); */
static int parsereadmats(char ***mat,int n,int m,FILE *fp,char *str);
/* static int readacomplex(char **p, double *rp, double *ip); */
static int parsereadmatx(int ***mat,int n,int m,FILE *fp,char *str,int
			 startvar);

/**************************************************************************/

int asread(char *fname, char *varname,void *p, int type) 
{
   char str[MAXSTR];
   int retval;
   char name[MAXSTR];
   int startvar;
   FILE *fp;
   int tint;
   char *charptr;
   
   fp = fopen(fname,"r");
   if(fp == NULL) {
      fprintf(stderr,"Error: input file not found\n");
      return(0);
   }
   while(fgets(str,MAXSTR,fp)) {
      retval = parsestr(str,name,&startvar);
      if(retval == 0) {
	 continue;
      }
      if(!strcmp(name,varname)) { /* the variable has been found */
	 switch(type) {
	 case TKUCHAR: case TKSCHAR:
            retval = sscanf(&str[startvar]," %d ",&tint);
	    charptr = (char *)p;
            *charptr = tint;
            break;
         case TKSHINT:
            retval = sscanf(&str[startvar]," %hd ",p);
            break;
	 case TKINT:
	    retval = sscanf(&str[startvar]," %d ",p);
	    break;
	 case TKFLOAT:
	    retval = sscanf(&str[startvar]," %f",p);
	    break;
	 case TKDOUBLE:
	    retval = sscanf(&str[startvar]," %lf",p);
	    break;
	 case TKLINT:
	    retval = sscanf(&str[startvar]," %ld",p);
	    break;
	 case TKHEX:
	    retval = sscanf(&str[startvar]," %x",p);
	    break;
	 case TKCHAR:	/* character ascii */
	    if((str[startvar] == '\\' & (iscommentstart(str[startvar+1]) ||
	 isquote(str[startvar+1]) || isspace(str[startvar+1]))  ) ||
               isquote(str[startvar]))
	       startvar++;  /* skip escape or quote */
	    retval = sscanf(&str[startvar],"%c",p);
	    break;
/* #ifdef CPP
	 case 12: case 13: 
	    {
	       double rp,ip;
	       char *strp;
	       Complex *compptr;
	       
	       strp = &str[startvar];
	       retval = readacomplex(&strp,&rp,&ip);
	       compptr = (Complex *)p;
	       *compptr = Complex(rp,ip);
	    }
	    break;
#else
	 case 12:
	    {
	       double rp,ip;
	       complx *fcompptr;
	       char *strp;
	       strp = &str[startvar];
	       if(retval = readacomplex(&strp,&rp,&ip)) {
		  fcompptr = (complx *)p;
		  fcompptr->re = (float)rp;
		  fcompptr->im = (float)ip;
	       }
	    }
	    break;
	 case 13:
	    {
	       double rp,ip;
	       dcomplx *dcompptr;
	       char *strp;
	       strp = &str[startvar];
	       if(retval = readacomplex(&strp,&rp,&ip)) {
		  dcompptr = (dcomplx *)p;
		  dcompptr->re = rp;
		  dcompptr->im = ip;
	       }
	    }
	    break;
#endif
*/
	 default: 
	    fprintf(stderr,"Type %d unknown\n",type);
	    fclose(fp);
	    return(0);
	 }
	 fclose(fp);
	 if(retval == 0) {	/* no input was read */
	    fprintf(stderr,"Error in reading %s\n",varname);
	    return(0);
	 }
	 return(1);
      }
   }
   fclose(fp);
   fprintf(stderr,"Variable %s not found in data file\n",varname);
   return(0);
}

/**************************************************************************/

/* read a vector from the data file */
int avread(char *fname, char *varname,void **p, int *n, int type) 
{
   char str[MAXSTR];
   int retval;
   char name[MAXSTR];
   int startvar;
   FILE *fp;
   int nr,nr1;

   
   fp = fopen(fname,"r");
   if(fp == NULL) {
      fprintf(stderr,"Error: input file not found\n");
      return(0);
   }
   while(fgets(str,MAXSTR,fp)) {
      retval = parsestr(str,name,&startvar);
      if(retval == 0) {
	 continue;
      }
      
      if(!strcmp(name,varname)) { /* the variable has been found */
	 nr = countargs(&str[startvar],type);
	 if(nr) {		/* valid input */
	    *p = (void *)malloc(nr*pbytes(type)); /* allocate memory */
	    switch(type) {
	    case TKUCHAR:
	       nr1 = readnumsc(&str[startvar],(char *)(*p),nr); 
               /* unsigned char */
	       break;
	    case TKSHINT:
	       nr1 = readnumshd(&str[startvar],(short int *)(*p),nr); 
               /* short int */
	       break;
	    case TKINT:
	       nr1 = readnumsd(&str[startvar],(int *)(*p),nr); /*  int */
	       break;
	    case TKFLOAT:
	       nr1 = readnumsf(&str[startvar],(float *)(*p),nr); /* float */
	       break;
	    case TKDOUBLE:
	       nr1 = readnumslf(&str[startvar],(double *)(*p),nr); /* double */
	       break;
	    case TKLINT:
	       nr1 = readnumsld(&str[startvar],(long int*)(*p),nr); /* long int */
	       break;
	    case TKSCHAR:
	       nr1 = readnumsc(&str[startvar],(*p),nr); /* signed char */
	       break;
	       /* case 12:
		  nr1 = readnumscompf(&str[startvar],(*p),nr);
		  break; */
#ifdef CPP
	       /*case 13:	
		*p = new Complex[nr];
		nr1 = readnumscomp(&str[startvar],*p,nr);
		break; */
#else
	       /* case 13:
		  nr1 = readnumscomplf(&str[startvar],*p,nr);
		  break; */
#endif
	    case TKCHAR: /* make sure there is enough space for the \0 at the */
	       /* end */
	       free(*p);
	       *p = (char *)malloc(nr+1); /* allocate memory */
	       nr1 = readnumss(&str[startvar],*p,nr);
	       nr1 = nr;
	       break;
	    case TKHEX:
	       nr1 = readnumsx(&str[startvar],(int *)(*p),nr); /*  int */
	       break;
	    default: 
	       fprintf(stderr,"Type %d unknown\n",type);
	       fclose(fp);
	       return(0);
	    }
	    if(nr != nr1) {
	       fprintf(stderr,"Error in reading %s\n",varname);
	       return(0);
	    }
	    *n = nr;
	    return(1);
	 } /* end if valid set of inputs */
      } /* end if variable found */
   } /* end while loop */
   fclose(fp);
   fprintf(stderr,"Variable %s not found in data file\n",varname);
   *n = 0;
   return(0);
}

/**************************************************************************/
/* read a matrix from the data file */
int amread(char *fname, char *varname,void ***p, int *m, int *n,int type) 
{
   char str[MAXSTR];
   int retval;
   char name[MAXSTR];
   int startvar;
   FILE *fp;
   fp = fopen(fname,"r");
   if(fp == NULL) {
      fprintf(stderr,"Error: input file not found\n");
      return(0);
   }
   while(fgets(str,MAXSTR,fp)) {
      retval = parsestr(str,name,&startvar);
      if(retval == 0) {
	 continue;
      }

      if(!strcmp(name,varname)) { /* the variable has been found */
	 retval = countmargs(&str[startvar],type,m,n,fp);
	 if(retval == 0) { fclose(fp); return(0); }
	 switch(type) {
	 case TKUCHAR: case TKSCHAR:
	    retval = parsereadmatc((char ***)p,*m,*n,fp,str,startvar);
	    break;
	 case TKSHINT:
	    retval = parsereadmathd((short int ***)p,*m,*n,fp,str,startvar);
	    break;
	 case TKINT:
	    retval = parsereadmatd((int ***)p,*m,*n,fp,str,startvar);
	    break;
	 case TKFLOAT:
	    retval = parsereadmatf((float ***)p,*m,*n,fp,str,startvar);
	    break;
	 case TKDOUBLE:
	    retval = parsereadmatlf((double ***)p,*m,*n,fp,str,startvar);
	    break;
	 case TKLINT:
	    retval = parsereadmatld((long int ***)p,*m,*n,fp,str,startvar);
	    break;
/*	 case 12:
		retval = parsereadmatCf((complx ***)p,*m,*n,fp,str,startvar);
		break;
	 case 13:
		retval = parsereadmatClf((dcomplx ***)p,*m,*n,fp,str,startvar);
		break;
*/
	 case TKCHAR:   /* array of character strings */
	    retval = parsereadmats((char ***)p,*m,*n,fp,&str[startvar]);
	    break;
	 case TKHEX:   /* integer array expressed in hex */
	    retval = parsereadmatx((int ***)p,*m,*n,fp,str,startvar);
	    break;
	 default: 
	    fprintf(stderr,"Type %d not implemented\n",type);
	    retval = 0;
	 }			/* end switch */
	 fclose(fp);
	 if(retval != 0) {
	    return(1);
	 }
	 return(1);
      } /* end if valid set of inputs */
   } /* end while loop */
   fclose(fp);
   fprintf(stderr,"Variable %s not found in data file\n",varname);
   *m = *n = 0;
   return(0);
}




/**************************************************************************/
/* read a (3-dimensional) tensor from the data file */
int atread(char *fname, char *varname, void ****p, int m, int *n, int
		   *q, int type) 
{
   char str[MAXSTR];
   int retval;
   char name[MAXSTR];
   int startvar;
   FILE *fp;
   int numrows;

/* the tensor is assumed to be stored in the form
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

   
   
   fp = fopen(fname,"r");
   if(fp == NULL) {
      fprintf(stderr,"Error: input file not found\n");
      return(0);
   }
   while(fgets(str,MAXSTR,fp)) {
      retval = parsestr(str,name,&startvar);
      if(retval == 0) {
	 continue;
      }
      if(!strcmp(name,varname)) { /* the variable has been found */
	 
	 retval = countmargs(&str[startvar],type,&numrows,q,fp);
	 if(retval == 0) return(0);
	 if(numrows % m) {	       /* not an even divisor */
	    fprintf(stderr,"Warning: cannot determine tensor\n");
	    return(0);
	 }
	 *n = numrows/m;
	 retval = parsereadtensor(p,m,*n,*q,fp,&str[startvar],type);
	 fclose(fp);
	 return(1);
      }
   }
   fprintf(stderr,"Variable %s not found in data file\n",varname);
   *n = *q = 0;
   return(0);
}



/**************************************************************************/
static int countmargs(char *s,int type,int *m,int *n,FILE *fp)
{
   long fpos;
   char numstr[MAXSTR];
   int nr,nc,nc1;

   nc = countargs(s,type);	/*  get the number of columns */
   if(!nc) return(0);
   fpos = ftell(fp);		/* reposition the file here */
   nr = 1;			/* we have read the first row already */
   while(fgets(numstr,MAXSTR,fp) && (nc1=plchecknumstr(numstr,type))) {
      nr++; /* search */
      if(nc1 > nc) nc = nc1;
   }

   fseek(fp,fpos,0);		/* reset file position */
   *m = nr;
   *n = nc;
   return(1);
}


/**************************************************************************/
/* count the number of arguments on this line */
static int countargs(char *s,int type)
{

   char *p = s;
   int numn = 0;
   char c;
   
   
   switch(type) {
   case TKUCHAR: case TKSHINT: case TKINT: case TKLINT: case TKSCHAR:
      /* integer type arguments */
      while(1) {
	 while((c = *p) && !isintnum(c)) p++; /* skip white space */
	 if(!c) return(numn);
	 while((c = *p) && isintnum(c)) p++; /* move past number */
	 numn++;
	 if(!c) return(numn);
      }
   case TKFLOAT: case TKDOUBLE:		/* floating type arguments */
      while(1) {
	 while((c = *p) && !isfnum(c)) p++; /* skip white space */
	 if(!c) return(numn);
	 while((c = *p) && isfnum(c)) p++; /* move past number */
	 numn++;
	 if(!c) return(numn);
      }
      /* case 12: case 13:
      {
	 double rp,ip;
	 while(1) {
	    if(readacomplex(&p,&rp,&ip)) {
	       numn++;
	       continue;
	    }
	    else break;
	 }
	 return(numn);
	 } */
   case TKCHAR:	/* string characters */
      if(s[(int)strlen(s)-1] == '\n') s[(int)strlen(s)-1] = 0;
      numn = 0;   
      while(*p) {
	 p++;
	 numn++;
      }
      return(numn);  /* add 1 for 0 terminator */
      break;
   case TKHEX:  /* hex numbers */
      while(1) {
	 while((c = *p) && !ishexnum(c)) p++; /* skip white space */
	 if(!c) return(numn);
	 while((c = *p) && ishexnum(c)) p++; /* move past number */
	 numn++;
	 if(!c) return(numn);
      }
   default:
      fprintf(stderr,"Type %d not supported for countargs\n");
      return(0);
   }
}



/**************************************************************************/
/* take an input string, and break it into the name, the address of where proper data
   begins and ends, and check for format.  Proper formats are the following:
   comments are delimited by ; or by % or by # or by // or 
   by using tradtional C comments, provided that the C comments begin
   and end on the SAME line

name 1   ; a single scalar value
name: 1  ; a single scalar value with colon after name
name 1 2 3  ; vector
name: 1 2 3 ; vector with colon
name 1,2,3  ; commas do not matter
name [1 2 3]   ; brackets are ignored 
name: 1,2,3
name: 1 2 3
name: 1 2 3
      4 5 6       ; matrix here
strname: This is my string  ; a string variable
strname This is my string   ; another string variable
strname "The quotes are preserved"  ; only comments are ignored

Matrices should not be broken by comment lines:
 name: 1 2 3
; No place to put a comment!
       4 5 6

White space is, of course, ignored.

*/

static int parsestr(char *str,char *name, int *startvar)
{
   char *p;			/* point to current position in input string */
   
   *startvar = 0;
   *name = 0;  

   plclearcomments(str);	/* get rid of any comments */
   p = str;
   if(plgetname(&p,name) == 0) {
      return(0);
   }
   if(plskwhcolcom(&p) == 0)
      return(0);
   *startvar = p - str;
   return(1);
}




/**************************************************************************/
/* get rid of any comments */
void plclearcomments(char *str)
{
   int i,j;
   int quoted = 0;
   char quotechar = 0;
   int ccomment = 0;

   for(i = 0; i < (int)strlen(str)-1; i++) {
      if(str[i] == '/' && str[i+1] == '*') {  /* beginning of comment found */
	 for(j = i+2; j < (int)strlen(str)-1; j++) {
	    if(str[j] == '*' && str[j+1] == '/') { /* end of comment found */
	       for(j = j+2; j <= (int)strlen(str); j++,i++) {
		  str[i] = str[j];
	       }
	       ccomment = 1;
	    }
	 }
      }
   }
   if(ccomment)
      plclearcomments(str);

   for(i = 0; i < strlen(str); i++) {
      if(isquote(str[i])) {
	 if(str[i] == quotechar) {  /* closed quote */
	    quoted--;
	 }
         else if(str[i] != quotechar) {  /* quote within a quote */
	    quoted++;
	    quotechar = str[i];
	 }
      }
      if(iscommentstart(str[i]) && !quoted) {
	 if(i && str[i-1] == '\\') continue;  /* escape char */
	 str[i] = 0;
	 break;
      }
      if(i < ((int)strlen(str)-1) && str[i] == '/' && str[i+1] == '/') {
	 str[i] = 0;
	 break;
      }
   }

}


/* skip any white space, determine of the resulting string is a valid */
 /* C identifier, and return 0.  Otherwise, return 1 to indicate error */

/**************************************************************************/
static int plgetname(char **p,char *name)
{
   plskwh(p);		/* skip all white space, return with */
				/* *p on first non-white character */
   if(!isCidstart(**p)) {		/* invalid identifier */
      return(0);
   }
   /* otherwise, copy in the string as long as it is */
   while(isCid(**p)) {
      *name = **p;
      name++;
      (*p)++;
   }
   *name = 0;
   return(1);
}
    
/**************************************************************************/
static void plskwh(char **p)
{

   while(**p && isspace(**p)) (*p)++;
}

/**************************************************************************/
static int plskwhcolcom(char **p) {
   char c;
   while((c = **p) && (isspace(c) || (c==':') || (c == ',') || (c == '[') ||
        (c == ']') || (c == '='))) (*p)++;
   if(!**p) return(0);
   else return(1);

}


/**************************************************************************/
static int plskiptoint(char **p) {
   char c;
   while((c = **p) && !isintnum(c)) (*p)++;
   if(!**p) return(0);
   else return(1);

}


/**************************************************************************/
static int plskiptohex(char **p) {
   char c;
   while((c = **p) && !ishexnum(c)) (*p)++;
   if(!**p) return(0);
   else return(1);

}

/**************************************************************************/
static int plskiptofloat(char **p) {
   char c;
   while((c = **p) && !isfnum(c)) (*p)++;
   if(!**p) return(0);
   else return(1);

}



/**************************************************************************/
static int plskipwhite(char **p) {
   char c;
   while((c = **p) && isspace(c)) (*p)++;
   if(!**p) return(0);
   else return(1);

}

/**************************************************************************/
static int plskiptocomp(char **p) {
   char c;
   while((c = **p) && !isfnum(c) && (c != '(')) (*p)++;
   if(!**p) return(0);
   else return(1);

}


/**************************************************************************/
/* checkstr makes sure that the current line contains numbers and/or */
 /* comments only, not new identifiers */
static int plchecknumstr(char *str,int type)
{
   char *p;			/* point to current position in input string */
   char name[MAXSTR];
   char c;

   plclearcomments(str);	/* get rid of any comments */
   p = str;
   if(plgetname(&p,name) == 1) {/* oops! a name was found, this is not */
				/* a line of numbers */
      return(0);
   }

   if(*p == '^' && type == TKCHAR) { /* continuation for arrays of strings */
	  p++;
	  while((c = *p) && isspace(c)) p++;
	  if(c) return(countargs(p,type));

	  return(0);   /* empty string */
   }

   return(countargs(p,type));
}


/**************************************************************************/
static int readnums(int type,char *str,char *a,int n)
{

   switch(type) {
   case TKUCHAR: case TKSCHAR:
      return(readnumsc(str,(char *)a,n));
   case TKSHINT:
      return(readnumshd(str,(short *)a,n));
   case TKINT:
      return(readnumsd(str,(int *)a,n));
   case TKFLOAT:
      return(readnumsf(str,(float *)a,n));
   case TKDOUBLE:
      return(readnumslf(str,(double *)a,n));
   case TKLINT:
      return(readnumsld(str,(long *)a,n));
   case TKHEX:
      return(readnumsx(str,(int *)a,n));
   case TKCHAR:
      return(readnumss(str,(char *)a,n));
   default:
      fprintf(stderr,"Type %d not defined in readnums\n",type);
      break;
   }
}


/**************************************************************************/
/* read from the string str up to n integers into a */
static int readnumsd(char *str,int *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptoint(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isintnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = atoi(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
/* read from the string str up to n integers into a */
static int readnumsx(char *str,int *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptohex(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && ishexnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      sscanf(temp,"%x",&(a[i]));
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
/* read from the string str up to n integers into a */
static int readnumsld(char *str,long int *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptoint(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isintnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = atol(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
/* read from the string str up to n integers into a */
static int readnumss(char *str,char *a,int n)
{
   char c;
   char *p;
   int j;

   plclearcomments(str);	/* get rid of any comments */
   p = str;
   if(plskwhcolcom(&p) == 0) {	/* in skipping, reached end of line */
      return(0);
   }

   j = 0;
   while((c = *p) && (j < n)) {
      a[j++] = c;
      p++;
   }
   a[j] = 0;
   return(strlen(a));
}


/**************************************************************************/
/* read from the string str up to n integers into a */
static int readnumshd(char *str,short int *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptoint(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isintnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = (short)atoi(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
static int readnumsf(char *str,float *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptofloat(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isfnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = (float)atof(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
static int readnumslf(char *str,double *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptofloat(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isfnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = (double)atof(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}


/**************************************************************************/
/* read from the string str up to n character integers into a */
static int readnumsc(char *str,char *a,int n)
{
   int i;
   char temp[MAXSTR];
   char c;
   char *p;
   int j;
   int numn;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)	/* read n numbers */
   {  
      if(plskiptoint(&p) == 0) {	/* in skipping, reached end of line */
	 return(numn);
      }
      
      j = 0;
      while((c = *p) && isintnum(c))	/* while numeric */
      {  temp[j++] = c;
	 p++;
      }
      temp[j] = 0;
      a[i] = (char)atoi(temp);
      numn++;
      if(!c) return(numn);
   }
   return(numn);
}




/**************************************************************************/
/*static int readnumscompf(char *str,complx *a,int n)
{
   int i;
   char *p;
   int numn;
   double rp,ip;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)
   {  
	  if(readacomplex(&p,&rp,&ip)) {
		 a[i].re = (float)rp;
		 a[i].im = (float)ip;
		 numn++;
	  }
   }
   return(numn);
}
*/


/**************************************************************************/
/* static int readnumscomplf(char *str,dcomplx *a,int n)
{
   int i;
   char *p;
   int numn;
   double rp,ip;

   p = str;
   numn = 0;
   for(i = 0; i < n; i++)
   {  
	  if(readacomplex(&p,&rp,&ip)) {
		 a[i].re = rp;
		 a[i].im = ip;
		 numn++;
	  }
   }
   return(numn);
}
*/


/**************************************************************************/
static int parsereadmatc(char ***mat,int n,int m,FILE *fp,char *str,
   int startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (char **)pcalloc_matrix(TKUCHAR,n,m,"parsemreadmatc");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumsc(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}


/**************************************************************************/
static int parsereadmathd(short int ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (short int **)pcalloc_matrix(TKSHINT,n,m,"parsemreadmathd");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumshd(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}

/**************************************************************************/
static int parsereadmatd(int ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (int **)pcalloc_matrix(TKINT,n,m,"parsemreadmatd");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumsd(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}


/**************************************************************************/
static int parsereadmatx(int ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (int **)pcalloc_matrix(TKINT,n,m,"parsemreadmatx");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumsx(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}


/**************************************************************************/
static int parsereadmatld(long int ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (long int **)pcalloc_matrix(TKDOUBLE,n,m,"parsemreadmatld");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumsld(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}


/**************************************************************************/
static int parsereadmatf(float ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (float **)pcalloc_matrix(TKFLOAT,n,m,"parsemreadmatf");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumsf(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}



/**************************************************************************/
static int parsereadmatlf(double ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (double **)pcalloc_matrix(TKDOUBLE,n,m,"parsemreadmatlf");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {			/* first time through, read from input */
				/* string */
	 strcpy(numstr,&str[startvar]);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
      }
      readnumslf(numstr,(*mat)[i],m);
   } /* end for */
   return(1);
}


/**************************************************************************/
/*
static int parsereadmatClf(dcomplx ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;

   *mat = (dcomplx **)pcalloc_matrix(DCOMPLX,n,m,"parsemreadmatClf");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++) {
	 (*mat)[i][j].re = 0;
	 (*mat)[i][j].im = 0;
      }
   for(i = 0; i < n; i++) {
      if(!i) {
	 strcpy(numstr,&str[startvar]);
      }
      else {
//	 fpos = ftell(fp);
	 retval = fgets(numstr,MAXSTR,fp);
//	 if(!retval || (plchecknumstr(numstr)==0)) {
//	    fseek(fp,fpos,0);
//	    break;
//	 }
//	 else {
//	    strcpy(numstr,str);
//	 }
//      }
      readnumscomplf(numstr,(*mat)[i],m);
   }
   return(1);
}
*/


/**************************************************************************/
/*
static int parsereadmatCf(complx ***mat,int n,int m,FILE *fp,char *str,int
	      startvar)
{
   int i,j;
   char numstr[MAXSTR];
   long fpos;
   char *retval;
   
   *mat = (complx **)pcalloc_matrix(COMPLX,n,m,"parsemreadmatCf");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++) {
	 (*mat)[i][j].re = 0;
	 (*mat)[i][j].im = 0;
      }
   for(i = 0; i < n; i++) {
      if(!i) {
	 strcpy(numstr,&str[startvar]);
      }
      else {
//	 fpos = ftell(fp);
	 retval = fgets(numstr,MAXSTR,fp);
//	 if(!retval || (plchecknumstr(numstr)==0)) {
//	    fseek(fp,fpos,0);
//	    break;
//	 }
//	 else {
//	    strcpy(numstr,str);
//	 }
//      }
      readnumscompf(numstr,(*mat)[i],m);
   }
   return(1);
}
*/


/**************************************************************************/
static int parsereadmats(char ***mat,int n,int m,FILE *fp,char *str)
{
   int i,j,j1;
   char numstr[MAXSTR];
   long fpos;
   char *retval;


   *mat = (char **)pcalloc_matrix(TKCHAR,n,m+1,"parsemreadmats");
   for(i = 0; i < n; i++)
      for(j = 0; j < m; j++)
	 (*mat)[i][j] = 0;
   for(i = 0; i < n; i++) {
      if(!i) {	/* first time through, read from input  string */
	 strcpy(numstr,str);
      }
      else {			/* after that, read from file */
	 retval = fgets(numstr,MAXSTR,fp);
	 if(numstr[(int)strlen(numstr)-1] == '\n')
	    numstr[(int)strlen(numstr)-1] = 0;
	 for(j = 0; j < strlen(numstr); j++) {
	    if(isspace(numstr[j]) || numstr[j] == '^') continue;
	    break;
	 }
	 for(j1 = 0; numstr[j]; j1++,j++)
	    numstr[j1] = numstr[j];
	 numstr[j1] = 0;
      }
      strncpy((*mat)[i],numstr,m);
      (*mat)[i][m] = 0;
   } /* end for */
   return(1);
}


/**************************************************************************/
static int parsereadtensor(void ****ten,int m,int n,int q,FILE *fp,
			   char *str,int type)
{
   int i,j,k,q1,j1;
   char numstr[MAXSTR];
   long fpos;
   char *retval;
   
   q1 = q;
   if(type == TKCHAR) {  /* allow one more byte for terminating \0 */
      q1 = q+1;
   }
   *ten = (void ***)pcalloc_tensor(type,m,n,q1,"parsemreadtenc");

   for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
	 if(!i && !j) {	/* first time through, read from input string */
	    strcpy(numstr,str);
	 }
	 else {			/* after that, read from file */
	    retval = fgets(numstr,MAXSTR,fp);
	    if(type == TKCHAR) {
	       if(numstr[(int)strlen(numstr)-1] == '\n')
		  numstr[(int)strlen(numstr)-1] = 0;
	       for(k = 0; k < strlen(numstr); k++) {
		  if(isspace(numstr[k]) || numstr[k] == '^') continue;
		  break;
	       }
	       for(j1 = 0; numstr[k]; j1++,k++)
		  numstr[j1] = numstr[k];
	       numstr[j1] = 0;
	    }
	 } /* end else read from file */
	 readnums(type,numstr,(*ten)[i][j],q);
      } /* end for j */
   } /* end for i */
   return(1);
}



/**************************************************************************/
/*
static int readacomplex(char **p, double *rp, double *ip)
{

   int j;
   char c;
   char temp[200];


   *rp = 0;
   *ip = 0;
   if(plskiptocomp(p) == 0) {
	  return(0);
   }
   if(**p == LPR)	{	
	  (*p)++;		
	  if(plskiptofloat(p) == 0) { 
		 return(0);
	  }
	  j = 0;
	  while((c = **p) && isfnum(c))  { 
		 temp[j++] = c;
		 (*p)++;
	  }
	  temp[j] = 0;
	  *rp = (double)atof(temp);
	  if(plskipwhite(p) == 0) {
		 return(1);
	  }
	  if(**p == ',') {
		 (*p)++;	
		 if(plskiptofloat(p) == 0) {
			goto skiptorpr;
		 }
		 j = 0;
		 while((c = **p) && isfnum(c))
			{  temp[j++] = c;
			   (*p)++;
			}
		 temp[j] = 0;
		 *ip = (double)atof(temp);
	  }
   skiptorpr:
	  if(plskipwhite(p) == 0) {
		 return(1);
	  }
	  if(**p != RPR) {
		 fprintf(stderr,"Syntax error\n");
			return(1);
	  }
	  (*p)++;
	  return(1);
   }
   else {		
	  j = 0;	
	  while((c = **p) && isfnum(c)) {
		 temp[j++] = c;
		 (*p)++;
	  }
	  temp[j] = 0;
	  *rp = atof(temp);
	  return(1);
   }
}
*/


/**********************************************************************
 *
 *  sort2.c -- a quicksort routine
 *
 **********************************************************************/

static double *alf;
static int *a2d;
static float *af;

/**********************************************************************/

static void quicksort2lfd(int left, int right)
{
    long i, j;
    double ref;
    int ref2;
    i = left;
    j = right;
    ref = alf[i];
    ref2 = a2d[i];
    while (i < j)
    {
	while (i < j && ref -  alf[j] < 0)
	    j--;
	if (i != j) {
	   a2d[i] = a2d[j];
	    alf[i++] = alf[j];
	}
	while (i < j && ref - alf[i] > 0)
	    i++;
	if (i != j) {
	   a2d[j] = a2d[i];
	    alf[j--] = alf[i];
	}
    }
    a2d[j] = ref2;
    alf[j] = ref;
    if (left < --j)
	quicksort2lfd(left, j);
    if (++i < right)
	quicksort2lfd(i, right);
}
/**********************************************************************/

/* sort array into INCREASING ORDER, and shuffle array2 at the same time
   if array2 = 0,1,...,N and A represents the original _unsorted_ array,
   then after the sort, A[array2[0]], A[array2[1]], ...
   represents the sorted data in array.
*/
void sort2lfd(int num_elements, double *array, int *array2)
{
    if (num_elements < 2)
	return;
    alf = array;
    a2d = array2;
    quicksort2lfd(0, num_elements - 1);
}

/**********************************************************************/
static void quicksort1lf(int left, int right)
{
    long i, j;
    double ref;
    int ref2;
    i = left;
    j = right;
    ref = alf[i];
    while (i < j)
    {
	while (i < j && ref -  alf[j] < 0)
	    j--;
	if (i != j) {
	    alf[i++] = alf[j];
	}
	while (i < j && ref - alf[i] > 0)
	    i++;
	if (i != j) {
	    alf[j--] = alf[i];
	}
    }
    alf[j] = ref;
    if (left < --j)
	quicksort1lf(left, j);
    if (++i < right)
	quicksort1lf(i, right);
}
/**********************************************************************/

void sort1lf(int num_elements, double *array)
{
    if (num_elements < 2)
	return;
    alf = array;
    quicksort1lf(0, num_elements - 1);
}

/**********************************************************************/

static void quicksort2fd(int left, int right)
{
    long i, j;
    float ref;
    int ref2;
    i = left;
    j = right;
    ref = af[i];
    ref2 = a2d[i];
    while (i < j)
    {
	while (i < j && ref -  af[j] < 0)
	    j--;
	if (i != j) {
	   a2d[i] = a2d[j];
	    af[i++] = af[j];
	}
	while (i < j && ref - af[i] > 0)
	    i++;
	if (i != j) {
	   a2d[j] = a2d[i];
	    af[j--] = af[i];
	}
    }
    a2d[j] = ref2;
    af[j] = ref;
    if (left < --j)
	quicksort2fd(left, j);
    if (++i < right)
	quicksort2fd(i, right);
}
/**********************************************************************/

/* sort array into INCREASING ORDER, and shuffle array2 at the same time
   if array2 = 0,1,...,N and A represents the original _unsorted_ array,
   then after the sort, A[array2[0]], A[array2[1]], ...
   represents the sorted data in array.
*/
void sort2fd(int num_elements, float *array, int *array2)
{
    if (num_elements < 2)
	return;
    af = array;
    a2d = array2;
    quicksort2fd(0, num_elements - 1);
}

/**********************************************************************/
static void quicksort1f(int left, int right)
{
    long i, j;
    float ref;
    int ref2;
    i = left;
    j = right;
    ref = af[i];
    while (i < j)
    {
	while (i < j && ref -  af[j] < 0)
	    j--;
	if (i != j) {
	    af[i++] = af[j];
	}
	while (i < j && ref - af[i] > 0)
	    i++;
	if (i != j) {
	    af[j--] = af[i];
	}
    }
    af[j] = ref;
    if (left < --j)
	quicksort1f(left, j);
    if (++i < right)
	quicksort1f(i, right);
}
/**********************************************************************/

void sort1f(int num_elements, float *array)
{
    if (num_elements < 2)
	return;
    af = array;
    quicksort1f(0, num_elements - 1);
}

/**********************************************************************/
/* Generate a unit-variance Gaussian random number.
   This function calls rand(), so you can control the see using
   void srand(unsigned int seed);

   Following Numerical Recipes
*/
static int graniset = 0;
static double grangset;

double gran(void)
{
   double rsq, v1, v2, fac;

   if(!graniset) {
	  graniset = 1;
	  do {
		 v1 = 2*(rand()/(double)RAND_MAX) - 1;
		 v2 = 2*(rand()/(double)RAND_MAX) - 1;
		 rsq = v1*v1 + v2*v2;
	  } while(rsq > 1 || rsq == 0);
	  fac = sqrt(-2*log(rsq)/rsq);
	  grangset = v1*fac;
	  return v2*fac;
   }
   else {
	  graniset = 0;
	  return grangset;
   }
}

double uran(void)
/* Generate a uniform random number between 0 and 1.  This 
   function calls rand(), so you can control the seed with srand().
*/
{
   return rand()/(double)RAND_MAX;
}

/*
Local Variables:
compile-command: "gcc -o testpmatlib -g testpmatlib.c pmatlib.c -lm"
End:
*/
