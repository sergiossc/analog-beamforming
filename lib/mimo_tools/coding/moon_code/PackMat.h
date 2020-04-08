#ifndef PACKMAT_H
#define PACKMAT_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <stdlib.h>
#include<cmath>
#include<math.h>
#include"matalloc.h"

using namespace std;

enum PackMode {HORIZPACK=0,VERTPACK=1};
enum StackMode {Default,Exact,Top,Truncate, Bottom, Center,Left,Right};

#ifndef MIN
#define MIN(a,b)  (a < b? a: b)
#endif
#ifndef MAX
#define MAX(a,b)  (a > b? a: b)
#endif

typedef unsigned int PKSIZE;	// integer size things are packed into
//typedef unsigned char PKSIZE;	// integer size things are packed into
typedef unsigned char CTSIZE;			// size of stuff used to count weights

class PackVec;
class PackMat;

class PackMat
{
   int numcols;
   int numrows;
   PKSIZE **Matrix;
   int numpackcols;				// number of packed columns
   int numpackcols1;			// number of entire packed columns
   int numpackcols2;			// number of partial packed columns
   int numpackrows;				// number of packed rows
   int numpackrows1;			// number of entire packed rows
   int numpackrows2;			// number of partial packed columns
   PackMode packmode;			// row packed or column packed
   int fracbits;				// number of bits left over

   static int *wtct;			// array used to count bit weight
   static int PkperCt;			// number of bytes per CTSIZE
   static int shifter;			// amount to shift by to divide
   static PKSIZE masker;		// bit mask of entire blob
   static StackMode sidestackmode;	// way of stacking side-by-side
   static StackMode topstackmode;	// way of stacking side-by-side

   friend class PackVec;		// allow PackVec to access the private stuff

   friend void stackside(PackMat & C,PackMat &A, PackMat & B,StackMode ssm,
						 PackMode pkmode);
   friend void stacktop(PackMat & C,PackMat &A, PackMat & B,StackMode tsm,
						 PackMode pkmode);
   friend void setval(PKSIZE **Matrix,int row, int col, int val, 
			PackMode packmode);



public:

   static int bitsper;			// number of bits in each element
   static void initPackMat();
   PackMat(void);
   PackMat(PackMat &D);
   PackMat(int rows, int cols, PackMode rowcolpack=HORIZPACK);
   PackMat(char *fname,PackMode rowcolpack = HORIZPACK);
   PackMat(string desc,PackMode rowcolpack = HORIZPACK);
   PackMat(istream &Densfile,PackMode rowcolpack = HORIZPACK);
   void PackMat_is(istream &Densfile,PackMode rowcolpack = HORIZPACK);
   int ReadFile(char *fname,PackMode rowcolpack = HORIZPACK);
   void fill(string fillstring);
   void fill(istream &Densfile);

   PackMat(PackVec &,PackMode rowcolpack = HORIZPACK);

   static void setSideStackMode(StackMode ssm) { sidestackmode = ssm;}
   static void setTopStackMode(StackMode tsm) { topstackmode = tsm;}

   void Size(int rows, int cols, PackMode pkmode=HORIZPACK);
   void Size(PackMat &A);
   void  resize(PackMat &A,int rows, int cols, PackMode pkmode=HORIZPACK);
   // resize, but keep the old stuff

   void PackMat::MakeRandom(double p=0);
   void PackMat::SetToZero();
   void PackMat::SetToOne();
   void PackMat::Identity(int n, PackMode pkmode=HORIZPACK);
   int PackMat::IsIdentity(void);
   inline void checkandset(PackMat &B,PackMode bpackmode,int brows, int bcols);
   // make sure there is exactly brows and bcols space
   void checkandset2(PackMat &B,PackMode bpackmode,int brows, int bcols);
   // make sure there is at least enough space
   void checkandresize2(PackMat &B,PackMode bpackmode,int brows, int bcols);
   // make sure there is enough space, and keep the old stuff

   // make sure there is at least brows and bcols space
   inline void setsize(PackMat &A, int rows, int cols, PackMode pkmode);
   inline void PackMat::callocmat(PackMat &A);
   void PackMat::clearfringe();
   void Submatrix(PackMat &submat,int row1, int row2, int col1, int col2,
						int subrow=0, int subcol=0,PackMode pkmode=HORIZPACK);
   void Subvector(PackVec &subvec,int row1, int row2,int col1, int col2,
				  int subrow=0);



   ~PackMat(void) ;

   PackMat &operator=(PackMat &DM2);

   int cols(void) 
   {return numcols;}

   int rows(void)
   {return numrows;}

   int mode(void)
   {return packmode;}

   int GetVal(int row, int col);
   void SetVal(int row, int col, int val);

   void SetMode(PackMode newmode, PackMat &DM2);

   void Mult(PackMat &DM, PackMat &Prod);
   void Multnocheck(PackMat &DM, PackMat &Prod);
   void Mult(PackVec &DM, PackVec &Prod);
   void Multnocheck(PackVec &DM, PackVec &Prod);
   void Add(PackMat &DM, PackMat &Sum);
   void Addnocheck(PackMat &DM, PackMat &Sum);

   void transpose(PackMat &T, PackMode newpackmode=HORIZPACK);
   void transposenocheck(PackMat &T, PackMode newpackmode=HORIZPACK);

   int getWeight(void);

   inline static int getWeight(PKSIZE in);

   int isIdentity(void);

   void PrintInfo(void);
   ostream& PackMat:: printpackmat(ostream &os) const ;


   void matinv(PackMat &Inv);
   int ludcmp(int &numind);
   void lubacksub(PackVec &b, PackVec &y);
   void lubacksub(PKSIZE *b, PKSIZE *y);
   int innerprod(PKSIZE *row, PKSIZE *col, int end);
   int innerprod2(PKSIZE *row, PKSIZE *col, int start, int end);

   void reducetosystematic(PackMat &G, int &numind);

   void PackMat::deletecolumn(int colno);
   int PackMat::getrowweight(int rowno);

};


// output stream
inline ostream&
operator<<(ostream &os, const PackMat &r) 
{
	r.printpackmat(os);
	return os;
}


class PackVec
{
   int numrows;
   PKSIZE *Vec;
   int numpackrows;				// number of packed columns
   int numpackrows1;			// number of entire packed columns
   int numpackrows2;			// number of partial packed columns
   int fracbits;				// number of bits left over

public:
   friend class PackMat;  // allow PackMat to access the private stuff

   static int bitsper;			// number of bits in each element

   static void initPackVec();


   PackVec(void);
   PackVec(int rows);
   PackVec(char *fname);
   PackVec(string desc);
   PackVec(istream &Densfile);
   void PackVec_is(istream &Densfile);
   PackVec(PackMat &mat);

   void Size(int rows);
   void Size(PackVec &A);
   void PackVec::MakeRandom(double p=0);
   void PackVec::SetToZero();
   void PackVec::SetToOne();
   inline void checkandset(PackVec &B,int brows);
   inline void checkandset2(PackVec &B,int brows);
   inline void setsize(PackVec &A, int rows);
   void PackVec::clearfringe();

   ~PackVec(void) ;

   PackVec &operator=(PackVec &DM2);

   int rows(void) 
   {return numrows;}

   int GetVal(int row);
   void SetVal(int row, int val);

   void Mult(PackMat &DM, PackVec &Prod);
   void Multnocheck(PackMat &DM, PackVec &Prod);
   void Add(PackVec &DM, PackVec &Sum);
   void Addnocheck(PackVec &DM, PackVec &Sum);

   int getWeight(void);

   inline int getWeight(PKSIZE in);

   void PrintInfo(void);
   ostream& PackVec:: printpackvec(ostream &os) const ;

   int innerprod(PackVec &row1, PackVec &row2, int end);
   int innerprod2(PackVec &row1, PackVec & row2, int start, int end);
   void Subvector(PackVec &subvec,int row1, int row2,int subrow=0);

   inline void PackVec::callocvec(PackVec &A);
   void PackVec::resize(PackVec &A,int rows);


};


// output stream
inline ostream&
operator<<(ostream &os, const PackVec &r) 
{
	r.printpackvec(os);
	return os;
}


void
stackside(PackMat & C,PackMat &A, PackMat & B, StackMode = Top,
		  PackMode pkmode=HORIZPACK);
void
stacktop(PackMat & C,PackMat &A, PackMat & B, StackMode = Left,
		 PackMode pkmode=HORIZPACK);
void setval(PKSIZE **Matrix,int row, int col, int val, 
			PackMode packmode=HORIZPACK);


#endif

   
/*
Local Variables:
compile-command: "g++ -c -g PackMat.cc"
End:
*/
