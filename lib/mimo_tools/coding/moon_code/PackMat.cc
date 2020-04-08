#include "PackMat.h"

using namespace std;

int PackMat::bitsper = 8*sizeof(PKSIZE);
PKSIZE  PackMat::masker = (1<<bitsper)-1;

int PackMat::shifter;				// amount to shift (e.g., >> 5 to /32)
int *PackMat::wtct = NULL;
int PackMat::PkperCt = sizeof(PKSIZE)/sizeof(CTSIZE);
StackMode PackMat::sidestackmode = Top;
StackMode PackMat::topstackmode = Left;

static PackMat U;				// the upper matrix for the LU decomposition
static int numind;				// number of linearly independent rows
static PackMat L;				// the lower matrix for the LU decomposition
static int *luindex= NULL;
static char printspace[5] = " ";

void
PackMat::initPackMat(void )
{
   int nbits = 8*sizeof(CTSIZE);
   int i,j;

   if(wtct) {					// if already allocated
	  delete[] wtct;
   }
   int N = (1<< (8*sizeof(CTSIZE)));// number of different weights to consider
   int wt;
   wtct = new int[N];
   for(i = 0; i < N; i++) {		// for each of these, find the weight
	  wt = 0;
	  for(j = 0; j < nbits; j++) {
		 if(i & (1<<j)) {
			wt++;
		 }
	  }
	  wtct[i] = wt;
   }
   shifter = 0;
   j = bitsper;
   while(j != 1) {
	  j /= 2;
	  shifter++;
   }
}


PackMat::PackMat(void)
{
   numcols = 0;
   numrows = 0;
   numpackcols = 0;
   numpackrows = 0;
   numpackcols1 = numpackcols2 = 0;
   numpackrows1 = numpackrows2 = 0;
   fracbits = 0;
   packmode = HORIZPACK;
   Matrix = 0;
}


PackMat::PackMat(PackMat &D)
{
   int i,j;
   setsize(*this,D.numrows,D.numcols,D.packmode);
   callocmat(*this);
   for(i = 0; i < D.numpackrows; i++) {
	  for(j = 0; j < D.numpackcols; j++) {
		 Matrix[i][j] = D.Matrix[i][j];
	  }
   }
}

inline void
PackMat::setsize(PackMat &A,int rows, int cols, PackMode pkmode)
{
   int i,j;
   A.numcols = cols;
   A.numrows = rows;
   A.packmode = pkmode;
   if(A.packmode == HORIZPACK) {
	  A.numpackrows = A.numpackrows1 = A.numrows;
	  A.numpackrows2 = 0;

	  A.numpackcols1 = A.numpackcols = A.numcols/bitsper;
	  A.numpackcols2 = 0;
	  A.fracbits = A.numcols%bitsper;
	  if(A.fracbits) {
		 A.numpackcols++;			// increment to hold the fractional int
		 A.numpackcols2 = 1;
	  }
   }
   else if(A.packmode == VERTPACK) {
	  A.numpackcols1 = A.numpackcols = A.numcols;
	  A.numpackcols2 = 0;

	  A.numpackrows1 = A.numpackrows = A.numrows/bitsper;
	  A.numpackrows2 = 0;
	  A.fracbits = A.numrows%bitsper;
	  if(A.fracbits) {
		 A.numpackrows++;
		 A.numpackrows2 = 1;
	  }
   }
   else {
	  cerr<<"Error: Invalid packmode"<<endl;
	  exit(-1);
   }
}

inline void
PackMat::callocmat(PackMat &A)
{
   int i,j;
   CALLOCMATRIX(A.Matrix,PKSIZE,A.numpackrows,A.numpackcols);
   // clear out the matrix
   for(i = 0; i < A.numpackrows; i++) {
	  for(j = 0; j < A.numpackcols; j++) {
		 A.Matrix[i][j] = 0;
	  }
   }
}

void
PackMat::resize(PackMat &A,int rows, int cols, PackMode pkmode)
// set A to the right size, but keep all the old stuff that was there before
{
   PKSIZE **NewMatrix;
   int oldnumpackrows, oldnumpackcols;
   int i,j;

   oldnumpackrows = A.numpackrows;
   oldnumpackcols = A.numpackcols;
   setsize(A,rows,cols,pkmode);

   CALLOCMATRIX(NewMatrix,PKSIZE,A.numpackrows,A.numpackcols);
   // clear out the matrix
   for(i = 0; i < oldnumpackrows; i++) {
	  for(j = 0; j < oldnumpackcols; j++) {
		 NewMatrix[i][j] = A.Matrix[i][j];
	  }
   }
   if(A.Matrix) FREEMATRIX(A.Matrix);
   A.Matrix = NewMatrix;
}



PackMat::PackMat(int rows, int cols, PackMode pkmode)
{
   int i,j;

   double temp1, temp2;

   setsize(*this,rows, cols, pkmode);
   callocmat(*this);
}

int
PackMat::ReadFile(char *fname,PackMode rowcolpack)
{
   ifstream Densefile(fname);
   if(!Densefile) {
	  cerr<<"Error: Unable to open input file "<<fname<<endl;
	  return 0;
   }
   PackMat_is(Densefile,rowcolpack);
   Densefile.close();
   return 1;
}

PackMat::PackMat(char *fname,PackMode rowcolpack)
{
   Matrix = 0;
   ifstream Densefile(fname);
   if(!Densefile) {
	  cerr<<"Error: Unable to open input file"<<fname<<endl;
	  exit(-1);
   }
   PackMat_is(Densefile,rowcolpack);
   Densefile.close();
}

PackMat::PackMat(string desc,PackMode rowcolpack)
{

   Matrix = 0;
   istringstream Densefile(desc);
   PackMat_is(Densefile,rowcolpack);
}


void
PackMat::fill(string desc)
{

   istringstream Densefile(desc);
   fill(Densefile);
}

PackMat::PackMat(istream &Densefile,PackMode rowcolpack)
{
   Matrix = 0;
   PackMat_is(Densefile,rowcolpack);
}

void
PackMat::PackMat_is(istream &Densefile,PackMode rowcolpack)
{

   Densefile>>numrows;
   Densefile>>numcols;
   setsize(*this,numrows, numcols,rowcolpack);
   if(Matrix) FREEMATRIX(Matrix);
   callocmat(*this);

   fill(Densefile);
}

void
PackMat::fill(istream &Densefile)
{
   int i, j, bit, k1;
   PKSIZE blob;

   if(packmode == HORIZPACK) {
	  for(i=0; i<numrows; i++) {
		 for(j=0;j< numpackcols1;j++) {
			blob = 0;
			for(k1 = 0; k1 < bitsper; k1++) {
			   Densefile >> bit;
//cout << bit << " ";
				blob = (bit << k1) | blob;
			}
			Matrix[i][j] = blob;
		 }
		 if(fracbits) {
			j = numpackcols1;
			blob = 0;
			for(k1 = 0; k1 < fracbits; k1++) {
			   Densefile >> bit;
//cout << bit << " ";
			   blob = (bit << k1) | blob;
			}
			Matrix[i][j] = blob;
		 }
	  }
   }
   else if(packmode == VERTPACK) {
	  for(i=0;i< numpackrows1;i++) {
		 for(j = 0; j < numcols; j++) Matrix[i][j] = 0;
		 for(k1 = 0; k1 < bitsper; k1++) {
			for(j = 0; j < numcols; j++) {
			   Densefile >> bit;
			   Matrix[i][j] = (bit << k1) | Matrix[i][j];
			}
		 }
	  }
	  if(fracbits) {
		 i = numpackrows1;
		 for(j = 0; j < numcols; j++) Matrix[i][j] = 0;
		 for(k1 = 0; k1 < fracbits; k1++) {
			for(j = 0; j < numcols; j++) {
			   Densefile >> bit;
			   Matrix[i][j] = (bit << k1) | Matrix[i][j];
			}
		 }
	  }
   }
}


PackMat::PackMat(PackVec &vec, PackMode rowcolpack)
{
   int i;
   if(rowcolpack == HORIZPACK) { // pack horizontally; make row vector
	  setsize(*this,1, vec.numrows,HORIZPACK);
	  callocmat(*this);
	  for(i = 0; i < numpackcols; i++) {
		 Matrix[0][i] = vec.Vec[i];
	  }
   }
   else {						// pack vertically; make a column vector
	  setsize(*this,vec.numrows,1,VERTPACK);
	  callocmat(*this);
	  for(i = 0; i < numpackrows; i++) {
		 Matrix[i][0] = vec.Vec[i];
	  }
   }
}


void
PackMat::Size(int rows, int cols, PackMode pkmode)
// resize, but do not attempt to preserve the contents
{
   int i,j;

   checkandset(*this,pkmode,rows,cols);
   for(i = 0; i < numpackcols; i++) {
	  for(j = 0; j < numpackrows; j++) {
		 Matrix[i][j] = 0;
	  }
   }
}

void
PackMat::Size(PackMat &A)
// resize to the same size as A, but do not attempt to preserve the contents
{
   Size(A.numrows,A.numcols,A.packmode);
}

PackMat &PackMat::operator=(PackMat &DM2)
{
   int i,j;

   checkandset(*this,DM2.packmode,DM2.numrows, DM2.numcols);
   
   for(i=0;i<numpackrows;i++)
	  for(j=0;j<numpackcols;j++)
		 Matrix[i][j] = DM2.Matrix[i][j];

   return *this;
}

PackMat::~PackMat(void)
{
   if(Matrix) FREEMATRIX(Matrix);
}

int PackMat::GetVal(int row, int col)
{
   int i,j,k;
   PKSIZE value;
   int col1, row1;
   int vb;
   int bitcol, bitrow;

   if(row >= numrows || col >= numcols) {
	  cerr<<"Error: Index exceeds matrix dimensions"<<endl;
	  exit(-1);
   }

   value = 0;

   if(packmode == HORIZPACK) {
	  col1 = col >> shifter;	// determine the integer number
	  value = Matrix[row][col1];

	  bitcol = col % bitsper;
	  vb = (value >> bitcol) & 1;
   }
   else if(packmode == VERTPACK) {
	  row1 = row >> shifter;
	  value = Matrix[row1][col];
	  bitrow = row % bitsper;
	  vb =  (value >> bitrow) & 1;
   }
   return vb;
}

void PackMat::SetVal(int row, int col, int val)
{
   int i,j,k, col1, row1;
   int bitcol, bitrow;
   PKSIZE value;
   
   if(row >= numrows || col >= numcols) {
	  cerr<<"Error: Index exceeds matrix dimensions"<<endl;
	  exit(-1);
   }
   
//    if(val != 0 && val != 1) {
// 	  cerr<<"Error: Illegal input value"<<endl;
// 	  cout<<"val = "<<val<<endl;
// 	  exit(-1);
//    }
   
   if(packmode == HORIZPACK) {
	  col1 = col >> shifter;
	  bitcol = col % bitsper;
	  if(!val) {					// set value to 0
		 Matrix[row][col1] &= ~(1<<bitcol);
	  }
	  else {						// set value to 1
		 Matrix[row][col1] |= (1<<bitcol);
	  }
   }
   
   else if(packmode == VERTPACK) {
	  row1 = row >> shifter;
	  bitrow = row % bitsper;
	  if(!val) {				// set value to 0
		 Matrix[row1][col] &= ~(1<<bitrow);
	  }
	  else {
		 Matrix[row1][col] |= (1<<bitrow);
	  }
   }
}



void PackMat::deletecolumn(int colno)
{
   int i,j;
   if(colno < 0 || colno >= numcols) return;
   if(colno == numcols-1) {
	  for(i = 0; i < numrows; i++) {
		 SetVal(i,colno,0);
	  }
   }
   else {
	  for(j = colno; j < numcols-1; j++) {
		 for(i = 0; i < numrows; i++) {
			SetVal(i,j,GetVal(i,j+1));
		 }
	  }
   }
   setsize(*this,numrows,numcols-1,packmode);
}


int PackMat::getrowweight(int rowno)
{
   int i;
   int wt = 0;
   if(rowno < 0 || rowno >= numrows) return 0;

   if(packmode == HORIZPACK) {	// the easy way
	  for(i = 0; i < numpackcols; i++)
		 wt += getWeight(Matrix[rowno][i]);
	  return wt;
   }
   else {
	  for(i = 0; i < numcols; i++) {
		 wt += GetVal(rowno,i);
	  }
   }
   return wt;
}

void PackMat::SetMode(PackMode newmode,PackMat &DM2)
{

   int i,j,k;
   int numr, numc;
   PKSIZE **NewMat;

   if(packmode == newmode) { 	// keep the mode
	  DM2 = *this;
	  return;
   }

   CALLOCMATRIX(NewMat,PKSIZE,numrows,numcols);

   for(i=0;i<numrows;i++) {
	  for(j=0;j<numcols;j++) {
		 k = GetVal(i,j);
		 setval(NewMat,i,j,k,newmode);
	  }
   }
   setsize(DM2,numrows,numcols,newmode);
   if(DM2.Matrix) FREEMATRIX(DM2.Matrix);
   DM2.Matrix = NewMat;
}

void PackMat::Submatrix(PackMat &submat,int row1, int row2, int col1, int col2,
						int subrow, int subcol,PackMode pkmode)
// extract the rectangular matrix, placing it starting at (subrow,subcol)
// in the submatrix, packing in the pkmode order
// if out of range in any direction, fill in with zeros
{
   int i,i1, j,j1, nrows, ncols;
   if(row1 < 0) row1 = 0;
   if(row1 >= numrows) row1 = numrows-1;
   if(col1 < 0) col1 = 0;
   if(col2 >= numcols) col2 = numcols-1;
   nrows = (row2-row1+1) + subrow;
   ncols = (col2-col1+1) + subcol;
   checkandset2(submat,pkmode,nrows,ncols);
   for(i=0, i1=subrow; i < (row2-row1+1); i++,i1++) {
	  for(j = 0, j1 = subcol; j < (col2-col1+1); j++,j1++) {
		 submat.SetVal(i1,j1,GetVal(i+row1,j+col1));
	  }
   }
}

void PackMat::Subvector(PackVec &subvec,int row1, int row2,int col1, int col2,
						int subrow) 
// extract the row or vector and place it (starting at subrow)
// in the vector subvec
{
   int i,i1, j,j1, nrows, ncols;
   if(row1 < 0) row1 = 0;
   if(row1 >= numrows) row1 = numrows-1;
   if(col1 < 0) col1 = 0;
   if(col2 >= numcols) col2 = numcols-1;

   if(row1==row2) {				// extract a row
	  ncols = (col2-col1+1) + subrow;
	  subvec.checkandset2(subvec,ncols);
	  for(i=0, i1=subrow; i < (col2-col1+1); i++,i1++) {
		 subvec.SetVal(i1,GetVal(row1,i+col1));
	  }
   }
   else if(col1==col2) {		// extract a column
	  nrows = (row2-row1+1) + subrow;
	  subvec.checkandset2(subvec,nrows);
	  for(i=0, i1=subrow; i < (row2-row1+1); i++,i1++) {
		 subvec.SetVal(i1,GetVal(i+row1,col1));
	  }
   }
   else {
	  cerr << "Either row or column must be constant\n";
	  exit(-1);
   }
}


void PackMat::Mult(PackMat &DM, PackMat &Prod)
{
// This function tests to make sure things are the right size
// and then calls the multiplication function

   if(packmode == HORIZPACK && DM.packmode == VERTPACK) {
	  if(numcols != DM.numrows) {  // correct size to multiply
		 cerr << "Error: Matrices not conformable in product\n";
		 exit(-1);
	  }
	  checkandset(Prod,HORIZPACK,numrows,DM.numcols);
	  Multnocheck(DM,Prod);
   }
   else {
	  cerr << "Multiplication not implemented for these matrix modes\n";
	  exit(-1);
   }
}


void PackMat::Multnocheck(PackMat &DM, PackMat &Prod)
{
// This function does not testing --- it assumes that
// everything is the right size and shape for multiplication

   int i,j,k,nj,j1,k1;
   int prodwt;

   nj = DM.numcols/bitsper;
   for(i = 0; i < numrows; i++) {
// cout << i << " " << flush;
	  // do all the ones that have the full 32 bits
	  j1 = 0;
	  for(j = 0; j < nj; j++) {
		 Prod.Matrix[i][j] = 0;
		 for(k1 = 0; k1 < bitsper; k1++) {
			prodwt = 0;
			for(k = 0; k < numpackcols;k ++) {
			   prodwt += getWeight(Matrix[i][k] & DM.Matrix[k][j1]);
			}
			Prod.Matrix[i][j] = ((prodwt % 2)<<k1) + Prod.Matrix[i][j];
			j1++;			   
		 }
	  }
	  // do the last fractbits
	  if(Prod.fracbits) {
		 j = nj;
		 Prod.Matrix[i][j] = 0;
		 for(k1 = 0; k1 < Prod.fracbits; k1++) {
			prodwt = 0;
			for(k = 0; k < numpackcols;k ++) {
			   prodwt += getWeight(Matrix[i][k] & DM.Matrix[k][j1]);
			}
			Prod.Matrix[i][j] = ((prodwt % 2)<<k1) + Prod.Matrix[i][j]; 
			j1++;
		 }
	  }
   }
}


void PackMat::Mult(PackVec &DM, PackVec &Prod)
// Multiply a matrix times a vector
{
// This function tests to make sure things are the right size
// and then calls the multiplication function

   if(packmode == HORIZPACK) {
	  if(numcols != numrows) {  // correct size to multiply
		 cerr << "Error: Matrices not conformable in product\n";
		 exit(-1);
	  }
	  Prod.checkandset(Prod,numrows);
	  Multnocheck(DM,Prod);
   }
   else {
	  cerr << "Multiplication not implemented for these matrix modes\n";
	  exit(-1);
   }
}


void PackMat::Multnocheck(PackVec &DM, PackVec &Prod)
{
// This function does not testing --- it assumes that
// everything is the right size and shape for multiplication

   int i,k,ibig,ilit;
   int prodwt;

   for(i = 0; i < DM.numpackrows; i++) {Prod.Vec[i] = 0;}
   for(i = 0; i < numrows; i++) {
	  ibig = i/bitsper;
	  ilit = i % bitsper;
	  prodwt = 0;
	  for(k = 0; k < numpackcols;k ++) {
		 prodwt += getWeight(Matrix[i][k] & DM.Vec[k]);
	  }
	  Prod.Vec[ibig] = ((prodwt % 2)<<ilit) + Prod.Vec[ibig];
   }
}


inline void
PackMat::checkandset(PackMat &B,PackMode bpackmode,int brows, int bcols)
// set up exactly brows and bcols in B, if it is not already there
{
   if(!(B.packmode == bpackmode && B.numrows == brows && B.numcols == bcols)) {
	  if(B.Matrix) {FREEMATRIX(B.Matrix);}
	  setsize(B,brows, bcols, bpackmode);
	  callocmat(B);
   }	  
}

inline void
PackMat::checkandset2(PackMat &B,PackMode bpackmode,int brows, int bcols)
// make sure there is at least brows and bcols space,
// but possibly more
{
   if(B.packmode != bpackmode || brows > B.numrows || bcols > B.numcols) {
	  // wrong mode, or not enough room already
	  if(B.Matrix) {FREEMATRIX(B.Matrix);}
	  setsize(B,MAX(brows,B.numrows), MAX(bcols,B.numcols), bpackmode);
	  callocmat(B);
   }	  
}

void
PackMat::checkandresize2(PackMat &B,PackMode bpackmode,int brows, int bcols)
// make sure there is at least brows and bcols space,
// but possibly more
// If necessary, resize, but keep the old stuff
{
   if(B.packmode != bpackmode || brows > B.numrows || bcols > B.numcols) {
	  // wrong mode, or not enough room already
//	  if(B.Matrix) {FREEMATRIX(B.Matrix);}
	  resize(B,MAX(brows,B.numrows), MAX(bcols,B.numcols), bpackmode);
   }	  
}




void PackMat::transpose(PackMat &T,PackMode newpackmode)
{
   int numc, numr;
   PKSIZE **NewMat;

   numc = numcols;
   numr = numrows;

   CALLOCMATRIX(NewMat,PKSIZE,numcols, numrows);

   int i, j;
   if(packmode != newpackmode) { // switching mode
	  for(i = 0; i < numpackrows; i++) {
		 for(j = 0; j < numpackcols; j++) {
			NewMat[j][i] = Matrix[i][j];
		 }
	  }
   }
   else {						// same mode
	  for(i = 0; i < numrows; i++) {
		 for(j = 0;j < numcols; j++) {
			setval(NewMat,j,i,GetVal(i,j));
		 }
	  }
   }
   setsize(T,numcols, numrows, newpackmode);
   if(T.Matrix) FREEMATRIX(T.Matrix);
   T.Matrix = NewMat;
}




void PackMat::clearfringe()
{
   int i,j;
   PKSIZE mask;
   if(packmode == HORIZPACK) {
	  mask = (1<<fracbits)-1;
	  for(i = 0; i < numrows; i++) {
		 Matrix[i][numpackcols1] &= mask;
	  }
   }
   else if(packmode == VERTPACK && fracbits) {
	  mask = (1<<fracbits)-1;
	  for(j = 0; j < numcols; j++) {
		 Matrix[numpackrows1][j] &= mask;
	  }
   }
}

int PackMat::getWeight() 		// get the weight of the entire matrix
{
   int i,j, wt;

   wt = 0;
   if(fracbits) {				// make sure there is no extra stuff on edges
	  clearfringe();
   }
   for(i = 0; i < numpackrows; i++) {
	  for(j = 0; j < numpackcols; j++) {
		 wt += getWeight(Matrix[i][j]);
	  }
   }
   return wt;
}

inline
int PackMat::getWeight(PKSIZE blob)
{
   int wt;
   int i;
   CTSIZE *ptr;					// point to the pieces of blob

   ptr = (CTSIZE *)&blob;
   wt = 0;
   for(i = 0; i < PkperCt; i++) {
	  wt += wtct[*ptr];
	  ptr++;					// point to the next bit
   }
   return wt;
}


// print function
ostream&
PackMat:: printpackmat(ostream &os) const 
{
   int i,j,k,i1;


   if(packmode == HORIZPACK) {
	  for(i = 0; i< numrows; i++) {
		 for(j = 0; j < numpackcols1; j++) {
			for(k = 0; k < bitsper; k++) {
			   if(Matrix[i][j] & (1<<k))
				  os << "1";
			   else
				  os << "0";
			   os << printspace;
			}
		 }
		 for(j = 0; j < numpackcols2; j++) {
			for(k = 0; k < fracbits; k++) {
			   if(Matrix[i][j+numpackcols1] & (1<<k))
				  os << "1";
			   else
				  os << "0";
			   os << printspace;
			}
		 }
		 os << endl;
	  }
   }
   else if(packmode == VERTPACK) {
	  for(i = 0; i < numpackrows1; i++) {
		 for(i1 = 0; i1 < bitsper; i1++) {
			for(j = 0; j < numcols; j++) {
			   if(Matrix[i][j] & (1<<i1)) 
				  os << "1";
			   else
				  os << "0";
			   os << printspace;
			}
			os << endl;
		 }
	  }
	  for(i = 0; i < numpackrows2; i++) {
		 for(i1 = 0; i1 < fracbits; i1++) {
			for(j = 0; j < numcols; j++) {
			   if(Matrix[i+numpackrows1][j] & (1<<i1)) 
				  os << "1";
			   else
				  os << "0";
			   os << printspace;
			}
			os << endl;
		 }
	  }
   }
   return os;
}


void PackMat::PrintInfo(void)
{
   
   if(packmode == HORIZPACK)
	  cout<<"Pack Mode = Column Pack"<<endl;
   else if(packmode == VERTPACK)
	  cout<<"PackMode = Row Pack"<<endl;
   cout<<"Rows = "<<numrows<<",  Cols = "<<numcols<<","<<endl;
   cout<<"Packrows: "<< numpackrows<< "  " << numpackrows1 << "  " <<
	  numpackrows2 << endl;
   cout<<"Packcols: "<< numpackcols<< "  " << numpackcols1 << "  " <<
	  numpackcols2 << endl;
   cout << "fracbits: " << fracbits << endl;

   cout << "bitsper=" << bitsper << "  PkperCt=" << PkperCt << 
	  "  shifter=" << shifter << endl;
}

void
PackMat::MakeRandom(double p)
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;
   PKSIZE mask;

   if(p == 0 || p == 0.5) {	// if equal prob of 0 or 1
	  for(i = 0; i < numpackrows; i++) {
		 for(j = 0; j < numpackcols; j++) {
			blob = 0;
			for(k = 0; k < sizeof(PKSIZE); k++) {
			   rno = 0xFF&(rand()>>8); // pick off a byte
			   blob = (blob << 8) + rno;
			}
			Matrix[i][j] = blob;
		 }
	  }
	  // now make sure there is nothing left over 
	  if(fracbits) {
		 clearfringe();
	  }
   }
   else {						// else go by the probabilities
	  double rnod;
	  for(i = 0; i < numrows; i++) {
		 for(j = 0; j < numcols;  j++) {
			rnod = double(rand())/double(RAND_MAX);
			if(rnod < p)
			   SetVal(i,j,1);
			else
			   SetVal(i,j,0);
		 }
	  }
   }
}


void
PackMat::SetToZero()
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;

   for(i = 0; i < numpackrows; i++) {
	  for(j = 0; j < numpackcols; j++) {
		 Matrix[i][j] = 0;
	  }
   }
}

void
PackMat::SetToOne()
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;

   for(i = 0; i < numpackrows; i++) {
	  for(j = 0; j < numpackcols; j++) {
		 Matrix[i][j] = masker;
	  }
   }
   // now make sure there is nothing left over 
   clearfringe();
}



void PackMat::Add(PackMat &DM, PackMat &Sum)
{
// This function tests to make sure things are the right size
// and then calls the addition function

   if(packmode == DM.packmode) {
	  if(!(numcols == DM.numcols && numrows == DM.numrows)) {  
          // correct size to add
		 cerr << "Error: Matrices not same size in sum\n";
		 exit(-1);
	  }
	  checkandset(Sum,packmode,numrows,numcols);
	  Addnocheck(DM,Sum);
   }
   else {
	  cerr << "Addition not implemented for these matrix modes\n";
	  exit(-1);
   }
}


void PackMat::Addnocheck(PackMat &DM, PackMat &Sum)
{
// This function does no testing --- it assumes that
// everything is the right size and shape for addition

   int i,j;
   PKSIZE temp;
   int prodwt;

   for(i = 0; i < numpackrows; i++) {
	  for(j = 0; j < numpackcols; j++) {
		 Sum.Matrix[i][j] = Matrix[i][j] ^ DM.Matrix[i][j];
	  }
   }
}

void
PackMat::Identity(int n, PackMode pkmode)
{
   int i,j;

   checkandset(*this,pkmode,n,n);
   for(i = 0; i < numpackrows; i++) {
	  for(j = 0; j< numpackcols; j++) {
		 Matrix[i][j] = 0;
	  }
   }
   for(i = 0; i < n; i++) {
	  SetVal(i,i,1);
   }
}

int PackMat::IsIdentity(void)
{
   int i;
   if(numrows != numcols) return 0;
   int wt = getWeight();
   if(wt != numrows) return 0;
   for(i = 0; i < numrows; i++) {
	  if(GetVal(i,i) != 1) return 0;
   }
   return 1;
}


int
PackMat::ludcmp(int &numindout)
// returns 0 if matrix is singular, and 1 if not singular
// numindout is the number of linearly independent columns
{
   int i, j, k,k1,j1, jbig, jlit;
   PKSIZE tmp2,jlitmask;
   int tmp1;
   int pivotfound;

   if(packmode != HORIZPACK) {
	  cerr << "LU decomposition for horizontally packed data only\n";
	  exit(-1);
   }
   if(luindex) delete[] luindex;
   luindex = new int[numrows];
   for(j = 0; j < numrows; j++) {
	  luindex[j] = j;
   }

   U = *this;
   L.Size(numrows,numrows,HORIZPACK);
   L.SetToZero();

//   for(j = 0; j < numcols-1; j++) {		// loop over columns
   for(j = 0; j < numcols; j++) {		// loop over columns
	  jbig = j/bitsper;
	  jlit = j % bitsper;
	  jlitmask = (1<<jlit);
	  // determine a pivot row
	  pivotfound = 0;
	  for(k1 = j; k1 < numrows; k1++) {
		 if(U.Matrix[k1][jbig]&jlitmask) {// nonzero bit --- pivot on this
			pivotfound = 1;
			break;
		 }
	  }
	  if(!pivotfound) {
//		 numind = MIN(numrows,numcols);
		 numind = j;
		 numindout = numind;
// cout <<"Bailing: numind(1)=" << numind << "\n";
// cout << "L\n" << L << endl;
// cout << "U\n" << U << endl;
// for(j = 0; j < numrows; j++) {
//  L.SetVal(j,j,1);
//  }
//  PackMat Uc;
//  U.SetMode(VERTPACK,Uc);
//  PackMat LU;
//  L.Mult(Uc,LU);
// cout << "LU\n" << LU << endl;
		 return 0;
	  }
//cout << "j=" << j << "  pivotrow=" << k1 << endl;
	  if(j != k1) {				// necessary to swap rows j and k1
		 for(j1 = jbig; j1 < numpackcols; j1++) {
			tmp2 = U.Matrix[k1][j1];
			U.Matrix[k1][j1] = U.Matrix[j][j1];
			U.Matrix[j][j1] = tmp2;
		 }
		 for(j1 = 0; j1 <= jbig; j1++) {
			tmp2 = L.Matrix[k1][j1];
			L.Matrix[k1][j1] = L.Matrix[j][j1];
			L.Matrix[j][j1] = tmp2;
		 }

//		 cout << "swapped rows\n" << U << endl;
		 tmp1 = luindex[j];
		 luindex[j] = luindex[k1];
		 luindex[k1] = tmp1;
	  }
	  for(i = j+1; i < numrows; i++) {
		 if(U.Matrix[i][jbig]&jlitmask) {
//			cout << "l("<<i<<","<<j<<")=" << 1 <<endl;
			L.Matrix[i][jbig] |= (jlitmask | L.Matrix[i][jbig]);
			// row i <- row i - mult * row j
			for(k = jbig/bitsper; k < numpackcols; k++) {
			   U.Matrix[i][k] = U.Matrix[i][k] ^ U.Matrix[j][k];
			}
		 }
	  }
// cout << "j=" << j << endl << "U=\n" << U << endl;
// cout << "L=\n" << L << endl;

   }
//    for(j = 0; j < numrows; j++) {
// 	  L.SetVal(j,j,1);
//    }
//  PackMat Uc;
//  U.SetMode(VERTPACK,Uc);
//  PackMat LU;
//  L.Mult(Uc,LU);
// cout << "LU\n" << LU << endl;

   numindout = numind = MIN(numrows,numcols);
   return 1;
}	 


void
PackMat::lubacksub(PackVec &b,PackVec &y)
{  
   // make sure there is enough room for y
   y.checkandset(y,numrows);
   lubacksub(b.Vec,y.Vec);
}


void PackMat::lubacksub(PKSIZE *b,PKSIZE *y)
// b is a packed one-dimensional array representing a RHS
{
   int i,ibig,ilit,sum;
   int indx,indxbig, indxlit;
   int bi;
   for(i = 0; i < numpackcols; i++) y[i] = 0;
   // do the backsubstitution to solve Ly = Pb
   for(i = 0; i < numrows; i++) { 
	  ibig = i/bitsper;
	  ilit = i % bitsper;
	  indx = luindex[i];
	  indxbig = indx/bitsper;
	  indxlit = indx % bitsper;
	  bi = ((b[indxbig]&(1<<indxlit)))>>indxlit;
	  sum = (bi + innerprod(L.Matrix[i],y,i)) % 2;  
//cout << "i=" << i << "  bi=" << bi << "  y=" << sum << endl;
	  if(sum) {
		 y[ibig] |= (1<<ilit);
	  }
   }
//cout << endl;
//cout << endl;
   // do forward subsitution to solve Uy = x
   i = numrows-1; 				// do the last row by itself
   ibig = i/bitsper;
   ibig = i % bitsper;
   bi = (y[ibig] & (1<<ilit)) >> ilit;
   if(bi)
	  y[ibig] |= (1<<ilit);
   else
	  y[ibig] &= ~(1<<ilit);
//cout << "i=" << numrows-1 << "  bi=" << bi << "  y=" << int(y[ibig]) << endl;

   for(i = numrows-2; i>= 0; i--) {	// now do all the rest of the rows
//cout << "i=" << i << " ";
	  ibig = i/bitsper;
	  ilit = i % bitsper;
	  bi = (y[ibig] & (1<<ilit)) >> ilit;
	  sum = (bi + innerprod2(U.Matrix[i],y,i+1,numrows-1)) % 2;
//cout << "i=" << i << "  yi=" << bi << "  newy=" << sum << endl;
	  if(sum)
		 y[ibig] |= (1<<ilit);
	  else
		 y[ibig] &= ~(1<<ilit);
   }

}

int PackMat::innerprod2(PKSIZE *row, PKSIZE *col, int start, int end)
// do the inner product between two row-packed quantities, 
// ending at the last bit
{
   int stbig, stlit, ebig, elit, i, nl;
   int sum=0;
   PKSIZE mask,mask1;
   int bitsdone=0;

   stbig = start/bitsper;
   stlit = start % bitsper;
   ebig = end/bitsper;
   elit = end % bitsper;
   mask1 = ((1<<(bitsper-stlit))-1) << stlit;
   mask = mask1;

//   cout << "start=" << start << "  end=" << end << endl;

   for(i = stbig; i < ebig; i++) {
//cout<<"a: i=" << i << "   row=" << int(row[i]) << "  col=" << int(col[i]) << 
//  "  mask=" << int(mask) << "  ip=" << getWeight(row[i]&col[i]&mask) << endl;
	  sum += getWeight(row[i] & col[i] & mask);
	  mask = masker;
	  stlit = 0;
   }
   nl = elit-stlit+1;
   mask = ((1<<nl)-1)<<stlit;
//cout<<"b: i=" << i << "  mask=" << int(mask) << "   row=" << int(row[i]) 
//<<"  col="<<int(col[i]) <<  "  ip=" << getWeight(row[i]&col[i]&mask) << endl;
 sum += getWeight(row[i] & col[i] & mask);
   return sum % 2;
}

int PackMat::innerprod(PKSIZE *row, PKSIZE *col, int end)
// do the inner product between two row-packed quantities, starting
// at the first bit and covering end bits.
{
   int i, ebig, elit;
   int sum= 0;

   ebig = end/bitsper; elit = end % bitsper;
   for(i = 0; i < ebig; i++) {
	  sum += getWeight(row[i] & col[i]);
   }
   if(elit) {
	  sum += getWeight((row[i]&col[i]) & ((1<<elit)-1));
   }
//cout << "ip:  row=" << int(row[0]) << " col=" << int(col[0]) << "  sum=" << 
//   int(sum) << "  elit=" << elit << "  mask=" << ((1<<elit)-1) << endl;
   return sum % 2;
}



void 
PackMat::matinv(PackMat &Inv)
{
   int i,ibig,ilit,j;
   int numind;

   setsize(Inv,numrows,numcols,VERTPACK);
   callocmat(Inv);

   if(packmode != HORIZPACK) {
	  cerr << "Matrix must be horizontal packed";
	  exit(-1);
   }
   PKSIZE *b = new PKSIZE[numpackcols+1];
   PKSIZE *y = new PKSIZE[numpackcols+1];
   if(!ludcmp(numind)) {
	  cerr << "Cannot invert matrix\n";
	  exit(-1);
   }
   for(i = 0; i < numpackcols+1; i++) {
	  b[i] = y[i] = 0;
   }
   for(i = 0; i < numrows; i++) {
	  ibig = i/bitsper;
	  ilit = i % bitsper;
	  b[ibig] = (1<<ilit);
	  lubacksub(b,y);
	  // copy into the columns of the inverse
	  for(j = 0; j < numpackcols; j++) {
		 Inv.Matrix[j][i] = y[j];
	  }
	  b[ibig] = 0;
   }
}

void 
PackMat::reducetosystematic(PackMat &G, int &numind)
// Given a matrix, reduce to a systematic form.  The matrix
// returned may have row permutations from the original
// The number of linearly independent rows (the dimension of the code)
// is returned in numind
{
   int i,j,k;
   PKSIZE bi;
   int jbig, jlit;

   G = *this;
   G.ludcmp(numind);
//   cout << "numind=" << numind << endl;

   // now work with the upper matrix
   G = U;
//cout << "G=\n" << G << endl;
   for(j = numind-1; j > 0; j--) {	// column
	  jbig = j/bitsper;
	  jlit = j % bitsper;
	  for(i = j-1; i >= 0; i--) { // row
		 bi = G.Matrix[i][jbig] & (1<<jlit);
//cout << "i=" << i << "  j=" << j << "  " << int(bi) << 
//			"  G[i][jbig]=" << int(G.Matrix[i][jbig]) << endl;
		 if(bi) {				// row[i] <- row[i] + row[j]
			for(k = jbig; k < numpackcols; k++) {
//cout << "k=" << k << " ";
			   G.Matrix[i][k] ^= G.Matrix[j][k];
			}
//cout << endl;
		 }
	  }
//	  cout << "j=" << j << "G=\n" << G << endl;
   }
}

//**********************************************************************


void
PackVec::initPackVec(void )
{
   PackMat::initPackMat();
}

PackVec::PackVec(void)
{
   numrows = 0;
   numpackrows = 0;
   numpackrows1 = numpackrows2 = 0;
   fracbits = 0;
   Vec = 0;
}

inline void
PackVec::setsize(PackVec &A,int rows)
{
   int i,j;
   A.numrows = rows;

   A.numpackrows1 = A.numpackrows = A.numrows/PackMat::bitsper;
   A.numpackrows2 = 0;
   A.fracbits = A.numrows%PackMat::bitsper;
   if(A.fracbits) {
	  A.numpackrows++;			// increment to hold the fractional int
	  A.numpackrows2 = 1;
   }
}

inline void
PackVec::callocvec(PackVec &A)
{
   int i;
   A.Vec = new PKSIZE[A.numpackrows];
   // clear out the matrix
   for(i = 0; i < A.numpackrows; i++) {
	  A.Vec[i] = 0;
   }
}



PackVec::PackVec(int rows)
{
   int i,j;

   double temp1, temp2;

   setsize(*this,rows);
   callocvec(*this);
}


PackVec::PackVec(PackMat &mat)
{
   int i;

   if(mat.numrows==1) {			// copy the first row into the vector
	  setsize(*this,mat.numcols);
	  Vec = new PKSIZE[numpackrows];
	  for(i = 0; i < numpackrows; i++) {
		 Vec[i] = mat.Matrix[0][i];
	  }
   }
   else if(mat.numcols == 1) {	// copy the first column into the vector
	  setsize(*this, mat.numrows);
	  Vec = new PKSIZE[numpackrows];
	  for(i = 0; i < numpackrows; i++) {
		 Vec[i] = mat.Matrix[i][0];
	  }
   }
   else {						// ambiguous
	  cerr << "Cannot uniquely build vector from matrix\n";
	  exit(-1);
   }

}


PackVec::PackVec(char *fname)
{

   ifstream Densefile(fname);
   if(!Densefile) {
	  cerr<<"Error: Unable to open input file"<<fname<<endl;
	  exit(-1);
   }
   PackVec_is(Densefile);
   Densefile.close();
}

PackVec::PackVec(string desc)
{
   istringstream Densefile(desc);
   PackVec_is(Densefile);
}

PackVec::PackVec(istream &Densefile)
{
   PackVec_is(Densefile);
}

void
PackVec::PackVec_is(istream &Densefile)
{
   int i,j,k1;
   PKSIZE blob;
   int bit;

   Densefile>>numrows;

   setsize(*this,numrows);
   callocvec(*this);

   for(i=0;i< numpackrows1;i++) {
	  for(k1 = 0; k1 < PackMat::bitsper; k1++) {
		 Densefile >> bit;
		 Vec[i] = (bit << k1) | Vec[i];
	  }
   }
   if(fracbits) {
	  i = numpackrows1;
	  for(k1 = 0; k1 < fracbits; k1++) {
		 Densefile >> bit;
		 Vec[i] = (bit << k1) | Vec[i];
	  }
   }
}

void
PackVec::Size(int rows)
// resize, but do not attempt to preserve the contents
{
   int i,j;

   checkandset(*this,rows);
   for(j = 0; j < numpackrows; j++) {
	  Vec[i] = 0;
   }
}

void
PackVec::Size(PackVec &A)
// resize to the same size as A, but do not attempt to preserve the contents
{
   Size(A.numrows);
}

PackVec &PackVec::operator=(PackVec &DM2)
{
   int i,j;

   checkandset(*this,DM2.numrows);
   
   for(i=0;i<numpackrows;i++)
	  Vec[i] = DM2.Vec[i];
   return *this;
}

PackVec::~PackVec(void)
{
   numrows = numpackrows = numpackrows1 = numpackrows2 = 0;
   
   if(Vec) delete[] Vec;
   Vec = 0;
}

int PackVec::GetVal(int row)
{
   int vb, row1, row2;

   row1 = row >> PackMat::shifter;
   row2 = row % PackMat::bitsper;

   vb = (Vec[row1] & (1<<row2)) >> row2;
   return vb;
}

void PackVec::SetVal(int row, int val)
{
   int i,j,k, col1, row1;
   int bitcol, bitrow;
   PKSIZE value;
   

   row1 = row >> PackMat::shifter;
   bitrow = row % PackMat::bitsper;
   if(!val) {				// set value to 0
	  Vec[row1] &= ~(1<<bitrow);
   }
   else {
	  Vec[row1] |= (1<<bitrow);
   }
}


void PackVec::Mult(PackMat &DM, PackVec &Prod)
{
// This function tests to make sure things are the right size
// and then calls the multiplication function

   if(DM.packmode == VERTPACK) {
	  if(numrows != DM.numrows) {  // correct size to multiply
		 cerr << "Error: Matrices not conformable in product\n";
		 exit(-1);
	  }
	  checkandset(Prod,DM.numcols);
	  Multnocheck(DM,Prod);
   }
   else {
	  cerr << "Multiplication not implemented for these matrix modes\n";
	  exit(-1);
   }
}


void PackVec::Multnocheck(PackMat &DM, PackVec &Prod)
{
// This function does not testing --- it assumes that
// everything is the right size and shape for multiplication

   int i,j,k,nj,j1,k1;
   PKSIZE temp;
   int prodwt;

   nj = DM.numcols/PackMat::bitsper;
   j1 = 0;
   for(j = 0; j < nj; j++) {
	  Prod.Vec[j] = 0;
	  for(k1 = 0; k1 < PackMat::bitsper; k1++) {
		 prodwt = 0;
		 for(k = 0; k < numpackrows;k ++) {
			prodwt += PackMat::getWeight(Vec[k] & DM.Matrix[k][j1]);
		 }
		 Prod.Vec[j] = ((prodwt % 2)<<k1) + Prod.Vec[j];
		 j1++;			   
	  }
   }
   // do the last fractbits
   if(Prod.fracbits) {
	  j = nj;
	  Prod.Vec[j] = 0;
	  for(k1 = 0; k1 < Prod.fracbits; k1++) {
		 prodwt = 0;
		 for(k = 0; k < numpackrows;k ++) {
			prodwt += PackMat::getWeight(Vec[k] & DM.Matrix[k][j1]);
		 }
		 Prod.Vec[j] = ((prodwt % 2)<<k1) + Prod.Vec[j]; 
		 j1++;
	  }
   }
}


inline void
PackVec::checkandset(PackVec &B,int brows)
{
   if(!(B.numrows == brows)) {
	  if(B.Vec) { delete[] B.Vec;}
	  setsize(B,brows);
	  callocvec(B);
   }	  
}


inline void
PackVec::checkandset2(PackVec &B,int brows)
// make sure there is at least brows, but possibly more
{
   if(B.numrows < brows) {
	  if(B.Vec) { delete[] B.Vec;}
	  setsize(B,brows);
	  callocvec(B);
   }  
}


void
PackVec::resize(PackVec &A,int rows)
// set A to the right size, but keep all the old stuff that was there before
{
   PKSIZE *NewVec;
   int oldnumpackrows;
   int i;

   oldnumpackrows = A.numpackrows;
   setsize(A,rows);
   NewVec = new PKSIZE[A.numpackrows];
   // clear out the matrix
   for(i = 0; i < oldnumpackrows; i++) {
		 NewVec[i] = A.Vec[i];
   }
   if(A.Vec) delete[] Vec;
   A.Vec = NewVec;
}


void PackVec::clearfringe()
{
   int i,j;
   PKSIZE mask;
   
   mask = (1<<fracbits)-1;
   Vec[numpackrows1] &= mask;
}


int PackVec::getWeight() 		// get the weight of the entire vector
{
   int i,j, wt;

   wt = 0;
   if(fracbits) {				// make sure there is no extra stuff on edges
	  clearfringe();
   }
   for(i = 0; i < numpackrows; i++) {
	  wt += PackMat::getWeight(Vec[i]);
   }
   return wt;

}

// print function
ostream&
PackVec:: printpackvec(ostream &os) const 
{
   int i,j,k,i1;

   for(i = 0; i < numpackrows1; i++) {
	  for(i1 = 0; i1 < PackMat::bitsper; i1++) {
		 if(Vec[i] & (1<<i1)) 
			os << "1";
		 else
			os << "0";
		 os << printspace;
	  }
   }
   for(i = 0; i < numpackrows2; i++) {
	  for(i1 = 0; i1 < fracbits; i1++) {
		 if(Vec[i+numpackrows1] & (1<<i1)) 
			os << "1";
		 else
			os << "0";
		 os << printspace;
	  }
   }
   return os;
}


void PackVec::PrintInfo(void)
{
   
   cout<<"Rows = "<<numrows <<"Packrows: "<< numpackrows<< "  " << 
	  numpackrows1 << "  " << numpackrows2 << endl;
   cout << "fracbits: " << fracbits << endl;
}

void
PackVec::MakeRandom(double p)
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;
   PKSIZE mask;

   if(p == 0 || p == 0.5) {	// if equal prob of 0 or 1
	  for(i = 0; i < numpackrows; i++) {
		 blob = 0;
		 for(k = 0; k < sizeof(PKSIZE); k++) {
			rno = 0xFFFF&(rand()>>8); // pick off two bytes
			blob = (blob << 16) + rno;
		 }
		 Vec[i] = blob;
	  }
	  // now make sure there is nothing left over 
	  if(fracbits) {
		 clearfringe();
	  }
   }
   else {						// else go by the probabilities
	  double rnod;
	  for(i = 0; i < numrows; i++) {
		 rnod = double(rand())/double(RAND_MAX);
		 if(rnod < p)
			SetVal(i,1);
		 else
			SetVal(i,0);
	  }
   }
}


void
PackVec::SetToZero()
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;

   for(i = 0; i < numpackrows; i++) {
	  Vec[i] = 0;
   }
}

void
PackVec::SetToOne()
{
   int i,j,k;
   unsigned int rno;
   PKSIZE blob;

   for(i = 0; i < numpackrows; i++) {
	  Vec[i] = PackMat::masker;
   }
   // now make sure there is nothing left over 
   clearfringe();
}

void PackVec::Add(PackVec &DM, PackVec &Sum)
{
// This function tests to make sure things are the right size
// and then calls the addition function

   if(!(numrows == DM.numrows)) {
	  // correct size to add
	  cerr << "Error: Matrices not same size in sum\n";
	  exit(-1);
   }
   checkandset(Sum,numrows);
   Addnocheck(DM,Sum);
}


void PackVec::Addnocheck(PackVec &DM, PackVec &Sum)
{
// This function does no testing --- it assumes that
// everything is the right size and shape for addition

   int i;

   for(i = 0; i < numpackrows; i++) {
	  Sum.Vec[i] = Vec[i] ^ DM.Vec[i];
   }
}


int PackVec::innerprod2(PackVec &row1, PackVec &row2, int start, int end)
// do the inner product between two row-packed quantities
{
   int stbig, stlit, ebig, elit, i, nl;
   int sum=0;
   PKSIZE mask,mask1;
   int bitsdone=0;

   stbig = start/PackMat::bitsper;
   stlit = start % PackMat::bitsper;
   ebig = end/PackMat::bitsper;
   elit = end % PackMat::bitsper;
   mask1 = ((1<<(PackMat::bitsper-stlit))-1) << stlit;
   mask = mask1;

   for(i = stbig; i < ebig; i++) {
	  sum += PackMat::getWeight(row1.Vec[i] & row2.Vec[i] & mask);
	  mask = PackMat::masker;
	  stlit = 0;
   }
   nl = elit-stlit+1;
   mask = ((1<<nl)-1)<<stlit;
   sum += PackMat::getWeight(row1.Vec[i] & row2.Vec[i] & mask);
   return sum % 2;
}

int PackVec::innerprod(PackVec &row1, PackVec &row2, int end)
// do the inner product between two row-packed quantities, starting
// at the first bit and covering end bits.
{
   int i, ebig, elit;
   int sum= 0;

   ebig = end/PackMat::bitsper; elit = end % PackMat::bitsper;
   for(i = 0; i < ebig; i++) {
	  sum += PackMat::getWeight(row1.Vec[i] & row2.Vec[i]);
   }
   if(elit) {
	  sum += PackMat::getWeight((row1.Vec[i]&row2.Vec[i]) & ((1<<elit)-1));
   }
   return sum % 2;
}


void PackVec::Subvector(PackVec &subvec,int row1, int row2,int subrow) 
// extract the row or vector and place it (starting at subrow)
// in the vector subvec
{
   int i,i1, nrows;

   if(row1 < 0) row1 = 0;
   if(row2 >= numrows) row2 = numrows-1;

   nrows = (row2-row1+1) + subrow;
   subvec.checkandset2(subvec,nrows);
   for(i=0, i1=subrow; i < (row2-row1+1); i++,i1++) {
	  subvec.SetVal(i1,GetVal(i+row1));
   }
}


void
stackside(PackMat & C,PackMat &A, PackMat & B, StackMode ssm,
		  PackMode pkmode)
// C = [A B]
// Stack side by size.  If not the same height, then use
// the stacksidemode value (exact, top, bottom, center1st, center2nd, truncate)
{
   int numzaboveA, numzaboveB, numzbelowA, numzbelowB, numArows, numBrows,
	  numCrows, numCcols;
   int i,j;
   int numcolA,numcolB, numrowA,numrowB;

   StackMode ss;

   if(ssm == Default) {
	  ss = PackMat::sidestackmode;
   }
   else {
	  ss = ssm;
   }
	  
   if(A.numrows && B.numrows && A.numrows != B.numrows && 
	  ss==Exact) { 
	  cerr << "Error in side stacking: not the same height\n";
	  exit(-1);
   }
   numcolA = A.numcols;
   numcolB = B.numcols;
   numrowA = A.numrows;
   numrowB = B.numrows;
   numCcols = numcolA + numcolB;


   numzaboveA = 0;				// default values
   numzaboveB = 0;
   numzbelowA = 0;
   numzbelowB = 0;
   numArows = numrowA;
   numBrows = numrowB;
   numCrows = MAX(numrowA,numrowB);
   
   switch(ss) {
   case Exact:
	  break;
   case Top:
	  numzaboveA = 0;
	  numzaboveB = 0;
	  numzbelowA = MAX(numrowB - numrowA,0);
	  numzbelowB = MAX(numrowA - numrowB,0);
	  break;
   case Truncate:
	  numArows = MIN(numrowA,numrowB);
	  numBrows = numArows;
	  numCrows = numArows;
	  break;
   case Bottom:
	  numzaboveA = MAX(numrowB - numrowA,0);
	  numzaboveB = MAX(numrowA - numrowB,0);
	  break;
   case Center:
	  if(A.numrows > B.numrows) {
		 numzaboveB = (numrowA - numrowB)/2;
		 numzbelowB = (numrowA - numrowB - numzaboveB);
	  }
	  else {
		 numzaboveA = (numrowB - numrowB)/2;
		 numzbelowA = (numrowB - numrowA - numzaboveA);
	  }
	  break;
   default:
	  cerr << "Error: invalid side stack mode\n";
	  exit(-1);
	  break;
   }
   C.checkandresize2(C,pkmode, numCrows, numCcols);
// cout << "C(resize)=\n" << C << endl;
//  cout << "numCcol=" << numCcols << endl;
//  cout << "numcolA=" << numcolA << endl;

   // The A column
   for(i = 0; i < numzaboveA; i++) {
	  for(j = 0; j < numcolA; j++) {
		 C.SetVal(i,j,0);
	  }
   }
   for(i = 0; i < numArows; i++) {
	  for(j = 0; j < numcolA; j++) {
		 C.SetVal(i+numzaboveA,j,A.GetVal(i,j));
	  }
   }
   for(i = 0; i < numzbelowA; i++) {
	  for(j = 0; j < numcolA; j++) {
		 C.SetVal(i+numzaboveA+numArows,j,0);
	  }
   }

   // The B column
   for(i = 0; i < numzaboveB; i++) {
	  for(j = 0; j < numcolB; j++) {
		 C.SetVal(i,j+numcolA,0);
	  }
   }
   for(i = 0; i < numBrows; i++) {
	  for(j = 0; j < numcolB; j++) {
		 C.SetVal(i+numzaboveB,j+numcolA,B.GetVal(i,j));
	  }
   }
   for(i = 0; i < numzbelowB; i++) {
	  for(j = 0; j < numcolB; j++) {
		 C.SetVal(i+numzaboveB+numBrows,j+numcolA,0);
	  }
   }
}


void
stacktop(PackMat & C,PackMat &A, PackMat & B, StackMode tsm, PackMode pkmode)
// C = [A;B]
// Stack A over B.  If not the same width, then use
// the stacktopmode value (exact, left, right, center, truncate)
{
   int numzleftA, numzleftB, numzrightA, numzrightB, numAcols, numBcols,
	  numCrows, numCcols;
   int i,j;
   int numcolA,numcolB, numrowA,numrowB;

   StackMode ts;

   if(tsm == Default) {
	  ts = PackMat::topstackmode;
   }
   else {
	  ts = tsm;
   }
	  
   if(A.numcols && B.numcols && A.numcols != B.numcols && 
	  ts==Exact) { 
	  cerr << "Error in top stacking: not the same width\n";
	  exit(-1);
   }
   numcolA = A.numcols;
   numcolB = B.numcols;
   numrowA = A.numrows;
   numrowB = B.numrows;
   numCrows = numrowA + numrowB;


   numzleftA = 0;				// default values
   numzleftB = 0;
   numzrightA = 0;
   numzrightB = 0;
   numAcols = numcolA;
   numBcols = numcolB;
   numCcols = MAX(numcolA,numcolB);
   
   switch(ts) {
   case Exact:
	  break;
   case Left:
	  numzleftA = 0;
	  numzleftB = 0;
	  numzrightA = MAX(numcolB - numcolA,0);
	  numzrightB = MAX(numcolA - numcolB,0);
	  break;
   case Truncate:
	  numAcols = MIN(numcolA,numcolB);
	  numBcols = numAcols;
	  numCcols = numAcols;
	  break;
   case Right:
	  numzleftA = MAX(numcolB - numcolA,0);
	  numzleftB = MAX(numcolA - numcolB,0);
	  break;
   case Center:
	  if(numcolA > numcolB) {
		 numzleftB = (numcolA - numcolB)/2;
		 numzrightB = (numcolA - numcolB - numzleftB);
	  }
	  else {
		 numzleftA = (numcolB - numcolB)/2;
		 numzrightA = (numcolB - numcolA - numzleftA);
	  }
	  break;
   default:
	  cerr << "Invalid top stack mode\n";
	  exit(-1);
	  break;
   }
   C.checkandresize2(C,pkmode, numCrows, numCcols);
//   cout << "C(resize)=\n" << C << endl;

   // The A row
   for(i = 0; i < numrowA; i++) {
	  for(j = 0; j < numzleftA; j++) {
		 C.SetVal(i,j,0);
	  }
   }
   for(i = 0; i < numrowA; i++) {
	  for(j = 0; j < numAcols; j++) {
		 C.SetVal(i,j+numzleftA,A.GetVal(i,j));
	  }
   }
   for(i = 0; i < numrowA; i++) {
	  for(j = 0; j < numzrightA; j++) {
		 C.SetVal(i,j+numzleftA+numAcols,0);
	  }
   }

   // The B row
   for(i = 0; i < numrowB; i++) {
	  for(j = 0; j < numzleftB; j++) {
		 C.SetVal(i+numrowA,j,0);
	  }
   }
   for(i = 0; i < numrowB; i++) {
	  for(j = 0; j < numBcols; j++) {
		 C.SetVal(i+numrowA,j+numzleftB,B.GetVal(i,j));
	  }
   }
   for(i = 0; i < numrowB; i++) {
	  for(j = 0; j < numzrightB; j++) {
		 C.SetVal(i+numrowA,j+numBcols+numzleftB,0);
	  }
   }
}

void setval(PKSIZE **Matrix,int row, int col, int val, 
					 PackMode packmode)
{
   int i,j,k, col1, row1;
   int bitcol, bitrow;
   PKSIZE value;
   
   if(packmode == HORIZPACK) {
	  col1 = col >> PackMat::shifter;
	  bitcol = col % PackMat::bitsper;
	  if(!val) {					// set value to 0
		 Matrix[row][col1] &= ~(1<<bitcol);
	  }
	  else {						// set value to 1
		 Matrix[row][col1] |= (1<<bitcol);
	  }
   }
   
   else if(packmode == VERTPACK) {
	  row1 = row >> PackMat::shifter;
	  bitrow = row % PackMat::bitsper;
	  if(!val) {				// set value to 0
		 Matrix[row1][col] &= ~(1<<bitrow);
	  }
	  else {
		 Matrix[row1][col] |= (1<<bitrow);
	  }
   }
}

/*
Local Variables:
compile-command: "g++ -c -g PackMat.cc"
End:
*/
