#include"PackMat.h"

using namespace std;

int main()
{
   int i,j,k;

   PackMat::initPackMat();

   PackMat b23(string("2 3 1 1 0 0 1 1"));

   cout << "b23" << b23 << endl;
   b23.PrintInfo();

   ifstream Densefile("A1.txt");
   if(!Densefile) {
	  cerr << "Cannote open file\n";
	  exit(-1);
   }
   PackMat A(Densefile);
   cout<<"A"<<endl;
   A.PrintInfo();
   cout << A << endl;
   
   PackMat B("B1.txt",VERTPACK);
   cout<<"B"<<endl;
   B.PrintInfo();
   cout << B << endl;


   PackMat Prod1;

   A.Mult(B,Prod1);

   cout<<"Prod1"<<endl;
   cout << Prod1 << endl;

   cout << "A (getval)" << endl;
   for(i = 0; i < A.rows(); i++) {
	  for(j = 0; j < A.cols(); j++) {
		 cout << A.GetVal(i,j) << " ";
	  }
	  cout << endl;
   }

   cout << "B (getval)" << endl;
   for(i = 0; i < B.rows(); i++) {
	  for(j = 0; j < B.cols(); j++) {
		 cout << B.GetVal(i,j) << " ";
	  }
	  cout << endl;
   }

   PackMat A2;
   A2 = B;
   for(i = 0; i < A2.rows(); i++) {
	  for(j = 0; j < A2.cols(); j++){
		 A2.SetVal(i,j,(j+1)%2);
	  }
   }
   cout << "A2 (setval)" << endl;
   cout << A2;

   A.SetMode(VERTPACK,A2);
   cout << "A2" << endl;
   A2.PrintInfo();
   cout << A2 << endl;
   
   srand(time(0));
   PackMat C;
   C.Size(10,10);
   C.MakeRandom(0.9);
   cout << "C" << endl;
   cout << C;
   C.SetToZero();
   cout << C;
   C.SetToOne();
   cout << C;
   cout << "Weight=" << C.getWeight() << endl;

   PackMat T;
   A.transpose(T,VERTPACK);
   cout << "T";
   T.PrintInfo();
   cout << T << endl;

   A2.PrintInfo();
   A2.transpose(T,HORIZPACK);
   cout << "A2'";
   cout << T << endl;

   PackMat D(A);
   cout << "D\n" << D << endl;
   PackMat S;
   A.Add(D,S);
   cout << "S\n" << S << endl;

   PackMat E = A;
   cout << "E\n" << E << endl;

   cout << "Square A:\n" << A << endl;
   int ni;
   A.ludcmp(ni);
//   int testbits[] = {1,0,1,0,0,0,0,1,1,1};   // 0
//   int testbits[] = {0,1,1,1,0,0,0,0,1,1};   // 1
//   int testbits[] = {1,0,1,1,1,0,0,0,0,0};   // 2
//   int testbits[] = {0,1,0,1,1,1,0,0,1,1};   // 3
//   int testbits[] = {1,0,1,0,1,1,1,0,0,0};   // 4
//   int testbits[] = {0,1,1,1,0,1,0,1,0,0};   // 5
//   int testbits[] = {1,0,1,1,0,0,0,0,1,1};   // 6
//   int testbits[] = {0,1,0,1,1,0,0,0,0,0};   // 7
//   int testbits[] = {1,0,1,0,1,0,0,0,1,1};   // 8
   int testbits[] = {0,1,1,1,0,1,1,1,0,1};   // 9
   int n = sizeof(testbits)/sizeof(int);
   cout << "nbits in n=" << n;
   int bitsper = PackMat::bitsper;
   int nbig = n/bitsper;
   int nlit = n % bitsper;
   PKSIZE *b = new PKSIZE[nbig+1];
   PKSIZE *y = new PKSIZE[nbig+1];
   for(i = 0; i < nbig+1; i++) {b[i] = 0;}


   int j1 = 0;
   for(i = 0; i < nbig; i++) {
	  for(j = 0; j < bitsper; j++) {
		 if(testbits[j1]) b[i] |= (1<<j);
		 j1++;
	  }
   }
   if(nlit) {
	  for(j= 0; j < nlit; j++) {
		 if(testbits[j1])
			b[i] |= (1<<j);
		 j1++;
	  }
   }
   cout <<" b=" << int(b[0]) << " " << int(b[1]) << endl;
   A.lubacksub(b,y);
   cout << "y=" << int(y[0]) << " " << int(y[1]) << endl;

//   PackVec bcol(string("10   1 0 1 0 0 0 0 1 1 1"));  // 0 
   PackVec bcol(string("10   0 1 1 1 0 0 0 0 1 1"));  // 1
//    PackVec bcol(string("10   1 0 1 1 1 0 0 0 0 0"));  // 2
//    PackVec bcol(string("10   0 1 0 1 1 1 0 0 1 1"));  // 3
//    PackVec bcol(string("10   1 0 1 0 1 1 1 0 0 0"));  // 4
//    PackVec bcol(string("10   0 1 1 1 0 1 0 1 0 0"));  // 5
//    PackVec bcol(string("10   1 0 1 1 0 0 0 0 1 1"));  // 6
//    PackVec bcol(string("10   0 1 0 1 1 0 0 0 0 0"));  // 7
//    PackVec bcol(string("10   1 0 1 0 1 0 0 0 1 1"));  // 8
//    PackVec bcol(string("10   0 1 1 1 0 1 1 1 0 1"));  // 9
   PackVec ycol;
   cout << "bcol=" << bcol << endl;
   A.lubacksub(bcol,ycol);
   cout << "ycol=" << ycol << endl;
				
   PackMat Ainv;
   A.matinv(Ainv);
   cout << "Ainv=" << Ainv << endl;
   PackMat Aprod;
   A.Mult(Ainv,Aprod);
   cout << "Aprod=" << Aprod << endl;

   PackMat G("A2.txt");
   cout << "G=\n" << G << endl;
   G.reducetosystematic(G,ni);

   PackMat A3("A3.txt");
   cout << "A3=\n" << A3 << endl;
   A3.ludcmp(ni);
   cout <<"ni=" << ni << endl;

   PackMat A4("A4.txt");
   cout << "A4=\n" << A4 << endl;
   A4.ludcmp(ni);
   A4.reducetosystematic(A4,ni);

   PackMat A5("A5.txt");
   cout << "A5= (tall, lin dep)\n" << A5 << endl;
   A5.ludcmp(ni);
   cout << "numind=" << ni << endl;

   A.SetMode(VERTPACK,A2);

   cout << "b=" << bcol << "  A=\n" << A2 << endl;
   PackVec bA;
   bcol.Mult(A2,bA);
   cout << "bA=" << bA << endl;

   cout << "b=" << bcol << "  A=\n" << A << endl;
   PackVec Ab;
   A.Mult(bcol,Ab);
   cout << "Ab=" << Ab << endl;
   for(i = 0; i < Ab.rows(); i++) {
	  cout << Ab.GetVal(i) << " ";
   }
   cout << endl;

   for(i = 0; i < Ab.rows(); i++) {
	  Ab.SetVal(i,(i+1) % 2);
   }
   cout << "Abnew=" << Ab << endl;

   cout << "Random:\n";
   for(i = 0; i < 10; i++) {
	  Ab.MakeRandom();
	  cout << Ab << "  weight=" << Ab.getWeight() << endl;
   }
   Ab.SetToZero();
   cout << "AB(0)=" << Ab << endl;
   Ab.SetToOne();
   cout << "AB(1)=" << Ab << "  weight=" << Ab.getWeight() << endl;
   cout << "b=" << bcol << endl;
   PackVec Svec;
   Ab.Add(bcol,Svec);
   cout << "Ab + b=" << Svec << endl;
   PackMat Asub(7,7);
   A.Submatrix(Asub,0,2,0,3,1,4);
   cout << "Asub\n" << Asub << endl;
   
   PackVec asub;
   A.Subvector(asub,2,2,1,4,3);
   cout << "asub=" << asub << endl;

   asub.SetToZero();
   A.Subvector(asub,2,4,3,3,3);
   cout << "asub=" << asub << endl;
   PackVec bsub(10);
   asub.Subvector(bsub,2,4,5);
   cout << "bsub=" << bsub << endl;

   PackMat Bsub1(bsub);
   cout << "Bsub1\n" << Bsub1 << endl;
   PackMat Bsub2(bsub,VERTPACK);
   cout << "Bsub2\n" << Bsub2 << endl;

   PackVec bsub1(Bsub1);
   cout << "bsub1=" << bsub1 << endl;
   PackVec bsub2(Bsub2);
   cout << "bsub2=" << bsub2 << endl;

   PackMat Bsub;
//   Asub.Size(3,4);
   A.Submatrix(Asub,0,2,0,3);
   A.Submatrix(Bsub,0,2,0,3);
   PackMat Cmat;
   A.setSideStackMode(Bottom);
   stackside(Cmat,Asub,Bsub,Truncate);
   cout << "Asub=\n" << Asub << endl;
   cout << "Bsub=\n" << Bsub << endl;
   cout << "Cmat=\n" << Cmat << endl;

   PackMat tmp(bsub1,VERTPACK);
   stackside(Cmat,Cmat,tmp,Truncate);
   cout << "Cmat(bsub1)=\n" << Cmat << endl;

   PackMat D1(string("3 2   0 1   1 0   1 1"));
   cout << "D1=\n" << D1 << endl;
   PackMat D2;
   for(i = 0; i < 3; i++) {
	  stackside(D2,D2,D1);
	  cout << "i=" << i << "  D2=\n" << D2 << endl;
   }

   Cmat.Size(0,0);
   stacktop(Cmat,Asub,Bsub,Truncate);
   cout << "Asub=\n" << Asub << endl;
   cout << "Bsub=\n" << Bsub << endl;
   cout << "Cmat=\n" << Cmat << endl;

   PackMat D3(string("2 3  0 1 0  1 1 0"));
   PackMat D4;
   for(i = 0; i < 3; i++) {
	  cout << "i=" << i << endl;
	  stacktop(D4,D3,D4);
	  cout << "i=" << i << "  D4=\n" << D4 << endl;
   }


   cout <<"Building big randoms\n";
   int Nbig = 100;
   PackMat Big1(Nbig,Nbig);
   PackMat Big2(Nbig,Nbig,VERTPACK);
   PackMat Big3(Nbig,Nbig);
   Big1.MakeRandom();
   Big2.MakeRandom();
   cout << "Calling Product\n";
   Big1.Mult(Big2,Big3);
   cout << "Product done\n";

   Nbig = 10;
   PackMat IT(Nbig,Nbig);
   IT.MakeRandom();
   PackMat ITinv;
   PackMat ITProd;
   IT.matinv(ITinv);
   IT.Mult(ITinv,ITProd);
   cout << "Isinverse? " << ITProd.IsIdentity() << endl;
   
}

/*
Local Variables:
compile-command: "g++ -o testPackMat -g testPackMat.cc PackMat.cc"
End:
*/
