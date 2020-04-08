/**************************************
*
*  Program: makerand --- make a file of random bytes
*
*
*  Todd K. Moon
*  Date: Dec. 19, 2002
*
***************************************/

#include <stdio.h>
#include <stdlib.h>

main(int argc, char *argv[])
{
   int nbytes;
   char *fname;
   FILE *fout;
   int i;
   unsigned char b;

   if(argc==1) {
	  printf("Usage: %s size fname\n",argv[0]);
	  exit(-1);
   }
   nbytes = atoi(argv[1]);
   fname = argv[2];
   fout = fopen(fname,"wb");
   for(i = 0; i < nbytes; i++) {
	  b = 256*((double)rand()/(double)RAND_MAX);
	  fwrite(&b,1,1,fout);
   }
   fclose(fout);
}
   

/*
Local Variables:
compile-command:"gcc -o makerand makerand.c"
End:
*/
