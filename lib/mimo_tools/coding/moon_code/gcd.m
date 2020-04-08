int gcd(int a, int b)
/* Compute g = (a,b) */
{  int g;
   while(b) {
      g = b;
      b = a % b;  /* compute remainder of a/b */
      a = g;
   }
   return a;
}
