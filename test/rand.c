#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int main(int argc, char **argv)
{
   int M, N;
   int i,j;
   M = atoi(argv[1]);
   N = atoi(argv[2]);

   srand((unsigned int)time(0));
   for(i=0; i<M; i++)
   {
	 for(j=0; j<N; j++)
	   printf("%lf ",rand()/((double)(RAND_MAX)));
	 printf("\n");
   }
   return 0;
}
