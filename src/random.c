#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<f2c.h>

int main(int argc, char** argv)
{
  int i,j; 
  int M,N;
  double max;
  M = atoi(argv[1]);
  N = atoi(argv[2]);
  max = atof(argv[3]);

  srand((unsigned int)time(0));
  for(i=0; i<M; i++)
  {
    for(j=0; j<N; j++)
	  printf("%f ", rand()/((double)(RAND_MAX)/max));
	printf("\n");
  }
  return 0;
}

