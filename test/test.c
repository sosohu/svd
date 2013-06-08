#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<sys/time.h>
#include<math.h>
#include<time.h>
#include<assert.h>

#define NUM_THREADS 4

/*
double serial(int M,int N,double* A)
{
  int i,j,k;
  double tmp;
  int count;
  count = 0;
  tmp = 0;
  while(count++<20*N)
  {
  for(i=0; i<M; i++)
  {
    for(k=0; k<N; k++)
	{
	  tmp = tmp + (*(A+i*N+k)*(*(A+i*N+k)));
	}
  }
  }
  return tmp;
}
*/

double parallel(int M,int N,double* A)
{
  int i,j,k,count,a;
  double tmp,q;
  int id;
  assert(M%4==0);
  tmp = 0;
  omp_set_num_threads(NUM_THREADS);
  count = 0;
  q = 1;
  while(count++ < 20)
  {
  a = 0;
  #pragma omp parallel default(shared) reduction(+:tmp) private(i,j,k,id,q)
  {	  
  id = omp_get_thread_num();
  while(a<N)
  {
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
    for(k=0; k<N; k++)
	{
	  tmp  = tmp + (*(A+i*N+k))*(*(A+i*N+k));
	}
    for(k=0; k<N; k++)
	{
	q  = (1+q)*(*(A+i*N+k));
	}	
  }
  #pragma omp master
  {
    a++;
  }
  #pragma omp barrier
  }//end for while
  }//end for pragma
  }
  return tmp;
}


double serial(int M,int N,double* A)
{
  int i,j,k,count,a;
  double tmp,q;
  int id;
  assert(M%4==0);
  tmp = 0;
  omp_set_num_threads(NUM_THREADS);
  count = 0;
  q = 1;
  while(count++<20)
  {
  a = 0;
  #pragma omp parallel default(shared) reduction(+:tmp) private(i,j,k,id,q)
  {
while(a<4)
{
  j = 0;
  id = omp_get_thread_num();
  while(j++<N/4)
  {
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
    for(k=0; k<N; k++)
	{
	  tmp  = tmp + (*(A+i*N+k))*(*(A+i*N+k));
	}
    for(k=0; k<N; k++)
	{
	  q  =  (1+q)*(*(A+i*N+k));
	}	
  }
  }
  #pragma omp master
  {
   a++;
  }
  #pragma omp barrier
}
  }
  }
  return tmp;
}


double parallel_asm(int M,int N,double* A)
{
  int i,j,k,k1,count;
  double tmp;
  int id;
  unsigned long mask = 16;
  tmp = 0;
  assert(M%8==0);
  omp_set_num_threads(NUM_THREADS);
  count = 0;
  while(count++<20*N)
  {
  #pragma omp parallel default(shared) reduction(+:tmp) private(i,j,k,k1,id)
  {
  id = omp_get_thread_num();
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
	if(((unsigned long)(A+i*N) % mask) == 0)
       k = 0;
	else
	{
	  k = 1;
	  tmp = tmp + (*(A+i*N))*(*(A+i*N));
	  tmp = tmp + (*(A+i*N+N-1))*(*(A+i*N+N-1));
	}
    asm 
	(
	 "mov.d %0,%4\n\t"
	 "L1:gsLQC1 $f2,$f0,0(%1)\n\t"
	 "add %2,%2,2\n\t"
	 "madd.d %0,%0,$f2,$f2\n\t"
	 "madd.d %0,%0,$f0,$f0\n\t"
	 "add %1,%1,16\n\t"
	 "blt %2,%3,L1\n\t"
     :"=f"(tmp)
	 :"r"(A+i*N+k),"r"(k),"r"(N-1),"f"(tmp)
	 :"$f0","$f2"
	);
	//
  }
  }
  }
  return tmp;
}

void small_test(double* A)
{
  double i,j;
  int k;
  i = 2.2;
  j = 2.2;
  if((unsigned long)(A)%16 == 0)
	k = 0;
  else
	k = 1;
  asm volatile
  (
   "gsSQC1 %0,%1,0(%2)\n\t"   
   :
   :"f"(i),"f"(j),"r"(A+k)
  );
}

void test()
{
  double* A;
  int i,k;
  unsigned long mask = 16;
  double fq,ft;
  A = (double*)malloc(10*sizeof(double));
  for(i=0; i<10; i++)
	*(A+i) = i;
  small_test(A);
  for(i=0; i<10; i++)
    printf("%f ",*(A+i));
  printf("\n");
}

void _move_data(int ac_p,int in_p,double *A,int M,int N)
{
  int p;
  double tmp;
  tmp = *(A+ac_p);
  if(ac_p >= in_p)
   for(p=ac_p; p>in_p; p--)
   {
    *(A+p) = *(A+p-1);
   }
  else
   for(p=ac_p; p<in_p; p++)
   {
    *(A+p) = *(A+p+1);
   }
  *(A+in_p) = tmp;
}

void TransposeInplaceQuick(int M, int N,double *A)
{
    int m =  N < M ? N : M;
    int i, j;
	double tmp;
	int p,ac_p,in_p;
    for (i = 0; i < m; ++i) 
	{
        for (j = i+1; j < m; ++j) 
		{
            tmp = *(A+i*N+j);
			*(A+i*N+j) = *(A+j*N+i);
			*(A+j*N+i) = tmp;
        }
    }

    if (m == N ) 
	{
        for (i = m; i < M; ++i) 
		{
            for (j = 0; j < N; ++j) 
			{
			  ac_p = i*N+j;
			  in_p = j*(N+i-m+1)+i;
              tmp = *(A+ac_p);
              for(p=ac_p; p>in_p; p--)
                 *(A+p) = *(A+p-1);
              *(A+in_p) = tmp;
            }
        }
    }
    else 
	{
        int last = (M-1)*N + N-1;
        int n = 0; 
        for (j = N - 1; j >= m; --j) 
		{
            int step = n;
            for (i = M - 1; i >= 0; --i, --last) 
			{
				ac_p = i*N+j-step;
				in_p = last;
                tmp = *(A+ac_p);
                for(p=ac_p; p<in_p; p++)
                  *(A+p) = *(A+p+1);
                *(A+in_p) = tmp;
                step -= (N - 1 - j);
            }
            n += M - 1;
        }
    }
    N = N + M;
	M = N - M;
	N = N - M;
}

void malloc_time(int M,int N)
{
  double* A = (double*)malloc(M*N*sizeof(double));
  int a = 10;
  int i,j;
  for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	  *(A+i*N+j) = i;
  free(A);
}

int main(int argc, char** argv)
{
  int M = atoi(argv[1]);
  int N = atoi(argv[2]);
  int i,j;
  double* A;
  double* B;
  double* C;
  struct timeval tpstart,tpend;
  unsigned long timeuse;
  double tmp;
  FILE *fp;

  A = (double*)malloc(M*N*sizeof(double));  
  B = (double*)malloc(M*N*sizeof(double));  
  C = (double*)malloc(M*N*sizeof(double));  
 
  fp = fopen("rand_out","r");
  if(fp == NULL)
  {
    printf("can not read the file rand_out \n");
	return 0;
  }
  for(i=0; i<M; i++)
    for(j=0; j<N; j++)
	{
	  fscanf(fp,"%lf ",A+i*N+j);
	  *(B+i*N+j) = *(A+i*N+j);
	  *(C+i*N+j) = *(A+i*N+j);
	}
  
  /*
  printf("A:\n");
  for(i=0; i<M; i++)
  {
    printf("\n");
	for(j=0; j<N; j++)
	  printf("%lf ",*(A+i*N+j));
  }
  */
  /*
  gettimeofday(&tpstart,NULL);
  //TransposeInplaceQuick(M,N,A);
  malloc_time(M,N);
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  printf("M:%d, N:%d, time use:%20.2fus\n",M,N,timeuse);
  */
  /*
  printf("\nAT:\n");
  for(i=0; i<N; i++)
  {
    printf("\n");
	for(j=0; j<M; j++)
	  printf("%lf ",*(A+i*M+j));
  }
  */

  /*
  for(i=0; i<M; i++)
  {
	printf("\n");
	for(j=0; j<N; j++)
	  printf(" %f",*(C+i*N+j));
  }*/
  /*
  gettimeofday(&tpstart,NULL);
  serial(M,N,A);
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  printf("serial time use:%20.2fus\n",timeuse);
  */
  test();
  printf("M: %d N: %d\n",M,N);
  gettimeofday(&tpstart,NULL);
  tmp=serial(M,N,A);
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  printf("serial time use:%20ld us %20.10f\n",timeuse,tmp);
  gettimeofday(&tpstart,NULL);
  tmp=parallel(M,N,B);
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  printf("paralle time use:%20ld us %20.10f\n",timeuse,tmp);
  gettimeofday(&tpstart,NULL);
  tmp=parallel_asm(M,N,C);
  gettimeofday(&tpend,NULL);
  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
  printf("paralle_asm time use:%20ld us %20.10f\n",timeuse,tmp);
  //test();
  /*
  for(i=0; i<M; i++)
  {
	printf("\n");
	for(j=0; j<N; j++)
	  printf(" %20.10f",*(C+i*N+j));
  }
  */
  free(A);
  free(B);
  free(C);

  return 0;
}
