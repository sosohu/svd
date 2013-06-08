/*
   The openmp program in c of the svd for a general double M-by-N matrix A
   The code is provided by Mr. Hu at the time of 2013.3 . 
   His email address is ustc.sosohu@gmail.com
 */

/*
    Purpose   
    ====================================================================   
    GEBRD reduces a general double M-by-N matrix A to upper or lower        
    bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.   
    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.     
    Arguments   
    =========   
    M       (input) INTEGER  假设矩阵按行存储（这个可以优化）  
            The number of rows in the matrix A.  M >= 0.   
    N       (input) INTEGER   
            The number of columns in the matrix A.  N >= 0.   
    A       (input/output) REAL array 
            On entry, the M-by-N general matrix
    D       (output) REAL array, N   
    Q       (output) REAL array M*M
    P       (output) REAL array N*N 
    INFO    (output) INTEGER   
            = 0:  successful exit   
   =========================================================================
 */


#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<math.h>
#include<omp.h>
#include<assert.h>

#define ITERATION 20
#define NUM_THREADS 4

extern double is_zero;
extern double tol;
//extern inline int mysgn(double s);
//extern inline double myabs(double s);
inline int mysgn(double s)
{
  if((s<is_zero)&&(s>-is_zero)) return 0;
  if(s>0.0)  return 1;
  return -1;
}

inline double myabs(double s)
{
  if(s<0.0) return -s;
  return s;
}


extern void TransposeQuick(int M, int N,double *A);

extern int Partion(int s,int t,double* A,int* index);

extern void QuickSort(int s,int t,double* A,int* index);

extern void  ChangeOrder(int M,int N,double* A,double* Q,double* P,double* D,int MlsN);

extern void orth_all(int start,int end,int M,int N,double* A,double* P,double* Q,int MlsN,
			  int* cond,double* tol_orth,int iteration,double last_ave_orth);

extern void orth_two(int st_1,int en_1,int st_2,int en_2,int M,int N,double* A,double* P, double* Q, int MlsN,
			  int* cond,double* tol_orth,int iteration,double last_ave_orth);

int myorder_gesvd_op_tran(int M,int N,double tol,double *A,double *D,
		  double *Q,double *P,int *INFO)
{
  int i,j,k,w,q;
  double a,b,c;
  int cond = 0;
  int iteration = 0;
  int id;
  int st_1,en_1,st_2,en_2;
  double tol_orth = 0;
  double last_ave_orth = 0;
  int MlsN = 0;
  
  *INFO = 0;
  if(M < 0)  *INFO = -1;
  else if(N < 0) *INFO = -2;

  if(*INFO < 0)
  {
    printf("error:%d\n",*INFO);
	exit(0);
  }
  
  omp_set_num_threads(NUM_THREADS);

  
  if(M<N) //转置矩阵A
  {
    MlsN = 1;
    TransposeQuick(M,N,A);
    N = N + M;
	M = N - M;
	N = N - M;
  }// 转置后 M N交换值  M > N
  
  assert(M%8==0);

  #pragma omp parallel default(shared) private(id,i)
  {
  id = omp_get_thread_num();
  if(MlsN==0)
   {
    //初始化Q为单位阵
    for(i=id*M/4; i<(id+1)*M/4; i++)
    {
      *(Q+i*M+i) = 1.0;
    }
   }
  else
   {
    //初始化P为单位阵
    for(i=id*M/4; i<(id+1)*M/4; i++)
    {
      *(P+i*M+i) = 1.0;
    }
   }
  }// end for pragma

  while(1)
  {
   cond = 0;
   last_ave_orth = 2*tol_orth/(M*(M-1));
   tol_orth = 0;
   #pragma omp parallel reduction(||:cond) reduction(+:tol_orth) default(shared) private(id)
   {
    id = omp_get_thread_num();
    // void orth_all(int start,int end,int M,int N,double* A,double* P,double* Q,int MlsN
	//		int* cond,double* tol_orth,int iteration,double last_ave_orth) 
    orth_all(id*M/4,(id+1)*M/4-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);

    #pragma omp barrier

	// (p1,p2) (p3,p4)
    #pragma omp sections
	{
    // void orth_two(int st_1,int en_1,int st_2,int en_2,int M,int N,double* A,double* P, double* Q, int MlsN,
	//	    int* cond,double* tol_orth,int iteration,double last_ave_orth)
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,2*M/8,3*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,3*M/8,4*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(4*M/8,5*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(5*M/8,6*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,3*M/8,4*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,2*M/8,3*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(4*M/8,5*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(5*M/8,6*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}

	// (p1,p3) (p2,p4)
    #pragma omp sections
	{
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(2*M/8,3*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(3*M/8,4*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(2*M/8,3*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(3*M/8,4*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}

	// (p1,p4) (p2,p3)
    #pragma omp sections
	{
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(2*M/8,3*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(3*M/8,4*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two(0*M/8,1*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(1*M/8,2*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(2*M/8,3*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	  #pragma omp section 
      orth_two(3*M/8,4*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth);
	}
   }// end for pragma
  
   if(cond == 0) break;

   if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}
  }// end for while
 
  ChangeOrder(M,N,A,Q,P,D,MlsN);  // 按行向量模的大小进行排序 

   
  return 0;
}
