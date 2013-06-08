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


void TransposeQuick(int M, int N,double *A)
{
  double* ATMP = (double*)malloc(M*N*sizeof(double));
  int i,j,id;
  assert(M%4==0);
  assert(N%4==0);
  #pragma omp parallel  default(shared) private(i,j,id)
  {
  id = omp_get_thread_num();
  for(i=id*M/4; i<(id+1)*M/4; i++)
	for(j=0; j<N; j++)
	  *(ATMP+i*N+j) = *(A+i*N+j);
  for(i=id*N/4; i<(id+1)*N/4; i++)
	for(j=0; j<M; j++)
	  *(A+i*M+j) = *(ATMP+j*N+i);
  }// end for pragma
  free(ATMP);
}

int Partion(int s,int t,double* A,int* index)
{
  double last = *(A+t);
  double tmp;
  int i,j,in;
  for(i=s,j=s-1; i<t; i++)
  {
	if((*(A+i)) >= last)  
	{
	 if(i>j+1)
	 {
	  tmp = *(A+i);
	  *(A+i) = *(A+j+1);
	  *(A+j+1) = tmp;
	  in = *(index+i);
	  *(index+i) = *(index+j+1);
	  *(index+j+1) = in;
	 }
	  j++;
	}      
  }
  tmp = last;
  *(A+t) = *(A+j+1);
  *(A+j+1) = tmp;
  in = *(index+t);
  *(index+t) = *(index+j+1);
  *(index+j+1) = in;
  return j+1;
}

void QuickSort(int s,int t,double* A,int* index)
{
  int q;
  if(s>=t)
	return;
  q = Partion(s,t,A,index);
  QuickSort(s,q-1,A,index);
  QuickSort(q+1,t,A,index);
}

void  ChangeOrder(int M,int N,double* A,double* Q,double* P,double* D,int MlsN)  // 按行向量模的大小进行排序 
{
  double* DTMP;
  double*ATMP;
  double* QP;
  int i,j,k,id;
  int* index;
  DTMP = (double*)malloc(M*sizeof(double));
  index = (int*)malloc(M*sizeof(int));
  ATMP = (double*)malloc(M*N*sizeof(double));
  QP = (double*)malloc(M*M*sizeof(double));
  #pragma omp parallel default(shared) private(id,i,j)
  {
   id = omp_get_thread_num();
   for(i=id*M/4; i<(id+1)*M/4; i++)
   {
	 *(DTMP+i) = 0;  
	 *(index+i) = i;
	 for(j=0; j<N; j++)
	 {
	   *(DTMP+i) = *(DTMP+i) + (*(A+i*N+j))*(*(A+i*N+j));
	   *(ATMP+i*N+j) = *(A+i*N+j);
	 }
	 if(MlsN==0)
	   for(j=0; j<M; j++)
		 *(QP+i*M+j) = *(Q+i*M+j);
	 else
	   for(j=0; j<M; j++)
		 *(QP+i*M+j) = *(P+i*M+j);
   }
  }// end for pragma

  QuickSort(0,M-1,DTMP,index);

  assert(N<=M);
  #pragma omp parallel default(shared) private(id,k,i,j)
  {
   id = omp_get_thread_num();
   for(k=id*M/4; k<(id+1)*M/4; k++)
   {
      i = *(index+k);
	  if(k<N)
	   *(D+k) = sqrt(*(DTMP+k));
      for(j=0; j<N; j++)
	  {
         *(A+k*N+j) = *(ATMP+i*N+j);
	  }
	  if(MlsN==0)
	   for(j=0; j<M; j++)
		 *(Q+k*M+j) = *(QP+k*M+j);
	  else
	   for(j=0; j<M; j++)
		 *(P+k*M+j) = *(QP+k*M+j);
   }
  }//end for pragma
  
  #pragma omp parallel default(shared) private(id,i,j)
  {
  id = omp_get_thread_num();
  if(MlsN==0)
  {
   for(i=id*N/4; i<(id+1)*N/4; i++)
	for(j=0; j<N; j++)
      *(P+i*N+j) = *(A+i*N+j)/(*(D+i));
  }
  else
  {
   for(i=id*N/4; i<(id+1)*N/4; i++)
    for(j=0; j<N; j++)
      *(Q+i*N+j) = *(A+i*N+j)/(*(D+i));
  }
  }//end for pragma

  free(DTMP);
  free(ATMP);
  free(QP);
}

int oe_gesvd_op_tran(int M,int N,double tol,double *A,double *D,
		  double *Q,double *P,int *INFO)
{
  int i,j,k,w,q;
  double a,b,c;
  double s,t,se;
  double tmp;
  double cs,sn;
  int cond = 0;
  double ele1, ele2;
  int iteration = 0;
  int* index; //索引
  int id;
  int st_index,en_index;
  int exe_count;
  int count_odd_even = 0;
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
  index = (int*)malloc(M*sizeof(int));
  if(index==NULL)
  {
    printf("can't alloc the space of index\n");
	return -2;
  }

  #pragma omp parallel default(shared) private(id,i)
  {
  id = omp_get_thread_num();
  if(MlsN==0)
   {
    //初始化Q为单位阵
    for(i=id*M/4; i<(id+1)*M/4; i++)
    {
      *(Q+i*M+i) = 1.0;
	  *(index+i) = i;
    }
   }
  else
   {
    //初始化P为单位阵
    for(i=id*M/4; i<(id+1)*M/4; i++)
    {
      *(P+i*M+i) = 1.0;
	  *(index+i) = i;
    }
   }
  }// end for pragma

  while(1)
  {
   cond = 0;
   exe_count = 0;
   count_odd_even = 0;
   last_ave_orth = 2*tol_orth/(M*(M-1));
   tol_orth = 0;
   #pragma omp parallel reduction(||:cond) reduction(+:tol_orth) default(shared) private(id,i,j,q,k,w,c,ele1,ele2,se,t,cs,sn,tmp,st_index,en_index)
   {
    id = omp_get_thread_num();
   while(exe_count < M)
  {
	if(count_odd_even == 0)	//处于偶数次
	{
	  st_index = id*M/4;
	  en_index = (id+1)*M/4;
	}
	else
	{
	  st_index = id*M/4+1;
	  en_index = (id+1)*M/4;
	}
    for(q=st_index; q<en_index; q=q+2)
	{
	  if(q+1>=M) break;
	  i = *(index+q);	 
	  j = *(index+q+1);
	  *(index+q) = j;    // 交换index
	  *(index+q+1) = i;
	    
	  c = 0.0;
	  for(k=0; k<N; k++)
         c = c + (*(A+i*N+k))*(*(A+j*N+k));
 

	  if(myabs(c) > tol) 
		cond = 1;  //计算能否跳出循环 
	  else
	  {
	  	  tol_orth = tol_orth + myabs(c);
		  continue;	
	  }
	  if(iteration==0)
	  {
	    tol_orth = tol_orth + myabs(c);
		continue;
	  }
	  if(myabs(c) < last_ave_orth) 
	  {
	    tol_orth = tol_orth + myabs(c);
		continue;
	  }

	  //每次旋转之前判断第i列和第j列的范数，把较大的调在前面，以保证奇异值从大到小排列
	  ele1 = 0.0;
	  ele2 = 0.0;
	  for(w=0; w<N; w++)
	  {
	     ele1 = ele1 + (*(A+i*N+w))*(*(A+i*N+w));
	     ele2 = ele2 + (*(A+j*N+w))*(*(A+j*N+w));
	  }

      se = (ele1 - ele2)/(2*c);
	  t = mysgn(se)/(myabs(se)+sqrt(1+se*se));
	  cs = 1/sqrt(1+t*t);
	  sn = cs * t;
      //更新AT的第i,j列
	  for(k=0; k<N; k++)
	  {
	    tmp = *(A+i*N+k);
		*(A+i*N+k) = cs*tmp + sn*(*(A+j*N+k));
		*(A+j*N+k) = -sn*tmp + cs*(*(A+j*N+k));
	  }
      if(MlsN==0)
	  {
	    //更新右奇异矩阵Q
        for(k=0; k<M; k++)
	    {
	     tmp = *(Q+i*M+k);
		 *(Q+i*M+k) = cs*tmp + sn*(*(Q+j*M+k));
		 *(Q+j*M+k) = -sn*tmp + cs*(*(Q+j*M+k));
	    }
	  }
	  else
	  {
	    //更新右奇异矩阵P
        for(k=0; k<M; k++)
	    {
	     tmp = *(P+i*M+k);
		 *(P+i*M+k) = cs*tmp + sn*(*(P+j*M+k));
		 *(P+j*M+k) = -sn*tmp + cs*(*(P+j*M+k));
        }
	  }
    }//end for for
    #pragma omp master
	{
     exe_count++;
	 if(count_odd_even == 0) count_odd_even = 1;
	 else count_odd_even = 0;
	}
    #pragma omp barrier 

  }// end for while

   }// end for pragma
  
   if(cond == 0) break;

   if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}
  }// end for while
 
  ChangeOrder(M,N,A,Q,P,D,MlsN);  // 按行向量模的大小进行排序 

  free(index);
   
  return 0;
}
