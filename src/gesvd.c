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
#include "omp.h"

#define ITERATION 20

double is_zero = 0.000001;

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

int gesvd(int M,int N,double tol,double *A,double *D,
		  double *Q,double *P,int *INFO)
{
  int i,j,k,w;
  double a,b,c;
  double s,t,se;
  double tmp;
  double cs,sn;
  int cond = 0;
  double ele1, ele2;
  int iteration = 0;
  
  *INFO = 0;
  if(M < 0)  *INFO = -1;
  else if(N < 0) *INFO = -2;

  if(*INFO < 0)
  {
    printf("error:%d\n",*INFO);
	exit(0);
  }
  
  
  if(M<=N)
{
  //初始化P为单位阵
  for(i=0; i<N; i++)
      *(P+i*N+i) = 1.0;

  //主计算过程
  while(cond==0)
  {
  cond = 1;
  //对于每个i<j迭代计算
  for(i=0; i<N-1; i++)
	for(j=i+1; j<N; j++)
	{
	  c = 0.0;
	  for(k=0; k<M; k++)
         c = c + (*(A+k*N+i))*(*(A+k*N+j));

	  if(myabs(c) > tol) cond = 0;  //计算能否跳出循环 
	  else continue;
	  //每次旋转之前判断第i列和第j列的范数，把较大的调在前面，以保证奇异值从大到小排列
	  ele1 = 0.0;
	  ele2 = 0.0;
	  for(w=0; w<M; w++)
	  {
	     ele1 = ele1 + (*(A+w*N+i))*(*(A+w*N+i));
	     ele2 = ele2 + (*(A+w*N+j))*(*(A+w*N+j));
	  }
	  if(ele1<ele2) //交换第i,j列
	  {
	    double exech;
		for(w=0; w<M; w++)
		{
		  exech = *(A+w*N+i);
		  *(A+w*N+i) = *(A+w*N+j);
		  *(A+w*N+j) = exech;
		}
		for(w=0; w<N; w++)
		{
		  exech = *(P+w*N+i);
		  *(P+w*N+i) = *(P+w*N+j);
		  *(P+w*N+j) = exech;
		}
		exech = ele1;
		ele1 = ele2;
		ele2 = exech;
	  }

      se = (ele1 - ele2)/(2*c);
	  t = mysgn(se)/(myabs(se)+sqrt(1+se*se));
	  cs = 1/sqrt(1+t*t);
	  sn = cs * t;
      //更新A的第i,j列
	  for(k=0; k<M; k++)
	  {
	    tmp = *(A+k*N+i);
		*(A+k*N+i) = cs*tmp + sn*(*(A+k*N+j));
		*(A+k*N+j) = -sn*tmp + cs*(*(A+k*N+j));
	  }
	  //更新右奇异矩阵P
      for(k=0; k<N; k++)
	  {
	    tmp = *(P+k*N+i);
		*(P+k*N+i) = cs*tmp + sn*(*(P+k*N+j));
		*(P+k*N+j) = -sn*tmp + cs*(*(P+k*N+j));
	  }
	}
    if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}
  }// end while
  //计算A每列的F范数，得到奇异值:该奇异值可能为0不是最终奇异值
  for(i=0; i<M; i++)
  {
	*(D+i) = 0.0;
	for(j=0; j<M; j++)
      *(D+i) = *(D+i) + (*(A+j*N+i))*(*(A+j*N+i));
	*(D+i) = sqrt(*(D+i));
  }
  //A的各列为相应奇异值的左奇异变量
  for(i=0; i<M; i++)
  {
    if((*(D+i)>-is_zero)&&(*(D+i)<is_zero))
	     ;
	else
	  for(j=0;j<M;j++)
        *(Q+j*M+i) = *(A+j*N+i)/(*(D+i));
  }
 }//end if(M<=N)===================================================================================
  //=================================================================================================
 else
 {
  //初始化Q为单位阵
  for(i=0; i<M; i++)
      *(Q+i*M+i) = 1.0;

  //主计算过程
  while(cond==0)
  {
  cond = 1;
  //对于每个i<j迭代计算
  for(i=0; i<M-1; i++)
	for(j=i+1; j<M; j++)
	{
	  
	  c = 0.0;
	  for(k=0; k<N; k++)
         c = c + (*(A+i*N+k))*(*(A+j*N+k));

	  if(myabs(c) > tol) cond = 0;  //计算能否跳出循环 
	  else continue;
	  //每次旋转之前判断第i列和第j列的范数，把较大的调在前面，以保证奇异值从大到小排列
	  ele1 = 0.0;
	  ele2 = 0.0;
	  for(w=0; w<N; w++)
	  {
	     ele1 = ele1 + (*(A+i*N+w))*(*(A+i*N+w));
	     ele2 = ele2 + (*(A+j*N+w))*(*(A+j*N+w));
	  }
	  if(ele1<ele2) //交换AT第i,j列
	  {
	    double exech;
		for(w=0; w<N; w++)
		{
		  exech = *(A+i*N+w);
		  *(A+i*N+w) = *(A+j*N+w);
		  *(A+j*N+w) = exech;
		}
		for(w=0; w<M; w++)
		{
		  exech = *(Q+i*M+w);
		  *(Q+i*M+w) = *(Q+j*M+w);
		  *(Q+j*M+w) = exech;
		}
		exech = ele1;
		ele1 = ele2;
		ele2 = exech;
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
	  //更新右奇异矩阵Q
      for(k=0; k<M; k++)
	  {
	    tmp = *(Q+i*M+k);
		*(Q+i*M+k) = cs*tmp + sn*(*(Q+j*M+k));
		*(Q+j*M+k) = -sn*tmp + cs*(*(Q+j*M+k));
	  }
	}
    if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}
  }// end while
  //计算A每列的F范数，得到奇异值:该奇异值可能为0不是最终奇异值
  for(i=0; i<N; i++)
  {
	*(D+i) = 0.0;
	for(j=0; j<N; j++)
      *(D+i) = *(D+i) + (*(A+i*N+j))*(*(A+i*N+j));
	*(D+i) = sqrt(*(D+i));
  }
  //A的各列为相应奇异值的左奇异变量
  for(i=0; i<N; i++)
  {
    if((*(D+i)>-is_zero)&&(*(D+i)<is_zero))
	  ;
	else
	  for(j=0;j<N;j++)
        *(P+i*N+j) = *(A+i*N+j)/(*(D+i));
  }
 }//end else
   
  return 0;
}
