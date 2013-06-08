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
extern int mysgn(double s);
extern double myabs(double s);
extern void TransposeQuick(int M, int N,double *A);
extern int Partion(int s,int t,double* A,int* index);
extern void QuickSort(int s,int t,double* A,int* index);
extern void  ChangeOrder(int M,int N,double* A,double* Q,double* P,double* D,int MlsN);  // 按行向量模的大小进行排序 

/*
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
*/

int odd_even_gesvd_op_a(int M,int N,double tol,double *A,double *D,
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
  unsigned long mask = 16;
  unsigned long off1,off2;
  
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
   #pragma omp parallel reduction(||:cond) reduction(+:tol_orth) default(shared) private(id,i,j,q,k,w,c,ele1,ele2,se,t,cs,sn,tmp,st_index,en_index,off1,off2)
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
	  /*
	  for(k=0; k<N; k++)
         c = c + (*(A+i*N+k))*(*(A+j*N+k)); 
      */
       //start----------------------------------------------------------
      off1 = (unsigned long)(A+i*N) % mask;
      off2 = (unsigned long)(A+j*N) % mask;
	  if(off1 == off2)
	  {
	    if(off1 == 0)
		  k = 0;
		else
		{
		  k = 1;
          c = c + (*(A+i*N))*(*(A+j*N)); 
          c = c + (*(A+i*N+N-1))*(*(A+j*N+N-1)); 
		}
       asm volatile
    	(
	    "mov.d %0,%4\n\t"
	    "L1: bge %2,%3,L2\n\t"
	    "gsLQC1 $f2,$f0,0(%1)\n\t"
	    "gsLQC1 $f6,$f4,0(%5)\n\t"
	    "madd.d %0,%0,$f6,$f2\n\t"
	    "madd.d %0,%0,$f4,$f0\n\t"
	    "add %2,%2,2\n\t"
	    "add %1,%1,16\n\t"
	    "add %5,%5,16\n\t"
	    "j L1\n\t"
	    "L2: nop\n\t" 
        :"=f"(c)
	    :"r"(A+i*N+k),"r"(k),"r"(N-1),"f"(c),"r"(A+j*N+k)
	    :"$f0","$f2","$f4","$f6"
    	);
	 }
	 else
	 {
	     if(off1 == 0)
		 {
		   k = 0;
           c = c + (*(A+i*N))*(*(A+j*N)); 
           asm volatile
           (
	       "mov.d %0,%4\n\t"
		   "l.d $f8,0\n\t"
	       "L3: bge %2,%3,L4\n\t"
	       "gsLQC1 $f2,$f0,0(%1)\n\t"
	       "gsLQC1 $f6,$f4,0(%5)\n\t"
	       "madd.d %0,%0,$f4,$f2\n\t"
	       "madd.d %0,%0,$f8,$f0\n\t"
		   "mov.d $f8,$f6\n\t"
	       "add %2,%2,2\n\t"
	       "add %1,%1,16\n\t"
	       "add %5,%5,16\n\t"
	       "j L3\n\t"
	       "L4: nop\n\t" 
           :"=f"(c)
	       :"r"(A+i*N+k),"r"(k),"r"(N-3),"f"(c),"r"(A+j*N+k+1)
	       :"$f0","$f2","$f4","$f6","$f8"
           );
           c = c + (*(A+i*N+N-2))*(*(A+j*N+N-2)); 
           c = c + (*(A+i*N+N-1))*(*(A+j*N+N-1)); 
		 }
		 else
		 {
		   k = 0;
           c = c + (*(A+i*N))*(*(A+j*N)); 
           asm volatile
           (
	       "mov.d %0,%4\n\t"
		   "l.d $f8,0\n\t"
	       "L5: bge %2,%3,L6\n\t"
	       "gsLQC1 $f2,$f0,0(%1)\n\t"
	       "gsLQC1 $f6,$f4,0(%5)\n\t"
	       "madd.d %0,%0,$f4,$f2\n\t"
	       "madd.d %0,%0,$f8,$f0\n\t"
		   "mov.d $f8,$f6\n\t"
	       "add %2,%2,2\n\t"
	       "add %1,%1,16\n\t"
	       "add %5,%5,16\n\t"
	       "j L5\n\t"
	       "L6: nop\n\t" 
           :"=f"(c)
	       :"r"(A+j*N+k),"r"(k),"r"(N-3),"f"(c),"r"(A+i*N+k+1)
	       :"$f0","$f2","$f4","$f6","$f8"
           );
           c = c + (*(A+i*N+N-2))*(*(A+j*N+N-2)); 
           c = c + (*(A+i*N+N-1))*(*(A+j*N+N-1)); 
		 }
	   }
         

	   //end----------------------------------------------------------------


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
	  /*
	  for(w=0; w<N; w++)
	  {
	     ele1 = ele1 + (*(A+i*N+w))*(*(A+i*N+w));
	     ele2 = ele2 + (*(A+j*N+w))*(*(A+j*N+w));
	  }
	  */


	  //start------------------------------------------------------------------------------
      if(off1 == 0)
		k = 0;
	  else
	  {
	    k = 1;
	    ele1 = ele1 + (*(A+i*N))*(*(A+i*N));
	    ele1 = ele1 + (*(A+i*N+N-1))*(*(A+i*N+N-1));
	  }
	  
       asm volatile
    	(
	    "mov.d %0,%4\n\t"
	    "L7: bge %2,%3,L8\n\t"
	    "gsLQC1 $f2,$f0,0(%1)\n\t"
	    "madd.d %0,%0,$f2,$f2\n\t"
	    "madd.d %0,%0,$f0,$f0\n\t"
	    "add %2,%2,2\n\t"
	    "add %1,%1,16\n\t"
	    "j L7\n\t"
	    "L8: nop\n\t" 
        :"=f"(ele1)
	    :"r"(A+i*N+k),"r"(k),"r"(N-1),"f"(ele1)
	    :"$f0","$f2"
    	);
	

      if(off2 == 0)
		k = 0;
	  else
	  {
	    k = 1;
	    ele2 = ele2 + (*(A+j*N))*(*(A+j*N));
	    ele2 = ele2 + (*(A+j*N+N-1))*(*(A+j*N+N-1));
	  }
       asm volatile
    	(
	    "mov.d %0,%4\n\t"
	    "L9: bge %2,%3,L10\n\t"
	    "gsLQC1 $f2,$f0,0(%1)\n\t"
	    "madd.d %0,%0,$f2,$f2\n\t"
	    "madd.d %0,%0,$f0,$f0\n\t"
	    "add %2,%2,2\n\t"
	    "add %1,%1,16\n\t"
	    "j L9\n\t"
	    "L10: nop\n\t" 
        :"=f"(ele2)
	    :"r"(A+j*N+k),"r"(k),"r"(N-1),"f"(ele2)
	    :"$f0","$f2"
    	);

	  //end------------------------------------------------------------------

      se = (ele1 - ele2)/(2*c);
	  t = mysgn(se)/(myabs(se)+sqrt(1+se*se));
	  cs = 1/sqrt(1+t*t);
	  sn = cs * t;
      //更新AT的第i,j列
	  /*
	  for(k=0; k<N; k++)
	  {
	    tmp = *(A+i*N+k);
		*(A+i*N+k) = cs*tmp + sn*(*(A+j*N+k));
		*(A+j*N+k) = -sn*tmp + cs*(*(A+j*N+k));
	  }
	  */
      
	  //start---------------------------------------------------------
       if(off1 == off2)
	   {
	     if(off1 == 0)
		   k = 0;
		 else
		 {
		   k = 1;
	       tmp = *(A+i*N);
		   *(A+i*N) = cs*tmp + sn*(*(A+j*N));
		   *(A+j*N) = -sn*tmp + cs*(*(A+j*N));
	       tmp = *(A+i*N+N-1);
		   *(A+i*N+N-1) = cs*tmp + sn*(*(A+j*N+N-1));
		   *(A+j*N+N-1) = -sn*tmp + cs*(*(A+j*N+N-1));
		 }
         asm volatile
      	 (
	     "L11: bge %1,%2,L12\n\t"
	     "gsLQC1 $f2,$f0,0(%0)\n\t"
	     "gsLQC1 $f6,$f4,0(%3)\n\t"
         "mov.d %7,$f0\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f0,$f8,%5,$f4\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f4,$f8,%4,$f4\n\t"
         "mov.d %7,$f2\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f2,$f8,%5,$f6\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f6,$f8,%4,$f6\n\t"
	     "gsSQC1 $f2,$f0,0(%0)\n\t"
	     "gsSQC1 $f6,$f4,0(%3)\n\t"
	     "add %1,%1,2\n\t"
	     "add %0,%0,16\n\t"
	     "add %3,%3,16\n\t"
	     "j L11\n\t"
	     "L12: nop\n\t" 
         :
	     :"r"(A+i*N+k),"r"(k),"r"(N-1),"r"(A+j*N+k),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp)
	     :"$f0","$f2","$f4","$f6","$f8"
    	 );
	   }
	   else
	   {
	     if(off1 == 0)
		 {
		   k = 0;
	       tmp = *(A+i*N);
		   *(A+i*N) = cs*tmp + sn*(*(A+j*N));
		   *(A+j*N) = -sn*tmp + cs*(*(A+j*N));
	       tmp = *(A+i*N+N-2);
		   *(A+i*N+N-2) = cs*tmp + sn*(*(A+j*N+N-2));
		   *(A+j*N+N-2) = -sn*tmp + cs*(*(A+j*N+N-2));
	       tmp = *(A+i*N+N-1);
		   *(A+i*N+N-1) = cs*tmp + sn*(*(A+j*N+N-1));
		   *(A+j*N+N-1) = -sn*tmp + cs*(*(A+j*N+N-1));
           asm volatile
      	   (
	       "L13: bge %1,%2,L14\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L13\n\t"
	       "L14: nop\n\t" 
           :
	       :"r"(A+i*N+k),"r"(k),"r"(N-3),"r"(A+j*N+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(A+j*N+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 }
		 else
		 {
		   k = 0;
	       tmp = *(A+i*N);
		   *(A+i*N) = cs*tmp + sn*(*(A+j*N));
		   *(A+j*N) = -sn*tmp + cs*(*(A+j*N));
	       tmp = *(A+i*N+N-2);
		   *(A+i*N+N-2) = cs*tmp + sn*(*(A+j*N+N-2));
		   *(A+j*N+N-2) = -sn*tmp + cs*(*(A+j*N+N-2));
	       tmp = *(A+i*N+N-1);
		   *(A+i*N+N-1) = cs*tmp + sn*(*(A+j*N+N-1));
		   *(A+j*N+N-1) = -sn*tmp + cs*(*(A+j*N+N-1));
           asm volatile
      	   (
	       "L15: bge %1,%2,L16\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L15\n\t"
	       "L16: nop\n\t" 
           :
	       :"r"(A+j*N+k),"r"(k),"r"(N-3),"r"(A+i*N+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(A+i*N+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 
		 }
	   }
	  //end---------------------------------------------------------


      if(MlsN==0)
	  {
	    //更新右奇异矩阵Q
		/*
        for(k=0; k<M; k++)
	    {
	     tmp = *(Q+i*M+k);
		 *(Q+i*M+k) = cs*tmp + sn*(*(Q+j*M+k));
		 *(Q+j*M+k) = -sn*tmp + cs*(*(Q+j*M+k));
	    }
		*/
      
	  off1 = (unsigned long)(Q+i*M) % mask;
      off2 = (unsigned long)(Q+j*M) % mask;

	  //start---------------------------------------------------------
       if(off1 == off2)
	   {
	     if(off1 == 0)
		   k = 0;
		 else
		 {
		   k = 1;
	       tmp = *(Q+i*M);
		   *(Q+i*M) = cs*tmp + sn*(*(Q+j*M));
		   *(Q+j*M) = -sn*tmp + cs*(*(Q+j*M));
	       tmp = *(Q+i*M+M-1);
		   *(Q+i*M+M-1) = cs*tmp + sn*(*(Q+j*M+M-1));
		   *(Q+j*M+M-1) = -sn*tmp + cs*(*(Q+j*M+M-1));
		 }
         asm volatile
      	 (
	     "L17: bge %1,%2,L18\n\t"
	     "gsLQC1 $f2,$f0,0(%0)\n\t"
	     "gsLQC1 $f6,$f4,0(%3)\n\t"
         "mov.d %7,$f0\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f0,$f8,%5,$f4\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f4,$f8,%4,$f4\n\t"
         "mov.d %7,$f2\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f2,$f8,%5,$f6\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f6,$f8,%4,$f6\n\t"
	     "gsSQC1 $f2,$f0,0(%0)\n\t"
	     "gsSQC1 $f6,$f4,0(%3)\n\t"
	     "add %1,%1,2\n\t"
	     "add %0,%0,16\n\t"
	     "add %3,%3,16\n\t"
	     "j L17\n\t"
	     "L18: nop\n\t" 
         :
	     :"r"(Q+i*M+k),"r"(k),"r"(M-1),"r"(Q+j*M+k),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp)
	     :"$f0","$f2","$f4","$f6","$f8"
    	 );
	   }
	   else
	   {
	     if(off1 == 0)
		 {
		   k = 0;
	       tmp = *(Q+i*M);
		   *(Q+i*M) = cs*tmp + sn*(*(Q+j*M));
		   *(Q+j*M) = -sn*tmp + cs*(*(Q+j*M));
	       tmp = *(Q+i*M+M-2);
		   *(Q+i*M+M-2) = cs*tmp + sn*(*(Q+j*M+M-2));
		   *(Q+j*M+M-2) = -sn*tmp + cs*(*(Q+j*M+M-2));
	       tmp = *(Q+i*M+M-1);
		   *(Q+i*M+M-1) = cs*tmp + sn*(*(Q+j*M+M-1));
		   *(Q+j*M+M-1) = -sn*tmp + cs*(*(Q+j*M+M-1));
           asm volatile
      	   (
	       "L19: bge %1,%2,L20\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L19\n\t"
	       "L20: nop\n\t" 
           :
	       :"r"(Q+i*M+k),"r"(k),"r"(M-3),"r"(Q+j*M+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(Q+j*M+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 }
		 else
		 {
		   k = 0;
	       tmp = *(Q+i*M);
		   *(Q+i*M) = cs*tmp + sn*(*(Q+j*M));
		   *(Q+j*M) = -sn*tmp + cs*(*(Q+j*M));
	       tmp = *(Q+i*M+M-2);
		   *(Q+i*M+M-2) = cs*tmp + sn*(*(Q+j*M+M-2));
		   *(Q+j*M+M-2) = -sn*tmp + cs*(*(Q+j*M+M-2));
	       tmp = *(Q+i*M+M-1);
		   *(Q+i*M+M-1) = cs*tmp + sn*(*(Q+j*M+M-1));
		   *(Q+j*M+M-1) = -sn*tmp + cs*(*(Q+j*M+M-1));
           asm volatile
      	   (
	       "L21: bge %1,%2,L22\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L21\n\t"
	       "L22: nop\n\t" 
           :
	       :"r"(Q+j*M+k),"r"(k),"r"(M-3),"r"(Q+i*M+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(Q+i*M+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 
		 }
	   }
	  //end---------------------------------------------------------
	  }
	  else
	  {
	    //更新右奇异矩阵P
		/*
        for(k=0; k<M; k++)
	    {
	     tmp = *(P+i*M+k);
		 *(P+i*M+k) = cs*tmp + sn*(*(P+j*M+k));
		 *(P+j*M+k) = -sn*tmp + cs*(*(P+j*M+k));
        }
		*/

	  off1 = (unsigned long)(P+i*M) % mask;
      off2 = (unsigned long)(P+j*M) % mask;

	  //start---------------------------------------------------------
       if(off1 == off2)
	   {
	     if(off1 == 0)
		   k = 0;
		 else
		 {
		   k = 1;
	       tmp = *(P+i*M);
		   *(P+i*M) = cs*tmp + sn*(*(P+j*M));
		   *(P+j*M) = -sn*tmp + cs*(*(P+j*M));
	       tmp = *(P+i*M+M-1);
		   *(P+i*M+M-1) = cs*tmp + sn*(*(P+j*M+M-1));
		   *(P+j*M+M-1) = -sn*tmp + cs*(*(P+j*M+M-1));
		 }
         asm volatile
      	 (
	     "L23: bge %1,%2,L24\n\t"
	     "gsLQC1 $f2,$f0,0(%0)\n\t"
	     "gsLQC1 $f6,$f4,0(%3)\n\t"
         "mov.d %7,$f0\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f0,$f8,%5,$f4\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f4,$f8,%4,$f4\n\t"
         "mov.d %7,$f2\n\t"
		 "mul.d $f8,%4,%7\n\t"
		 "madd.d $f2,$f8,%5,$f6\n\t"
		 "mul.d $f8,%6,%7\n\t"
		 "madd.d $f6,$f8,%4,$f6\n\t"
	     "gsSQC1 $f2,$f0,0(%0)\n\t"
	     "gsSQC1 $f6,$f4,0(%3)\n\t"
	     "add %1,%1,2\n\t"
	     "add %0,%0,16\n\t"
	     "add %3,%3,16\n\t"
	     "j L23\n\t"
	     "L24: nop\n\t" 
         :
	     :"r"(P+i*M+k),"r"(k),"r"(M-1),"r"(P+j*M+k),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp)
	     :"$f0","$f2","$f4","$f6","$f8"
    	 );
	   }
	   else
	   {
	     if(off1 == 0)
		 {
		   k = 0;
	       tmp = *(P+i*M);
		   *(P+i*M) = cs*tmp + sn*(*(P+j*M));
		   *(P+j*M) = -sn*tmp + cs*(*(P+j*M));
	       tmp = *(P+i*M+M-2);
		   *(P+i*M+M-2) = cs*tmp + sn*(*(P+j*M+M-2));
		   *(P+j*M+M-2) = -sn*tmp + cs*(*(P+j*M+M-2));
	       tmp = *(P+i*M+M-1);
		   *(P+i*M+M-1) = cs*tmp + sn*(*(P+j*M+M-1));
		   *(P+j*M+M-1) = -sn*tmp + cs*(*(P+j*M+M-1));
           asm volatile
      	   (
	       "L25: bge %1,%2,L26\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L25\n\t"
	       "L26: nop\n\t" 
           :
	       :"r"(P+i*M+k),"r"(k),"r"(M-3),"r"(P+j*M+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(P+j*M+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 }
		 else
		 {
		   k = 0;
	       tmp = *(P+i*M);
		   *(P+i*M) = cs*tmp + sn*(*(P+j*M));
		   *(P+j*M) = -sn*tmp + cs*(*(P+j*M));
	       tmp = *(P+i*M+M-2);
		   *(P+i*M+M-2) = cs*tmp + sn*(*(P+j*M+M-2));
		   *(P+j*M+M-2) = -sn*tmp + cs*(*(P+j*M+M-2));
	       tmp = *(P+i*M+M-1);
		   *(P+i*M+M-1) = cs*tmp + sn*(*(P+j*M+M-1));
		   *(P+j*M+M-1) = -sn*tmp + cs*(*(P+j*M+M-1));
           asm volatile
      	   (
	       "L27: bge %1,%2,L28\n\t"
	       "gsLQC1 $f2,$f0,0(%0)\n\t"
	       "gsLQC1 $f6,$f4,0(%3)\n\t"
           "mov.d %7,$f0\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f0,$f8,%5,%8\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d %8,$f8,%4,%8\n\t"
           "mov.d %7,$f2\n\t"
		   "mul.d $f8,%4,%7\n\t"
		   "madd.d $f2,$f8,%5,$f4\n\t"
		   "mul.d $f8,%6,%7\n\t"
		   "madd.d $f4,$f8,%4,$f4\n\t"
	       "gsSQC1 $f2,$f0,0(%0)\n\t"
	       "gsSQC1 $f4,%8,0(%3)\n\t"
		   "mov.d %8,$f6\n\t"
	       "add %1,%1,2\n\t"
	       "add %0,%0,16\n\t"
	       "add %3,%3,16\n\t"
	       "j L27\n\t"
	       "L28: nop\n\t" 
           :
	       :"r"(P+j*M+k),"r"(k),"r"(M-3),"r"(P+i*M+k+1),"f"(cs),"f"(sn),"f"(-sn),"f"(tmp),"f"(P+i*M+k)
	       :"$f0","$f2","$f4","$f6","$f8"
    	   );
		 
		 }
	   }
	  //end---------------------------------------------------------


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
