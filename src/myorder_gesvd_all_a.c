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

/*
void orth_row(int i, int j,int M,int N, double* A,double* P,double* Q,int MlsN,int* cond, double* tol_orth,int iteration,double last_ave_orth,double* ele)
{
  int k,w;
  double c,ele1,ele2,se,cs,sn,t,tmp;
	  c = 0.0;
	  for(k=0; k<N; k++)
         c = c + (*(A+i*N+k))*(*(A+j*N+k));
	  if(myabs(c) > tol) 
		*cond = 1;  //计算能否跳出循环 
	  else
	  {
	  	*tol_orth = *tol_orth + myabs(c);
		return;	
	  }
	  if(iteration==0)
	  {
	    *tol_orth = *tol_orth + myabs(c);
		return;
	  }
	  if(myabs(c) < last_ave_orth) 
	  {
	    *tol_orth = *tol_orth + myabs(c);
		return;
	  }
	  /*
	  ele1 = 0.0;
	  ele2 = 0.0;
	  for(w=0; w<N; w++)
	  {
	     ele1 = ele1 + (*(A+i*N+w))*(*(A+i*N+w));
	     ele2 = ele2 + (*(A+j*N+w))*(*(A+j*N+w));
	  }
	  *
	  ele1 = *(ele+i);
	  ele2 = *(ele+j);
      se = (ele1 - ele2)/(2*c);
	  t = mysgn(se)/(myabs(se)+sqrt(1+se*se));
	  cs = 1/sqrt(1+t*t);
	  sn = cs * t;
	  //更新ele_i ele_j
	  *(ele+i) = cs*cs*ele1 + sn*sn*ele2 + 2*cs*sn*c;
	  *(ele+j) = sn*sn*ele1 + cs*cs*ele2 - 2*cs*sn*c;
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
}
*/


inline void orth_row(int i, int j,int M,int N, double* A,double* P,double* Q,int MlsN,int* cond, double* tol_orth,int iteration,double last_ave_orth,double* ele)
{
  int k,w;
  double c,ele1,ele2,se,cs,sn,t,tmp;
  unsigned long off1,off2;
  unsigned long mask = 16;

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
	    "add %1,%1,16\n\t"
	    "add %2,%2,2\n\t"
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
		*cond = 1;  //计算能否跳出循环 
	  else
	  {
	  	  *tol_orth = *tol_orth + myabs(c);
		  return;	
	  }
	  if(iteration==0)
	  {
	    *tol_orth = *tol_orth + myabs(c);
		return;
	  }
	  if(myabs(c) < last_ave_orth) 
	  {
	    *tol_orth = *tol_orth + myabs(c);
		return;
	  }

	  //每次旋转之前判断第i列和第j列的范数，把较大的调在前面，以保证奇异值从大到小排列
	  /*
	  ele1 = 0.0;
	  ele2 = 0.0;
	  */
	  ele1 = *(ele+i);
	  ele2 = *(ele+j);
	  
	  /*
	  for(w=0; w<N; w++)
	  {
	     ele1 = ele1 + (*(A+i*N+w))*(*(A+i*N+w));
	     //ele2 = ele2 + (*(A+j*N+w))*(*(A+j*N+w));
	  }
	  */
	  

      /*
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
	    //"L7: bge %2,%3,L8\n\t"
		"L7: nop\n\t"
	    "gsLQC1 $f2,$f0,0(%1)\n\t"
	    "madd.d %0,%0,$f2,$f2\n\t"
	    "madd.d %0,%0,$f0,$f0\n\t"
	    "add %2,%2,2\n\t"
	    "add %1,%1,16\n\t"
	    "blt %2,%3,L7\n\t"
	    //"L8:nop\n\t" 
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
*/
      se = (ele1 - ele2)/(2*c);
	  t = mysgn(se)/(myabs(se)+sqrt(1+se*se));
	  cs = 1/sqrt(1+t*t);
	  sn = cs * t;
	  //更新ele_i ele_j
	  *(ele+i) = cs*cs*ele1 + sn*sn*ele2 + 2*cs*sn*c;
	  *(ele+j) = sn*sn*ele1 + cs*cs*ele2 - 2*cs*sn*c;
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
}


void orth_all_a(int start,int end,int M,int N,double* A,double* P,double* Q,int MlsN,
			  int* cond,double* tol_orth,int iteration,double last_ave_orth,double* ele)
{
  int i,j,k,w;
  
  if(start >= end) return;
  for(i=start; i<end; i++)
   for(j=i+1; j<=end; j++)	 
   {
     orth_row(i,j,M,N,A,P,Q,MlsN,cond,tol_orth,iteration,last_ave_orth,ele);
   }//end for for
}

void orth_two_a(int st_1,int en_1,int st_2,int en_2,int M,int N,double* A,double* P, double* Q, int MlsN,
			  int* cond,double* tol_orth,int iteration,double last_ave_orth,double* ele)
{
  int i,j,k,w;
  if((st_1>en_1)||(st_2>en_2))  return;
  if(!((en_1<st_2)||(en_2<st_1))) return;
  for(i=st_1; i<=en_1; i++)
   for(j=st_2; j<=en_2; j++)	 
   {
     orth_row(i,j,M,N,A,P,Q,MlsN,cond,tol_orth,iteration,last_ave_orth,ele);
   }//end for for
}

int myorder_gesvd_all_a(int M,int N,double tol,double *A,double *D,
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
  double* ele;
  
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
  
  ele = (double*)calloc(M,sizeof(double));

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
	  *(ele+i) = 0;
	  for(j=0; j<N; j++)
		*(ele+i) = *(ele+i) + (*(A+i*N+j))*(*(A+i*N+j));
    }
   }
  else
   {
    //初始化P为单位阵
    for(i=id*M/4; i<(id+1)*M/4; i++)
    {
      *(P+i*M+i) = 1.0;
	  *(ele+i) = 0;
	  for(j=0; j<N; j++)
		*(ele+i) = *(ele+i) + (*(A+i*N+j))*(*(A+i*N+j));
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
    orth_all_a(id*M/4,(id+1)*M/4-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);

    #pragma omp barrier

	// (p1,p2) (p3,p4)
    #pragma omp sections
	{
    // void orth_two(int st_1,int en_1,int st_2,int en_2,int M,int N,double* A,double* P, double* Q, int MlsN,
	//	    int* cond,double* tol_orth,int iteration,double last_ave_orth)
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,2*M/8,3*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,3*M/8,4*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(4*M/8,5*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(5*M/8,6*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,3*M/8,4*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,2*M/8,3*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(4*M/8,5*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(5*M/8,6*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}

	// (p1,p3) (p2,p4)
    #pragma omp sections
	{
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(2*M/8,3*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(3*M/8,4*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(2*M/8,3*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(3*M/8,4*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}

	// (p1,p4) (p2,p3)
    #pragma omp sections
	{
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(2*M/8,3*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(3*M/8,4*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}
    #pragma omp sections
	{
	  #pragma omp section
      orth_two_a(0*M/8,1*M/8-1,7*M/8,8*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(1*M/8,2*M/8-1,6*M/8,7*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(2*M/8,3*M/8-1,5*M/8,6*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	  #pragma omp section 
      orth_two_a(3*M/8,4*M/8-1,4*M/8,5*M/8-1,M,N,A,P,Q,MlsN,&cond,&tol_orth,iteration,last_ave_orth,ele);
	}
   }// end for pragma
  
   if(cond == 0) break;

   if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}
  }// end for while
 
  ChangeOrder(M,N,A,Q,P,D,MlsN);  // 按行向量模的大小进行排序 
   
  return 0;
}
