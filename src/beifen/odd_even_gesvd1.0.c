/*
   The openmp program in c of the svd for a general real M-by-N matrix A
   The code is provided by Mr. Hu at the time of 2013.3 . 
   His email address is ustc.sosohu@gmail.com
 */

/*
    Purpose   
    ====================================================================   
    GEBRD reduces a general real M-by-N matrix A to upper or lower        
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

extern real is_zero;
extern int mysgn(real s);
extern real myabs(real s);

//real is_zero = 0.000001;

/*
int mysgn(real s)
{
  if((s<is_zero)&&(s>-is_zero)) return 0;
  if(s>0.0)  return 1;
  return -1;
}

real myabs(real s)
{
  if(s<0.0) return -s;
  return s;
}
*/

int odd_even_gesvd(int M,int N,real tol,real *A,real *D,
		  real *Q,real *P,int *INFO)
{
  int i,j,k,w,q;
  real a,b,c;
  real s,t,se;
  real tmp;
  real cs,sn;
  int cond = 0;
  real ele1, ele2;
  int iteration = 0;
  int* index; //索引
  int id;
  int st_index,en_index;
  int exe_count;
  int count_odd_even = 0;
  
  *INFO = 0;
  if(M < 0)  *INFO = -1;
  else if(N < 0) *INFO = -2;

  if(*INFO < 0)
  {
    printf("error:%d\n",*INFO);
	exit(0);
  }
  
  omp_set_num_threads(NUM_THREADS);
  
  if(M<=N)
{
  assert(N%8==0);
  index = (int*)malloc(N*sizeof(int));
  if(index==NULL)
  {
    printf("can't alloc the space of index\n");
	return -2;
  }
  
  #pragma omp parallel default(shared) private(id,i)
  {
  id = omp_get_thread_num();
  //初始化P为单位阵
  for(i=id*N/4; i<(id+1)*N/4; i++)
  {
     *(P+i*N+i) = 1.0;
	 *(index+i) = i;
  }
  }// end for pragma

  while(1)
  {
   cond = 0;
   exe_count = 0;
   count_odd_even = 0;
   #pragma omp parallel reduction(||:cond) default(shared) private(id,i,j,q,k,w,c,ele1,ele2,se,t,cs,sn,tmp,st_index,en_index)
   {
    id = omp_get_thread_num();
   while(exe_count < N)
  {
	if(count_odd_even == 0)	//处于偶数次
	{
	  st_index = id*N/4;
	  en_index = (id+1)*N/4;
	}
	else
	{
	  st_index = id*N/4+1;
	  en_index = (id+1)*N/4;
	}
    for(q=st_index; q<en_index; q=q+2)
	{
	  if(q+1>=N) break;
	  i = *(index+q);	 
	  j = *(index+q+1);
	  *(index+q) = j;    // 交换index
	  *(index+q+1) = i;
      //对A的第i,j列做正交化
	  c = 0.0;
	  for(k=0; k<M; k++)
         c = c + (*(A+k*N+i))*(*(A+k*N+j));
	  if(myabs(c) > tol)
		  cond = 1;  //不正交
	  else  continue; //正交

	  //每次旋转之前判断第i列和第j列的范数，把较大的调在前面，以保证奇异值从大到小排列
	  ele1 = 0.0;
	  ele2 = 0.0;
	  for(w=0; w<M; w++)
	  {
	     ele1 = ele1 + (*(A+w*N+i))*(*(A+w*N+i));
	     ele2 = ele2 + (*(A+w*N+j))*(*(A+w*N+j));
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
      
	}//end for for
    #pragma omp master
	{
     exe_count++;
	 if(count_odd_even == 0) count_odd_even = 1;
	 else count_odd_even = 0;
	}
    #pragma omp barrier 
   }//end for while
   }//end for pragma
   
   if(cond == 0) break;

   if((++iteration) > ITERATION) {printf("the iteration is over\n"); break;}

  }//end for while

  assert(M%4==0);
  //计算A每列的F范数，得到奇异值:该奇异值可能为0不是最终奇异值
  #pragma omp parallel default(shared) private(id,i,j)
  {
   id = omp_get_thread_num();
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
	*(D+i) = 0.0;
	for(j=0; j<M; j++)
      *(D+i) = *(D+i) + (*(A+j*N+i))*(*(A+j*N+i));
	*(D+i) = sqrt(*(D+i));
  }
  //A的各列为相应奇异值的左奇异变量
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
    if((*(D+i)>-is_zero)&&(*(D+i)<is_zero))
	     ;
	else
	  for(j=0;j<M;j++)
        *(Q+j*M+i) = *(A+j*N+i)/(*(D+i));
  }
  }// end for pragma
  free(index);
 }//end if(M<=N)===================================================================================
  //=================================================================================================
 else
 {
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
  //初始化Q为单位阵
  for(i=id*M/4; i<(id+1)*M/4; i++)
  {
      *(Q+i*M+i) = 1.0;
	  *(index+i) = i;
  }
  }// end for pragma

  while(1)
  {
   cond = 0;
   exe_count = 0;
   count_odd_even = 0;
   #pragma omp parallel reduction(||:cond) default(shared) private(id,i,j,q,k,w,c,ele1,ele2,se,t,cs,sn,tmp,st_index,en_index)
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

	  if(myabs(c) > tol) cond = 1;  //计算能否跳出循环 
	  else continue;
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
	  //更新右奇异矩阵Q
      for(k=0; k<M; k++)
	  {
	    tmp = *(Q+i*M+k);
		*(Q+i*M+k) = cs*tmp + sn*(*(Q+j*M+k));
		*(Q+j*M+k) = -sn*tmp + cs*(*(Q+j*M+k));
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

  assert(N%4==0);
  //计算A每列的F范数，得到奇异值:该奇异值可能为0不是最终奇异值
  #pragma omp parallel default(shared) private(id,i,j)
  {
   id = omp_get_thread_num();
  for(i=id*N/4; i<(id+1)*N/4; i++)
  {
	*(D+i) = 0.0;
	for(j=0; j<N; j++)
      *(D+i) = *(D+i) + (*(A+i*N+j))*(*(A+i*N+j));
	*(D+i) = sqrt(*(D+i));
  }
  //A的各列为相应奇异值的左奇异变量
  for(i=id*N/4; i<(id+1)*N/4; i++)
  {
    if((*(D+i)>-is_zero)&&(*(D+i)<is_zero))
	  ;
	else
	  for(j=0;j<N;j++)
        *(P+i*N+j) = *(A+i*N+j)/(*(D+i));
  }
  }//end for pragma
  free(index);
 }//end else
   
  return 0;
}
