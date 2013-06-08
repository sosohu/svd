#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<f2c.h>
#include<omp.h>
extern int gesvd(int M,int N,double tol,double *A,
			     double *D,double *Q,double *P,int *INFO);  // (1): 普通单边jacobi串行算法
extern int odd_even_gesvd(int M,int N,double tol,double *A, // (2): 普通单边jacobi并行算法
			     double *D,double *Q,double *P,int *INFO); 
extern int myorder_gesvd(int M,int N,double tol,double *A, // (4):  (3)+每轮迭代剔除内积小的计算
			     double *D,double *Q,double *P,int *INFO); 
extern int myorder_gesvd_tran(int M,int N,double tol,double *A, // (4):  (3)+每轮迭代剔除内积小的计算
			     double *D,double *Q,double *P,int *INFO); 
extern int myorder_gesvd_op_tran(int M,int N,double tol,double *A, // (5):   (4)+在M<N时，转置矩阵，使得按列存储
			     double *D,double *Q,double *P,int *INFO); 
extern int myorder_gesvd_all_a(int M,int N,double tol,double *A,
			     double *D,double *Q,double *P,int *INFO); 
extern int myorder_gesvd_op_tran_ele(int M,int N,double tol,double *A, //（6）：  （5）+嵌入汇编代码
			     double *D,double *Q,double *P,int *INFO); 

int output_orthogonal(int size,double *A) 
{ // A*AT
  double *ATA = (double*)malloc(size*size*sizeof(double));
  int i,j,k;
  for(i=0; i<size; i++)
  {
	for(j=0; j<size; j++)
	{
	  *(ATA+i*size+j) = 0.0;
	  for(k=0;k<size;k++)
	     *(ATA+i*size+j) = *(ATA+i*size+j) + (*(A+i*size+k))*(*(A+j*size+k));
	  printf("%lf ",*(ATA+i*size+j));
	}
    printf("\n");
  }
}

int output_orthogonal_con(int size,double *A)
{ // AT*A
  double *ATA = (double*)malloc(size*size*sizeof(double));
  int i,j,k;
  for(i=0; i<size; i++)
  {
	for(j=0; j<size; j++)
	{
	  *(ATA+i*size+j) = 0.0;
	  for(k=0;k<size;k++)
	     *(ATA+i*size+j) = *(ATA+i*size+j) + (*(A+k*size+i))*(*(A+k*size+j));
	  printf("%lf ",*(ATA+i*size+j));
	}
    printf("\n");
  }
}

int main(int argc, char** argv)
{
	int M,N;
	int *INFO;
	int i,j;
	FILE *fp;
	double tol;
	double *A,*ATMP,*ATMP0,*ATMP1,*ATMP2,*ATMP3,*ATMP4,*D,*Q,*P,*PTMP,*QTMP,*PTMP2,*QTMP2;
	int minMN;
	struct timeval tpstart,tpend;
	unsigned long timeuse;
    
	M = atoi(argv[1]);
	N = atoi(argv[2]);
	minMN = N<M?N:M;
	tol = atof(argv[3]);
	A = (double*)malloc(M*N*sizeof(double));
	ATMP = (double*)malloc(M*N*sizeof(double));
	ATMP0 = (double*)malloc(M*N*sizeof(double));
	ATMP1 = (double*)malloc(M*N*sizeof(double));
	ATMP2 = (double*)malloc(M*N*sizeof(double));
	ATMP3 = (double*)malloc(M*N*sizeof(double));
	ATMP4 = (double*)malloc(M*N*sizeof(double));
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
	INFO = (int*)calloc(1,sizeof(int));
	if((NULL==A)||(NULL==D)||(NULL==Q)||(NULL==P)||(NULL==INFO))
    {
	  printf("error: can not alloc the space\n");
	  return -2;
	}
	//从rand_out文件读取初始A矩阵
    fp = fopen(argv[4],"r");
	if(NULL == fp)
	{
	  printf("error: no input file! \n");
	  return -1;
	}
	for(i=0;i<M;i++)
	  for(j=0;j<N;j++)
	  {
		fscanf(fp,"%lf",A+i*N+j);
		*(ATMP+i*N+j) = *(A+i*N+j);
		*(ATMP0+i*N+j) = *(A+i*N+j);
		*(ATMP1+i*N+j) = *(A+i*N+j);
		*(ATMP2+i*N+j) = *(A+i*N+j);
		*(ATMP3+i*N+j) = *(A+i*N+j);
		*(ATMP4+i*N+j) = *(A+i*N+j);
	  }
	
    gettimeofday(&tpstart,NULL);
    //调用gesvd函数
    gesvd(M,N,tol,A,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("gesvd time use:%20ld us\n",timeuse);


    free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    //调用odd_even_gesvd函数
    odd_even_gesvd(M,N,tol,ATMP,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("odd_even_gesvd time use:%20ld us\n",timeuse);
	
    free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    //调用odd_even_gesvd函数
    myorder_gesvd(M,N,tol,ATMP0,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("myorder_gesvd time use:%20ld us\n",timeuse);
	
    free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    //调用odd_even_gesvd_optim函数
    myorder_gesvd_tran(M,N,tol,ATMP1,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("myorder_gesvd_tran time use:%20ld us\n",timeuse);
    
    
    free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    myorder_gesvd_op_tran(M,N,tol,ATMP2,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("myorder_gesvd_op_tran time use:%20ld us\n",timeuse);
   
	
    free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    //调用oe_gesvd_op_tran函数
    myorder_gesvd_op_tran_ele(M,N,tol,ATMP3,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("myorder_gesvd_op_tran_ele time use:%20ld us\n",timeuse);
    
	
	free(D); free(Q); free(P);
	D = (double*)calloc(minMN,sizeof(double));
	Q = (double*)calloc(M*M,sizeof(double));
	P = (double*)calloc(N*N,sizeof(double));
    gettimeofday(&tpstart,NULL);
    //调用odd_even_gesvd_op_a函数
    myorder_gesvd_all_a(M,N,tol,ATMP4,D,Q,P,INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("myorder_gesvd_all_a time use:%20ld us\n",timeuse);

    /*
	//将INFO,D,Q,P输出
	printf("the INFO: %d\n\n",*INFO);
    printf("the D:\n");
	for(i=0; i<minMN; i++)
	  printf("%lf ",*(D+i));
	printf("\n\nthe Q:\n");
	for(i=0; i<M; i++)
	{
	  for(j=0; j<M; j++)
         printf("%lf ",*(Q+i*M+j));
	  printf("\n");
	}
	printf("\nthe P:\n");
	for(i=0; i<N; i++)
	{
	  for(j=0; j<N; j++)
         printf("%lf ",*(P+i*N+j));
	  printf("\n");
	}
	printf("the Q*QT\n");
	output_orthogonal(M,Q);
	printf("the P*PT\n");
	output_orthogonal(N,P);
	*/
	free(A);
	free(ATMP);
	free(ATMP0);
	free(ATMP1);
	free(ATMP2);
	free(ATMP3);
	free(ATMP4);
	free(D);
	free(Q);
	free(P);
	free(INFO);
	fclose(fp);
	return 0;
}

