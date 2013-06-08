#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include<omp.h>
#include<INCLUDE/f2c.h>
#include<INCLUDE/clapack.h>


extern int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, 
		   doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);

integer MAX(integer m, integer n)
{
  integer max;
  max = m>n? m:n;
  return max;
}

integer MIN(integer m, integer n)
{
  integer max;
  max = m<n? m:n;
  return max;
}

int main(int argc, char** argv)
{
	integer  *M,*N;
	integer  *INFO,*LDA,*LDU,*LDVT;
	integer  i,j;
	FILE *fp;
	doublereal *A,*S,*U,*VT,*WORK,*E;
	integer  minMN;
	integer  *LWORK;
	struct timeval tpstart,tpend;
	unsigned long timeuse;
	char c = 'A';
	char* jobu = &c;
	char* jobvt = &c;
    

	*M = atoi(argv[1]);
	*N = atoi(argv[2]);
	minMN = (*N)<(*M)?(*N):(*M);
	*LDA = *M;
	*LDU = *M;
	*LDVT = *N;
	LWORK = MAX(3*MIN(*M,*N)+MAX(*M,*N),5*MIN(*M,*N));
	A = (doublereal*)malloc((*LDA)*(*N)*sizeof(doublereal));
	S = (doublereal*)malloc(minMN*sizeof(doublereal));
	U = (doublereal*)malloc((*LDU)*(*M)*sizeof(doublereal));
	VT = (doublereal*)malloc((*LDVT)*(*N)*sizeof(doublereal));
	WORK =  (doublereal*)calloc(*LWORK,sizeof(doublereal));
	INFO = (integer*)calloc(1,sizeof(integer));
	if((NULL==A)||(NULL==S)||(NULL==U)||(NULL==VT)||(NULL==WORK))
    {
	  printf("error: can not alloc the space\n");
	  return -2;
	}
	//从rand_out文件读取初始A矩阵
    fp = fopen(argv[3],"r");
	if(NULL == fp)
	{
	  printf("error: no input file! \n");
	  return -1;
	}
	for(i=1;i<*M;i++)
	  for(j=1;j<*N;j++)
	  {
		fscanf(fp,"%lf",A+i*(*N)+j);
	  }
	
	
	printf("M:%ld, N: %ld, LDA: %ld\n",*M,*N,*LDA);
    gettimeofday(&tpstart,NULL);
    //调用gesvd函数
    dgesvd_(jobu, jobvt, M, N, A, LDA, S, U, LDU, VT, &LDVT, WORK, LWORK, INFO);
    gettimeofday(&tpend,NULL);
	//计算函数时间
	timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;
	printf("gesvd time use:%20ld us\n",timeuse);
    
    

	//将INFO,D,Q,P输出
	printf("the INFO: %d\n\n",*INFO);
    printf("the D:\n");
	for(i=0; i<minMN; i++)
	  printf("%lf ",*(S+i));
	free(A);
	free(S);
	free(VT);
	free(U);
	free(WORK);
	free(INFO);
	fclose(fp);
	return 0;
}

