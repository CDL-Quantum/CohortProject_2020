//* matrix and vector classes */
// include diagmat*vector and the complex version
#include "matvec.h"
#include <sys/types.h>
#include <sys/times.h>
matrix::matrix(int rows,int columns)
{
/* constructor for matrix */
  if (rows*columns<=0) exit(1); // test for small size
  Rows=rows;
  Columns=columns;
  TheMatrix = new double[rows*columns];
  b=new double *[rows];
  b[0]=TheMatrix;
  for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=0.;
  if (!TheMatrix) exit(1); // test for excessive memory
}
matrix::matrix(const matrix& m)
{
/*  copy constructor for matrix */
  Rows=m.Rows;
  Columns=m.Columns;
  TheMatrix = new double[Rows*Columns];
  if (!TheMatrix) exit(1); // test for excessive memory
  for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=m.TheMatrix[i];
}
matrix& matrix::operator=(const matrix& arg)
{
#ifdef DEBUG
  if (Rows != arg.Rows || Columns != arg.Columns)  {
    std::cerr<<"matrix sizes not matching"<<std::endl;
    exit(1);
  }
#endif
  Rows=arg.Rows;
  Columns=arg.Columns;
  for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=arg.TheMatrix[i];
  return *this;
}
matrix::~matrix()
{
  delete[] TheMatrix;
}
matrix operator+(const matrix& arg1, const matrix& arg2)
{
  matrix sum(arg1.Rows,arg2.Columns);
  int Rows=arg1.Rows;
  int Columns=arg2.Columns;
  for (int i=0;i<Rows*Columns;i++) 
    sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2.TheMatrix[i];
  if (!sum.TheMatrix) exit(1);
  return sum;
}
matrix operator-(const matrix& arg1, const matrix& arg2)
{
  matrix dif(arg1.Rows,arg2.Columns);
  int Rows=arg1.Rows;
  int Columns=arg2.Columns;
  for (int i=0;i<Rows*Columns;i++) 
    dif.TheMatrix[i]=arg1.TheMatrix[i]-arg2.TheMatrix[i];
  if (!dif.TheMatrix) exit(1);
  return dif;
}
const matrix operator*(const matrix& arg1, const matrix& arg2)
{
// matrix product
  matrix prod(arg1.Rows,arg2.Columns);
  int Rows=arg1.Rows;
  int Columns=arg2.Columns;
#ifdef BLAS
  char transa = 'N';
  char transb = 'N';
  int m=arg1.Rows;
  int n=arg2.Columns;
  int k=arg1.Columns;
  double alpha=1.;
  double beta=0.;
  int lda=arg1.Rows;
  int ldb=arg2.Rows;
  int ldc=prod.Rows;
  FORTRAN(dgemm)(&transa,&transb,&m,&n,&k,&alpha,arg1.TheMatrix,&lda,
	 arg2.TheMatrix,&ldb,&beta,prod.TheMatrix,&ldc);
  return prod;
#else
  for (int i=0;i<Rows;i++) {
	  for (int j=0;j<Columns;j++) {
			prod.TheMatrix[i+j*Rows]=0.;
			for (int k=0;k<arg1.Columns;k++)
				prod.TheMatrix[i+j*Rows]+=arg1.TheMatrix[i+k*Rows]
				*arg2.TheMatrix[k+j*arg2.Rows];
		}
	}
	if (!prod.TheMatrix) exit(1);
	return prod;
#endif
}
const matrix transpose(const matrix& arg1)
{
  int Rows=arg1.Rows;
  int Columns=arg1.Columns;
  matrix trans(Columns,Rows);
  for (int i=0;i<Rows;i++) for (int j=0;j<Columns;j++) 
    trans.TheMatrix[j+i*Columns]=arg1.TheMatrix[i+j*Rows];
  if (!trans.TheMatrix) exit(1);
  return trans;
}
const matrix operator*(const double arg1, const matrix& arg2)
{
// multiplication by a scalar
	matrix prod(arg2.Rows,arg2.Columns);
	int Rows=arg2.Rows;
    int Columns=arg2.Columns;
    for (int i=0;i<Rows;i++) for (int j=0;j<Columns;j++) 
        prod.TheMatrix[i+j*Rows]=arg1*arg2.TheMatrix[i+j*Rows];
	return prod;
}
const vector operator*(const double arg1, const vector& arg2)
{
// multiplication by a scalar
	vector prod(arg2.Sz);
	int Sz=arg2.Sz;
	for (int i=0;i<Sz;i++) prod.TheVector[i]=arg1*arg2.TheVector[i];
	return prod;
}
const vector operator*(const matrix& arg1, const vector& arg2)
{
// matrix vector product
	vector prod(arg1.Rows);
	int Rows=arg1.Rows;
	int Sz=arg2.Sz;
#ifdef BLAS
	char trans = 'N';
	int m=arg1.Rows;
	int n=arg1.Columns;
	double alpha=1.;
	double beta=0.;
	int lda=arg1.Rows;
	int incx=1;
	int incy=1;
	FORTRAN(dgemv)(&trans,&m,&n,&alpha,arg1.TheMatrix,&lda,
			arg2.TheVector,&incx,&beta,prod.TheVector,&incy);
	return prod;
#endif
	for (int i=0;i<Rows;i++) {
		prod.TheVector[i]=0.;
		for (int k=0;k<arg1.Columns;k++)
		prod.TheVector[i]+=arg1.TheMatrix[i+k*Rows]*arg2.TheVector[k];
	}
	return prod;
}
//vector operator*(const symmatrix& arg1, const vector& arg2)
//{
// matrix vector product
//	vector prod(arg2.Sz);
//	int Rows=arg1.Rows;
//	int Sz=arg2.Sz;
// #ifdef BLAS
// #endif
//	for (int i=0;i<Rows;i++) {
//		prod.TheVector[i]=0.;
//		for (int k=0;k<arg1.Columns;k++)
//		prod.TheVector[i]+=arg1.TheMatrix[i+k*Rows]*arg2.TheVector[k];
//	}
//	return prod;
//}

const matrix operator*(const matrix& arg1, const diagmat& arg2)
{
	int Rows=arg1.Rows;
	int Columns=arg1.Columns;
	matrix prod(Rows,Columns);
	int Sz=arg2.Sz;
	for (int i=0;i<Rows;i++)
		for (int j=0;j<Columns;j++)
		  prod.TheMatrix[i+j*Rows]=arg1.TheMatrix[i+j*Rows]*arg2.TheMatrix[j];
	return prod;
}
const cmatrix operator*(const cmatrix& arg1, const cdiagmat& arg2)
{
	int Rows=arg1.Rows;
	int Columns=arg1.Columns;
	cmatrix prod(Rows,Columns);
	int Sz=arg2.Sz;
	for (int i=0;i<Rows;i++)
		for (int j=0;j<Columns;j++)
		  prod.TheMatrix[i+j*Rows]=arg1.TheMatrix[i+j*Rows]*arg2.TheMatrix[j];
	return prod;
}
const matrix  operator*(const diagmat& arg1, const matrix& arg2)
{
	int Rows=arg2.Rows;
	int Columns=arg2.Columns;
	matrix prod(Rows,Columns);
	int Sz=arg1.Sz;
	for (int i=0;i<Rows;i++)
		for (int j=0;j<Columns;j++)
		  prod.TheMatrix[i+j*Rows]=arg2.TheMatrix[i+j*Rows]*arg1.TheMatrix[i];
	return prod;
}
const vector  operator*(const diagmat& arg1, const vector& arg2)
{
	int Sz=arg1.Sz;
	vector prod(Sz);
	//for (int i=0;i<Sz;i++)
	//  prod.TheVector[i]=arg1.TheMatrix[i]*arg2.TheVector[i];
	FORTRAN(dmatvecprod)(prod.TheVector,arg1.TheMatrix,arg2.TheVector,&Sz);
	return prod;
}
const cvector  operator*(const cdiagmat& arg1, const cvector& arg2)
{
	int Sz=arg1.Sz;
	cvector prod(Sz);
	for (int i=0;i<Sz;i++)
		  prod.TheVector[i]=arg1.TheMatrix[i]*arg2.TheVector[i];
	return prod;
}
const cmatrix  operator*(const cdiagmat& arg1, const cmatrix& arg2)
{
	int Rows=arg2.Rows;
	int Columns=arg2.Columns;
	cmatrix prod(Rows,Columns);
	int Sz=arg1.Sz;
	for (int i=0;i<Rows;i++)
		for (int j=0;j<Columns;j++)
		  prod.TheMatrix[i+j*Rows]=arg2.TheMatrix[i+j*Rows]*arg1.TheMatrix[i];
	return prod;
}
matrix  operator+(const matrix& arg1, const diagmat& arg2)
{
    matrix sum=arg1;
    int sz=arg2.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg2.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
cmatrix  operator+(const cmatrix& arg1, const cdiagmat& arg2)
{
    cmatrix sum=arg1;
    int sz=arg2.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg2.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
matrix  operator+(const diagmat& arg1, const matrix& arg2)
{
    matrix sum=arg2;
    int sz=arg1.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg1.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
cmatrix  operator+(const cdiagmat& arg1, const cmatrix& arg2)
{
    cmatrix sum=arg2;
    int sz=arg1.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg1.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
matrix  operator-(const matrix& arg1, const diagmat& arg2)
{
    matrix sum=arg1;
    int sz=arg2.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]-=arg2.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
cmatrix  operator-(const cmatrix& arg1, const cdiagmat& arg2)
{
    cmatrix sum=arg1;
    int sz=arg2.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]-=arg2.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
matrix  operator-(const diagmat& arg1, const matrix& arg2)
{
    matrix sum=(-1)*arg2;
    int sz=arg1.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg1.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
cmatrix  operator-(const cdiagmat& arg1, const cmatrix& arg2)
{
    cmatrix sum=complex(1.,0)*arg2;
    int sz=arg1.Sz;
    for (int i=0;i<sz;i++)
        sum.TheMatrix[i+i*sz]+=arg1.TheMatrix[i];
    if (!sum.TheMatrix) exit(1);
    return sum;
}
const diagmat operator*(const diagmat& arg1, const diagmat& arg2)
{
	int Sz=arg1.Sz;
	diagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1.TheMatrix[i]*arg2.TheMatrix[i];
	return prod;
}
const cdiagmat operator*(const cdiagmat& arg1, const cdiagmat& arg2)
{
	int Sz=arg1.Sz;
	cdiagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1.TheMatrix[i]*arg2.TheMatrix[i];
	return prod;
}
const diagmat operator*(const double arg1, const diagmat& arg2)
{
	int Sz=arg2.Sz;
	diagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1*arg2.TheMatrix[i];
	return prod;
}
const cdiagmat operator*(const complex arg1, const cdiagmat& arg2)
{
	int Sz=arg2.Sz;
	cdiagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1*arg2.TheMatrix[i];
	return prod;
}
const diagmat operator*(const diagmat& arg1, const double arg2)
{
	int Sz=arg1.Sz;
	diagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1.TheMatrix[i]*arg2;
	return prod;
}
const cdiagmat operator*(const cdiagmat& arg1, const complex arg2)
{
	int Sz=arg1.Sz;
	cdiagmat prod(Sz);
	for (int i=0;i<Sz;i++)
		prod.TheMatrix[i]=arg1.TheMatrix[i]*arg2;
	return prod;
}
diagmat operator+(const diagmat& arg1, const diagmat& arg2)
{
	int Sz=arg1.Sz;
	diagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2.TheMatrix[i];
	return sum;
}
cdiagmat operator+(const cdiagmat& arg1, const cdiagmat& arg2)
{
	int Sz=arg1.Sz;
	cdiagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2.TheMatrix[i];
	return sum;
}
diagmat operator+(const double arg1, const diagmat& arg2)
{
	int Sz=arg2.Sz;
	diagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1+arg2.TheMatrix[i];
	return sum;
}
cdiagmat operator+(const complex arg1, const cdiagmat& arg2)
{
	int Sz=arg2.Sz;
	cdiagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1+arg2.TheMatrix[i];
	return sum;
}
diagmat operator+(const diagmat& arg1, const double arg2)
{
	int Sz=arg1.Sz;
	diagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2;
	return sum;
}
cdiagmat operator+(const cdiagmat& arg1, const complex arg2)
{
	int Sz=arg1.Sz;
	cdiagmat sum(Sz);
	for (int i=0;i<Sz;i++)
		sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2;
	return sum;
}
diagmat operator-(const diagmat& arg1, const diagmat& arg2)
{
    int Sz=arg1.Sz;
    diagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1.TheMatrix[i]-arg2.TheMatrix[i];
    return sum;
}
cdiagmat operator-(const cdiagmat& arg1, const cdiagmat& arg2)
{
    int Sz=arg1.Sz;
    cdiagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1.TheMatrix[i]-arg2.TheMatrix[i];
    return sum;
}
diagmat operator-(const double arg1, const diagmat& arg2)
{
    int Sz=arg2.Sz;
    diagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1-arg2.TheMatrix[i];
    return sum;
}
cdiagmat operator-(const complex arg1, const cdiagmat& arg2)
{
    int Sz=arg2.Sz;
    cdiagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1-arg2.TheMatrix[i];
    return sum;
}
diagmat operator-(const diagmat& arg1, const double arg2)
{
    int Sz=arg1.Sz;
    diagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1.TheMatrix[i]-arg2;
    return sum;
}
cdiagmat operator-(const cdiagmat& arg1, const complex arg2)
{
    int Sz=arg1.Sz;
    cdiagmat sum(Sz);
    for (int i=0;i<Sz;i++)
        sum.TheMatrix[i]=arg1.TheMatrix[i]-arg2;
    return sum;
}
diagmat inverse(const diagmat& v)
{
	int Sz=v.Sz;
    diagmat invr(Sz);
    for (int i=0;i<Sz;i++)
        invr.TheMatrix[i]=1./v.TheMatrix[i];
    return invr;
}
matrix inverse(const matrix& a)
{
// for symmetric positive definite matrix
    	matrix inv=a;
	char UPLO='U';
	int INFO=0;
	int N=a.Rows;
	int LDA=N;
	FORTRAN(dpotrf)(&UPLO,&N,inv.TheMatrix,&LDA,&INFO);
	FORTRAN(dpotri)(&UPLO,&N,inv.TheMatrix,&LDA,&INFO);
	if (INFO != 0) {
		if (INFO < 0) {std::cerr<<"illegal value for argument "<<INFO;
				std::cerr<<"in dpotri (inverse)"<<std::endl;}
		if (INFO > 0) std::cerr<<"inverse could not be computed"<<std::endl;
	}
	// fill the lower triangle
	for (int i=0;i<N;i++) 
		for (int j=0;j<i;j++)
			inv(i,j)=inv(j,i);
	return inv;
}
double det(const matrix& a)
{
  int n=a.Rows;
  int IFAIL=0;
  matrix A=transpose(a)*a;
  /* vector WKSPCE(n);
     f03aaf_(A.TheMatrix,&n,&n,&DET,WKSPCE.TheVector,&IFAIL);
     std::cerr<<"det to be implemented"<<std::endl;
     exit(1); */
  vector v=diag(A);
  double DET=1.;
  for (int i=0;i<n;i++) DET*=v(i);
  DET=sqrt(DET);
  std::cout<<DET<<std::endl;
  return DET;
}
cdiagmat inverse(const cdiagmat& v)
{
	int Sz=v.Sz;
    cdiagmat invr(Sz);
    for (int i=0;i<Sz;i++)
        invr.TheMatrix[i]=complex(1.,0.)/v.TheMatrix[i];
    return invr;
}
diagmat sqrtdiagmat(const diagmat& v)
{
	int Sz=v.Sz;
    diagmat sq(Sz);
    for (int i=0;i<Sz;i++)
        sq.TheMatrix[i]=sqrt(v.TheMatrix[i]);
    return sq;
}
cdiagmat sqrtdiagmat(const cdiagmat& v)
{
	int Sz=v.Sz;
    cdiagmat sq(Sz);
    for (int i=0;i<Sz;i++)
        sq.TheMatrix[i]=sqrt(v.TheMatrix[i]);
    return sq;
}

const double operator*(const vector& arg1, const vector& arg2)
{
	int Sz=arg2.Sz;
	double dot=0.;
	//for (int i=0;i<Sz;i++) 
	//  dot+=arg1.TheVector[i]*arg2.TheVector[i];
	FORTRAN(dotprod)(arg1.TheVector,arg2.TheVector,&Sz,&dot);
	return dot;
}
vector operator+(const vector& arg1, const vector& arg2)
{
	vector sum(arg1.Sz);
	int Sz=arg1.Sz;
	for (int i=0;i<Sz;i++) 
		sum.TheVector[i]=arg1.TheVector[i]+arg2.TheVector[i];
	if (!sum.TheVector) exit(1);
	return sum;
}
vector operator+(const double& arg1, const vector& arg2)
{
	vector sum(arg2.Sz);
	int Sz=arg2.Sz;
	for (int i=0;i<Sz;i++)
        sum.TheVector[i]=arg1+arg2.TheVector[i];
    if (!sum.TheVector) exit(1);
    return sum;
}
cvector operator+(const complex& arg1, const cvector& arg2)
{
	cvector sum(arg2.Sz);
	int Sz=arg2.Sz;
	for (int i=0;i<Sz;i++)
        sum.TheVector[i]=arg1+arg2.TheVector[i];
    if (!sum.TheVector) exit(1);
    return sum;
}
vector operator+(const vector& arg1, const double& arg2)
{
    vector sum(arg1.Sz);
    int Sz=arg1.Sz;
    for (int i=0;i<Sz;i++)
        sum.TheVector[i]=arg2+arg1.TheVector[i];
    if (!sum.TheVector) exit(1);
    return sum;
}
cvector operator+(const cvector& arg1, const complex& arg2)
{
    cvector sum(arg1.Sz);
    int Sz=arg1.Sz;
    for (int i=0;i<Sz;i++)
        sum.TheVector[i]=arg2+arg1.TheVector[i];
    if (!sum.TheVector) exit(1);
    return sum;
}
vector operator-(const vector& arg1, const vector& arg2)
{
	vector diff(arg1.Sz);
	int Sz=arg1.Sz;
	for (int i=0;i<Sz;i++) 
		diff.TheVector[i]=arg1.TheVector[i]-arg2.TheVector[i];
	if (!diff.TheVector) exit(1);
	return diff;
}
vector::vector(int sz)
{
	if (sz<=0) exit(1); // test for small size
	Sz=sz;
	TheVector= new double[sz];
	if (!TheVector) exit(1); // test for excessive memory
	for (int i=0;i<Sz;i++) TheVector[i]=0.;
}
diagmat::diagmat(int sz)
{
	if (sz<=0) exit(1); // test for small size
	Sz=sz;
	TheMatrix= new double[sz];
	for (int i=0;i<Sz;i++) TheMatrix[i]=0.;
	if (!TheMatrix) exit(1); // test for excessive memory
}
cdiagmat::cdiagmat(int sz)
{
	if (sz<=0) exit(1); // test for small size
	Sz=sz;
	TheMatrix= new complex[sz];
	for (int i=0;i<Sz;i++) TheMatrix[i]=complex(0.,0.);
	if (!TheMatrix) exit(1); // test for excessive memory
}
vector::vector(const vector& v)
{
/* copy constructor for vector */
	Sz=v.Sz;
	TheVector= new double[Sz];
	if (!TheVector) exit(1); // test for excessive memory
	for (int i=0;i<Sz;i++) TheVector[i]=v.TheVector[i];
}
diagmat::diagmat(const diagmat& v)
{
/* copy constructor for diagmat */
	Sz=v.Sz;
	TheMatrix= new double[Sz];
	if (!TheMatrix) exit(1); // test for excessive memory
	for (int i=0;i<Sz;i++) TheMatrix[i]=v.TheMatrix[i];
}
cdiagmat::cdiagmat(const cdiagmat& v)
{
/* copy constructor for diagmat */
	Sz=v.Sz;
	TheMatrix= new complex[Sz];
	if (!TheMatrix) exit(1); // test for excessive memory
	for (int i=0;i<Sz;i++) TheMatrix[i]=v.TheMatrix[i];
}
vector& vector::operator=(const vector& arg)
{
#ifdef DEBUG
	if (Sz != arg.Sz)  {
		std::cerr<<"vector sizes not matching"<<std::endl;
		exit(1);
	}
#endif
    Sz=arg.Sz;
    for (int i=0;i<Sz;i++) TheVector[i]=arg.TheVector[i];
    return *this;
}
diagmat& diagmat::operator=(const diagmat& arg)
{
#ifdef DEBUG
	if (Sz != arg.Sz)  {
		std::cerr<<"diagmat sizes not matching"<<std::endl;
		exit(1);
	}
#endif
    Sz=arg.Sz;
    for (int i=0;i<Sz;i++) TheMatrix[i]=arg.TheMatrix[i];
    return *this;
}
cdiagmat& cdiagmat::operator=(const cdiagmat& arg)
{
#ifdef DEBUG
	if (Sz != arg.Sz)  {
		std::cerr<<"cdiagmat sizes not matching"<<std::endl;
		exit(1);
	}
#endif
    Sz=arg.Sz;
    for (int i=0;i<Sz;i++) TheMatrix[i]=arg.TheMatrix[i];
    return *this;
}
vector::~vector()
{
	delete[] TheVector;
}
diagmat::~diagmat()
{
	delete[] TheMatrix;
}
cdiagmat::~cdiagmat()
{
	delete[] TheMatrix;
}
vector inverse(const vector& v)
{
	int Sz=v.Sz;
	vector inv(Sz);
	for (int i=0;i<Sz;i++) inv.TheVector[i]=1./v.TheVector[i];
	if (!inv.TheVector) exit(1);
	return inv;
}
vector sqrtvec(const vector& v)
{
	int Sz=v.Sz;
	vector sqrtv(Sz);
	for (int i=0;i<Sz;i++) sqrtv.TheVector[i]=sqrt(v.TheVector[i]);
	if (!sqrtv.TheVector) exit(1);
	return sqrtv;
}
vector compdiag(const matrix& PRe,const matrix& PIm)
{
  int nsize=PIm.Rows;
  vector eval(nsize);
  vector WK1(nsize);
  vector WK2(nsize);
  vector WK3(nsize);
  int IFAIL=0;
  
  /*f02awf_(PRe.TheMatrix,&nsize,PIm.TheMatrix,&nsize,&nsize,eval.TheVector,
    WK1.TheVector,WK2.TheVector,WK3.TheVector,&IFAIL);*/
  std::cerr<<"compdiag to be implemented"<<std::endl;
  exit(1);
  
  
  if (!eval.TheVector) exit(1);
  return eval;
}
vector compgene(cmatrix& A,cmatrix& B)
{
    int i;
	int N=A.Rows;
	int IFAIL=0;
	int MATV=0;
	int*  ITER=new int[N];
	vector ALFR(N);
	vector ALFI(N);
	vector BETA(N);
	double EPS1=0.;
	matrix AR(N,N);
	matrix AI(N,N);
	matrix BR(N,N);
	matrix BI(N,N);
	for (i=0;i<N;i++) for (int j=0;j<N;j++) {
		AR(i,j)=real(A(i,j));
		AI(i,j)=imag(A(i,j));
		BR(i,j)=real(B(i,j));
		BI(i,j)=imag(B(i,j));
	}
	matrix VR(N,N);
	matrix VI(N,N);

	/* f02gjf_(&N,AR.TheMatrix,&N,AI.TheMatrix,&N,BR.TheMatrix,&N,
			BI.TheMatrix,&N,&EPS1,ALFR.TheVector,ALFI.TheVector,
			BETA.TheVector,&MATV,VR.TheMatrix,&N, 
           VI.TheMatrix,&N,ITER,&IFAIL); */
	std::cerr<<"compgene to be implemented"<<std::endl;
        exit(1);

	if (IFAIL != 0) std::cerr <<"error in compgene "<<std::endl;
	if (!ALFR.TheVector) exit(1);
	for (i=0;i<N;i++) ALFR(i)=ALFR(i)/BETA(i);
// sort !
	char ORDER='A';
	int M1=1;
	int M2=N;
	/* m01caf_(ALFR.TheVector,&M1,&M2,&ORDER,&IFAIL);*/
	 std::cerr<<"sort to be implemented"<<std::endl;
        exit(1);

	if (IFAIL != 0) std::cerr <<"error in m01caf_ "<<std::endl;
	return ALFR;
}
vector hermdiag(matrix& PRe,matrix& PIm)
{
    int nsize=PIm.Rows;
    vector eval(nsize);
    vector WK1(nsize);
    vector WK2(nsize);
    vector WK3(nsize);
	matrix TR(nsize,nsize);
	matrix TI(nsize,nsize);
    int IFAIL=0;
	/* f02axf_(PRe.TheMatrix,&nsize,PIm.TheMatrix,&nsize,&nsize,
		eval.TheVector,TR.TheMatrix,&nsize,TI.TheMatrix,&nsize,WK1.TheVector,
		WK2.TheVector, WK3.TheVector,&IFAIL);*/
	std::cerr<<"hermdiag to be implemented"<<std::endl;
        exit(1);

	PRe=TR;
	PIm=TI;
	if (IFAIL != 0) std::cerr <<"error in hermdiag "<<std::endl;
	if (!eval.TheVector) exit(1);
	return eval;
}
cmatrix complexm(const matrix &Re,const matrix &Im)
{
	cmatrix cm(Re.Rows,Re.Columns);
	for (int i=0;i<Re.Rows;i++) {
		for (int j=0;j<Re.Columns;j++) {
			cm(i,j)=complex(Re(i,j),Im(i,j));
		}
	}
	return cm;
}
cmatrix complexm(const matrix &Re)
{
	cmatrix cm(Re.Rows,Re.Columns);
	for (int i=0;i<Re.Rows;i++) {
		for (int j=0;j<Re.Columns;j++) {
			cm(i,j)=complex(Re(i,j),0.);
		}
	}
	return cm;
}
matrix realm(cmatrix &a)
{
	matrix m(a.Rows,a.Columns);
	for (int i=0;i<a.Rows;i++) {
		for (int j=0;j<a.Columns;j++) {
			m(i,j)=real(a(i,j));
		}
	}
	return m;
}
matrix imagm(cmatrix &a)
{
	matrix m(a.Rows,a.Columns);
	for (int i=0;i<a.Rows;i++) {
		for (int j=0;j<a.Columns;j++) {
			m(i,j)=imag(a(i,j));
		}
	}
	return m;
}
cdiagmat complexm(diagmat &Re,diagmat &Im)
{
	cdiagmat cm(Re.Sz);
	for (int i=0;i<Re.Sz;i++) cm(i)=complex(Re(i),Im(i));
	return cm;
}
cdiagmat complexm(diagmat &Re)
{
	cdiagmat cm(Re.Sz);
	for (int i=0;i<Re.Sz;i++) cm(i)=complex(Re(i));
	return cm;
}
diagmat realm(cdiagmat &a)
{
	diagmat m(a.Sz);
	for (int i=0;i<a.Sz;i++) m(i)=real(a(i));
	return m;
}
diagmat imagm(cdiagmat &a)
{
	diagmat m(a.Sz);
	for (int i=0;i<a.Sz;i++) m(i)=imag(a(i));
	return m;
}
vector hermdl(cmatrix& a)
{
	char UPLO='U';
	char JOBZ='V';
	int INFO=0;
	int N=a.Rows;
	int LDA=N;
	int LWORK=2*N;
	cvector WORK(LWORK);
	vector W(N);
	vector RWORK(3*N);
	FORTRAN(cheev)(&JOBZ,&UPLO,&N,a.TheMatrix,&LDA,W.TheVector,WORK.TheVector,
		&LWORK,RWORK.TheVector,&INFO);
	if (INFO != 0) std::cerr<<"diagonalization failed"<<std::endl;
	return W;
}
// construct identity matrix
vector diaggen(matrix& a)
{
	int N=a.Rows;
	char JOBVL='V';
	char JOBVR='V';
	int INFO=0;
	int LDA=N;
	vector WR(N);
	vector WI(N);
	int LDVL=N;
	int LDVR=N;
	matrix VL(LDVL,N);
	matrix VR(LDVR,N);
	int LWORK=4*N;
	vector WORK(LWORK);
	vector W(N);
	FORTRAN(dgeev)(&JOBVL,&JOBVR,&N,a.TheMatrix,&LDA,WR.TheVector, 
	WI.TheVector,VL.TheMatrix,&LDVL,VR.TheMatrix,&LDVR,WORK.TheVector, 
	&LWORK,&INFO);
	if (INFO != 0) std::cerr<<"diagonalization failed"<<std::endl; 
	return WR; 
}
vector diag(matrix& a)
  // returns eigenvalues in a vector and transforms the matrix argument
  // for symmetric matrices
{
  int N=a.Rows;
  char UPLO='U';
  char JOBZ='V';
  int INFO=0;
  int LDA=N;
  int LWORK=3*N;
  vector WORK(LWORK);
  vector W(N);
  FORTRAN(dsyev)(&JOBZ,&UPLO,&N,a.TheMatrix,&LDA,W.TheVector,WORK.TheVector,
	 &LWORK,&INFO);
  if (INFO != 0) std::cerr<<"diagonalization failed"<<std::endl;
  return W;
}
cvector diag(cmatrix& a,cmatrix &VL,cmatrix &VR)
  // returns eigenvalues in a vector and transforms the matrix argument
  // for symmetric matrices
{
  int N=a.Rows;
  char JOBVL='V';
  char JOBVR='V';
  int INFO=0;
  int LDA=N;
  int LDVL=N;
  int LDVR=N;
  int LWORK=2*N;
  cvector WORK(LWORK);
  cvector W(N);
  vector RWORK(2*N);

  FORTRAN(zgeev)(&JOBVL,&JOBVR,&N,a.TheMatrix,&LDA,W.TheVector,VL.TheMatrix,&LDVL,
        VR.TheMatrix,&LDVR,WORK.TheVector,&LWORK,
         RWORK.TheVector,&INFO);
  if (INFO != 0) std::cerr<<"diagonalization failed"<<std::endl;
  return W;
}
cmatrix inverse(const cmatrix &A)
{
  int N=A.Rows;
  int INFO=0;
  int LDA=N;
  int LDB=N;
  int NRHS=N;
  int*  IPIV=new int[N];
  cmatrix B(N,N);
  cmatrix AP=A;
  for (int i=0;i<N;i++) B(i,i)=complex(1.,0.);
  FORTRAN(zgesv)(&N,&NRHS,AP.TheMatrix,&LDA,IPIV,B.TheMatrix,&LDB,&INFO);
  return B;
}
// implement plain ql
matrix ql(vector &D,vector &E)
{
  int N=D.Sz;
  char COMPZ='I';
  matrix Z(N,N); for (int i=0;i<N;i++) Z(i,i)=1.;
  int INFO=0;
  int LDZ=N;
  vector WORK(2*N-2);
  FORTRAN(dsteqr)(&COMPZ,&N,D.TheVector,E.TheVector,Z.TheMatrix,&LDZ,
	  WORK.TheVector,&INFO);
  return Z;
}
// complex stuff
cmatrix::cmatrix(int rows,int columns)
{
/* constructor for cmatrix */
	if (rows*columns<=0) exit(1); // test for small size
	Rows=rows;
	Columns=columns;
	TheMatrix = new complex[rows*columns];
	for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=complex(0.,0.);
	if (!TheMatrix) exit(1); // test for excessive memory
}
cmatrix::cmatrix(const cmatrix& m)
{
/*  copy constructor for cmatrix */
	Rows=m.Rows;
	Columns=m.Columns;
	TheMatrix = new complex[Rows*Columns];
	if (!TheMatrix) exit(1); // test for excessive memory
	for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=m.TheMatrix[i];
}
cmatrix& cmatrix::operator=(const cmatrix& arg)
{
#ifdef DEBUG
	if (Rows != arg.Rows || Columns != arg.Columns)  {
		std::cerr<<"matrix sizes not matching"<<std::endl;
		exit(1);
	}
#endif
	Rows=arg.Rows;
	Columns=arg.Columns;
	for (int i=0;i<Rows*Columns;i++) TheMatrix[i]=arg.TheMatrix[i];
	return *this;
}
cmatrix::~cmatrix()
{
	delete[] TheMatrix;
}
complex& cmatrix::operator() (int rindex,int cindex) const
{
#ifdef DEBUG
    if (rindex < 0 || rindex >=Rows || cindex< 0 || cindex >=Columns) {
        std::cerr<<"Indices exceed cmatrix size:\n"<<"rindex="<<rindex<<",cindex="<<cindex<<"Actual cmatrix size is: "<<Rows<<" by "<<Columns<<std::endl;
        exit(1);
    }
#endif
    return TheMatrix[cindex*Rows+rindex];
}
cmatrix operator+(const cmatrix& arg1, const cmatrix& arg2)
{
	cmatrix sum(arg1.Rows,arg2.Columns);
	int Rows=arg1.Rows;
	int Columns=arg2.Columns;
	for (int i=0;i<Rows*Columns;i++) 
		sum.TheMatrix[i]=arg1.TheMatrix[i]+arg2.TheMatrix[i];
	if (!sum.TheMatrix) exit(1);
	return sum;
}
cmatrix operator-(const cmatrix& arg1, const cmatrix& arg2)
{
	cmatrix dif(arg1.Rows,arg2.Columns);
	int Rows=arg1.Rows;
	int Columns=arg2.Columns;
	for (int i=0;i<Rows*Columns;i++) 
		dif.TheMatrix[i]=arg1.TheMatrix[i]-arg2.TheMatrix[i];
	if (!dif.TheMatrix) exit(1);
	return dif;
}
const cmatrix operator*(const cmatrix& arg1, const cmatrix& arg2)
{
// cmatrix product
	cmatrix prod(arg1.Rows,arg2.Columns);
	int Rows=arg1.Rows;
	int Columns=arg2.Columns;
#ifdef BLAS
	char transa = 'N';
	char transb = 'N';
	int m=arg1.Rows;
	int n=arg2.Columns;
	int k=arg1.Columns;
	complex alpha=complex(1.,0.);
	complex beta=complex(0.,0.);
	int lda=arg1.Rows;
	int ldb=arg2.Rows;
	int ldc=prod.Rows;
	FORTRAN(zgemm)(&transa,&transb,&m,&n,&k,&alpha,arg1.TheMatrix,&lda,
			arg2.TheMatrix,&ldb,&beta,prod.TheMatrix,&ldc);
	return prod;
#endif
	for (int i=0;i<Rows;i++) {
		for (int j=0;j<Columns;j++) {
			prod.TheMatrix[i+j*Rows]=complex(0.,0.);
			for (int k=0;k<arg1.Columns;k++)
				prod.TheMatrix[i+j*Rows]+=arg1.TheMatrix[i+k*Rows]
				*arg2.TheMatrix[k+j*arg2.Rows];
		}
	}
	if (!prod.TheMatrix) exit(1);
	return prod;
}
cmatrix transpose(const cmatrix& arg1)
{
  int Rows=arg1.Rows;
  int Columns=arg1.Columns;
  cmatrix trans(Columns,Rows);
  for (int i=0;i<Rows;i++) for (int j=0;j<Columns;j++) 
    trans.TheMatrix[j+i*Columns]=conj(arg1.TheMatrix[i+j*Rows]);
  if (!trans.TheMatrix) exit(1);
  return trans;
}
cdiagmat transpose(const cdiagmat& arg1)
{
  int Rows=arg1.Sz;
  cdiagmat trans(Rows);
  for (int i=0;i<Rows;i++) trans.TheMatrix[i]=conj(arg1.TheMatrix[i]);
  if (!trans.TheMatrix) exit(1);
  return trans;
}
const cmatrix operator*(const complex arg1, const cmatrix& arg2)
{
// multiplication by a scalar
	cmatrix prod(arg2.Rows,arg2.Columns);
	int Rows=arg2.Rows;
    int Columns=arg2.Columns;
    for (int i=0;i<Rows;i++) for (int j=0;j<Columns;j++) 
        prod.TheMatrix[i+j*Rows]=arg1*arg2.TheMatrix[i+j*Rows];
	if (!prod.TheMatrix) exit(1);
	return prod;
}
const cvector operator*(const complex arg1, const cvector& arg2)
{
// multiplication by a scalar
	cvector prod(arg2.Sz);
	int Sz=arg2.Sz;
	for (int i=0;i<Sz;i++) prod.TheVector[i]=arg1*arg2.TheVector[i];
	if (!prod.TheVector) exit(1);
	return prod;
}
const cvector operator*(const cmatrix& arg1, const cvector& arg2)
{
// cmatrix cvector product
	cvector prod(arg1.Rows);
	int Rows=arg1.Rows;
	int Sz=arg2.Sz;
#ifdef BLAS
	char trans = 'N';
	int m=arg1.Rows;
	int n=arg1.Columns;
	complex alpha=complex(1.,0.);
	complex beta=complex(0.,0.);
	int lda=arg1.Rows;
	int incx=1;
	int incy=1;
	FORTRAN(zgemv)(&trans,&m,&n,&alpha,arg1.TheMatrix,&lda,
			arg2.TheVector,&incx,&beta,prod.TheVector,&incy);
	return prod;
#else
	for (int i=0;i<Rows;i++) {
		prod.TheVector[i]=complex(0.,0.);
		for (int k=0;k<arg1.Columns;k++)
		prod.TheVector[i]+=arg1.TheMatrix[i+k*Rows]*arg2.TheVector[k];
	}
	if (!prod.TheVector) exit(1);
	return prod;
#endif
}
cvector zmatvec(const cmatrix& arg1, const cvector& arg2,complex scalar)
{
// cmatrix cvector product scaled by alpha
	cvector prod(arg1.Rows);
	int Rows=arg1.Rows;
	int Sz=arg2.Sz;
	char trans = 'N';
	int m=arg1.Rows;
	int n=arg1.Columns;
	complex alpha=scalar;
	complex beta=complex(0.,0.);
	int lda=arg1.Rows;
	int incx=1;
	int incy=1;
	FORTRAN(zgemv)(&trans,&m,&n,&alpha,arg1.TheMatrix,&lda,
			arg2.TheVector,&incx,&beta,prod.TheVector,&incy);
	return prod;
}
const complex operator*(const cvector& arg1, const cvector& arg2)
{
	int Sz=arg2.Sz;
	complex dot=complex(0.,0.);
	for (int i=0;i<Sz;i++) 
	  dot+=conj(arg1.TheVector[i])*arg2.TheVector[i];
	//cdotprod_(arg1.TheVector,arg2.TheVector,&Sz,&dot);
	return dot;
}
cvector cconj(const cvector& arg1)
{
	int Sz=arg1.Sz;
	cvector cv(Sz);
	for (int i=0;i<Sz;i++) cv(i)=conj(arg1.TheVector[i]);
	return cv;
}
vector realpart(const cvector& arg1)
{
	int Sz=arg1.Sz;
	vector rp(Sz);
	for (int i=0;i<Sz;i++) rp(i)=real(arg1.TheVector[i]);
	return rp;
}
vector imagpart(const cvector& arg1)
{
	int Sz=arg1.Sz;
	vector rp(Sz);
	for (int i=0;i<Sz;i++) rp(i)=imag(arg1.TheVector[i]);
	return rp;
}
cvector operator+(const cvector& arg1, const cvector& arg2)
{
	cvector sum(arg1.Sz);
	int Sz=arg1.Sz;
	for (int i=0;i<Sz;i++) 
		sum.TheVector[i]=arg1.TheVector[i]+arg2.TheVector[i];
	if (!sum.TheVector) exit(1);
	return sum;
}
cvector operator-(const cvector& arg1, const cvector& arg2)
{
	cvector diff(arg1.Sz);
	int Sz=arg1.Sz;
	for (int i=0;i<Sz;i++) 
		diff.TheVector[i]=arg1.TheVector[i]-arg2.TheVector[i];
	if (!diff.TheVector) exit(1);
	return diff;
}

cvector::cvector(int sz)
{
	if (sz<=0) exit(1); // test for small size
	Sz=sz;
	TheVector= new complex[sz];
	for (int i=0;i<Sz;i++) TheVector[i]=complex(0.,0.);
	if (!TheVector) exit(1); // test for excessive memory
}
cvector::cvector(const cvector& v)
{
/* copy constructor for cvector */
	Sz=v.Sz;
	TheVector= new complex[Sz];
	if (!TheVector) exit(1); // test for excessive memory
	for (int i=0;i<Sz;i++) TheVector[i]=v.TheVector[i];
}
cvector& cvector::operator=(const cvector& arg)
{
#ifdef DEBUG
	if (Sz != arg.Sz)  {
		std::cerr<<"vector sizes not matching"<<std::endl;
		exit(1);
	}
#endif
    Sz=arg.Sz;
    for (int i=0;i<Sz;i++) TheVector[i]=arg.TheVector[i];
    return *this;
}
cvector::~cvector()
{
	delete[] TheVector;
}
void cnormalise(cvector& arg1)
{
  int i;
  int Sz=arg1.Sz;
  complex dot=complex(0.,0.);
  //for (int i=0;i<Sz;i++) 
  //	dot+=conj(arg1.TheVector[i])*arg1.TheVector[i];
  dot=arg1*arg1;
  for (i=0;i<Sz;i++) arg1.TheVector[i]/=sqrt(dot);
}
void normalise(vector& arg1)
{
  int i;
  int Sz=arg1.Sz;
  double dot=0.;
  //for (int i=0;i<Sz;i++) 
  //	dot+=arg1.TheVector[i]*arg1.TheVector[i];
  dot=arg1*arg1;
  for (i=0;i<Sz;i++) arg1.TheVector[i]/=sqrt(dot);
}
matrix HCSCE(matrix &a,matrix &B)
{
  matrix A=a;
  char UPLO='U';
  int INFO=0;
  int ITYPE=1;
  int LDA=A.Rows;
  int LDB=LDA;
  int N=LDA;
  FORTRAN(spotrf)(&UPLO,&N,B.TheMatrix,&LDA,&INFO);
  std::cout<<INFO<<std::endl;
  FORTRAN(ssygst)(&ITYPE,&UPLO,&N,A.TheMatrix,&LDA,B.TheMatrix,&LDB,&INFO);
  std::cout<<INFO<<std::endl;
  // fill the lower triangle
  for (int i=0;i<N;i++) 
    for (int j=0;j<i;j++)
      A(i,j)=A(j,i);
  return A;
}

void DumpVector(const diagmat &pot4d,int nsize4d,FILE* filename)
{
  fwrite(pot4d.TheMatrix,sizeof(pot4d.TheMatrix),nsize4d,filename);
  return;
}
void ReadVector(diagmat &pot4d,int nsize4d,FILE* filename)
{
  fread(pot4d.TheMatrix,sizeof(pot4d.TheMatrix),nsize4d,filename);
  return;
}
void DumpVector(const vector &v,int n,FILE* filename)
{
  fwrite(v.TheVector,sizeof(v.TheVector),n,filename);
  return;
}
void ReadVector(vector &v,int n,FILE* filename)
{
  fread(v.TheVector,sizeof(v.TheVector),n,filename);
  return;
}
