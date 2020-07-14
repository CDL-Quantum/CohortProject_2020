using namespace std; 
#include <stddef.h>
#include <stdlib.h>
#include <new>
//#include <new.h>
#include <complex>
//#include <complex.h>
#include <iostream>
//#include <iostream.h>
#include <fstream>
//#include <fstream.h>
#include <strstream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define EXTERN extern "C"
#if HP
#define FORTRAN(x) (x)
#define complex complex<double>
#elif AIX
#define FORTRAN(x) (x)
#define EXTERNC extern "C"
#elif KCC
#define FORTRAN(x) (x)
#define complex complex<double>
#define EXTERNC 
#elif SGI 
#define FORTRAN(x) (x ## _)
#elif LINUX
#define FORTRAN(x) (x ## _)
#define EXTERNC extern "C"
#define complex std::complex<double>
#elif SUN
#define FORTRAN(x) (x ## _)
#else
#error fortran naming convention must be defined!
#endif
class cmatrix;
class matrix;
class cvector;
class diagmat;
class cdiagmat;
// class symmatrix;
// class hermmatrix;
class vector {
  int Sz;
	double* TheVector;
 public:
	vector(int);
	vector(const vector& v);
	~vector();
	double& operator() (int Index) const;
	vector& operator=(const vector& arg);

/* operations to code (friends) */
	friend void DumpVector(const vector &v,int n,FILE* filename);
	friend void ReadVector(vector &v,int n,FILE* filename);
	friend const vector operator*(const matrix&, const vector&);
	friend vector operator+(const vector& arg1, const vector& arg2);
	friend vector operator-(const vector& arg1, const vector& arg2);
	friend const vector operator*(const double arg1, const vector& arg2);
	friend vector inverse(const vector& v);
	friend matrix inverse(const matrix& a);
	friend double det(const matrix& a);
	friend vector sqrtvec(const vector& v);
	friend vector compdiag(const matrix& PRe,const matrix& PIm);
	friend vector hermdiag(matrix& PRe,matrix& PIm);
	friend const double operator*(const vector& arg1, const vector& arg2);
//	friend vector operator*(const symmatrix& arg1, const vector& arg2);
	friend const vector  operator*(const diagmat& arg1, const vector& arg2);
	friend vector diag(matrix& a);
	friend vector Hpsifort(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
			       diagmat &I1,diagmat &I2,matrix &Tphi,
			       diagmat &pot4d,vector &v,int Rpoints, 
			       int nsize1,int nsize2);
	friend cvector Hpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,diagmat &pot3d,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort3dfast(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		       diagmat &I2,matrix &Tphi,diagmat &pot3d,
		       cvector &v,int nsize1,int nsize2,vector &vRe,
		       vector &vIm,vector &uRe,vector &wRe,vector &uIm,
		       vector &wIm,cvector &u);

	friend int arpackinterface(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
				   diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,
				   int niter,double emin,
				   double emax,int maxnev,double tolerance,int syml,
				   vector &v,vector &d);
	friend void arpackinterface4d(matrix &tmatR,matrix &Ttheta1 ,matrix &Ttheta2,
		       matrix &Ttheta1R,matrix &Ttheta2R,diagmat &Inertia,
		       matrix &Tphi,matrix &sinphippluscos,diagmat &cosphi,
		       matrix &cosppmsinmcos,matrix &cotmat1,matrix &cotmat2,
		       matrix &delcot1,matrix &delcot2,
		       diagmat &InertiaMR, diagmat &InertiaMR2,
		       diagmat &pot4d,vector &initialv,int Rpoints,
		       int nsize1,int nsize2,int niter,double emin,
		       double emax,int nenergies,double tolerance);
	friend void arpackinterface4dPI(matrix &tmatR,diagmat &Ir0R,
					diagmat &IR,
					matrix &Hrot,diagmat &inertia,
					matrix &Tphi,
					matrix &delcot,matrix &cotmat,
					diagmat &cosphi,
					matrix &phiA,
					matrix &phiB,diagmat &pot4dsym,
					vector &initialv,
					int Rpoints,int nsize1,int nsize2,
					int niter,
					double emin,double emax,int maxnev,
					double tolerance);
	friend vector Hpsifort3dfastreal(matrix &Ttheta1,
					 matrix &Ttheta2,diagmat &I1,
					 diagmat &I2,matrix &Tphi,
					 diagmat &pot3d,
					 vector &v,int nsize1,int nsize2,
					 vector &u,vector &w);
	friend cvector Tpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort4d(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
				  diagmat &I1,
				  diagmat &I2,matrix &Tphi,diagmat &pot4d,
				  cvector &v,int nsize1,int nsize2,int sizeR);
	friend matrix ql(vector &D,vector &E);
	friend vector diaggen(matrix& a);
	friend vector operator+(const double& arg1, const vector& arg2);
	friend vector operator+(const vector& arg1, const double& arg2);
	friend void normalise(vector& arg1);
	friend vector hermdl(cmatrix& a);
	friend vector compgene(cmatrix& A,cmatrix& B);
	friend vector realpart(const cvector& arg1);
	friend vector imagpart(const cvector& arg1);
	friend void FT(cvector &v,int flag);
	friend vector Hpsifort4dreal(matrix &TR,matrix &Ttheta1,
				      matrix &Ttheta2,diagmat &I1,
				      diagmat &I2,matrix &Tphi,
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR,vector &uRe,vector &wRe);
	friend vector Hpsifort4drealE(matrix &TR,
				      matrix &Ttheta1,matrix &Ttheta2,
				      matrix &Ttheta1R,matrix &Ttheta2R,
				      diagmat &Inertia,
				      matrix &Tphi,matrix &cosphipp,
				      diagmat &cosphi,matrix &sinphip,
				      matrix &cotmat1,matrix &cotmat2,
				      matrix &delcot1,
				      matrix &delcot2,diagmat &InertiaMR,
				      diagmat &InertiaMR2,
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR);
	friend void LanczosPI(matrix &tmatR,diagmat &IroR,diagmat &IR,
			      matrix &Hrot,diagmat &inertia,
	       matrix &Tphi,matrix &delcot,matrix &cotmat,diagmat &cosphi,
	       matrix &phiA,matrix &phiB,diagmat &pot4dsym,vector &v,
	       int Rpoints,int nsize1,int nsize2,int niter,double emin,
	       double emax);
	friend cvector diag(cmatrix& a,cmatrix &VL,cmatrix &VR);
	friend cmatrix inverse(const cmatrix &a);
};
class cmatrix;
class vector;
class cvector {
	int Sz;
	complex* TheVector;
public:
	cvector(int);
	cvector(const cvector& v);
	~cvector();
	complex& operator() (int Index) const;
	cvector& operator=(const cvector& arg);

	friend const complex operator*(const cvector& arg1, const cvector& arg2);
	friend cvector operator+(const cvector& arg1, const cvector& arg2);
	friend cvector operator+(const complex& arg1, const cvector& arg2);
	friend cvector operator+(const cvector& arg1, const complex& arg2);
	friend cvector operator-(const cvector& arg1, const cvector& arg2);
	friend const cvector operator*(const complex arg1, const cvector& arg2);
	friend const cvector operator*(const cmatrix&, const cvector&);
//	friend cvector operator*(const hermmatrix& arg1, const cvector& arg2);
	friend const cvector operator*(const cdiagmat& arg1, const cvector& arg2);
	friend void cnormalise(cvector& arg1);
	friend vector hermdl(cmatrix& a);	
	friend cvector cconj(const cvector& arg1);
	friend vector realpart(const cvector& arg1);
	friend vector imagpart(const cvector& arg1);
	friend cvector zmatvec(const cmatrix& arg1, const cvector& arg2,complex scalar);
	friend cmatrix ql(cvector &D,cvector &E);
	friend  void FT(cvector &v,int flag);
	friend cvector Hpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,diagmat &pot3d,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort3dfast(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		       diagmat &I2,matrix &Tphi,diagmat &pot3d,
		       cvector &v,int nsize1,int nsize2,vector &vRe,
		       vector &vIm,vector &uRe,vector &wRe,vector &uIm,
		       vector &wIm,cvector &u);
	friend vector Hpsifort3dfastreal(matrix &Ttheta1,
					 matrix &Ttheta2,diagmat &I1,
					 diagmat &I2,matrix &Tphi,
					 diagmat &pot3d,
					 vector &v,int nsize1,int nsize2,
					 vector &u,vector &w);

	friend cvector Tpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort4d(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
				  diagmat &I1,
				  diagmat &I2,matrix &Tphi,diagmat &pot4d,
				  cvector &v,int nsize1,int nsize2,int sizeR);
	friend cvector diag(cmatrix& a,cmatrix &VL,cmatrix &VR);
	friend cmatrix inverse(const cmatrix &a);

};
class matrix {
	int Rows;
	int Columns;
	double* TheMatrix;
	double **b;
public:
	matrix(int Rows,int Columns);
 	matrix(const matrix& m); 
	~matrix();
	//double& operator() (int rindex,int cindex) const;
	matrix& operator=(const matrix& arg);
/* indexing */
	double& operator() (int rindex,int cindex) const;
/* operations to code (friends) */
	friend matrix operator+(const matrix&, const matrix&);
	friend matrix operator-(const matrix& arg1, const matrix& arg2);
	friend const matrix operator*(const matrix&, const matrix&);
	friend const vector operator*(const matrix&, const vector&); 
	friend matrix HCSCE(matrix &a,matrix &B);
	friend const matrix operator*(const double arg1, const matrix& arg2);
	friend const matrix transpose(const matrix& arg1);
	friend vector compdiag(const matrix& PRe,const matrix& PIm);
	friend vector hermdiag(matrix& PRe,matrix& PIm);
	friend vector diag(matrix& a);
	friend vector Hpsifort(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
			       diagmat &I1,diagmat &I2,matrix &Tphi,
			       diagmat &pot4d,vector &v,int Rpoints, 
			       int nsize1,int nsize2);
	friend cvector Hpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,diagmat &pot3d,
				  cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort3dfast(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		       diagmat &I2,matrix &Tphi,diagmat &pot3d,
		       cvector &v,int nsize1,int nsize2,vector &vRe,
		       vector &vIm,vector &uRe,vector &wRe,vector &uIm,
		       vector &wIm,cvector &u);
	friend vector Hpsifort3dfastreal(matrix &Ttheta1,
					 matrix &Ttheta2,diagmat &I1,
					 diagmat &I2,matrix &Tphi,
					 diagmat &pot3d,
					 vector &v,int nsize1,int nsize2,
					 vector &u,vector &w);
	
	friend int arpackinterface(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		       diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,
		       int niter,double emin,
				    double emax,int maxnev,double tolerance,int syml,
				    vector &v,vector &d);
	friend void arpackinterface4d(matrix &tmatR,matrix &Ttheta1 ,matrix &Ttheta2,
		       matrix &Ttheta1R,matrix &Ttheta2R,diagmat &Inertia,
		       matrix &Tphi,matrix &sinphippluscos,diagmat &cosphi,
		       matrix &cosppmsinmcos,matrix &cotmat1,matrix &cotmat2,
		       matrix &delcot1,matrix &delcot2,
		       diagmat &InertiaMR, diagmat &InertiaMR2,
		       diagmat &pot4d,vector &initialv,int Rpoints,
		       int nsize1,int nsize2,int niter,double emin,
		       double emax,int nenergies,double tolerance);
	friend void arpackinterface4dPI(matrix &tmatR,diagmat &Ir0R,
					diagmat &IR,
					matrix &Hrot,diagmat &inertia,
					matrix &Tphi,
					matrix &delcot,matrix &cotmat,
					diagmat &cosphi,
					matrix &phiA,
					matrix &phiB,diagmat &pot4dsym,
					vector &initialv,
					int Rpoints,int nsize1,int nsize2,
					int niter,
					double emin,double emax,int maxnev,
					double tolerance);
	friend cvector Tpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort4d(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
				  diagmat &I1,
				  diagmat &I2,matrix &Tphi,diagmat &pot4d,
				  cvector &v,int nsize1,int nsize2,int sizeR);
	friend matrix ql(vector &D,vector &E);
	friend vector diaggen(matrix& a);
	friend const matrix  operator*(const matrix&, const diagmat&);
	friend const matrix  operator*(const diagmat&, const matrix&);
	friend matrix  operator+(const matrix& arg1, const diagmat& arg2);
	friend matrix  operator+(const diagmat& arg1, const matrix& arg2);
	friend matrix  operator-(const matrix& arg1, const diagmat& arg2);
	friend matrix  operator-(const diagmat& arg1, const matrix& arg2);
	friend vector compgene(cmatrix& A,cmatrix& B);
	friend matrix inverse(const matrix& a);
	friend double det(const matrix& a);
	friend cmatrix complexm(const matrix &Re,const matrix &Im);
	friend cmatrix complexm(const matrix &Re);
	friend matrix realm(cmatrix &a);
	friend matrix imagm(cmatrix &a);
	friend vector Hpsifort4dreal(matrix &TR,matrix &Ttheta1,
				      matrix &Ttheta2,diagmat &I1,
				      diagmat &I2,matrix &Tphi,
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR,vector &uRe,vector &wRe);
	
	friend vector Hpsifort4drealE(matrix &TR,
				      matrix &Ttheta1,matrix &Ttheta2,
				      matrix &Ttheta1R,matrix &Ttheta2R,
				      diagmat &Inertia,
				      matrix &Tphi,matrix &cosphipp,
				      diagmat &cosphi,matrix &sinphip,
				      matrix &cotmat1,matrix &cotmat2,
				      matrix &delcot1,
				      matrix &delcot2,diagmat &InertiaMR,
				      diagmat &InertiaMR2,
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR);
	friend void LanczosPI(matrix &tmatR,diagmat &IroR,diagmat &IR,
			      matrix &Hrot,diagmat &inertia,
	       matrix &Tphi,matrix &delcot,matrix &cotmat,diagmat &cosphi,
	       matrix &phiA,matrix &phiB,diagmat &pot4dsym,vector &v,
	       int Rpoints,int nsize1,int nsize2,int niter,double emin,
	       double emax);
};
class diagmat {
	int Sz;
	double* TheMatrix;
public:
	diagmat(int);
	diagmat(const diagmat& v);
	~diagmat();
	double& operator() (int Index) const;
	diagmat& operator=(const diagmat& arg);

/* operations to code (friends) */
	friend void DumpVector(const diagmat &pot4d,int nsize4d,FILE* filename);
	friend void ReadVector(diagmat &pot4d,int nsize4d,FILE* filename);
	friend const matrix  operator*(const matrix&, const diagmat&);
	friend const matrix  operator*(const diagmat&, const matrix&);
	friend const diagmat operator*(const diagmat& arg1, const diagmat& arg2);
	friend const diagmat operator*(const double arg1, const diagmat& arg2);
	friend const diagmat operator*(const diagmat& arg1, const double arg2);
	friend const vector  operator*(const diagmat& arg1, const vector& arg2);
	friend diagmat operator+(const diagmat& arg1, const diagmat& arg2);
	friend diagmat operator+(const double arg1, const diagmat& arg2);
	friend diagmat operator+(const diagmat& arg1, const double arg2);
	friend matrix  operator+(const matrix& arg1, const diagmat& arg2);
	friend matrix  operator+(const diagmat& arg1, const matrix& arg2);
	friend matrix  operator-(const matrix& arg1, const diagmat& arg2);
	friend matrix  operator-(const diagmat& arg1, const matrix& arg2);
	friend diagmat operator-(const diagmat& arg1, const diagmat& arg2);
	friend diagmat operator-(const double arg1, const diagmat& arg2);
	friend diagmat operator-(const diagmat& arg1, const double arg2);
	friend diagmat inverse(const diagmat& v);
	friend vector Hpsifort(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
			       diagmat &I1,diagmat &I2,matrix &Tphi,
			       diagmat &pot4d,vector &v,int Rpoints, 
			       int nsize1,int nsize2);
	friend cvector Hpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,diagmat &pot3d,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort3dfast(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
		       diagmat &I2,matrix &Tphi,diagmat &pot3d,
		       cvector &v,int nsize1,int nsize2,vector &vRe,
		       vector &vIm,vector &uRe,vector &wRe,vector &uIm,
		       vector &wIm,cvector &u);
	friend vector Hpsifort3dfastreal(matrix &Ttheta1,
					 matrix &Ttheta2,diagmat &I1,
					 diagmat &I2,matrix &Tphi,
					 diagmat &pot3d,
					 vector &v,int nsize1,int nsize2,
					 vector &u,vector &w);
	
	friend int  arpackinterface(matrix &ddr,diagmat &G11,diagmat &G22,diagmat &G33,diagmat &G12,
		       diagmat &G13,diagmat &G23,diagmat &VaplusV,int size,
		       int niter,double emin,
				    double emax,int maxnev,double tolerance,int syml,
				    vector &v,vector &d);
	friend void arpackinterface4d(matrix &tmatR,matrix &Ttheta1 ,matrix &Ttheta2,
		       matrix &Ttheta1R,matrix &Ttheta2R,diagmat &Inertia,
		       matrix &Tphi,matrix &sinphippluscos,diagmat &cosphi,
		       matrix &cosppmsinmcos,matrix &cotmat1,matrix &cotmat2,
		       matrix &delcot1,matrix &delcot2,
		       diagmat &InertiaMR, diagmat &InertiaMR2,
		       diagmat &pot4d,vector &initialv,int Rpoints,
		       int nsize1,int nsize2,int niter,double emin,
		       double emax,int nenergies,double tolerance);
	friend void arpackinterface4dPI(matrix &tmatR,diagmat &Ir0R,
					diagmat &IR,
					matrix &Hrot,diagmat &inertia,
					matrix &Tphi,
					matrix &delcot,matrix &cotmat,
					diagmat &cosphi,
					matrix &phiA,
					matrix &phiB,diagmat &pot4dsym,
					vector &initialv,
					int Rpoints,int nsize1,int nsize2,
					int niter,
					double emin,double emax,int maxnev,
					double tolerance);
	friend cvector Tpsifort3d(matrix &Ttheta1,matrix &Ttheta2,diagmat &I1,
				 diagmat &I2,matrix &Tphi,
				 cvector &v,int nsize1,int nsize2);
	friend cvector Hpsifort4d(matrix &TR,matrix &Ttheta1,matrix &Ttheta2,
				  diagmat &I1,
				  diagmat &I2,matrix &Tphi,diagmat &pot4d,
				  cvector &v,int nsize1,int nsize2,int sizeR);
	friend diagmat sqrtdiagmat(const diagmat& v);
	friend cdiagmat complexm(diagmat &Re,diagmat &Im);
	friend cdiagmat complexm(diagmat &Re);
	friend diagmat realm(cdiagmat &a);
	friend diagmat imagm(cdiagmat &a);
	friend vector Hpsifort4dreal(matrix &TR,matrix &Ttheta1,
				      matrix &Ttheta2,diagmat &I1,
				      diagmat &I2,matrix &Tphi,
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR,vector &uRe,vector &wRe);
	
	friend vector Hpsifort4drealE(matrix &TR,
				      matrix &Ttheta1,matrix &Ttheta2,
				      matrix &Ttheta1R,matrix &Ttheta2R,
				      diagmat &Inertia,
				      matrix &Tphi,matrix &cosphipp,
				      diagmat &cosphi,matrix &sinphip,
				      matrix &cotmat1,matrix &cotmat2,
				      matrix &delcot1,
				      matrix &delcot2,diagmat &InertiaMR,
				      diagmat &InertiaMR2,	    
				      diagmat &pot4d,
				      vector &v,int nsize1,int nsize2,
				      int sizeR);
	friend void LanczosPI(matrix &tmatR,diagmat &IroR,diagmat &IR,
			      matrix &Hrot,diagmat &inertia,
	       matrix &Tphi,matrix &delcot,matrix &cotmat,diagmat &cosphi,
	       matrix &phiA,matrix &phiB,diagmat &pot4dsym,vector &v,
	       int Rpoints,int nsize1,int nsize2,int niter,double emin,
	       double emax);
};
class cdiagmat {
	int Sz;
	complex* TheMatrix;
public:
	cdiagmat(int);
	cdiagmat(const cdiagmat& v);
	~cdiagmat();
	complex& operator() (int Index) const;
	cdiagmat& operator=(const cdiagmat& arg);

/* operations to code (friends) */
	friend const cmatrix  operator*(const cmatrix&, const cdiagmat&);
	friend const cmatrix  operator*(const cdiagmat&, const cmatrix&);
	friend const cdiagmat operator*(const cdiagmat& arg1, const cdiagmat& arg2);
	friend const cdiagmat operator*(const complex arg1, const cdiagmat& arg2);
	friend const cdiagmat operator*(const cdiagmat& arg1, const complex arg2);
	friend const cvector  operator*(const cdiagmat& arg1, const cvector& arg2);
	friend cdiagmat operator+(const cdiagmat& arg1, const cdiagmat& arg2);
	friend cdiagmat operator+(const complex arg1, const cdiagmat& arg2);
	friend cdiagmat operator+(const cdiagmat& arg1, const complex arg2);
	friend cmatrix  operator+(const cmatrix& arg1, const cdiagmat& arg2);
	friend cmatrix  operator+(const cdiagmat& arg1, const cmatrix& arg2);
	friend cmatrix  operator-(const cmatrix& arg1, const cdiagmat& arg2);
	friend cmatrix  operator-(const cdiagmat& arg1, const cmatrix& arg2);
	friend cdiagmat operator-(const cdiagmat& arg1, const cdiagmat& arg2);
	friend cdiagmat operator-(const complex arg1, const cdiagmat& arg2);
	friend cdiagmat operator-(const cdiagmat& arg1, const complex arg2);
	friend cdiagmat inverse(const cdiagmat& v);
	friend cdiagmat sqrtdiagmat(const cdiagmat& v);
	friend cdiagmat complexm(diagmat &Re,diagmat &Im);
	friend cdiagmat complexm(diagmat &Re);
	friend diagmat realm(cdiagmat &a);
	friend diagmat imagm(cdiagmat &a);
	friend cdiagmat transpose(const cdiagmat& arg1);
};
class cmatrix {
	int Rows;
	int Columns;
	complex* TheMatrix;
public:
	cmatrix(int Rows,int Columns);
 	cmatrix(const cmatrix& m); 
	~cmatrix();
	cmatrix& operator=(const cmatrix& arg);
/* indexing */
	complex& operator() (int rindex,int cindex) const;
/* operations to code (friends) */
	friend cmatrix operator+(const cmatrix&, const cmatrix&);
	friend cmatrix operator-(const cmatrix& arg1, const cmatrix& arg2);
	friend const cmatrix operator*(const cmatrix&, const cmatrix&);
	friend const cvector operator*(const cmatrix&, const cvector&); 
	friend const cmatrix operator*(const complex arg1, const cmatrix& arg2);
	friend cmatrix transpose(const cmatrix& arg1);
	friend vector hermdl(cmatrix& a);
	friend const cmatrix  operator*(const cmatrix&, const cdiagmat&);
	friend const cmatrix  operator*(const cdiagmat&, const cmatrix&);
	friend cmatrix  operator+(const cmatrix& arg1, const cdiagmat& arg2);
	friend cmatrix  operator+(const cdiagmat& arg1, const cmatrix& arg2);
	friend cmatrix  operator-(const cmatrix& arg1, const cdiagmat& arg2);
	friend cmatrix  operator-(const cdiagmat& arg1, const cmatrix& arg2);
	friend vector compgene(cmatrix& A,cmatrix& B);
	friend cmatrix complexm(const matrix &Re,const matrix &Im);
	friend cmatrix complexm(const matrix &Re);
	friend matrix realm(cmatrix &a);
	friend matrix imagm(cmatrix &a);
	friend cmatrix ql(cvector &D,cvector &E);
	friend cvector zmatvec(const cmatrix& arg1, const cvector& arg2,complex scalar);
	friend cvector diag(cmatrix& a,cmatrix &VL,cmatrix &VR);
	friend cmatrix inverse(const cmatrix &a);
};
inline double& vector::operator() (int index) const
{
#ifdef DEBUG
    if (index < 0 || index >=Sz) {
		cerr<<"Index exceeds vector size:\n"<<"index="<<index<<"Actual vector size is: "<<Sz<<endl;
		exit(1);
	}
#endif
    return TheVector[index];
}
inline double& diagmat::operator() (int index) const
{
#ifdef DEBUG
    if (index < 0 || index >=Sz) {
		cerr<<"Index exceeds vector size:\n"<<"index="<<index<<"Actual vector size is: "<<Sz<<endl;
		exit(1);
	}
#endif
    return TheMatrix[index];
}
inline complex& cdiagmat::operator() (int index) const
{
#ifdef DEBUG
    if (index < 0 || index >=Sz) {
		cerr<<"Index exceeds vector size:\n"<<"index="<<index<<"Actual vector size is: "<<Sz<<endl;
		exit(1);
	}
#endif
    return TheMatrix[index];
}
inline complex& cvector::operator() (int index) const
{
#ifdef DEBUG
    if (index < 0 || index >=Sz) {
		cerr<<"Index exceeds cvector size:\n"<<"index="<<index<<"Actual cvector size is: "<<Sz<<endl;
		exit(1);
	}
#endif
    return TheVector[index];
}
inline double& matrix::operator() (int rindex,int cindex) const
{
#ifdef DEBUG
  if (rindex < 0 || rindex >=Rows || cindex< 0 || cindex >=Columns) {
    cerr<<"Indices exceed matrix size:\n"<<"rindex="<<rindex<<",cindex="<<cindex<<"Actual matrix size is: "<<Rows<<" by "<<Columns<<endl;
    exit(1);
  }
#endif
  return TheMatrix[cindex*Rows+rindex];
}
// class symmatrix: public matrix
// {
// 	friend vector operator*(const symmatrix& arg1, const vector& arg2);
// };
// class hermmatrix : public cmatrix
// {
// 	friend cvector operator*(const hermmatrix& arg1, const cvector& arg2);
// };
void gridsize(int &nsize,double &dx,double xmax,double mass,double energy);
void thegrid(int nsize,double dx,vector &gridvec);
void potential(int nsize,vector &gridvec,double V0,double a,
double xmax,double lambda,double eta,vector &potvec,vector &absreacvec,
vector &absprodvec,double epsilon0);
void kinetic(matrix &tmat,int nsize,double dx, double mass);
void brute(matrix& tmat,vector& potvec,vector& absreacvec,vector& absprodvec,
double &Ne,int nsize,double E,double Ntheo);
EXTERN void FORTRAN(dsyev)(char *JOBZ,char *UPLO,int *N,double *,int *LDA,
double *,double *, int *LWORK,int *INFO);
void lanc(int niter,const matrix H,const int nsize,const vector start);
EXTERN void FORTRAN(f02axf)(double*,int *,double*,int *,int *,double *,
        double*,int *,double*,int *,double*,double*,double*,int *IFAIL);
EXTERN void FORTRAN(cheev)(char* JOBZ,char* UPLO,int* N,complex*,
		int* LDA,double*,complex*,int*,double *,int *);
vector davidson(cmatrix& A,const int niter,cvector& v,
const int nsize);
EXTERN void FORTRAN(f02gjf)(int*,double *,int*,double *,int*,double *,int*,
           double *,int*,double *,double *,double *,
            double *,int *,double *,int *,
           double *,int*,int*,int*);
EXTERN void FORTRAN(m01caf)(double *,int*,int*,char*,int*);
EXTERN void FORTRAN(f01aaf)(double *,int*,int*,double *,int*,double *,int*);
EXTERN void FORTRAN(f03aaf)(double*,int*,int*,double*,double*,int*);
EXTERN void FORTRAN(f02abf)(double *,int *,int *,double *,double *,int *,double
*,int *);
EXTERN void FORTRAN(dgemm)(char *transa,char *transb,int *m,int *n,
int *k,double *alpha,const double *arg1,int *lda,const double *arg2,
int *ldb,double *beta,double *prod,int *ldc);
EXTERN void FORTRAN(zgemm)(char *transa,char *transb,int *m,int *n,
int *k,complex *alpha,const complex *arg1,int *lda,const complex *arg2,
int *ldb,complex *beta,complex *prod,int *ldc);
EXTERN void FORTRAN(dgemv)(char *trans,int *m,int *n,double *alpha,
const double * a,int *lda,const double *x,int *incx,double *beta,double *prod,
int *incy);
EXTERN void FORTRAN(zgemv)(char *trans,int *m,int *n,complex *alpha,
const complex *a,int *lda,const complex *x,int *incx,complex *beta,
complex *prod,int *incy);
EXTERN void FORTRAN(zgesvd)(char *JOBU,char *JOBVT,int *M,int *N,
				complex *A,int *LDA,double *S,complex *U,
				int *LDU,complex *VT,int *LDVT,complex *WORK,
	            int *LWORK,double *RWORK,int *INFO);
EXTERN void FORTRAN(dpotri)(char *UPLO,int *N,double *a,int *LDA,int *INFO);
EXTERN void FORTRAN(dpotrf)(char *UPLO,int *N,double *a,int *LDA,int *INFO);
EXTERN void FORTRAN(dgeev)(char *JOBVL,char *JOBVR,int *N,double *a,int *LDA,
		       double *WR,double *WI,double *VL,int *LDVL,double *VR,
		       int *LDVR,double *WORK,int *LWORK,int *INFO);
EXTERN void FORTRAN(dgetri)(int *N,double *inv,int *LDA,int *IPIV,double *WORK,
int *LWORK,int *INFO);
EXTERN void FORTRAN(four1)(double *v,int *N,int *flag);
EXTERN void FORTRAN(spotrf)(char *UPLO,int *N,double *B,int *LDA,int *INFO);
EXTERN void  FORTRAN(ssygst)(int * ITYPE,char * UPLO,int *N,double *A,int *LDA,
			 double *B,int *LDB,int *INFO);
EXTERN void FORTRAN(dsteqr)(char *COMPZ,int *N,double *D,double *E,double *Z,
			int *LDZ,double *WORK,int *INFO);
EXTERN void FORTRAN(zsteqr)(char *COMPZ,int *N,complex *D,complex *E,complex *Z,
			int *LDZ,complex *WORK,int *INFO);

EXTERN void  FORTRAN(dmatvecprod)(double *prod,double *arg1,double *arg2,int *Sz);
EXTERN void FORTRAN(dotprod)(double *arg1,double *arg2,int *Sz,double *dot);
EXTERN void FORTRAN(cdotprod)(complex *arg1,complex *arg2,int *Sz,complex *dot);
EXTERN void FORTRAN(zgeev)(char *JOBVL,char *JOBVR,int *N,complex *a,int *LDA,
                       complex *W,complex *VL,int *LDVL,
        complex *VR,int *LDVR,complex *WORK,int *LWORK,
        double *RWORK,int *INFO);
EXTERN void  FORTRAN(zgesv)(int *N,int *NRHS,complex *A,int *LDA,int *IPIV,
				complex *B,int *LDB,int *INFO);
