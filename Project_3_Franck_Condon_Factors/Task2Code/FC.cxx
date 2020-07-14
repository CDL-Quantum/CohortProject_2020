/*
 * FC is a program to calculate Franck-Condon factors
Copyright (C) 1999,2006, 2020 Pierre-Nicholas Roy

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#############                                                   ##############
Note: This is not the full version, but a modified version by Matthew Schmidt
with permission from Pierre-Nicholas Roy for use with the CDL Cohort Project 2020.
The full version can be found on: https://github.com/pnroy/FC
Modifications:
- added comments (including references to relevant papers)
- removed output files containing matrices and delta vector
- removed capability to input Force Constants instead of normal modes
#############                                                   ##############
*/


#include "FC.h"
double int00=1.;
int factrl(int n);
vector triplelmn(int ni,int nj,int nk,int nl,int nm,int nn,
		 int modei,int modej,int modek,int model,
		 int modem,int moden,
		 matrix& Pmat,matrix &Qmat,matrix &Rmat,
		 diagmat& II,vector& delta);
vector tripleij(int ni,int nj,int nk,int nl,int nm,
	       int modei,int modej,int modek,int model,int modem,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta);
vector triplei(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta);
vector triplelm(int ni,int nj,int nk,int nl,int nm,
	       int modei,int modej,int modek,int model,int modem,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta);
vector triple(int ni,int nj,int nk,int modei,int modej,int modek,
matrix& Pmat,diagmat& II,vector& delta);
complex HermiteH(int n,complex x);
complex sigma(int k,matrix &Qmat,matrix &Rmat,vector &delta);
complex tau(int k,matrix &Pmat,diagmat &II,vector &delta);
double oi0(int k,int n,double int00,matrix &Pmat,diagmat &II,vector &delta);
double o0i(int k,int n,double int00,matrix& Qmat,matrix& Rmat,vector&
delta);
double o11(int k,int i,double int00,matrix& Rmat,diagmat& II,
matrix& Pmat,matrix& Qmat,vector& delta);
double o21(int k,int i,double int00,matrix& Rmat,diagmat& II,
matrix& Pmat,matrix& Qmat,vector& delta);
double  boltzman(double EoverT);
double oki(int nkp1,int ni,int k,int i,double int00,matrix& Rmat,diagmat&
II,matrix& Pmat,matrix& Qmat,vector& delta);
double okl00(int nkp1,int nl,int k,int l,double int00,diagmat& II,
matrix& Pmat,vector& delta);
double lorentz(double omega,double w0,double gamma);
vector lor(double pos,double signal,int size,double gamma,vector grid);
double oklk0(int nk,int nl,int nkp1p,int k,int l,double int00,diagmat II,
matrix Pmat, matrix &Qmat, matrix& Rmat,vector& delta);
double oklkl(int nk,int nl,int nkp1p,int nlp,int k,int l,
double int00,diagmat II,matrix Pmat, matrix &Qmat, matrix& Rmat,
vector& delta);
double o00ij(int nip1,int nj,int i,int j,double int00,matrix& Qmat,
	     matrix& Rmat,vector& delta);
double position(int modei,int modej,int modek,
		int model,int modem,int moden,
		int cati,int catj,int catk,
		int neul,int neum,int neun,
		vector &omegaion,vector &omega);
vector triplehot(int ni,int nj,int nk,int modei,int modej,int modek,
	      matrix& Qmat,matrix &Rmat,diagmat &II,vector& delta);
vector singlecold(int ni,int modei,matrix& Pmat,diagmat& II,vector& delta);
vector singlehot(int ni,int modei,matrix& Qmat,matrix &Rmat,vector& delta);
vector doublecold(int ni,int nj,int modei,int modej,
		  matrix& Pmat,diagmat& II,vector& delta);
vector doublehot(int ni,int nj,int modei,int modej,
	      matrix& Qmat,matrix &Rmat,vector& delta);
vector ijhot(int ni,int nj,int modei,int modej,matrix &Pmat,matrix &Qmat,
	     matrix &Rmat,diagmat &II,vector &delta);
vector ijkhot(int ni,int nj,int nk,int modei,int modej,int modek,
	      matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	      vector &delta);
vector triplel(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta);
vector ijklhot(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	      matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	      vector &delta);
vector iklhot(int ni,int nk,int nl,int modei,int modek,int model,
	      matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	      vector &delta);
double const speedofl=137.0359895;
int main(int argc,char** argv)
{
  if (argc!=2) {
    std::cerr<<"usage: "<<argv[0]<<" <input filename>"<<std::endl;
    exit(0);
  }
  std::strstream filenamein;
  filenamein<<argv[1]<<std::ends;
  ifstream file(filenamein.str());


  int i,j,k,l,lp,m,n,natom;
  int modei;
  int modej;
  int modek;
  int model;
  int modem;
  int moden;
  double pos;
  char Units;
  file >> Units;
  file >> natom;
  int N;

  N=3*natom-6;

  vector massvec(natom);
  diagmat sqmass(3*natom);
  diagmat invsqmass(3*natom);
  for (i=0;i<natom;i++) file >> massvec(i);  
  /* get masses */
  massvec=(1./autoamu)*massvec;
  for (i=0;i<natom;i++) {
    for (j=0;j<3;j++) {
      sqmass(i*3+j)=sqrt(massvec(i));
      invsqmass(i*3+j)=1./sqrt(massvec(i));
    }
  }
  vector geon(3*natom);
  vector geoc(3*natom);
  vector omega(N);
  vector omegaion(N);

  //Lmatrices: eigenvectors of mass-weighted Hessian (unitless)
  //Represented by L and L' (Eqs. 1 & 2) in Roy, Carrington 1995
  //Represented by O_I (Eq. 37) in Quesada 2019 Paper

  //Omega values are sqrt of eigenvalues of mass-weighted Hessian
  //Represented by Omega (Eq. 37) in Quesada 2019 Paper

  matrix Lmat(3*natom,N);
  matrix Lmation(3*natom,N);
  
  for (i=0;i<2;i++) {				
    /*loop over species
      get geometries */
    if (i==0) {
      for (l=0;l<(3*natom);l++) file >> geon(l);
    }
    if (i==1) {
      for (l=0;l<(3*natom);l++) file >> geoc(l);
    }
    for (j=0;j<N;j++)  {       
      if (i==0) {
	file >> omega(j);
	omega(j)/=hatocm;
      }
      if (i==1) {
	file >> omegaion(j);
	omegaion(j)/=hatocm;
      }
      for (k=0;k<(3*natom);k++) {
	if (i==0) file >> Lmat(k,j);
	if (i==1) file >> Lmation(k,j);
      }
    }
  }

  double wmin,wmax,fwhm,gamma,v00;
  double T;
  double res;  
  /* resolution of output data file */
  int cationquanta;
  int neutralquanta;
  file >> wmin;
  file >> wmax;
  file >> T;
  file >> fwhm;
  file >> res;
  file >> v00;
  /* neutral quanta */
  file >>neutralquanta;
  /* cation quanta */
  file >>cationquanta;
  double labellimit;
  file >>labellimit;
  int numberofactivemodes=3;
  file>>numberofactivemodes;
  // calculate center of mass
  vector CMi(3);
  vector CMf(3);
  double totalmass=0.;
  for (i=0;i<natom;i++) {
    totalmass+=massvec(i);
    for (j=0;j<3;j++) {
      CMi(j)+=massvec(i)*geon(i*3+j);
      CMf(j)+=massvec(i)*geoc(i*3+j);    
    }
  }
  CMi=(1./totalmass)*CMi;
  CMf=(1./totalmass)*CMf;
  //subtract eq. geometries from center of mass
  for (i=0;i<natom;i++) {
    for (j=0;j<3;j++) {
      geon(i*3+j)-=CMi(j);
      geoc(i*3+j)-=CMf(j);    
    }
  }
  
  /* perform orthogonalisation or not? */
  char ortho_of_S='N';
  char ortho_of_L;
  ortho_of_L='Y';
  
  /* set limits */
  int ni,nj,nk,nl,nm,nn;
  /* definition of gamma */
  gamma=fwhm;
  vector Rvec=geoc-geon;

  /* Convert Units to Bohr */
  if (Units == 'A') Rvec=(1./.529177)*Rvec;  
  if (Units != 'A' && Units !='B') std::cerr<<"wrong Units flag"<<std::endl;

  vector Rvecmass=sqmass*Rvec;

  /* modified Gram-Schmidt orthogonalisation on both L's */

  matrix oLmat=Lmat; 
  matrix oLmation=Lmation;
  for (k=0;k<N;k++) { 
    double proj1=0.0;
    double proj2=0.0;
    for (i=0;i<(3*natom);i++) {
      proj1+=oLmat(i,k)*oLmat(i,k); 
      proj2+=oLmation(i,k)*oLmation(i,k);
    }
    proj1=sqrt(proj1);
    proj2=sqrt(proj2);
    for (i=0;i<(3*natom);i++) {
      oLmat(i,k)=oLmat(i,k)/proj1; oLmation(i,k)=oLmation(i,k)/proj2; }
    for (j=k+1;j<N;j++) { 
      proj1=0.0; proj2=0.0; 
      for (i=0;i<(3*natom);i++) { 
	proj1+=oLmat(i,k)*oLmat(i,j); 
	proj2+=oLmation(i,k)*oLmation(i,j);
      } 
      for (i=0;i<(3*natom);i++) { 
	oLmat(i,j)=oLmat(i,j)-oLmat(i,k)*proj1;
	oLmation(i,j)=oLmation(i,j)-oLmation(i,k)*proj2; 
      } 
    }
  } 

  /* matrix diff=Lmat-oLmat; */
  if (ortho_of_L == 'Y') {
    Lmat=oLmat;
    Lmation=oLmation;
  }

  diagmat lambda(N);
  diagmat lambdaion(N);
  for (j=0;j<N;j++) {
    lambda(j)=sqrt(omega(j)); lambdaion(j)=sqrt(omegaion(j)); 
  }


  //S matrix is Duschinsky transformation 
  //Represented by O_D in Quesada 2019 paper (see note after Eq. 42c)
  //Represented by S (Eq. 4) in Roy, Carrington 1995 paper 
  matrix OSmat=transpose(Lmation)*Lmat;
  matrix Smat=OSmat;

  /* orthogonalisation of S */ 
  for (k=0;k<N;k++) {
    double proj1=0.0; 
    for (i=0;i<N;i++) {
      proj1+=Smat(i,k)*Smat(i,k); 
    } proj1=sqrt(proj1); 
    for (i=0;i<N;i++) { 
      Smat(i,k)=Smat(i,k)/proj1; 
    } 
   for (j=k+1;j<N;j++) { 
      proj1=0.0; 
      for (i=0;i<N;i++) {
	proj1+=Smat(i,k)*Smat(i,j); 
      } for (i=0;i<N;i++) {
	Smat(i,j)=Smat(i,j)-Smat(i,k)*proj1; 
      } 
    } 
  }
  if (ortho_of_S == 'N') Smat=OSmat;

  matrix orthoS=transpose(Smat)*Smat;
  matrix orthoL=transpose(Lmat)*Lmat;
  matrix orthoLion=transpose(Lmation)*Lmation;
  
  vector dvec=transpose(Lmation)*Rvecmass;   //Q (Eq. 5) in Roy, Carrington 1995
  /* for sharp initial in terms of final */
  matrix Jmat=transpose(Lmat)*Lmation;
  vector K=transpose(Lmat)*sqmass*(geoc-geon);
  
  diagmat II(N);
  for (j=0;j<N;j++) II(j)=1.; //identity matrix
  matrix JJ=lambdaion*(Smat*inverse(lambda)); //A (Eq. 42b) in Quesada 2019 Paper
  matrix Qmat=inverse((II+transpose(JJ)*JJ));
  matrix Pmat=JJ*Qmat*transpose(JJ);
  matrix Rmat=Qmat*transpose(JJ);

  //Displacement Vector, delta
  //d (Eq. 42c) in Quesada 2019 Paper
  vector delta=lambdaion*dvec;

  /*  <0|0> transition  */
  for (j=0;j<N;j++) int00*=(omegaion(j)/omega(j));
  int00=pow(int00,.25);
  int00*=exp(log(sqrt(det(Qmat)))+((double) (N) /2.)*log(2.));
  int00*=exp(-.5*delta*((II-Pmat)*delta))*sqrt(det(Smat));
  std::cout<<"00 band "<<int00<<std::endl;

  /*  spectrum */
  /*  definition of the spectral window */
  double signal,fcfac;
  double o;
  double length=wmax-wmin;
  int size=(int) ceil(length/res);
  vector sumlor(size);
  vector grid(size);
  for (i=0;i<size;i++) grid(i)=((double) i)*res+wmin;
  int temperature=(int) T;

  std::strstream specname;
  std::strstream specfname;
  specname<<filenamein.str()<<".sticks.out"<<ends;
  specfname<<filenamein.str()<<".spec.out"<<ends;
  ofstream spec(specname.str());
  ofstream specf(specfname.str());
  std::strstream fcfilename;
  fcfilename<<filenamein.str()<<".fc.tex"<<ends;
  ofstream fcfile(fcfilename.str());
	fcfile<<"\\documentclass[12pt]{article}"<<std::endl;
	fcfile<<"\\begin{document}"<<std::endl;
	fcfile<<"\\include{epsf}"<<std::endl;
	fcfile<<"\\unitlength1in"<<std::endl;
	fcfile<<"\\renewcommand{\\baselinestretch}{1.5}"<<std::endl;
	fcfile<<"\\begin{table}"<<std::endl;
	fcfile<<"\\begin{tabular}{rcrrr}"<<std::endl;
	fcfile<<"Peak& Assignment&Position (cm$^{-1}$)";
	fcfile<<"& Rel. FC factor& Signal\\\\"<<std::endl;
	fcfile<<"\\hline"<<std::endl;
  fcfile.precision(3);
  fcfile.setf(ios::fixed,ios::floatfield);


  labellimit/=(int00*int00);
  
  /*  <0|0> band */
  int fccount=0;
  pos=v00;
  if (pos <wmax && pos >wmin) {
    signal=1.;
    //modification below and to all signal output to spec file
    spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
    fccount+=1;
    fcfile<<fccount<<"& $1^0_0$&"<<pos<<"&1&1\\\\"<<std::endl;
    sumlor=sumlor+lor(pos,signal,size,gamma,grid);
  } 
  std::cout<<"start recursive calculation"<<std::endl;  
  int istart=1;
  int jstart=1;
  int kstart=1;
  int lstart=1;
  int mstart=1;
  int nstart=1;
  for (modei=0;modei<(numberofactivemodes);modei++) {
    /* <i|0> */
    ni=cationquanta;
    vector t=singlecold(ni,
			modei,
			Pmat,II,delta);
    moden=0; n=0;
    modem=0; m=0;
    model=0; l=0;
    modek=0; k=0;
    modej=0; j=0;
    for (i=istart;i<ni;i++) {
       pos=position(
		   modei,modej,modek,
		   model,modem,moden,
		   i,j,k,l,m,n,
		   omegaion,omega);
       pos+=v00;
       if (pos <wmax && pos >wmin) {
	o=t(((((i)))));
	fcfac=o*o/(int00*int00);
	signal=fcfac*
	  boltzman((omega(model)*((double) l)/T
		    +omega(modem)*((double) m)/T
		    +omega(moden)*((double) n)/T));
	sumlor=sumlor+lor(pos,signal,size,gamma,grid);
	if (signal >= labellimit ) {
	  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
	  fccount+=1;
	  fcfile<<fccount<<"&$";
	  fcfile<<modei+1<<"_"<<l<<"^"<<i<<"$& ";
	  fcfile<<pos;
	  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
	}
       }
    }
    /* <0|i> */
    ni=neutralquanta;
    vector thot=singlehot(ni,
			  modei,
			  Qmat,Rmat,delta);
    moden=0; n=0;
    modem=0; m=0;
    model=0; l=0;
    modek=0; k=0;
    modej=0; j=0;
    for (i=istart;i<ni;i++) {
      pos=position(
		   modej,modej,modek,
		   modei,modem,moden,
		   j,j,k,i,m,n,
		   omegaion,omega);
      pos+=v00;
      if (pos <wmax && pos >wmin) {
	o=thot(((((i)))));
	fcfac=o*o/(int00*int00);
	signal=fcfac*
	  boltzman((omega(modei)*((double) i)/T
		    +omega(modem)*((double) m)/T
		    +omega(moden)*((double) n)/T));
	sumlor=sumlor+lor(pos,signal,size,gamma,grid);	
	if (signal >= labellimit ) {
	  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
	  fccount+=1;
	  fcfile<<fccount<<"&$";
	  fcfile<<modei+1<<"_"<<i<<"^"<<l<<"&$ ";
	  fcfile<<pos;
	  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
	}
      }
    }
    for (modej=0;modej<(numberofactivemodes);modej++) {
      /* <i|j> */
      ni=cationquanta;
      nj=neutralquanta;
      vector t=ijhot(ni,nj,modei,modej,Pmat,Qmat,Rmat,II,delta);
      moden=0; n=0;
      modem=0; m=0;
      model=0; l=0;
      modek=0; k=0;
      for (i=istart;i<ni;i++) {
	for (j=jstart;j<nj;j++) {
	  pos=position(
		       modei,modek,modek,
		       modej,modem,moden,
		       i,k,k,j,m,n,
		       omegaion,omega);
	  pos+=v00;
	  if (pos <wmax && pos >wmin) {
	    o=t(((((i*nj+j)))));
	    fcfac=o*o/(int00*int00);
	    signal=fcfac*
	      boltzman((omega(modej)*((double) j)/T
			+omega(modem)*((double) m)/T
			+omega(moden)*((double) n)/T));
	    sumlor=sumlor+lor(pos,signal,size,gamma,grid);
	    if (signal >= labellimit ) {
	      spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
	      fccount+=1;
	      fcfile<<fccount<<"&$";
	      fcfile<<modei+1<<"_"<<l<<"^"<<i<<" ";
	      fcfile<<modej+1<<"_"<<j<<"^"<<l<<"$& ";
	      fcfile<<pos;
	      fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
	    }
	  }
	}
      }
    }
    for (modej=(modei+1);modej<numberofactivemodes;modej++) {
      /* <ij|0> */
      ni=cationquanta/2+1;
      nj=cationquanta/2+1;
      vector t=doublecold(ni,nj,
			  modei,modej,
			  Pmat,II,delta);
      moden=0; n=0;
      modem=0; m=0;
      model=0; l=0;
      modek=0; k=0;
      for (i=istart;i<ni;i++) {
	for (j=jstart;j<nj;j++) {
	  pos=position(
		       modei,modej,modek,
		       model,modem,moden,
		       i,j,k,l,m,n,
		       omegaion,omega);
	  pos+=v00;
	  if (pos <wmax && pos >wmin) {
	    o=t(((((i*nj+j)))));
	    fcfac=o*o/(int00*int00);
	    signal=fcfac*
	      boltzman((omega(model)*((double) l)/T
			+omega(modem)*((double) m)/T
			+omega(moden)*((double) n)/T));
	    sumlor=sumlor+lor(pos,signal,size,gamma,grid);
	    if (signal >= labellimit ) {
	      spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
	      fccount+=1;
	      fcfile<<fccount<<"&$";
	      fcfile<<modei+1<<"_"<<l<<"^"<<i<<" ";
	      fcfile<<modej+1<<"_"<<m<<"^"<<j<<"$& ";
	      fcfile<<pos;
	      fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
	    }
	  }
	}
      }
      /* <0|ij> */
      ni=neutralquanta/2+1;
      nj=neutralquanta/2+1;
      vector t2=doublehot(ni,nj,
			  modei,modej,
			  Qmat,Rmat,delta);
      moden=0; n=0;
      modem=0; m=0;
      model=0; l=0;
      modek=0; k=0;
      for (i=istart;i<ni;i++) {
	for (j=jstart;j<nj;j++) {
	  pos=position(
		       modek,modek,modek,
		       modei,modej,moden,
		       k,k,k,i,j,n,
		       omegaion,omega);
	  pos+=v00;
	  if (pos <wmax && pos >wmin) {
	    o=t2(((((i*nj+j)))));
	    fcfac=o*o/(int00*int00);
	    signal=fcfac*
	      boltzman((omega(modei)*((double) i)/T
			+omega(modej)*((double) j)/T
			+omega(moden)*((double) n)/T));
	    sumlor=sumlor+lor(pos,signal,size,gamma,grid);
	    if (signal >= labellimit ) {
	 
	      spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
	      fccount+=1;
	      fcfile<<fccount<<"&$";
	      fcfile<<modei+1<<"_"<<i<<"^"<<l<<" ";
	      fcfile<<modej+1<<"_"<<j<<"^"<<l<<"$& ";
	      fcfile<<pos;
	      fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
	    }
	  }
	}
      }
      for (modek=0;modek<(numberofactivemodes);modek++) { 
	/* <ij|k> */
	ni=cationquanta/2 +1;
	nj=cationquanta/2 +1;
	nk=neutralquanta;
	vector t=ijkhot(ni,nj,nk,
			modei,modej,modek,Pmat,Qmat,Rmat,II,delta);
	moden=0; n=0;
	modem=0; m=0;
	model=0; l=0;
	for (i=istart;i<ni;i++) {
	  for (j=jstart;j<nj;j++) {
	    for (k=kstart;k<nk;k++) {
	      pos=position(
			   modei,modej,model,
			   modek,modem,moden,
			   i,j,l,k,m,n,
			   omegaion,omega);
	      pos+=v00;
	      if (pos <wmax && pos >wmin) {
		o=t(((((i*nj+j)*nk+k))));
		fcfac=o*o/(int00*int00);
		signal=fcfac*
		  boltzman((omega(modek)*((double) k)/T
			    +omega(modem)*((double) m)/T
			    +omega(moden)*((double) n)/T));
		sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		if (signal >= labellimit ) {
		  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
		  fccount+=1;
		  fcfile<<fccount<<"&$";
		  fcfile<<modei+1<<"_"<<l<<"^"<<i<<" ";
		  fcfile<<modej+1<<"_"<<m<<"^"<<j<<" ";
		  fcfile<<modek+1<<"_"<<k<<"^"<<l<<"$& ";
		  fcfile<<pos;
		  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		}
	      }
	    }
	  }
	}
	/* <k|ij> */
	nk=cationquanta;
	nj=neutralquanta/2+1;
	ni=neutralquanta/2+1;
	vector t3=iklhot(nk,ni,nj,modek,modei,modej,Pmat,Qmat,Rmat,II,delta);
	moden=0; n=0;
	modem=0; m=0;
	model=0; l=0;
	for (i=istart;i<ni;i++) {
	  for (j=jstart;j<nj;j++) {
	    for (k=kstart;k<nk;k++) {
	      pos=position(
			   modek,model,model,
			   modei,modej,moden,
			   k,l,l,i,j,n,
			   omegaion,omega);
	      pos+=v00;
	      if (pos <wmax && pos >wmin) {
		o=t3(((((k*ni+i)*nj+j))));
		fcfac=o*o/(int00*int00);
		signal=fcfac*
		  boltzman((omega(modei)*((double) i)/T
			    +omega(modej)*((double) j)/T
			    +omega(moden)*((double) n)/T));
		sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		if (signal >= labellimit ) {
		  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
		  fccount+=1;
		  fcfile<<fccount<<"&$";
		  fcfile<<modei+1<<"_"<<i<<"^"<<l<<" ";
		  fcfile<<modej+1<<"_"<<j<<"^"<<l<<" ";
		  fcfile<<modek+1<<"_"<<n<<"^"<<k<<"$& ";
		  fcfile<<pos;
		  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		}
	      }
	    }
	  }
	}
	for (model=(modek+1);model<(numberofactivemodes);model++) {
	  /* <ij|kl> */
	  ni=cationquanta/2 +1;
	  nj=cationquanta/2 +1;
	  nk=neutralquanta/2+1;
	  nl=neutralquanta/2+1;
	  vector t=ijklhot(ni,nj,nk,nl,modei,modej,modek,model,
			   Pmat,Qmat,Rmat,II,delta);
	  moden=0; n=0;
	  modem=0; m=0;
	  for (i=istart;i<ni;i++) {
	    for (j=jstart;j<nj;j++) {
	      for (k=kstart;k<nk;k++) {
		for (l=lstart;l<nl;l++) {
		  pos=position(
			       modei,modej,moden,
			       modek,model,moden,
			       i,j,n,k,l,n,
			       omegaion,omega);
		  pos+=v00;
		  if (pos <wmax && pos >wmin) {
		    o=t(((((i*nj+j)*nk+k)*nl+l)));
		    fcfac=o*o/(int00*int00);
		    signal=fcfac*
		      boltzman((omega(modek)*((double) k)/T
				+omega(model)*((double) l)/T
				+omega(moden)*((double) n)/T));
		    sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		    if (signal >= labellimit ) {
		      spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
		      fccount+=1;
		      fcfile<<fccount<<"&$";
		      fcfile<<modei+1<<"_"<<m<<"^"<<i<<" ";
		      fcfile<<modej+1<<"_"<<m<<"^"<<j<<" ";
		      fcfile<<modek+1<<"_"<<k<<"^"<<n<<" ";
		      fcfile<<model+1<<"_"<<l<<"^"<<n<<"$& ";
		      fcfile<<pos;
		      fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      for (modek=(modej+1);modek<(numberofactivemodes);modek++) { 
	/* <ijk|0> */
	ni=cationquanta/3+1;
	nj=cationquanta/3+1;
	nk=cationquanta/3+1;
	vector t=triple(ni,nj,nk,
			modei,modej,modek,
			Pmat,II,delta);
	moden=0; n=0;
	modem=0; m=0;
	model=0; l=0;
	for (i=istart;i<ni;i++) {
	  for (j=jstart;j<nj;j++) {
	    for (k=kstart;k<nk;k++) {
	      pos=position(
			   modei,modej,modek,
			   model,modem,moden,
			   i,j,k,l,m,n,
			   omegaion,omega);
	      pos+=v00;
	      if (pos <wmax && pos >wmin) {
		o=t(((((i*nj+j)*nk+k))));
		fcfac=o*o/(int00*int00);
		signal=fcfac*
		  boltzman((omega(model)*((double) l)/T
			    +omega(modem)*((double) m)/T
			    +omega(moden)*((double) n)/T));
		sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		if (signal >= labellimit ) {
		  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
		  fccount+=1;
		  fcfile<<fccount<<"&$";
		  fcfile<<modei+1<<"_"<<l<<"^"<<i<<" ";
		  fcfile<<modej+1<<"_"<<m<<"^"<<j<<" ";
		  fcfile<<modek+1<<"_"<<n<<"^"<<k<<"$& ";
		  fcfile<<pos;
		  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		}
	      }
	    }
	  }
	}
	/* <0|ijk> */
	ni=neutralquanta/3+1;
	nj=neutralquanta/3+1;
	nk=neutralquanta/3+1;
	vector t2=triplehot(ni,nj,nk,modei,modej,modek,Qmat,Rmat,II,delta);
	moden=0; n=0;
	modem=0; m=0;
	model=0; l=0;
	for (i=istart;i<ni;i++) {
	  for (j=jstart;j<nj;j++) {
	    for (k=kstart;k<nk;k++) {
	      pos=position(
			   moden,moden,moden,
			   modei,modej,modek,
			   n,n,n,i,j,k,
			   omegaion,omega);
	      pos+=v00;
	      if (pos <wmax && pos >wmin) {
		o=t2(((((i*nj+j)*nk+k))));
		fcfac=o*o/(int00*int00);
		signal=fcfac*
		  boltzman((omega(modei)*((double) i)/T
			    +omega(modej)*((double) j)/T
			    +omega(modek)*((double) k)/T));
		sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		if (signal >= labellimit ) {
		  spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
		  fccount+=1;
		  fcfile<<fccount<<"&$";
		  fcfile<<modei+1<<"_"<<i<<"^"<<l<<" ";
		  fcfile<<modej+1<<"_"<<j<<"^"<<l<<" ";
		  fcfile<<modek+1<<"_"<<k<<"^"<<l<<"$& ";
		  fcfile<<pos;
		  fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		}
	      }
	    }
	  }
	}

	for (model=0;model<(numberofactivemodes);model++) {
	  /* <ijk|l> */
#ifdef DUSHQUICK
	  if ((model==modei) || (model==modej) || (model==modek) )
#endif
	    {
	    ni=cationquanta/3+1;
	    nj=cationquanta/3+1;
	    nk=cationquanta/3+1;
	    nl=neutralquanta;
	    vector t=triplel(ni,nj,nk,nl,
			     modei,modej,modek,model,
			     Pmat,Qmat,Rmat,II,delta);
	    moden=0; n=0;
	    modem=0; m=0;
	    for (i=istart;i<ni;i++) {
	      for (j=jstart;j<nj;j++) {
		for (k=kstart;k<nk;k++) {
		  for (l=lstart;l<nl;l++) {
		    pos=position(
				 modei,modej,modek,
				 model,modem,moden,
				 i,j,k,l,m,n,
				 omegaion,omega);
		    pos+=v00;
		    if (pos <wmax && pos >wmin) {
		      o=t(((((i*nj+j)*nk+k)*nl+l)));
		      fcfac=o*o/(int00*int00);
		      signal=fcfac*
			boltzman((omega(model)*((double) l)/T
				  +omega(modem)*((double) m)/T
				  +omega(moden)*((double) n)/T));
		      sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		      if (signal >= labellimit ) {
			spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
			fccount+=1;
			fcfile<<fccount<<"&$";
			fcfile<<modei+1<<"_"<<n<<"^"<<i<<" ";
			fcfile<<modej+1<<"_"<<m<<"^"<<j<<" ";
			fcfile<<modek+1<<"_"<<n<<"^"<<k<<" ";
			fcfile<<model+1<<"_"<<l<<"^"<<n<<"$& ";
			fcfile<<pos;
			fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		      }
		    }
		  }
		}
	      }
	    }
	    /* <l|ijk> */
	    ni=neutralquanta/3+1;
	    nj=neutralquanta/3+1;
	    nk=neutralquanta/3+1;
	    nl=cationquanta;
	    vector t2=triplei(nl,ni,nj,nk,model,modei,modej,modek,
			      Pmat,Qmat,Rmat,II,delta);
	    moden=0; n=0;
	    modem=0; m=0;
	    for (i=istart;i<ni;i++) {
	      for (j=jstart;j<nj;j++) {
		for (k=kstart;k<nk;k++) {
		  for (l=lstart;l<nl;l++) {
		    pos=position(
				 model,moden,moden,
				 modei,modej,modek,
				 l,n,n,i,j,k,
				 omegaion,omega);
		    pos+=v00;
		    if (pos <wmax && pos >wmin) {
		      o=t2(((((l*ni+i)*nj+j)*nk+k)));
		      fcfac=o*o/(int00*int00);
		      signal=fcfac*
			boltzman((omega(modei)*((double) i)/T
				  +omega(modej)*((double) j)/T
				  +omega(modek)*((double) k)/T));
		      sumlor=sumlor+lor(pos,signal,size,gamma,grid);
		      if (signal >= labellimit ) {
			spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
			fccount+=1;
			fcfile<<fccount<<"&$";
			fcfile<<modei+1<<"_"<<i<<"^"<<n<<" ";
			fcfile<<modej+1<<"_"<<j<<"^"<<n<<" ";
			fcfile<<modek+1<<"_"<<k<<"^"<<n<<" ";
			fcfile<<model+1<<"_"<<n<<"^"<<l<<"$& ";
			fcfile<<pos;
			fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  for (modem=(model+1);modem<(numberofactivemodes);modem++) {
#ifdef DUSHQUICK
	    	    if ((model==modei && modem==modej) || (model==modei && modem==modek) || (model==modej && modem==modek) ) 
#endif
	      {
	      /* <ijk|lm> */
	      ni=cationquanta/3+1;
	      nj=cationquanta/3+1;
	      nk=cationquanta/3+1;
	      nl=neutralquanta/2+1;
	      nm=neutralquanta/2+1;
	      vector t=triplelm(ni,nj,nk,nl,nm,
				modei,modej,modek,model,
				modem,Pmat,Qmat,Rmat,II,delta);
	      moden=0; n=0;
	      for (i=istart;i<ni;i++) {
		for (j=jstart;j<nj;j++) {
		  for (k=kstart;k<nk;k++) {
		    for (l=lstart;l<nl;l++) {
		      for (m=mstart;m<nm;m++) {
			pos=position(
				     modei,modej,modek,
				     model,modem,moden,
				     i,j,k,l,m,n,
				     omegaion,omega);
			pos+=v00;
			if (pos <wmax && pos >wmin) {
			  o=t(((((i*nj+j)*nk+k)*nl+l)*nm+m));
			  fcfac=o*o/(int00*int00);
			  signal=fcfac*
			    boltzman((omega(model)*((double) l)/T
				      +omega(modem)*((double) m)/T
				      +omega(moden)*((double) n)/T));
			  sumlor=sumlor+lor(pos,signal,size,gamma,grid);
			  if (signal >= labellimit ) {
			    spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
			    fccount+=1;
			    fcfile<<fccount<<"&$";
			    fcfile<<modei+1<<"_"<<n<<"^"<<i<<" ";
			    fcfile<<modej+1<<"_"<<n<<"^"<<j<<" ";
			    fcfile<<modek+1<<"_"<<n<<"^"<<k<<" ";
			    fcfile<<model+1<<"_"<<l<<"^"<<n<<" ";
			    fcfile<<modem+1<<"_"<<m<<"^"<<n<<"$& ";
			    fcfile<<pos;
			    fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
			  }
			}
		      }
		    }
		  }
		}
	      }
	      /* <lm|ijk> */
	      ni=neutralquanta/3+1;
	      nj=neutralquanta/3+1;
	      nk=neutralquanta/3+1;
	      nl=cationquanta/2+1;
	      nm=cationquanta/2+1;
	      vector t2=tripleij(nl,nm,ni,nj,nk,
				 model,modem,modei,modej,modek,
				 Pmat,Qmat,Rmat,II,delta);
	      moden=0; n=0;
	      for (i=istart;i<ni;i++) {
		for (j=jstart;j<nj;j++) {
		  for (k=kstart;k<nk;k++) {
		    for (l=lstart;l<nl;l++) {
		      for (m=mstart;m<nm;m++) {
			pos=position(
				     model,modem,moden,
				     modei,modej,modek,
				     l,m,n,i,j,k,
				     omegaion,omega);
			pos+=v00;
			if (pos <wmax && pos >wmin) {
			  o=t2(((((l*nm+m)*ni+i)*nj+j)*nk+k));
			  fcfac=o*o/(int00*int00);
			  signal=fcfac*
			    boltzman((omega(modei)*((double) i)/T
				      +omega(modej)*((double) j)/T
				      +omega(modek)*((double) k)/T));
			  sumlor=sumlor+lor(pos,signal,size,gamma,grid);
			  if (signal >= labellimit ) {
			    spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
			    fccount+=1;
			    fcfile<<fccount<<"&$";
			    fcfile<<modei+1<<"_"<<i<<"^"<<n<<" ";
			    fcfile<<modej+1<<"_"<<j<<"^"<<n<<" ";
			    fcfile<<modek+1<<"_"<<k<<"^"<<n<<" ";
			    fcfile<<model+1<<"_"<<n<<"^"<<l<<" ";
			    fcfile<<modem+1<<"_"<<n<<"^"<<m<<"$& ";
			    fcfile<<pos;
			    fcfile<<"& "<<fcfac<<"& "<<signal<<"\\\\"<<std::endl;
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	    for (moden=(modem+1);moden<(numberofactivemodes);moden++) {
	      /* <ijk|lmn> */
	      ni=cationquanta/3+1;
	      nj=cationquanta/3+1;
	      nk=cationquanta/3+1;
	      nl=neutralquanta/3+1;
	      nm=neutralquanta/3+1;
	      nn=neutralquanta/3+1;

#ifdef DUSHQUICK
	      if (modei==model && modej==modem && modek==moden)
#endif
		{
		vector t=triplelmn(ni,nj,nk,nl,nm,nn,
				   modei,modej,modek,model,
				   modem,moden,Pmat,Qmat,Rmat,II,delta);
		for (i=istart;i<ni;i++) {
		  for (j=jstart;j<nj;j++) {
		    for (k=kstart;k<nk;k++) {
		      for (l=lstart;l<nl;l++) {
			for (m=mstart;m<nm;m++) {
			  for (n=nstart;n<nn;n++) {
			    pos=position(
					 modei,modej,modek,
					 model,modem,moden,
					 i,j,k,l,m,n,
					 omegaion,omega);
			    pos+=v00;
			    if (pos <wmax && pos >wmin) {
			      o=t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n);
			      fcfac=o*o/(int00*int00);
			      signal=fcfac*
				boltzman((omega(model)*((double) l)/T
					  +omega(modem)*((double) m)/T
					  +omega(moden)*((double) n)/T));
			      sumlor=sumlor+lor(pos,signal,size,gamma,grid);
			      if (signal >= labellimit ) {
				spec<<pos<<" "<<signal*(int00*int00)<<std::endl;
				fccount+=1;
				fcfile<<fccount<<"&$";
				fcfile<<modei+1<<"_"<<l<<"^"<<i<<" ";
				fcfile<<modej+1<<"_"<<m<<"^"<<j<<" ";
				fcfile<<modek+1<<"_"<<n<<"^"<<k<<"$&";
				fcfile<<pos;
				fcfile<<" &"<<fcfac<<" &"<<signal<<"\\\\"<<std::endl;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  fcfile<<"\\end{tabular}"<<std::endl;
  fcfile<<"\\end{table}"<<std::endl;
  fcfile<<"\\end{document}"<<std::endl;
  std::cout<<"end of calculation"<<std::endl;
  
  /* normalise */
  /*	normalise(sumlor); */
  /* output of fudge */
  /* S/N simulation */
  double snlimit=labellimit*(2./(M_PI*gamma));
  for (i=0;i<size;i++) {
    specf.precision(3);
    specf.setf(ios::fixed,ios::floatfield);
    specf<<grid(i)<<"  ";
    specf.precision(6);
    specf.setf(ios::scientific,ios::floatfield);
    if ((sumlor(i)) >= snlimit)
      specf<<sumlor(i)*(int00*int00)<<std::endl;
    else
      specf<<snlimit*(int00*int00)<<std::endl;
  }
  /* combinations */
  /* stat file */
  if (Units == 'A') std::cout<<"Geometries in Angstroms"<<std::endl;
  if (Units == 'B') std::cout<<"Geometries in Bohrs"<<std::endl;
  std::cout<<"frequency list"<<std::endl;
  std::cout<<"mode  neutral  ion"<<std::endl;
  std::cout<<" number of bands calculated: "<<fccount<<std::endl;
  double ZPEn=0.;
  double ZPEc=0.;
  for (i=0;i<N;i++) {
    std::cout<<i+1<<"      "<<omega(i)*hatocm<<"     "<<omegaion(i)*hatocm<<std::endl;
    ZPEn+=omega(i);
    ZPEc+=omegaion(i);
  }
  ZPEn*=(.5*hatocm);
  ZPEc*=(.5*hatocm);
  std::cout<<"ZPE neutral= "<<ZPEn<<std::endl;
  std::cout<<"ZPE cation= "<<ZPEc<<std::endl;
}

double position(int modei,int modej,int modek,
		int model,int modem,int moden,
		int cati,int catj,int catk,
		int neul,int neum,int neun,
		vector &omegaion,vector &omega)
{
  double pos=(
	      ( ((double) cati)*omegaion(modei)
		+((double) catj)*omegaion(modej)
		+((double) catk)*omegaion(modek)
		)
	      -(
		((double) neul)*omega(model)
		+((double) neum)*omega(modem)
		+((double) neun)*omega(moden)
		)
	      )*hatocm;
  return pos;
}
double oi0(int k,int n,double int00,matrix& Pmat,diagmat& II,vector& delta)
{
  complex t=tau(k,Pmat,II,delta);
  double o=int00/sqrt((double) factrl(n))
    *real(pow(complex(.5-Pmat(k,k),0.),(double) n/2.)
	  *HermiteH(n,t) );
  return o;
}
double o0i(int k,int n,double int00,matrix& Qmat,matrix& Rmat,vector& delta)
{
  complex s=sigma(k,Qmat,Rmat,delta);
  double o=pow(-1.,(double) n)*int00/sqrt((double) factrl(n))
    *real(pow(complex(.5-Qmat(k,k),0.),(double) n/2.)*HermiteH(n,s) );
  return o;
}
complex tau(int k,matrix& Pmat,diagmat& II,vector& delta)
{
  complex t= complex(((II-Pmat)*delta)(k),0.)/
    sqrt(complex(1.-2.*Pmat(k,k),0.));
  return t;
}
complex sigma(int k,matrix& Qmat,matrix& Rmat,vector& delta)
{
  complex s=complex((Rmat*delta)(k),0.)
    /sqrt(complex(1.-2.*Qmat(k,k),0.));
  return s;
}
int factrl(int n)
{
  int f=1;
  if (n==0) return 1;
  for (int i=1;i<=n;i++) f*=i;
  return f;
}
complex HermiteH(int n,complex x)
{
  if (n==0) return complex(1.,0.);
  if (n==1) return (complex(2.,0.)*x);
  
  complex val;
  val=complex(2.,0.)*x*HermiteH(n-1,x)-(complex(((double) (2*(n-1))),0.)
					*HermiteH(n-2,x));
  return val;
}
double  boltzman(double EoverT)
{
  double Kpercm=0.69509;
  EoverT*=hatocm; 
  /* au to 1/cm */
  return exp(-EoverT/Kpercm);
}
double lorentz(double w,double w0,double gamma)
{
  double PL=gamma/((2.*M_PI)*((w-w0)*(w-w0)+gamma*gamma/4.));
  return PL;
}
vector lor(double pos,double signal,int size,double gamma,vector grid)
{
  vector sp(size);
  for (int i=0;i<size;i++) sp(i)=signal*lorentz(grid(i),pos,gamma);
  return sp;
}
vector singlecold(int ni,int modei,matrix& Pmat,diagmat& II,vector& delta)
{
  vector t(ni);
  t(0)=int00;
  int i;
  /* single */
  for (i=1;i<ni;i++) {
    t(i)=oi0(modei,i,int00,Pmat,II,delta);
  }
  return t;
}
vector singlehot(int ni,int modei,matrix& Qmat,matrix &Rmat,vector& delta)
{
  vector t(ni);
  t(0)=int00;
  int i;
  /* single */
  for (i=1;i<ni;i++) {
    t(i)=o0i(modei,i,int00,Qmat,Rmat,delta);
  }
  return t;
}
vector doublecold(int ni,int nj,int modei,int modej,
	      matrix& Pmat,diagmat& II,vector& delta)
{
  vector t(ni*nj);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  t(0)=int00;
  int i,j;
  /* single */
  j=0;
  for (i=1;i<ni;i++) {
    t((i*nj+j))=oi0(modei,i,int00,Pmat,II,delta);
  }
  i=0;
  for (j=1;j<nj;j++) {
    t((i*nj+j))=oi0(modej,j,int00,Pmat,II,delta);
  }
  /* double */
  j=1;
  for (i=1;i<ni;i++)
    t((i*nj+j))=(m2Pmat(modej,modei))
      *sqrt( ((double) i) / ((double) j) )
      *t(((i-1)*nj+j-1))
      +sqrt(2./((double) j))*PD(modej)
      *t((i*nj+j-1));
  for (i=1;i<ni;i++) 
    for (j=2;j<nj;j++) {
      t((i*nj+j))=(m2Pmat(modej,modej))
	*sqrt( ((double) (j-1) ) / ((double) j) )
	*t((i*nj+j-2))
	+(m2Pmat(modej,modei))
	*sqrt( ((double) i)/ ((double) j) )
	*t(((i-1)*nj+j-1))
	+sqrt(2./((double) j))*PD(modej)
	*t((i*nj+j-1));
    }
  return t;
}
vector doublehot(int ni,int nj,int modei,int modej,
	      matrix& Qmat,matrix &Rmat,vector& delta)
{
  vector t(ni*nj);
  vector RD=((Rmat)*delta);
  t(0)=int00;
  int i,j;
  /* single */
  j=0;
  for (i=1;i<ni;i++) {
    t((i*nj+j))=o0i(modei,i,int00,Qmat,Rmat,delta);
  }
  i=0;
  for (j=1;j<nj;j++) {
    t((i*nj+j))=o0i(modej,j,int00,Qmat,Rmat,delta);
  }
  /* double */
  j=1;
  for (i=1;i<ni;i++) 
    t((i*nj+j))=(2.*Qmat(modej,modei))
      *sqrt( ((double) i) / ((double) j) )
      *t(((i-1)*nj+j-1))
      -sqrt(2./((double) j))*RD(modej)
      *t((i*nj+j-1));
  for (i=1;i<ni;i++) 
    for (j=2;j<nj;j++) {
      t((i*nj+j))=(2.*Qmat(modej,modej)-1.)
	*sqrt( ((double) (j-1) ) / ((double) j) )
	*t((i*nj+j-2))
	+(2.*Qmat(modej,modei))
	*sqrt( ((double) i)/ ((double) j) )
	*t(((i-1)*nj+j-1))
	-sqrt(2./((double) j))*RD(modej)
	*t((i*nj+j-1));
    }
  return t;
}
vector ijhot(int ni,int nj,int modei,int modej,matrix &Pmat,matrix &Qmat,
	     matrix &Rmat,diagmat &II,vector &delta)
{
  vector t(ni*nj);
  t(0)=int00;
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  int i,j;
  /* single */
  j=0;
  for (i=1;i<ni;i++) {
    t((i*nj+j))=oi0(modei,i,int00,Pmat,II,delta);
  }
  i=0;
  for (j=1;j<nj;j++) {
    t((i*nj+j))=o0i(modej,j,int00,Qmat,Rmat,delta);
  }
  i=1;
  for (j=1;j<nj;j++)  
    t((i*nj+j))=(m2Rmat(modej,modei))
      *sqrt( ((double) j)/ ((double) i) )
      *t(((i-1)*nj+j-1))
      +sqrt(2./((double) i))*PD(modei)
      *t(((i-1)*nj+j));
  for (j=1;j<nj;j++) 
    for (i=2;i<ni;i++) {
      t((i*nj+j))=(m2Pmat(modei,modei))
	*sqrt( ((double) (i-1) ) / ((double) i) )
	*t(((i-2)*nj+j))
	+(m2Rmat(modej,modei))
	*sqrt( ((double) j)/ ((double) i) )
	*t(((i-1)*nj+j-1))
	+sqrt(2./((double) i))*PD(modei)
	*t(((i-1)*nj+j));
    }
  return t;
}
vector ijkhot(int ni,int nj,int nk,int modei,int modej,int modek,
	      matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	      vector &delta)
{
  vector t(ni*nj*nk);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k;
  vector dcold=doublecold(ni,nj,modei,modej,Pmat,II,delta);
  k=0;
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) {
      t((i*nj+j)*nk+k)=dcold(i*nj+j);
    }
  vector sij=ijhot(ni,nk,modei,modek,Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++) 
    for (k=0;k<nk;k++) {
      t((i*nj+j)*nk+k)=sij(i*nk+k);
    }
  sij=ijhot(nj,nk,modej,modek,Pmat,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++) {
      t((i*nj+j)*nk+k)=sij(j*nk+k);
    }
  j=1;
  for (i=1;i<ni;i++)
    for (k=1;k<nk;k++) 
      t((i*nj+j)*nk+k)=
	//	(m2Rmat(modej,modek)) modified because typo in Ref.
	(m2Rmat(modek,modej))
	*sqrt( ((double) k) / ((double) j) )
	*t(((i)*nj+j-1)*nk+k-1)
	+(m2Pmat(modej,modei))
	*sqrt( ((double) i) / ((double) j) )
	*t(((i-1)*nj+j-1)*nk+k)
	+sqrt(2./((double) j))*PD(modej)
	*t((i*nj+j-1)*nk+k);
  for (i=1;i<ni;i++)
    for (k=1;k<nk;k++) 
      for (j=2;j<nj;j++) {
	t((i*nj+j)*nk+k)=
	  //	  (m2Rmat(modej,modek))  modified because typo in Ref.
	  (m2Rmat(modek,modej))
	  *sqrt( ((double) k) / ((double) j) )
	  *t(((i)*nj+j-1)*nk+k-1)
	  +(m2Pmat(modej,modei))
	  *sqrt( ((double) i) / ((double) j) )
	  *t(((i-1)*nj+j-1)*nk+k)
	  +(m2Pmat(modej,modej))
	  *sqrt( ((double) (j-1)) / ((double) j) )
	  *t(((i)*nj+j-2)*nk+k)
	  +sqrt(2./((double) j))*PD(modej)
	  *t((i*nj+j-1)*nk+k);
      }
  return t;
}
vector iklhot(int ni,int nk,int nl,int modei,int modek,int model,
	      matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	      vector &delta)
{
  /* type <i|kl> */
  vector t(ni*nl*nk);
  vector RD=((Rmat)*delta);
  matrix m2Qmat=(2.*Qmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,k,l;
  i=0;
  vector dhot=doublehot(nk,nl,modek,model,Qmat,Rmat,delta);
  for (k=0;k<nk;k++)
    for (l=0;l<nl;l++) {
      t((i*nk+k)*nl+l)=dhot(k*nl+l);
    }
  vector sij=ijhot(ni,nk,modei,modek,Pmat,Qmat,Rmat,II,delta);
  l=0;
  for (i=0;i<ni;i++) 
    for (k=0;k<nk;k++) {
      t((i*nk+k)*nl+l)=sij(i*nk+k);
    }
  sij=ijhot(ni,nl,modei,model,Pmat,Qmat,Rmat,II,delta);
  k=0;
  for (i=0;i<ni;i++) 
    for (l=0;l<nl;l++) {
      t((i*nk+k)*nl+l)=sij(i*nl+l);
    }
  l=1;
  for (i=1;i<ni;i++)
    for (k=1;k<nk;k++) 
      t((i*nk+k)*nl+l)=
	m2Rmat(model,modei)*sqrt( ((double) i)/((double) l))
	*t(((i-1)*nk+k)*nl+l-1)
	+m2Qmat(model,modek)*sqrt( ((double) k)/((double) l))
	*t(((i)*nk+k-1)*nl+l-1)
	-sqrt( 2./((double) l))*RD(model)
	*t(((i)*nk+k)*nl+l-1);
  for (i=1;i<ni;i++)
    for (k=1;k<nk;k++) 
      for (l=2;l<nl;l++) {
	t((i*nk+k)*nl+l)=
	  m2Rmat(model,modei)*sqrt( ((double) i)/((double) l))
	  *t(((i-1)*nk+k)*nl+l-1)
	  +m2Qmat(model,modek)*sqrt( ((double) k)/((double) l))
	  *t(((i)*nk+k-1)*nl+l-1)
	  +(m2Qmat(model,model))*sqrt( ((double) (l-1))/((double) l))
	  *t(((i)*nk+k)*nl+l-2)
	  -sqrt( 2./((double) l))*RD(model)
	  *t(((i)*nk+k)*nl+l-1);
      }
  return t;
}
vector ijklhot(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	       matrix &Pmat,matrix &Qmat,matrix &Rmat,diagmat &II,
	       vector &delta)
{
  /* type <ij|kl> */
  vector t(ni*nj*nk*nl);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l;

  vector sijk=ijkhot(ni,nj,nk,modei,modej,modek,Pmat,Qmat,Rmat,II,delta);
  l=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) {
	t(((i*nj+j)*nk+k)*nl+l)=sijk((i*nj+j)*nk+k);	
      }
  sijk=ijkhot(ni,nj,nl,modei,modej,model,Pmat,Qmat,Rmat,II,delta);
  k=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++)
      for (l=0;l<nl;l++) {
	  t(((i*nj+j)*nk+k)*nl+l)=sijk((i*nj+j)*nl+l);
	}
  vector sikl=iklhot(ni,nk,nl,modei,modek,model,Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++)
      for (l=0;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=sikl((i*nk+k)*nl+l);
      }
  sikl=iklhot(nj,nk,nl,modej,modek,model,Pmat,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++)
      for (l=0;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=sikl((j*nk+k)*nl+l);
      }  

  i=1;
  for (j=1;j<nj;j++) 
    for (k=1;k<nk;k++) 
      for (l=1;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=
	  m2Rmat(modei,modek)*sqrt(((double) k)/((double) i))*
	  t((((i-1)*nj+j)*nk+k-1)*nl+l)
	  +m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	  t((((i-1)*nj+j)*nk+k)*nl+l-1)
	  +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))*
	  t((((i-1)*nj+j-1)*nk+k)*nl+l)
	  +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
      }  
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++) 
	for (l=1;l<nl;l++) {
	  t(((i*nj+j)*nk+k)*nl+l)=
	    //	    m2Rmat(modei,modek)*sqrt(((double) k)/((double) i))*
	    m2Rmat(modek,modei)*sqrt(((double) k)/((double) i))* 
	    t((((i-1)*nj+j)*nk+k-1)*nl+l)
	    //	    +m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	    +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	    t((((i-1)*nj+j)*nk+k)*nl+l-1)
	    +m2Pmat(modei,modei)*sqrt(((double) (i-1))/((double) i))*
	    t((((i-2)*nj+j)*nk+k)*nl+l)
	    +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))*
	    t((((i-1)*nj+j-1)*nk+k)*nl+l)
	    +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
	}  
  return t;
}
vector triple(int ni,int nj,int nk,int modei,int modej,int modek,
	      matrix& Pmat,diagmat& II,vector& delta)
{
  vector t(ni*nj*nk);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  t(0)=int00;
  int i,j,k;
  vector dcold=doublecold(ni,nj,modei,modej,Pmat,II,delta);
  k=0;
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) {
      t((i*nj+j)*nk+k)=dcold(i*nj+j);
    }
  dcold=doublecold(ni,nk,modei,modek,Pmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++) {
      t((i*nj+j)*nk+k)=dcold(i*nk+k);
    }
  dcold=doublecold(nj,nk,modej,modek,Pmat,II,delta);
  i=0;
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++) {
      t((i*nj+j)*nk+k)=dcold(j*nk+k);
    }
  k=1;
  for (i=1;i<ni;i++)
    for (j=1;j<nj;j++)
      t((i*nj+j)*nk+k)=
	(m2Pmat(modek,modei))
	  *sqrt( ((double) i) / ((double) k) )
	  *t(((i-1)*nj+j)*nk+k-1)
	  +(m2Pmat(modek,modej))
	  *sqrt( ((double) j) / ((double) k) )
	  *t(((i)*nj+j-1)*nk+k-1)
	  +sqrt(2./((double) k))*PD(modek)
	  *t((i*nj+j)*nk+k-1);
  for (i=1;i<ni;i++)
    for (j=1;j<nj;j++)
      for (k=2;k<nk;k++) {
	t((i*nj+j)*nk+k)=
	  (m2Pmat(modek,modei))
	  *sqrt( ((double) i) / ((double) k) )
	  *t(((i-1)*nj+j)*nk+k-1)
	  +(m2Pmat(modek,modej))
	  *sqrt( ((double) j) / ((double) k) )
	  *t(((i)*nj+j-1)*nk+k-1)
	  +(m2Pmat(modek,modek))
	  *sqrt( ((double) k-1) / ((double) k) )
	  *t(((i)*nj+j)*nk+k-2)
	  +sqrt(2./((double) k))*PD(modek)
	  *t((i*nj+j)*nk+k-1);
      }
  return t;
}
vector triplehot(int ni,int nj,int nk,int modei,int modej,int modek,
	      matrix& Qmat,matrix &Rmat,diagmat &II,vector& delta)
{
  vector t(ni*nj*nk);
  vector RD=((Rmat)*delta);
  matrix m2Qmat=(2.*Qmat-II);
  t(0)=int00;
  int i,j,k;
  
  /* double */
  i=0;
  vector dhot=doublehot(nj,nk,modej,modek,Qmat,Rmat,delta);
  for (j=0;j<nj;j++)
    for (k=0;k<nk;k++)
      t((i*nj+j)*nk+k)=dhot(j*nk+k);
  j=0;
  dhot=doublehot(ni,nk,modei,modek,Qmat,Rmat,delta);
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++)
      t((i*nj+j)*nk+k)=dhot(i*nk+k);
  k=0;
  dhot=doublehot(ni,nj,modei,modej,Qmat,Rmat,delta);
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++)
      t((i*nj+j)*nk+k)=dhot(i*nj+j);
  i=1;
  for (j=1;j<nj;j++)
    for (k=1;k<nk;k++) {
      t((i*nj+j)*nk+k)=
	m2Qmat(modei,modej)*sqrt(((double) (j))/((double) i))
	*t(((i-1)*nj+j-1)*nk+k)
	+m2Qmat(modei,modek)*sqrt(((double) (k))/((double) i))
	*t(((i-1)*nj+j)*nk+k-1)
	-sqrt(2./((double) i))*RD(modei)*t(((i-1)*nj+j)*nk+k);
    }
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++)
      for (k=1;k<nk;k++) {
	t((i*nj+j)*nk+k)=
	  m2Qmat(modei,modei)*sqrt(((double) (i-1))/((double) i))
	  *t(((i-2)*nj+j)*nk+k)
	  +m2Qmat(modei,modej)*sqrt(((double) (j))/((double) i))
	  *t(((i-1)*nj+j-1)*nk+k)
	  +m2Qmat(modei,modek)*sqrt(((double) (k))/((double) i))
	  *t(((i-1)*nj+j)*nk+k-1)
	  -sqrt(2./((double) i))*RD(modei)*t(((i-1)*nj+j)*nk+k);
      }
  return t;
}
vector triplel(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta)
{
  vector t(ni*nj*nk*nl);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l;
  vector tcold=triple(ni,nj,nk,modei,modej,modek,Pmat,II,delta);
  l=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++)
      for (k=0;k<nk;k++) {
	t(((i*nj+j)*nk+k)*nl +l)=tcold((i*nj+j)*nk+k);
      }
  vector sijk=ijkhot(ni,nj,nl,modei,modej,model,Pmat,Qmat,Rmat,II,delta);
  k=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (l=0;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=sijk((i*nj+j)*nl+l);	
      }
  sijk=ijkhot(ni,nk,nl,modei,modek,model,Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=sijk((i*nk+k)*nl+l);	
      }
  sijk=ijkhot(nj,nk,nl,modej,modek,model,Pmat,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=sijk((j*nk+k)*nl+l);	
      }
  i=1;
  for (j=1;j<nj;j++) 
    for (k=1;k<nk;k++) 
      for (l=1;l<nl;l++) {
	t(((i*nj+j)*nk+k)*nl+l)=
	  //m2Rmat(modei,modek)*sqrt(((double) k)/((double) i))*
	  m2Pmat(modei,modek)*sqrt(((double) k)/((double) i))*
	  t((((i-1)*nj+j)*nk+k-1)*nl+l)
	  //	  +m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	  +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	  t((((i-1)*nj+j)*nk+k)*nl+l-1)
	  +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))*
	  t((((i-1)*nj+j-1)*nk+k)*nl+l)
	  +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
	}  
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++) 
	for (l=1;l<nl;l++) {
	  t(((i*nj+j)*nk+k)*nl+l)=
	    //	    m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	    m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	    t((((i-1)*nj+j)*nk+k)*nl+l-1)
	    +m2Pmat(modei,modei)*sqrt(((double) (i-1))/((double) i))*
	    t((((i-2)*nj+j)*nk+k)*nl+l)
	    +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))*
	    t((((i-1)*nj+j-1)*nk+k)*nl+l)
	    +m2Pmat(modei,modek)*sqrt(((double) k)/((double) i))*
	    t((((i-1)*nj+j)*nk+k-1)*nl+l)
	    +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
	}      
  return t;
}
vector triplelm(int ni,int nj,int nk,int nl,int nm,
	       int modei,int modej,int modek,int model,int modem,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta)
{
  vector t(ni*nj*nk*nl*nm);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l,m;
  /* depends on 
     <ijk|m>
     <ijk|l>
     <ij|lm> <ik|lm> <jk|lm>     */
  i=0; 
  vector oijkl=ijklhot(nj,nk,nl,nm,modej,modek,model,modem,
		       Pmat,Qmat,Rmat,II,delta);
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++)  
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    oijkl(((j*nk+k)*nl+l)*nm+m);  
  j=0; 
  oijkl=ijklhot(ni,nk,nl,nm,modei,modek,model,modem,
			       Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++) 
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++)  
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    oijkl(((i*nk+k)*nl+l)*nm+m);
  k=0; 
  oijkl=ijklhot(ni,nj,nl,nm,modei,modej,model,modem,
		Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) 
      for (l=0;l<nl;l++)  
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    oijkl(((i*nj+j)*nl+l)*nm+m);
  vector tijkl=triplel(ni,nj,nk,nl,modei,modej,modek,model,
			    Pmat,Qmat,Rmat,II,delta);
  m=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
	for (l=0;l<nl;l++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    tijkl(((i*nj+j)*nk+k)*nl+l);  
  tijkl=triplel(ni,nj,nk,nm,modei,modej,modek,modem,		    
		Pmat,Qmat,Rmat,II,delta);
  l=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    tijkl(((i*nj+j)*nk+k)*nm+m);  
  i=1;
  for (j=1;j<nj;j++) 
    for (k=1;k<nk;k++) 
      for (l=1;l<nl;l++) 
	for (m=1;m<nm;m++) {
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    //	    m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	    m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	    t(((((i-1)*nj+j)*nk+k)*nl+l-1)*nm+m)
	    //	    +m2Rmat(modei,modem)*sqrt(((double) m)/((double) i))*
	    +m2Rmat(modem,modei)*sqrt(((double) m)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m-1)
	    +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))
	    *t(((((i-1)*nj+j-1)*nk+k)*nl+l)*nm+m)
	    +m2Pmat(modei,modek)*sqrt(((double) k)/((double) i))
	    *t(((((i-1)*nj+j)*nk+k-1)*nl+l)*nm+m)
	    +sqrt(2./((double) i))*PD(modei)
	    *t((((i-1)*nj+j)*nk+k)*nl+l);
	}    
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++) 
	for (l=1;l<nl;l++) 
	  for (m=1;m<nm;m++) {
	    t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	      //	      m2Rmat(modei,model)*sqrt(((double) l)/((double) i))*
	      m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k)*nl+l-1)*nm+m)
	      //	      +m2Rmat(modei,modem)*sqrt(((double) m)/((double) i))*
	      +m2Rmat(modem,modei)*sqrt(((double) m)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m-1)
	      +m2Pmat(modei,modei)
	      *sqrt(((double) (i-1))/((double) i))
	      *t(((((i-2)*nj+j)*nk+k)*nl+l)*nm+m)
	      +m2Pmat(modei,modej)*sqrt(((double) j)/((double) i))
	      *t(((((i-1)*nj+j-1)*nk+k)*nl+l)*nm+m)
	      +m2Pmat(modei,modek)*sqrt(((double) k)/((double) i))
	      *t(((((i-1)*nj+j)*nk+k-1)*nl+l)*nm+m)
	      +sqrt(2./((double) i))*PD(modei)
	      *t((((i-1)*nj+j)*nk+k)*nl+l);
	  }    
  return t;
}
vector triplei(int ni,int nj,int nk,int nl,
	       int modei,int modej,int modek,int model,
	       matrix& Pmat,matrix &Qmat,matrix &Rmat,
	       diagmat& II,vector& delta)
{
  vector t(ni*nj*nk*nl);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l;
  /* <i|jkl> */
  /* depends on
     <0|jkl> i=0
     <i|jk>   l=0
     <i|jl>   k=0
     <i|kl>   j=0
     */
  vector thot=triplehot(nj,nk,nl,modej,modek,model,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++) 
	t(((i*nj+j)*nk+k)*nl +l)=thot((j*nk+k)*nl+l);
  vector sikl=iklhot(ni,nj,nk,modei,modej,modek,Pmat,Qmat,Rmat,II,delta);
  l=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++)
	t(((i*nj+j)*nk+k)*nl+l)=sikl((i*nj+j)*nk+k);	
  sikl=iklhot(ni,nk,nl,modei,modek,model,Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++) 
	t(((i*nj+j)*nk+k)*nl+l)=sikl((i*nk+k)*nl+l);	
  sikl=iklhot(ni,nj,nl,modei,modej,model,Pmat,Qmat,Rmat,II,delta);
  k=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (l=0;l<nl;l++)
	t(((i*nj+j)*nk+k)*nl+l)=sikl((i*nj+j)*nl+l);	
  i=1;
  for (j=1;j<nj;j++) 
    for (k=1;k<nk;k++) 
      for (l=1;l<nl;l++)
	t(((i*nj+j)*nk+k)*nl+l)=
	  m2Rmat(modej,modei)*sqrt(((double) j)/((double) i))* // modified
	  t((((i-1)*nj+j-1)*nk+k)*nl+l)
	  +m2Rmat(modek,modei)*sqrt(((double) k)/((double) i))*
	  t((((i-1)*nj+j)*nk+k-1)*nl+l)
	  +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	  t((((i-1)*nj+j)*nk+k)*nl+l-1)
	  +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++) 
	for (l=1;l<nl;l++)
	  t(((i*nj+j)*nk+k)*nl+l)=
	    m2Rmat(modej,modei)*sqrt(((double) j)/((double) i))*
	    t((((i-1)*nj+j-1)*nk+k)*nl+l)
	    +m2Rmat(modek,modei)*sqrt(((double) k)/((double) i))*
	    t((((i-1)*nj+j)*nk+k-1)*nl+l)
	    +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	    t((((i-1)*nj+j)*nk+k)*nl+l-1)
	    +m2Pmat(modei,modei)*sqrt(((double) (i-1))/((double) i))*
	    t((((i-2)*nj+j)*nk+k)*nl+l)
	    +sqrt(2./((double) i))*PD(modei)*t((((i-1)*nj+j)*nk+k)*nl+l);
  return t;
}
vector tripleij(int ni,int nj,int nk,int nl,int nm,
		int modei,int modej,int modek,int model,int modem,
		matrix& Pmat,matrix &Qmat,matrix &Rmat,
		diagmat& II,vector& delta)
{
  vector t(ni*nj*nk*nl*nm);
  vector PD=((II-Pmat)*delta);
  matrix m2Pmat=(2.*Pmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l,m;
  /* <ij|klm> */
  /* depends on 
     <ij|kl>   m=0
     <ij|km>   l=0
     <ij|lm>   k=0
     <i|klm>   j=0
     <j|klm>   i=0
     */
  m=0; 
  vector oijkl=ijklhot(ni,nj,nk,nl,modei,modej,modek,model,
			       Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
	for (l=0;l<nl;l++)  
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    oijkl(((i*nj+j)*nk+k)*nl+l);  
  l=0; 
  oijkl=ijklhot(ni,nj,nk,nm,modei,modej,modek,modem,
		Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	   oijkl(((i*nj+j)*nk+k)*nm+m);  
  k=0; 
  oijkl=ijklhot(ni,nj,nl,nm,modei,modej,model,modem,
		Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++) 
    for (j=0;j<nj;j++) 
      for (l=0;l<nl;l++)  
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    oijkl(((i*nj+j)*nl+l)*nm+m);
  vector tijkl=triplei(ni,nk,nl,nm,modei,modek,model,modem,
		       Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++) 
	for (m=0;m<nm;m++)
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    tijkl(((i*nk+k)*nl+l)*nm+m);

  tijkl=triplei(nj,nk,nl,nm,modej,modek,model,modem,		    
		Pmat,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++)
      for (l=0;l<nl;l++) 
	for (m=0;m<nm;m++) 
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    tijkl(((j*nk+k)*nl+l)*nm+m);  
  i=1;
  for (j=1;j<nj;j++) 
    for (k=1;k<nk;k++) 
      for (l=1;l<nl;l++) 
	for (m=1;m<nm;m++) {
	  t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	    m2Rmat(modek,modei)*sqrt(((double) k)/((double) i))*
	    t(((((i-1)*nj+j)*nk+k-1)*nl+l)*nm+m)
	    +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	    t(((((i-1)*nj+j)*nk+k)*nl+l-1)*nm+m)
	    +m2Rmat(modem,modei)*sqrt(((double) m)/((double) i))*
	    t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m-1)
	    +m2Pmat(modei,modej)*sqrt(((double) (j))/((double) i))*
	    t(((((i-1)*nj+j-1)*nk+k)*nl+l)*nm+m)
	    +sqrt(2./((double) i))*PD(modei)
	    *t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m);
	}    
  for (i=2;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++) 
	for (l=1;l<nl;l++) 
	  for (m=1;m<nm;m++) {
	    t((((i*nj+j)*nk+k)*nl+l)*nm+m)=
	      m2Rmat(modek,modei)*sqrt(((double) k)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k-1)*nl+l)*nm+m)
	      +m2Rmat(model,modei)*sqrt(((double) l)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k)*nl+l-1)*nm+m)
	      +m2Rmat(modem,modei)*sqrt(((double) m)/((double) i))*
	      t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m-1)
	      +m2Pmat(modei,modei)*sqrt(((double) (i-1))/((double) i))*
	      t(((((i-2)*nj+j)*nk+k)*nl+l)*nm+m)
	      +m2Pmat(modei,modej)*sqrt(((double) (j))/((double) i))*
	      t(((((i-1)*nj+j-1)*nk+k)*nl+l)*nm+m)
	      +sqrt(2./((double) i))*PD(modei)
	      *t(((((i-1)*nj+j)*nk+k)*nl+l)*nm+m);
	  }
  return t;
}
vector triplelmn(int ni,int nj,int nk,int nl,int nm,int nn,
		 int modei,int modej,int modek,int model,
		 int modem,int moden,
		 matrix& Pmat,matrix &Qmat,matrix &Rmat,
		 diagmat& II,vector& delta)
{
  vector t(ni*nj*nk*nl*nm*nn);
  vector RD=((Rmat)*delta);
  matrix m2Qmat=(2.*Qmat-II);
  matrix m2Rmat=2.*Rmat;
  t(0)=int00;
  int i,j,k,l,m,n;
  /* depends on
     <ijk|lm> n=0
     <ijk|nm> l=0
     <ijk|ln> m=0
     <ij|lmn> k=0
     <ik|lmn> j=0
     <jk|lmn> i=0
     */
  vector tijklm=triplelm(ni,nj,nk,nl,nm,
			modei,modej,modek,model,modem,
			Pmat,Qmat,Rmat,II,delta);
  n=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
	for (l=0;l<nl;l++) 
	  for (m=0;m<nm;m++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijklm((((i*nj+j)*nk+k)*nl+l)*nm+m);
  tijklm=triplelm(ni,nj,nk,nl,nn,
		  modei,modej,modek,model,moden,
		  Pmat,Qmat,Rmat,II,delta);
  m=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
        for (l=0;l<nl;l++) 
          for (n=0;n<nn;n++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijklm((((i*nj+j)*nk+k)*nl+l)*nn+n);
  l=0;
  tijklm=triplelm(ni,nj,nk,nm,nn,
		  modei,modej,modek,modem,moden,
		  Pmat,Qmat,Rmat,II,delta);
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++) 
      for (k=0;k<nk;k++) 
        for (m=0;m<nm;m++) 
          for (n=0;n<nn;n++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijklm((((i*nj+j)*nk+k)*nm+m)*nn+n);  
  vector tijlmn=tripleij(nj,nk,nl,nm,nn,
			 modej,modek,model,modem,moden,
			 Pmat,Qmat,Rmat,II,delta);
  i=0;
  for (j=0;j<nj;j++) 
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++)
	for (m=0;m<nm;m++) 
	  for (n=0;n<nn;n++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijlmn((((j*nk+k)*nl+l)*nm+m)*nn+n);
  tijlmn=tripleij(ni,nk,nl,nm,nn,
		  modei,modek,model,modem,moden,
		  Pmat,Qmat,Rmat,II,delta);
  j=0;
  for (i=0;i<ni;i++)
    for (k=0;k<nk;k++) 
      for (l=0;l<nl;l++)
	for (m=0;m<nm;m++) 
	  for (n=0;n<nn;n++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijlmn((((i*nk+k)*nl+l)*nm+m)*nn+n);
  tijlmn=tripleij(ni,nj,nl,nm,nn,
		  modei,modej,model,modem,moden,
		  Pmat,Qmat,Rmat,II,delta);
  k=0;
  for (i=0;i<ni;i++)
    for (j=0;j<nj;j++)
      for (l=0;l<nl;l++)
	for (m=0;m<nm;m++) 
	  for (n=0;n<nn;n++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      tijlmn((((i*nj+j)*nl+l)*nm+m)*nn+n);
  n=1;
  for (i=1;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++)
	for (l=1;l<nl;l++) 
	  for (m=1;m<nm;m++)
	    t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
	      m2Rmat(moden,modei)*sqrt(((double) (i))/((double) n))
	      *t((((((i-1)*nj+j)*nk+k)*nl+l)*nm+m)*nn+n-1)
	      +m2Rmat(moden,modej)*sqrt(((double) (j))/((double) n))
	      *t(((((i*nj+j-1)*nk+k)*nl+l)*nm+m)*nn+n-1)
	      +m2Rmat(moden,modek)*sqrt(((double) (k))/((double) n))
	      *t(((((i*nj+j)*nk+k-1)*nl+l)*nm+m)*nn+n-1)
	      +m2Qmat(moden,model)*sqrt(((double) (l))/((double) n))
	      *t(((((i*nj+j)*nk+k)*nl+l-1)*nm+m)*nn+n-1)
	      +m2Qmat(moden,modem)*sqrt(((double) (m))/((double) n))
	      *t(((((i*nj+j)*nk+k)*nl+l)*nm+m-1)*nn+n-1)
	      -sqrt(2./((double) n))*RD(moden)
	      *t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n-1);  
  for (i=1;i<ni;i++)
    for (j=1;j<nj;j++) 
      for (k=1;k<nk;k++)
	for (l=1;l<nl;l++) 
	  for (m=1;m<nm;m++) 	    
	    for (n=2;n<nn;n++)
	      t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n)=
		m2Rmat(moden,modei)*sqrt(((double) (i))/((double) n))
		*t((((((i-1)*nj+j)*nk+k)*nl+l)*nm+m)*nn+n-1)
		+m2Rmat(moden,modej)*sqrt(((double) (j))/((double) n))
		*t(((((i*nj+j-1)*nk+k)*nl+l)*nm+m)*nn+n-1)
		+m2Rmat(moden,modek)*sqrt(((double) (k))/((double) n))
		*t(((((i*nj+j)*nk+k-1)*nl+l)*nm+m)*nn+n-1)
		+m2Qmat(moden,model)*sqrt(((double) (l))/((double) n))
		*t(((((i*nj+j)*nk+k)*nl+l-1)*nm+m)*nn+n-1)
		+m2Qmat(moden,modem)*sqrt(((double) (m))/((double) n))
		*t(((((i*nj+j)*nk+k)*nl+l)*nm+m-1)*nn+n-1)
		+m2Qmat(moden,moden)
		*sqrt(((double) (n-1))/((double) n))
		*t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n-2)
		-sqrt(2./((double) n))*RD(moden)
		*t(((((i*nj+j)*nk+k)*nl+l)*nm+m)*nn+n-1);	      
  return t;
}
