#ifndef  _FAI_2DFSDESOLVER_H__
#define _FAI_2DFSDESOLVER_H__
#include<cmath>
#include <iostream>
#include"kiss_fft.h"
#define  complex kiss_fft_cpx 
typedef unsigned int uint_32;
/*..............................*/


/*declaration of functions*/
complex cpmult(complex,complex);
complex cpminus(complex,complex);
complex cpdivid(complex,complex);
complex cpplus(complex,complex);
complex cpdividr(complex,double);
complex cpmultr(complex,double);
complex rdividcp(double,complex);
double cpmodule(complex z);
double cpmodule_sq(complex z);
double postv_func(double x,double y);
void Time_frac_diffusion_2Deq_solver(uint_32 N,uint_32 level,double X_L,double X_R,double Y_Low,\
										  double Y_Upp,double* Frac_div_appro_coeff,complex** rhs,double eps,double iter_num_arr[1],double tol);
double** offlaplacex_generate(uint_32* M_array,uint_32 M_array_len,double X_L,double X_R,double Y_Low, double Y_Upp);
double** offlaplacey_generate(uint_32* M_array,uint_32 M_array_len,double X_L,double X_R,double Y_Low, double Y_Upp);
complex** mainlaplace_generate(uint_32* M_array,uint_32 M_array_len,double** offlaplacex,double** offlaplacey,double X_L,double X_R,double Y_Low, double Y_Upp);
uint_32 Vcycle_BLTDTDB(uint_32* M_array,uint_32 M_agrray_len,double** offlaplacex,double** offlaplacey,complex** mainlaplace,complex deltaFTFC,complex* rhs,complex* sol,double tol);
double Infnorm(complex* vec,uint_32 vec_len);
void get_residual(double* offlaplacex, double* offlaplacey,complex* mainlaplace,complex* right_hand_side,complex* sol,uint_32 M,complex* residual);
void In_Cplt_LU_Sol(double* offlaplacex,double*offlaplacey,complex* mainlaplace, complex* rhs,complex* sol);
void Restr_Operator(complex* rh,complex* f2h,uint_32 M);
void Interp_Operator(complex*u2h,complex*uh,uint_32 M);
void tracking(double* offlaplacex,complex* mainlaplace ,complex* x,complex* b,uint_32 size);
void AXLGS_smoother(double* offlaplacex,double* offlaplacey,complex* mainlaplace,complex* initialgs,complex* rhs,uint_32 M,uint_32 time,bool iszeroinitial);
double _Pow_int(double x, int n);
/*.....................................................*/


/*definition of functions*/
double _Pow_int(double x, int n) {
	double answer = 1.;
	for (int i=0; i<n; ++i) {
		answer *= x;
	}
	return answer;
}

void Time_frac_diffusion_2Deq_solver(uint_32 N,uint_32 level,double X_L,double X_R,double Y_Low,\
										  double Y_Upp,double* Frac_div_appro_coeff,complex** rhs,double eps,double iter_num_arr[1],double tol)
{
	if(level<2)
	{
		std::cout<<"FAIsolver warning: level="<<level<<"is too small "<<std::endl;
	}
	if (tol<=0)
	{
		std::cout<<"FAI warning: the tolerance must be positive number!!!"<<std::endl;
	}
	double oneoverN=1.0/(double)N;
	double delta=pow(eps,oneoverN);
	double* Ddelta=new double[N];
	complex* fin=new complex[N],*out=new complex[N];
	
	uint_32 i;
	for (i=0;i<N;i++)
	{
		Ddelta[i]=_Pow_int(delta,i);
	}
	
	complex *FTFDAC=new complex[N];//fourier transformation of discrete Fractional derivative coefficients
	for (i=0;i<N;i++)
	{
		FTFDAC[i].r=Frac_div_appro_coeff[i]*Ddelta[i];
		FTFDAC[i].i=0;
	}

	//set up cfg for ffts later
	kiss_fft_cfg cfg;
	cfg=kiss_fft_alloc(N,0,NULL,NULL);

	kiss_fft(cfg,FTFDAC,FTFDAC);//do fft for FTFDAC

	for (i=0;i<N;i++)
	{
		FTFDAC[i].i=-FTFDAC[i].i;
		fin[i].i=0;
	}
	
	uint_32 M_array_len=level-1;
	uint_32* M_array=new uint_32[M_array_len];
	for (i=level;i>1;i--)
	{
		M_array[level-i]=_Pow_int(2,i)-1;
	}
	
	uint_32 j,Msquare=M_array[0]*M_array[0];
	// do fft for the right hand side
	for (i=0;i<Msquare;i++)
	{
		for (j=0;j<N;j++)
		{
			fin[j].r=rhs[j][i].r*Ddelta[j];
		}
		kiss_fft(cfg,fin,out);
		for (j=0;j<N;j++)
		{
			rhs[j][i].r=oneoverN*out[j].r;
			rhs[j][i].i=-oneoverN*out[j].i;
		}
	}
	delete[]out;
	double** offlaplacex,**offlaplacey;
	complex** mainlaplace;
	
	offlaplacex=offlaplacex_generate(M_array,M_array_len,X_L,X_R,Y_Low,Y_Upp);
	offlaplacey=offlaplacey_generate(M_array,M_array_len,X_L,X_R,Y_Low,Y_Upp);
	mainlaplace=mainlaplace_generate(M_array,M_array_len,offlaplacex,offlaplacey,X_L,X_R,Y_Low,Y_Upp);
	complex deltaFTFC=FTFDAC[0];
	complex* initialsolver=new complex[Msquare];
	iter_num_arr[0]=0;
	iter_num_arr[0]+=Vcycle_BLTDTDB(M_array,M_array_len,offlaplacex,offlaplacey,mainlaplace,deltaFTFC,rhs[0],initialsolver,tol);
	delete[]rhs[0];rhs[0]=initialsolver;initialsolver=NULL;
	for (i=1;i<N;i++)
	{
		deltaFTFC.i=FTFDAC[i].i;
		deltaFTFC.r=FTFDAC[i].r-FTFDAC[i-1].r;
		initialsolver=new complex[Msquare];
		iter_num_arr[0]+=Vcycle_BLTDTDB(M_array,M_array_len,offlaplacex,offlaplacey,mainlaplace,deltaFTFC,rhs[i],initialsolver,tol);
		delete[]rhs[i];rhs[i]=initialsolver;initialsolver=NULL;
	}
	delete[]FTFDAC;
	for (i=0;i<M_array_len;i++)
	{
		delete[]offlaplacex[i];
		delete[]offlaplacey[i];
		delete[]mainlaplace[i];
	}
	delete[]offlaplacex;
	delete[]offlaplacey;
	delete[]mainlaplace;
	for (j=0;j<Msquare;j++)
	{
		for (i=0;i<N;i++)
		{
			fin[i]=rhs[i][j];
		}
		kiss_fft(cfg,fin,fin);
			
		for (i=0;i<N;i++)
		{
			rhs[i][j].r=fin[i].r/Ddelta[i];
		}
	}
	delete[]fin;
	kiss_fft_cleanup();
	kiss_fft_free(cfg);
	delete[]Ddelta;
	iter_num_arr[0]/=(double)N;
}
uint_32 Vcycle_BLTDTDB(uint_32* M_array,uint_32 M_agrray_len,double** offlaplacex,double** offlaplacey,complex** mainlaplace,complex deltaFTFC,complex* rhs,complex* sol,double tol)
{
	uint_32 i,j,tempMs,Ms;
	for (j=0;j<M_agrray_len;j++)
	{
		tempMs=M_array[j]*M_array[j];
		for (i=0;i<tempMs;i++)
		{
			mainlaplace[j][i].i=deltaFTFC.i;
			mainlaplace[j][i].r+=deltaFTFC.r;
		}
	}

	if (M_agrray_len<=1)
	{
		In_Cplt_LU_Sol(offlaplacex[0],offlaplacey[0],mainlaplace[0],rhs,sol);
		return 0;
	}
	Ms=M_array[0]*M_array[0];
	
	double normrhs=Infnorm(rhs,Ms),curerror;
	if(normrhs==0)
	{
		return 0;
	}
	
	complex* residual;
	complex** solutions=new complex*[M_agrray_len-1];
	complex** rhsides=new complex*[M_agrray_len-1];
	uint_32 timepre=3,timepost=3,Marrlenm1=M_agrray_len-1,itertime=0;
	bool iszeroinital=0;
	for (i=0;i<Marrlenm1;i++)
	{
		tempMs=M_array[i+1]*M_array[i+1];
		solutions[i]=new complex[tempMs];
		rhsides[i]=new complex[tempMs];
	}
	while(true)
	{	
		//pre-Vcycle
		iszeroinital=(itertime==0?1:0);
		AXLGS_smoother(offlaplacex[0],offlaplacey[0],mainlaplace[0],sol,rhs,M_array[0],timepre,iszeroinital);
		residual=new complex[M_array[0]*M_array[0]];
		get_residual(offlaplacex[0],offlaplacey[0],mainlaplace[0],rhs,sol,M_array[0],residual);
		curerror=Infnorm(residual,Ms)/normrhs;
		if (curerror<tol)
		{
			break;
		}
		
		for (j=1;j<Marrlenm1;j++)
		{
			//restriction
			Restr_Operator(residual,rhsides[j-1],M_array[j-1]);
			//presmoothing
			AXLGS_smoother(offlaplacex[j],offlaplacey[j],mainlaplace[j],solutions[j-1],rhsides[j-1],M_array[j],timepre,1);
			delete[]residual;
			residual=new complex[M_array[j]*M_array[j]];
			get_residual(offlaplacex[j],offlaplacey[j],mainlaplace[j],rhsides[j-1],solutions[j-1],M_array[j],residual);
			
		}
	
		j=Marrlenm1;
		Restr_Operator(residual,rhsides[j-1],M_array[j-1]);
		delete[]residual;

		//exactly solve at the coarsest grid
		In_Cplt_LU_Sol(offlaplacex[j],offlaplacey[j],mainlaplace[j],rhsides[j-1],solutions[j-1]);

		//post-Vcycle
		for (j=Marrlenm1-1;j>0;j--)
		{
			//interpolation
			Interp_Operator(solutions[j],solutions[j-1],M_array[j+1]);
			//smoothing
			AXLGS_smoother(offlaplacex[j],offlaplacey[j],mainlaplace[j],solutions[j-1],rhsides[j-1],M_array[j],timepost,0);
		}
	
		Interp_Operator(solutions[0],sol,M_array[1]);
		//smoothing on the finest grid
		
		AXLGS_smoother(offlaplacex[0],offlaplacey[0],mainlaplace[0],sol,rhs,M_array[0],timepost,0);
		itertime++;
	}
    delete[]residual;
	for (i=0;i<Marrlenm1;i++)
	{
		delete[]solutions[i];
		delete[]rhsides[i];
	}
	delete[]solutions;
	delete[]rhsides;
	return itertime;
}
double** offlaplacex_generate(uint_32* M_array,uint_32 M_array_len,double X_L,double X_R,double Y_Low, double Y_Upp)
{
	double** offlaplacex=new double*[M_array_len];
	double h1, oneovreh1s,h2,x,y;
	uint_32 m,i,Mp1,Mm1,j,inds;
	for (m=0;m<M_array_len;m++)
	{
		
		//parameter settings
		Mp1=M_array[m]+1;
		h1=(X_R-X_L)/(double)(Mp1);
		oneovreh1s=1.0/(h1*h1);
		h2=(Y_Upp-Y_Low)/(double)(Mp1);
		Mm1=M_array[m]-1;
		
		offlaplacex[m]=new double[Mm1*M_array[m]];

		x=X_L+0.5*h1;
		
		for (j=1;j<Mp1;j++)
		{
			inds=(j-1)*Mm1;
			y=Y_Low+j*h2;
			for (i=0;i<Mm1;i++)
			{
				offlaplacex[m][inds+i]=oneovreh1s*postv_func(x+(i+1)*h1,y);
			}
		}
		
	}
	return offlaplacex;
}

double** offlaplacey_generate(uint_32* M_array,uint_32 M_array_len,double X_L,double X_R,double Y_Low, double Y_Upp)
{
	double** offlaplacey=new double*[M_array_len];
	double h1, oneovreh2s,h2,y,ys;
	uint_32 m,i,Mp1,Mm1,j,M,inds;
	for (m=0;m<M_array_len;m++)
	{
		//parameter settings
		M=M_array[m];
		Mp1=M+1;
		h1=(X_R-X_L)/(double)(Mp1);
		h2=(Y_Upp-Y_Low)/(double)(Mp1);
		oneovreh2s=1.0/(h2*h2);
		Mm1=M-1;
		offlaplacey[m]=new double[Mm1*M];
		y=Y_Low+0.5*h2;
		for (j=0;j<Mm1;j++)
		{
			inds=j*M-1;
			ys=y+(j+1)*h2;
			for (i=1;i<Mp1;i++)
			{
				offlaplacey[m][inds+i]=oneovreh2s*postv_func(X_L+i*h1,ys);
			}
		}
	}

	return offlaplacey;
}

complex** mainlaplace_generate(uint_32* M_array,uint_32 M_array_len,double** offlaplacex,double** offlaplacey,double X_L,double X_R,double Y_Low, double Y_Upp)
{
	complex** mainlaplace=new complex*[M_array_len];
	double h1, oneovreh2s,h2,oneovreh1s,ypm,xp1,xp05,xpm,yp1,yp05,xpm05,ypm05;
	uint_32 m,i,Mp1,Mm1,j,M,indsy,Mm2,indsypre,indsx;
	for (m=0;m<M_array_len;m++)
	{
		//parameter settings
		M=M_array[m];
		Mp1=M+1;
		h1=(X_R-X_L)/(double)(Mp1);
		h2=(Y_Upp-Y_Low)/(double)(Mp1);
		oneovreh2s=1.0/(h2*h2);
		oneovreh1s=1.0/(h1*h1);
		Mm1=M-1;
		mainlaplace[m]=new complex[M*M];
		Mm2=M-2;
		xp1=X_L+h1;
		yp1=Y_Low+h2;
		xp05=X_L+0.5*h1;
		yp05=Y_Low+0.5*h2;
		xpm=X_L+M*h1;
		ypm=Y_Low+M*h2;
		xpm05=X_L+(M+0.5)*h1;
		ypm05=Y_Low+(M+0.5)*h2;
		mainlaplace[m][0].r=offlaplacex[m][0]+offlaplacey[m][0]+oneovreh1s*postv_func(xp05,yp1)+oneovreh2s*postv_func(xp1,yp05);
		for (i=1;i<Mm1;i++)
		{
			mainlaplace[m][i].r=offlaplacex[m][i-1]+offlaplacex[m][i]+offlaplacey[m][i]+oneovreh2s*postv_func(xp1+i*h1,yp05);
		}
		mainlaplace[m][Mm1].r=offlaplacex[m][Mm2]+offlaplacey[m][Mm1]+oneovreh1s*postv_func(xpm05,yp1)+oneovreh2s*postv_func(xpm,yp05);

		for (j=1;j<Mm1;j++)
		{
			indsy=j*M;
			indsypre=(j-1)*M;
			indsx=j*Mm1;
			mainlaplace[m][indsy].r=offlaplacey[m][indsy]+offlaplacey[m][indsypre]+offlaplacex[m][indsx]+oneovreh1s*postv_func(xp05,yp1+j*h2);
			for (i=1;i<Mm1;i++)
			{
				mainlaplace[m][indsy+i].r=offlaplacey[m][indsy+i]+offlaplacey[m][indsypre+i]+offlaplacex[m][indsx+i]+offlaplacex[m][indsx+i-1];
			}
			mainlaplace[m][indsy+Mm1].r=offlaplacey[m][indsy+Mm1]+offlaplacey[m][indsypre+Mm1]+offlaplacex[m][indsx+Mm2]+oneovreh1s*postv_func(xpm05,yp1+j*h2);
		}

		indsypre=Mm2*M;
		indsx=Mm1*Mm1;
		indsy=Mm1*M;
		mainlaplace[m][indsy].r=offlaplacex[m][indsx]+oneovreh1s*postv_func(xp05,ypm)+offlaplacey[m][indsypre]+oneovreh2s*postv_func(xp1,ypm05);
		for (i=1;i<Mm1;i++)
		{
			mainlaplace[m][indsy+i].r=offlaplacex[m][indsx+i]+offlaplacex[m][indsx+i-1]+offlaplacey[m][indsypre+i]+oneovreh2s*postv_func(xp1+i*h1,ypm05);
		}
		mainlaplace[m][indsy+Mm1].r=offlaplacex[m][indsx+Mm2]+offlaplacey[m][indsypre+Mm1]+oneovreh1s*postv_func(xpm05,ypm)+oneovreh2s*postv_func(xpm,ypm05);
	}
	return mainlaplace;
}
inline complex cpmult(complex z1,complex z2)
{
	complex z;
	z.r=z1.r*z2.r-z1.i*z2.i;
	z.i=z1.i*z2.r+z1.r*z2.i;
	return z;
}
inline complex cpminus(complex z1,complex z2)
{
	complex z;
	z.r=z1.r-z2.r;
	z.i=z1.i-z2.i;
	return z;
}
inline complex cpplus(complex z1,complex z2)
{
	complex z;
	z.r=z1.r+z2.r;
	z.i=z1.i+z2.i;
	return z;
}
inline complex cpdivid(complex z1,complex z2)
{
	complex z;
	double module;
	module=_Pow_int(z2.r,2)+_Pow_int(z2.i,2);
	z.r=(z1.r*z2.r+z1.i*z2.i)/module;
	z.i=(z1.i*z2.r-z1.r*z2.i)/module;
	return z;
}
inline complex cpdividr(complex z,double r)
{
	z.r=z.r/r;
	z.i=z.i/r;
	return z;
}
inline complex rdividcp(double r,complex z2)
{
	complex z;
	double module;
	module=_Pow_int(z2.r,2)+_Pow_int(z2.i,2);
	z.r=r*z2.r/module;
	z.i=-r*z2.i/module;
	return z;
}
inline complex cpmultr(complex z,double r)
{
	z.r*=r;
	z.i*=r;
	return z;
}
inline double cpmodule(complex z)
{
	double m;
	m=sqrt(z.r*z.r+z.i*z.i);
	return m;
}

inline double cpmodule_sq(complex z)
{
	return z.r*z.r+z.i*z.i;
}

inline double postv_func(double x,double y)
{
	return exp(x*y);
}
double Infnorm(complex* vec,uint_32 vec_len)
{
	double max_mod=cpmodule_sq(vec[0]),cur_mod;
	for (uint_32 i=1;i<vec_len;i++)
	{
		cur_mod=cpmodule_sq(vec[i]);
		if (max_mod<cur_mod)
		{
			max_mod=cur_mod;
		}
	}
	return sqrt(max_mod);
}

void get_residual(double* offlaplacex, double* offlaplacey,complex* mainlaplace,complex* rhs,complex* sol,uint_32 M,complex* residual)
{
	uint_32 i,j,indy,indypre,indx,indnexts,Mm1=M-1,Mm2=M-2;
	indnexts=M;
	residual[0]=cpplus(rhs[0],cpminus(cpplus(cpmultr(sol[1],offlaplacex[0]),\
		cpmultr(sol[indnexts],offlaplacey[0])),cpmult(sol[0],mainlaplace[0])));
	for (i=1;i<Mm1;i++)
	{
		residual[i]=cpplus(rhs[i],cpminus(cpplus(cpplus(cpmultr(sol[i-1],offlaplacex[i-1]),\
			cpmultr(sol[i+1],offlaplacex[i])),cpmultr(sol[indnexts+i],offlaplacey[i])),\
			cpmult(mainlaplace[i],sol[i])));
	}
	residual[Mm1]=cpplus(rhs[Mm1],cpminus(cpplus(cpmultr(sol[Mm2],offlaplacex[Mm2]),\
		cpmultr(sol[indnexts+Mm1],offlaplacey[Mm1])),cpmult(sol[Mm1],mainlaplace[Mm1])));
	for (j=1;j<Mm1;j++)
	{
		indypre=(j-1)*M;
		indy=j*M;
		indx=j*Mm1;
		indnexts=(j+1)*M;
		residual[indy]=cpplus(rhs[indy],cpminus(cpplus(cpplus(cpmultr(sol[indypre],\
			offlaplacey[indypre]),cpmultr(sol[indnexts],offlaplacey[indy])),\
			cpmultr(sol[indy+1],offlaplacex[indx])),cpmult(sol[indy],mainlaplace[indy])));
		for (i=1;i<Mm1;i++)
		{
			residual[indy+i]=cpplus(rhs[indy+i],cpminus(cpplus(cpmultr(sol[indnexts+i],\
				offlaplacey[indy+i]),cpplus(cpmultr(sol[indypre+i],offlaplacey[indypre+i]),\
				cpplus(cpmultr(sol[indy+i-1],offlaplacex[indx+i-1]),cpmultr(sol[indy+i+1],\
				offlaplacex[indx+i])))),cpmult(sol[indy+i],mainlaplace[indy+i])));
		}

		residual[indy+Mm1]=cpplus(rhs[indy+Mm1],cpminus(cpplus(cpmultr(sol[indy+Mm2],\
			offlaplacex[indx+Mm2]),cpplus(cpmultr(sol[indypre+Mm1],offlaplacey[indypre+Mm1]),\
			cpmultr(sol[indnexts+Mm1],offlaplacey[indy+Mm1]))),cpmult(sol[indy+Mm1],\
			mainlaplace[indy+Mm1])));
	}

	indypre=Mm2*M;
	indy=Mm1*M;
	indx=Mm1*Mm1;
	residual[indy]=cpplus(rhs[indy],cpminus(cpplus(cpmultr(sol[indy+1],offlaplacex[indx]),\
		cpmultr(sol[indypre],offlaplacey[indypre])),cpmult(sol[indy],mainlaplace[indy])));
	for (i=1;i<Mm1;i++)
	{
		residual[indy+i]=cpplus(rhs[indy+i],cpminus(cpplus(cpmultr(sol[indypre+i],offlaplacey[indypre+i]),\
			cpplus(cpmultr(sol[indy+i-1],offlaplacex[indx+i-1]),cpmultr(sol[indy+i+1],\
			offlaplacex[indx+i]))),cpmult(sol[indy+i],mainlaplace[indy+i])));
	}

	residual[indy+Mm1]=cpplus(rhs[indy+Mm1],cpminus(cpplus(cpmultr(sol[indy+Mm2],\
		offlaplacex[indx+Mm2]),cpmultr(sol[indypre+Mm1],offlaplacey[indypre+Mm1])),cpmult(sol[indy+Mm1],\
		mainlaplace[indy+Mm1])));
}
inline void In_Cplt_LU_Sol(double* offlaplacex,double*offlaplacey,complex* mainlaplace, complex* rhs,complex* sol)
{
	//use the direct solver "Incomplete LU factorization" to solve the linear system derived from the coarsest grid 
	int i,k,r,ub,cc,ind;
	complex**Ui_ipk=new complex*[4];
	complex**Lipk_i=new complex*[3];
	complex EofA[4],t;
	EofA[0]=mainlaplace[0];EofA[1].r=-offlaplacex[0];EofA[1].i=0;EofA[2].r=0;EofA[2].i=0;
	EofA[3].r=-offlaplacey[0];EofA[3].i=0;
	Ui_ipk[0]=new complex[9];
	for(i=1;i<4;i++)
	{
		Ui_ipk[i]=new complex[9-i];
		Lipk_i[i-1]=new complex[9-i];
	}

	//compute L_i_0 and U_o_i
	Ui_ipk[0][0]=EofA[0];
	for(i=1;i<4;i++)
	{
		Ui_ipk[i][0]=EofA[i];
		Lipk_i[i-1][0]=cpdivid(EofA[i],Ui_ipk[0][0]);
	}
	//compute Ur_i and Li_r
	for(r=1;r<9;r++)
	{
		ub=std::min(r+3,8)+1;
		EofA[0]=mainlaplace[r];
		if (r<6)
		{
			EofA[3].r=-offlaplacey[r];
		}
		cc=(r+1)%3;
		if(cc!=0)
		{
			ind=2*(r+1-cc)/3+cc-1;
			EofA[1].r=-offlaplacex[ind];
		}
		else
		{
			EofA[1].r=0;
		}
		for(i=r;i<ub;i++)
		{
			t.r=0;
			t.i=0;
			for(k=std::max(i-3,0);k<r;k++)
				t=cpminus(t,cpmult(Lipk_i[r-k-1][k],Ui_ipk[i-k][k]));
			Ui_ipk[i-r][r]=cpplus(EofA[i-r],t);
		}

		for(i=r+1;i<ub;i++)
		{
			t.r=0;
			t.i=0;
			for(k=std::max(i-3,0);k<r;k++)
				t=cpminus(t,cpmult(Lipk_i[i-k-1][k],Ui_ipk[r-k][k]));
			Lipk_i[i-r-1][r]=cpplus(EofA[i-r],t);
			Lipk_i[i-r-1][r]=cpdivid(Lipk_i[i-r-1][r],Ui_ipk[0][r]);
		}
	}

	if (sol==NULL)
	{
		sol=new complex[9];
	}
	//solve Ly=f
	sol[0]=rhs[0];
	for(i=1;i<9;i++)
	{
		t.r=0;
		t.i=0;
		for(k=std::max(i-3,0);k<i;k++)
			t=cpminus(t,cpmult(Lipk_i[i-k-1][k],sol[k]));
		sol[i]=cpplus(rhs[i],t);
	}
	//solve Ux=y
	sol[8]=cpdivid(sol[8],Ui_ipk[0][8]);
	for(i=7;i>-1;i--)
	{
		t.r=0;
		t.i=0;
		ub=std::min(i+4,9);
		for(k=std::min(i+1,8);k<ub;k++)
			t=cpminus(t,cpmult(Ui_ipk[k-i][i],sol[k]));
		sol[i]=cpplus(sol[i],t);
		sol[i]=cpdivid(sol[i],Ui_ipk[0][i]);
	}
	delete[]Ui_ipk[0];
	for(i=1;i<4;i++)
	{
		delete[]Ui_ipk[i];
		delete[]Lipk_i[i-1];
	}
	delete[]Ui_ipk;
	delete[]Lipk_i;
}

inline void Restr_Operator(complex* rh,complex* f2h,uint_32 M)
{
	uint_32 i,j,nM,ti,tip1,tip2,indh,ind2h,indhp1,indhp2;
	nM=(M-1)/2;
	for (j=0;j<nM;j++)
	{
		indh=2*j*M;
		indhp1=indh+M;
		indhp2=indh+2*M;
		ind2h=j*nM;
		for (i=0;i<nM;i++)
		{
			ti=2*i;tip1=ti+1;tip2=ti+2;
			f2h[ind2h+i]=cpplus(cpmultr(rh[indh+ti],0.0625),cpplus(cpmultr(rh[indh+tip1],0.125),cpplus(cpmultr(rh[indh+tip2],0.0625),cpplus(cpmultr(rh[indhp1+ti],0.125),\
				cpplus(cpmultr(rh[indhp1+tip1],0.25),cpplus(cpmultr(rh[indhp1+tip2],0.125),cpplus(cpmultr(rh[indhp2+ti],0.0625),cpplus(cpmultr(rh[indhp2+tip1],0.125),cpmultr(rh[indhp2+tip2],0.0625)))))))));
		}
	}
}



inline void Interp_Operator(complex*u2h,complex*uh,uint_32 M)
{
	uint_32 i,j,nM,Mm1=M-1,oddjmnm,evenjmnm,jp1m,jm,nMm2,nMm1;
	nM=2*M+1;
	nMm2=nM-2;
	nMm1=nM-1;
	//perform the interpolation and correction
	for(j=0;j<Mm1;j++)
	{
		oddjmnm=(2*j+1)*nM;
		evenjmnm=(2*j+2)*nM;
		jp1m=(j+1)*M;
		jm=j*M;
		for(i=0;i<Mm1;i++)
		{
			//compute the u_2i+1_2j+1 term in finer grid
			uh[2*i+1+oddjmnm]=cpplus(uh[2*i+1+oddjmnm],u2h[i+jm]);
			//compute the u_2i+2_2j+1 term in finer grid
			uh[2*i+2+oddjmnm]=cpplus(uh[2*i+2+oddjmnm],cpmultr(cpplus(u2h[i+jm],u2h[i+1+jm]),0.5));

			//compute the u_2i+1_2j+2 term in finer grid
			uh[2*i+1+evenjmnm]=cpplus(uh[2*i+1+evenjmnm],cpmultr(cpplus(u2h[i+jm],u2h[i+jp1m]),0.5));

			//compute the u_2i+2_2j+2 term in finer grid
			uh[2*i+2+evenjmnm]=cpplus(uh[2*i+2+evenjmnm],cpplus(cpmultr(cpplus(u2h[i+jm],u2h[i+1+jm]),0.25),cpplus(cpmultr(u2h[i+jp1m],0.25),cpmultr(u2h[i+1+jp1m],0.25))));
		}

		//compute the u_2i+1_2j+1 term in finer grid at i=M1-1
		uh[nMm2+oddjmnm]=cpplus(uh[nMm2+oddjmnm],u2h[Mm1+jm]);

		//compute the u_2i+1_2j+2 term in finer grid at i=M1-1
		uh[nMm2+evenjmnm]=cpplus(uh[nMm2+evenjmnm],cpmultr(cpplus(u2h[Mm1+jm],u2h[Mm1+jp1m]),0.5));

		//compute the u_i_2j+1 term in finer grid at i=0,nM1-1
		uh[oddjmnm]=cpplus(uh[oddjmnm],cpmultr(u2h[jm],0.5));
		uh[nMm1+oddjmnm]=cpplus(uh[nMm1+oddjmnm],cpmultr(u2h[Mm1+jm],0.5));

		//compute the u_i_2j+2 term in finer grid at i=0
		uh[evenjmnm]=cpplus(uh[evenjmnm],cpmultr(cpplus(u2h[jm],u2h[jp1m]),0.25));

		//compute the u_i_2j+2 term in finer grid at i=nM1-1
		uh[nMm1+evenjmnm]=cpplus(uh[nMm1+evenjmnm],cpmultr(cpplus(u2h[Mm1+jm],u2h[Mm1+jp1m]),0.25));
	}

	//perform the remaining interpolation and correction
	oddjmnm=(2*M-1)*nM;
	evenjmnm=(nM-1)*nM;
	jm=(M-1)*M;

	for(i=0;i<Mm1;i++)
	{
		//compute the u_2i+1_2j+1 term in finer grid at j=M2-1
		uh[2*i+1+oddjmnm]=cpplus(uh[2*i+1+oddjmnm],u2h[i+jm]);

		//compute the u_2i+2_2j+1 term in finer grid at j=M2-1
		uh[2*i+2+oddjmnm]=cpplus(uh[2*i+2+oddjmnm],cpmultr(cpplus(u2h[i+jm],u2h[i+1+jm]),0.5));

		//compute the u_2i+1_j term in finer grid at j=0,nM2-1
		uh[2*i+1]=cpplus(uh[2*i+1],cpmultr(u2h[i],0.5));
		uh[2*i+1+evenjmnm]=cpplus(uh[2*i+1+evenjmnm],cpmultr(u2h[i+jm],0.5));

		//compute the u_2i+2_0 term in finer grid 
		uh[2*i+2]=cpplus(uh[2*i+2],cpmultr(cpplus(u2h[i],u2h[i+1]),0.25));

		//compute the u_2i+2_nM2-1 term in finer grid 
		uh[2*i+2+evenjmnm]=cpplus(uh[2*i+2+evenjmnm],cpmultr(cpplus(u2h[i+jm],u2h[i+1+jm]),0.25));
	}
	//compute the u_2i+1_2j+1 term in finer grid at i=M1-1,j=M2-1
	uh[nMm2+oddjmnm]=cpplus(uh[nMm2+oddjmnm],u2h[i+jm]);

	//compute the u_2i+1_j term in finer grid at i=M1-1, j=0,nM2-1
	uh[nMm2]=cpplus(uh[nMm2],cpmultr(u2h[Mm1],0.5));
	uh[nMm2+evenjmnm]=cpplus(uh[nMm2+evenjmnm],cpmultr(u2h[Mm1+jm],0.5));

	//compute u_i_2j+1 at i=0,nM1-1, j=M2-1
	uh[oddjmnm]=cpplus(uh[oddjmnm],cpmultr(u2h[jm],0.5));
	uh[nMm1+oddjmnm]=cpplus(uh[nMm1+oddjmnm],cpmultr(u2h[Mm1+jm],0.5));

	//perform the remaining interpolation and correction at boundary points
	uh[0]=cpplus(uh[0],cpmultr(u2h[0],0.25));
	uh[evenjmnm]=cpplus(uh[evenjmnm],cpmultr(u2h[jm],0.25));
	uh[nMm1]=cpplus(uh[nMm1],cpmultr(u2h[Mm1],0.25));
	uh[nMm1+evenjmnm]=cpplus(uh[nMm1+evenjmnm],cpmultr(u2h[Mm1+jm],0.25));
}



inline void tracking(double* offlaplacex,complex* mainlaplace ,complex* x,complex* b,uint_32 size)
{

	uint_32 i;
	complex* bmab=new complex[size-1];
	complex* beta=new complex[size-1];
	beta[0]=rdividcp(-offlaplacex[0],mainlaplace[0]);
	for(i=0;i<size-2;i++)
	{
		bmab[i]=cpplus(mainlaplace[i+1],cpmultr(beta[i],offlaplacex[i]));
		beta[i+1]=rdividcp(-offlaplacex[i+1],bmab[i]);
	}
	bmab[size-2]=cpplus(mainlaplace[size-1],cpmultr(beta[size-2],offlaplacex[size-2]));
	x[0]=cpdivid(b[0],mainlaplace[0]);
	for(i=1;i<size;i++)
		x[i]=cpdivid(cpplus(b[i],cpmultr(x[i-1],offlaplacex[i-1])),bmab[i-1]);
	delete[]bmab;
	int j;
	for(j=size-2;j>-1;j--)
	{
		x[j]=cpminus(x[j],cpmult(beta[j],x[j+1]));
	}
	delete[]beta;
}

inline void AXLGS_smoother(double* offlaplacex,double* offlaplacey,complex* mainlaplace,complex* initialgs,complex* rhs,uint_32 M,uint_32 time,bool iszeroinitial)
{
	if (time<1)
	{
		std::cout<<"FAI warning: the iteration time must be larger that zero!!"<<std::endl;
		return;
	}
	uint_32 Mm1,i,j,Mm2,vecpreinds,vecnextinds,veccurinds,offxcurinds,k,ntime=time;
	Mm1=M-1;
	Mm2=M-2;
	complex* temprhs=new complex[M];
	if (iszeroinitial)
	{
		ntime-=1;
		for (j=1;j<Mm1;j+=2)
		{
			vecpreinds=(j-1)*M;
			veccurinds=j*M;
			vecnextinds=(j+1)*M;
			offxcurinds=j*Mm1;

			//update the red line
			tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,rhs+veccurinds,M);
		}

        

		//update black lines
		vecnextinds=M;
		for (i=0;i<M;i++)
		{
			temprhs[i]=cpplus(rhs[i],cpmultr(initialgs[vecnextinds+i],offlaplacey[i]));
		}
		//update the first black line
		tracking(offlaplacex,mainlaplace,initialgs,temprhs,M);


		//update the remaining black lines
		for (j=2;j<Mm1;j+=2)
		{
			vecpreinds=(j-1)*M;
			veccurinds=j*M;
			vecnextinds=(j+1)*M;
			offxcurinds=j*Mm1;

			for (i=0;i<M;i++)
			{
				temprhs[i]=cpplus(rhs[veccurinds+i],cpplus(cpmultr(initialgs[vecpreinds+i],\
					offlaplacey[vecpreinds+i]),cpmultr(initialgs[vecnextinds+i],\
					offlaplacey[veccurinds+i])));
			}
			//update the black line
			tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,temprhs,M);
		}

		//the last black line
		vecpreinds=Mm2*M;
		veccurinds=Mm1*M;
		offxcurinds=Mm1*Mm1;
		for (i=0;i<M;i++)
		{
			temprhs[i]=cpplus(rhs[veccurinds+i],cpmultr(initialgs[vecpreinds+i],offlaplacey[vecpreinds+i]));
		}
		tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,temprhs,M);
	}


	for (k=0;k<ntime;k++)
	{
		for (j=1;j<Mm1;j+=2)
		{
			vecpreinds=(j-1)*M;
			veccurinds=j*M;
			vecnextinds=(j+1)*M;
			offxcurinds=j*Mm1;
			//compute the temporary right hand side
			
			for (i=0;i<M;i++)
			{	
				temprhs[i]=cpplus(rhs[veccurinds+i],cpplus(cpmultr(initialgs[vecpreinds+i],\
					offlaplacey[vecpreinds+i]),cpmultr(initialgs[vecnextinds+i],\
					offlaplacey[veccurinds+i])));
			}
				
			//update the red line
			tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,temprhs,M);
		}
		
		//update black lines
		vecnextinds=M;
		for (i=0;i<M;i++)
		{
			temprhs[i]=cpplus(rhs[i],cpmultr(initialgs[vecnextinds+i],offlaplacey[i]));
		}

		//update the first black line
		tracking(offlaplacex,mainlaplace,initialgs,temprhs,M);

		//update the remaining black lines
		for (j=2;j<Mm1;j+=2)
		{
			vecpreinds=(j-1)*M;
			veccurinds=j*M;
			vecnextinds=(j+1)*M;
			offxcurinds=j*Mm1;

			for (i=0;i<M;i++)
			{
				temprhs[i]=cpplus(rhs[veccurinds+i],cpplus(cpmultr(initialgs[vecpreinds+i],\
					offlaplacey[vecpreinds+i]),cpmultr(initialgs[vecnextinds+i],\
					offlaplacey[veccurinds+i])));
			}

			//update the black line
			tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,temprhs,M);
		}

		//the last black line
		vecpreinds=Mm2*M;
		veccurinds=Mm1*M;
		offxcurinds=Mm1*Mm1;
		for (i=0;i<M;i++)
		{
			temprhs[i]=cpplus(rhs[veccurinds+i],cpmultr(initialgs[vecpreinds+i],offlaplacey[vecpreinds+i]));
		}
		tracking(offlaplacex+offxcurinds,mainlaplace+veccurinds,initialgs+veccurinds,temprhs,M);
	}

	delete[]temprhs;
}
#endif
