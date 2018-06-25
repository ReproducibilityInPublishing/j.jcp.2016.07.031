#include <iostream>
#include "FAI_2DFSDEsolver.h"
#include"gamma.h"
#include <time.h>
#include <iomanip>
#define PI  3.1415926535897932384626433832

// Some problem parameters
const uint_32 level=9;
const uint_32 N=32;
const double X_L=-1, X_R=1, Y_Low=-1, Y_Upp=1, T=0.5, alpha=1./6.;
const double eps=0.5e-8;
const double tol=1.0e-8;

double inf_norm_real(double* vec, uint_32 veclen);
double averageiter(double* vec, uint_32 veclen);

int main() {
	// Some variable initialization
	uint_32 i, j;
	clock_t t1, t2;
	double s;

	// Start clock
	t1=clock();

	// Set the inner spacial grid resolution.
	uint_32 M=_Pow_int(2, level)-1;
	uint_32 Ms=M*M;

	// Set the inner spacial grid resolution.
	double tau=T/(double)N;
	double h1 = (X_R-X_L)/(double)(M+1);
	double h2 = (Y_Upp-Y_Low)/(double)(M+1);

	// Computing fractional derivative coefficients
	double miu = pow(tau, alpha)*Gamma(2-alpha);
	double oneovermiu = 1.0/miu;
	double onema=1-alpha;
	double apkp1 = 1;
	double apk = 0;
	double* ak=new double[N];
	for(i=0; i<N; i++) {
		ak[i]=apkp1-apk;
		apk=apkp1;
		apkp1=pow(i+2, onema);
	}
	double* Frac_div_appro_coeff=new double[N];
	Frac_div_appro_coeff[0]=oneovermiu;
	for (i=1; i<N; i++) {
		Frac_div_appro_coeff[i]=oneovermiu*(ak[i]-ak[i-1]);
	}
	delete[]ak;

	complex** rhs=new complex*[N];

	double fovergama4ma=Gamma(4)/Gamma(4-alpha);
	double t;
	double* xarry=new double[M];
	double* yarray=new double[M];
	double* exactsol=new double[Ms];

	double mh1;
	double mh2;
	double tmalpha=3-alpha;
	double tq;
	double tp3ma;
	double oneoverh1s=1/(h1*h1);
	double oneoverh2s=1/(h2*h2);
	double x_s=X_L+0.5*h1;
	double x_e=X_L+(M+0.5)*h1;
	double y_s=Y_Low+0.5*h2;
	double y_e=Y_Low+(M+0.5)*h2;
	double oneoh1sxl=oneoverh1s*X_L;
	double oneoh1sxr=oneoverh1s*X_R;
	double oneoh2syl=oneoverh2s*Y_Low;
	double oneoh2syu=oneoverh2s*Y_Upp;
	uint_32 Mm1=M-1;
	uint_32 Mmmm1=Mm1*M;

	for (i=0; i<M; i++) {
		xarry[i]=X_L+(i+1)*h1;
		yarray[i]=Y_Low+(i+1)*h2;
	}

	uint_32 indy;
	uint_32 n;

	for (j=0; j<M; j++) {
		indy=j*M;

		for (i=0; i<M; i++) {
			exactsol[indy+i]=xarry[i]*yarray[j];
		}
	}

	for (i=0; i<M; i++) {
		xarry[i]*=xarry[i];
		yarray[i]*=yarray[i];
	}

	double* expxy=new double[Ms];

	for (i=0; i<Ms; i++) {
		expxy[i]=exp(exactsol[i]);
	}

	for (n=0; n<N; n++) {
		// t != 0 for first step??
		t=(n+1)*tau;
		tq=t*t*t;
		tp3ma=pow(t, tmalpha);

		rhs[n]=new complex[Ms];

		for (j=0; j<M; j++) {
			mh2=Y_Low+(j+1)*h2;
			indy=j*M;

			for (i=0; i<M; i++) {
				mh1=X_L+(i+1)*h1;
				rhs[n][indy+i].r=tp3ma*exactsol[indy+i]*fovergama4ma-expxy[indy+i]*tq*
				                 exp(xarry[i]+yarray[j]);
				rhs[n][indy+i].i=0;
			}

			rhs[n][indy].r+=(oneoh1sxl*postv_func(x_s, mh2)*tq*mh2);
			rhs[n][indy+Mm1].r+=(oneoh1sxr*postv_func(x_e, mh2)*tq*mh2);
			//rhs[n][indy].r+=(tq*x_s*mh2);
			//rhs[n][indy+Mm1].r+=(tq*mh2*x_e);
		}

		for (i=0; i<M; i++) {
			mh1=X_L+(i+1)*h1;
			rhs[n][i].r+=(oneoh2syl*postv_func(mh1, y_s)*tq*mh1);
			//rhs[n][i].r+=(oneoh2syl*postv_func(mh1, y_s)*tq*mh1);
			rhs[n][Mmmm1+i].r+=(oneoh2syu*postv_func(mh1, y_e)*tq*mh1);
			//rhs[n][Mmmm1+i].r+=(oneoh2syu*postv_func(mh1, y_e)*tq*mh1);
		}
	}

	double iter_num_arr[1];
	Time_frac_diffusion_2Deq_solver(N, level, X_L, X_R, Y_Low, Y_Upp,
	                                Frac_div_appro_coeff, rhs, eps, iter_num_arr, tol);
	t2=clock();
	delete[]Frac_div_appro_coeff;
	double norm1, error, *tempsolver=new double[Ms], max=0, max1;
	norm1=inf_norm_real(exactsol, Ms);
	norm1*=(T*T*T);

	for (n=0; n<N; n++) {
		t=(n+1)*tau;
		tq=t*t*t;

		for (i=0; i<Ms; i++) {
			tempsolver[i]=rhs[n][i].r-tq*exactsol[i];
		}

		max1=inf_norm_real(tempsolver, Ms);

		if(max<max1) {
			max=max1;
		}
	}

	std::cout<<"averageiter:"<<iter_num_arr[0]<<std::endl;
	error=max/norm1;
	std::cout<<"the relative error under infinite norm is:"<<error<<std::endl;
	std::cout.setf(std::ios::fixed);
	s=(double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout<<std::setprecision(7)<<"the running time is :"<<s<<std::endl;
	/*
	 for (i=0;i<Ms;i++)
	 {
		 std::cout<<exactsol[i]*pow(T,onepalpha)-rhs[N-1][i].r<<std::endl;
	 }
	*/
	delete[]exactsol;
	delete[]tempsolver;

	for (i=0; i<N; i++) {
		delete[]rhs[i];
	}
}
double inf_norm_real(double* vec, uint_32 veclen) {
	double max=fabs(vec[0]);

	for (uint_32 i=1; i<veclen; i++) {
		if (max<fabs(vec[i])) {
			max=fabs(vec[i]);
		}
	}

	return max;
}
double averageiter(double* vec, uint_32 veclen) {
	double average=vec[0];

	for (uint_32 i=1; i<veclen; i++) {
		average+=vec[i];
	}

	return average/(double)veclen;
}
