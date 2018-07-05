#include <iostream>
#include "FAI_2DFSDEsolver.h"
#include"gamma.h"
#include <time.h>
#include <iomanip>
#include <cmath>
#include "example2.h"
double inf_norm_real(double* vec, uint_32 veclen);

void example2(double& average_iter, double& cpu_time, double& infnorm_error, const uint_32 level, const uint_32 N, const double X_L, const double X_R, const double Y_Low, const double Y_Upp, const double T, const double alpha, const double eps, const double tol) {
	uint_32 M, i, j, Ms;
	clock_t t1, t2;
	double s;
	t1=clock();
	M=_Pow_int(2, level)-1;
	Ms=M*M;
	double tau=T/(double)N, miu, oneovermiu;
	miu=pow(tau, alpha)*Gamma(2-alpha);
	double h1, h2;
	h1=(X_R-X_L)/(double)(M+1);
	h2=(Y_Upp-Y_Low)/(double)(M+1);
	oneovermiu=1.0/miu;
	double onema=1-alpha, apkp1, apk;
	double* ak=new double[N];
	apk=0;
	apkp1=1;

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
	double t, *xarry=new double[M], *yarray=new double[Ms];
	double* exactsol=new double[Ms];
	double mh1, mh2, tmalpha=3-alpha, oneoverh1s=1/(h1*h1), oneoverh2s=1/(h2*h2),
	                 tq, tp3ma;
	double x_s=X_L+0.5*h1, x_e=X_L+(M+0.5)*h1, y_s=Y_Low+0.5*h2,
	       y_e=Y_Low+(M+0.5)*h2;
	double oneoh1sxl=oneoverh1s*X_L, oneoh1sxr=oneoverh1s*X_R,
	       oneoh2syl=oneoverh2s*Y_Low, oneoh2syu=oneoverh2s*Y_Upp;
	uint_32 Mm1=M-1, Mmmm1=Mm1*M;

	for (i=0; i<M; i++) {
		xarry[i]=X_L+(i+1)*h1;
		yarray[i]=Y_Low+(i+1)*h2;
	}

	uint_32 indy, n;

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
				                 (xarry[i]+yarray[j]);
				rhs[n][indy+i].i=0;
			}

			rhs[n][indy].r+=(oneoh1sxl*postv_func(x_s, mh2)*tq*mh2);
			rhs[n][indy+Mm1].r+=(oneoh1sxr*postv_func(x_e, mh2)*tq*mh2);
		}

		for (i=0; i<M; i++) {
			mh1=X_L+(i+1)*h1;
			rhs[n][i].r+=(oneoh2syl*postv_func(mh1, y_s)*tq*mh1);
			rhs[n][Mmmm1+i].r+=(oneoh2syu*postv_func(mh1, y_e)*tq*mh1);
		}
	}

	double iter_num_arr = 0.;
	Time_frac_diffusion_2Deq_AIMGM_solver(N, level, X_L, X_R, Y_Low, Y_Upp,
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

	average_iter = iter_num_arr;
	//std::cerr<<"averageiter: "<<iter_num_arr<<std::endl;
	error=max/norm1;
	infnorm_error = error;
	//std::cerr.setf(std::ios::scientific);
	//std::cerr<<"the relative error under infinite norm is: "<<error<<std::endl;
	//std::cerr.unsetf(std::ios::scientific);
	//std::cerr.setf(std::ios::fixed);
	s=(double)(t2-t1)/CLOCKS_PER_SEC;
	cpu_time = s;
	//std::cerr<<std::setprecision(7)<<"the running time is : "<<s<<std::endl;
	/*
	 for (i=0;i<Ms;i++)
	 {
		 std::cerr<<exactsol[i]*pow(T,onepalpha)-rhs[N-1][i].r<<std::endl;
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
inline double postv_func(double x, double y) {
	return exp(x*y);
}
