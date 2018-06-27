#include <iostream>
#include "FAI_2DFSDEsolver.h"
#include"gamma.h"
#include <time.h>
#include <iomanip>
#include <cmath>
double inf_norm_real(double* vec, uint_32 veclen);
double averageiter(double* vec, uint_32 veclen);

void example3(double& average_iter, double& cpu_time, double& infnorm_error, const uint_32 level, const uint_32 N, const double X_L, const double X_R, const double Y_Low, const double Y_Upp, const double T, const double alpha, const double eps, const double tol) {
	// Some variable initialization
	uint_32 i, j;
	clock_t t1, t2;
	double s;

	// Start clock
	t1=clock();

	// Set the inner spacial grid resolution.
	uint_32 M=_Pow_int(2, level)-1;
	uint_32 Ms=M*M;

	// Define grid spacings
	double tau=T/(double)N;
	double h1, h2;
	h1=(X_R-X_L)/(double)(M+1);
	h2=(Y_Upp-Y_Low)/(double)(M+1);

	//Not sure what this is for
	double miu = pow(tau, alpha)*Gamma(2-alpha);
	double oneovermiu = 1.0/miu;
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

	double fovergmfma=Gamma(4)/Gamma(4-alpha);
	double t;
	double* expx=new double[M];
	double* expy=new double[M];
	double* exactsol=new double[Ms];

	double mh1;
	double mh2;
	double tmalpha=3-alpha;
	double tq;
	double tptma;
	double oneoverh1s=1/(h1*h1);
	double oneoverh2s=1/(h2*h2);
	double x_s=X_L+0.5*h1;
	double x_e=X_L+(M+0.5)*h1;
	double y_s=Y_Low+0.5*h2;
	double y_e=Y_Low+(M+0.5)*h2;
	double expxl=exp(X_L);
	double expxr=exp(X_R);
	double expyl=exp(Y_Low);
	double expyu=exp(Y_Upp);
	double oneoh1sxl=oneoverh1s*expxl;
	double oneoh1sxr=oneoverh1s*expxr;
	double oneoh2syl=oneoverh2s*expyl;
	double oneoh2syu=oneoverh2s*expyu;
	double fovemtptma;
	uint_32 Mm1=M-1;
	uint_32 Mmmm1=Mm1*M;

	for (i=0; i<M; i++) {
		expx[i]=exp(X_L+(i+1)*h1);
		expy[i]=exp(Y_Low+(i+1)*h2);
	}

	uint_32 indy;
	uint_32 n;

	for (j=0; j<M; j++) {
		indy=j*M;

		for (i=0; i<M; i++) {
			exactsol[indy+i]=expx[i]*expy[j];
		}
	}

	delete[]expx;
	delete[]expy;

	double* expxy=new double[Ms];

	for (j=0; j<M; j++) {
		indy=j*M;
		mh2=Y_Low+(j+1)*h2;

		for (i=0; i<M; i++) {
			expxy[indy+i]=exp(mh2*(X_L+(i+1)*h1));
		}
	}

	for (n=0; n<N; n++) {
		t=(n+1)*tau;
		tptma=pow(t, tmalpha);
		tq=t*t*t;
		fovemtptma=fovergmfma*tptma;
		rhs[n]=new complex[Ms];

		for (j=0; j<M; j++) {
			mh2=Y_Low+(j+1)*h2;
			indy=j*M;

			for (i=0; i<M; i++) {
				mh1=X_L+(i+1)*h1;
				rhs[n][indy+i].r=exactsol[indy+i]*(fovemtptma-tq*expxy[indy+i]*(2+mh1+mh2));
				rhs[n][indy+i].i=0;
			}

			rhs[n][indy].r+=(oneoh1sxl*postv_func(x_s, mh2)*tq*exp(mh2));
			rhs[n][indy+Mm1].r+=(oneoh1sxr*postv_func(x_e, mh2)*tq*exp(mh2));
		}

		for (i=0; i<M; i++) {
			mh1=X_L+(i+1)*h1;
			rhs[n][i].r+=(oneoh2syl*postv_func(mh1, y_s)*tq*exp(mh1));
			rhs[n][Mmmm1+i].r+=(oneoh2syu*postv_func(mh1, y_e)*tq*exp(mh1));
		}
	}

	delete[]expxy;
	double* iter_num_arr=new double[N];
	Time_frac_diffusion_2Deq_solver(N, level, X_L, X_R, Y_Low, Y_Upp,
	                                Frac_div_appro_coeff, rhs, eps, iter_num_arr, tol);
	t2=clock();
	delete[]Frac_div_appro_coeff;
	double norm1, error, *tempsolver=new double[Ms], max=0, max1;
	norm1=inf_norm_real(exactsol, Ms);
	norm1*=pow(T, 3);

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

	average_iter = averageiter(iter_num_arr, N);
	//std::cout<<"averageiter:"<<averageiter(iter_num_arr, N)<<std::endl;
	error=max/norm1;
	infnorm_error = error;
	//std::cout<<"the relative error under infinite norm is:"<<error<<std::endl;
	s=(double)(t2-t1)/CLOCKS_PER_SEC;
	cpu_time = s;
	//std::cout.setf(std::ios::fixed);
	//std::cout<<std::setprecision(7)<<"the running time is :"<<s<<std::endl;
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
inline double postv_func(double x, double y) {
	return exp(x*y);
}
