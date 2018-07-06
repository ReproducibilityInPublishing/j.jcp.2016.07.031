#include <iostream>
#include"gamma.h"
#include <time.h>
#include <iomanip>
#include <cmath>
#include "example1.h"
double inf_norm_real(double* vec, uint_32 veclen);

// Solve the tridiagonal matrix with (1+2*r) on the diagonal, -r on the off diagonals.
// Solves the x system with constant j.
void solveA(double* u, double* b, uint_32 j,  double r, uint_32 M) {
	// Build c' factors
	double c[M-1];
	c[0] = -r/(1+2*r);
	for(uint_32 i=1; i<M-1; ++i) {
		c[i] = -r/(1+2*r+r*c[i-1]);
	}
	// Build d' factors
	uint_32 indy = j*(M+2);
	double d[M];
	d[0] = b[indy+1]/(1+2*r);
	for(uint_32 i=1; i<M; ++i) {
		d[i] = (b[indy+i+1]+r*d[i-1])/(1+2*r+r*c[i-1]);
	}
	u[indy+M] = d[M-1];
	for(uint_32 i=M-1; i>0; --i) {
		u[indy+i] = d[i-1]-c[i-1]*u[indy+i+1];
	}
}

// Solve the tridiagonal matrix with (1+2*r) on the diagonal, -r on the off diagonals.
// Solves the x system with constant j.
void solveB(double* u, double* b, uint_32 i,  double r, uint_32 M) {
	// Build c' factors
	double c[M-1];
	c[0] = -r/(1+2*r);
	for(uint_32 i=1; i<M-1; ++i) {
		c[i] = -r/(1+2*r+r*c[i-1]);

	}
	// Build d' factors
	double d[M];
	d[0] = b[(M+2)+i]/(1+2*r);
	for(uint_32 j=1; j<M; ++j) {
		uint_32 indy = (j+1)*(M+2);
		d[j] = (b[indy+i]+r*d[j-1])/(1+2*r+r*c[j-1]);
	}
	u[M*(M+2)+i] = d[M-1];
	for(uint_32 j=M-1; j>0; --j) {
		uint_32 indy = j*(M+2);
		u[indy+i] = d[j-1]-c[j-1]*u[(j+1)*(M+2)+i];
	}
}

uint_32 _Pow_int(uint_32 x, uint_32 n) {
	uint_32 answer = 1;
	for(uint_32 i=0; i<n; ++i) {
		answer *= x;
	}
	return answer;
};

void example1_BDADI(double& cpu_time, double& infnorm_error, const uint_32 level, const uint_32 N, const double X_L, const double X_R, const double Y_Low, const double Y_Upp, const double T, const double alpha, const double eps, const double tol) {
	clock_t t1, t2;
	double s;
	t1=clock();
	uint_32 M=_Pow_int(2, level)-1;
	uint_32 Mp2 = M+2;
	uint_32 Ms=Mp2*Mp2;
	double tau=T/(double)N, mu;
	mu=pow(tau, alpha)*Gamma(2-alpha);
	double h1, h2;
	h1=(X_R-X_L)/(double)(M+1);
	h2=(Y_Upp-Y_Low)/(double)(M+1);
	double onema=1-alpha, apkp1, apk;
	double* ak=new double[N];
	apk=0;
	apkp1=1;

	for(uint_32 i=0; i<N; i++) {
		ak[i]=apkp1-apk;
		apk=apkp1;
		apkp1=pow(i+2, onema);
	}

	double** u=new double*[N];
	double** f=new double*[N];
	double gammafac=1./Gamma(3-alpha);
	double t;
	double* exactsol=new double[Ms];
	double tmalpha=2-alpha, oneoverh1s=1/(h1*h1), oneoverh2s=1/(h2*h2),
	                 ts, tp2ma;
	double x_s=X_L+0.5*h1, x_e=X_L+(M+0.5)*h1, y_s=Y_Low+0.5*h2,
	       y_e=Y_Low+(M+0.5)*h2;
	uint_32 Mm1=M-1, Mmmm1=Mm1*M;

	// Initialize solver

	for (uint_32 j=0; j<Mp2; j++) {
		uint_32 indy=j*Mp2;
		double y = Y_Low+j*h2;

		for (uint_32 i=0; i<Mp2; i++) {
			double x = X_L+i*h1;
			exactsol[indy+i]=sin(x)*sin(y);
		}
	}

	for (uint_32 n=0; n<N; n++) {
		t=(n+1)*tau;
		ts=t*t;
		tp2ma=pow(t, tmalpha);
		f[n]=new double[Ms];
		u[n]=new double[Ms];

		for(uint_32 i=0; i<Ms; ++i) {
			f[n][i] = 0.;
			u[n][i] = 0.;
		}

		for (uint_32 j=0; j<Mp2; ++j) {
			uint_32 indy=j*Mp2;

			for (uint_32 i=0; i<Mp2; i++) {
				f[n][indy+i] = 2*exactsol[indy+i]*(tp2ma*gammafac+ts);
			}
		}
		// X boundary condition
		for(uint_32 j=0; j<Mp2; ++j) {
			u[n][j*Mp2] = exactsol[j*Mp2]*ts;
			u[n][j*Mp2+(Mp2-1)] = exactsol[j*Mp2+(Mp2-1)]*ts;
		}
		// Y boundary conditions
		for(uint_32 i=1; i<Mp2-1; ++i) {
			u[n][i] = exactsol[i]*ts;
			u[n][(Mp2-1)*Mp2+i] = exactsol[(Mp2-1)*Mp2+i]*ts;
		}
	}

	// Solver starting

	double* u_star = new double[Ms];
	double* b = new double[Ms];
	double* u_star_oj = new double[Mp2];
	double* u_star_mj = new double[Mp2];
	double mid_y_fac = 1+2*oneoverh2s;
	double off_y_fac = -oneoverh2s;
	double mid_x_fac = 1+2*oneoverh1s;
	double off_x_fac = -oneoverh1s;

	// Initial time step

	//u_star boundaries
	u_star_oj[0] = 0.;
	u_star_mj[0] = 0.;
	u_star_oj[Mp2-1] = 0.;
	u_star_mj[Mp2-1] = 0.;
	for(uint_32 j=1; j<Mp2-1; ++j) {
		u_star_oj[j] = mid_y_fac*u[0][j*Mp2]+off_y_fac*u[0][(j-1)*Mp2]+off_y_fac*u[0][(j+1)*Mp2];
		u_star_mj[j] = mid_y_fac*u[0][j*Mp2+(Mp2-1)]+off_y_fac*u[0][(j-1)*Mp2+(Mp2-1)]+off_y_fac*u[0][(j+1)*Mp2+(Mp2-1)];
	}

	for(uint_32 j=1; j<Mp2-1; ++j) {
		uint_32 indy = j*Mp2;
		for(uint_32 i=1;i<Mp2-1; ++i) {
			b[indy+i] = mu*f[0][indy+i];
		}
		//Boundary conditions for u_star
		b[indy+1] += mu*oneoverh1s*u_star_oj[j];
		b[indy+Mp2-2] += mu*oneoverh1s*u_star_mj[j];
	}

	// Solve the first M systems
	for(uint_32 j=1; j<Mp2-1; ++j) {
		solveA(u_star, b, j, mu*oneoverh1s, M);
	}
	// Handle intermediate boundary conditions
	for(uint_32 i=1; i<Mp2-1; ++i) {
		u_star[i] += mu*oneoverh2s*u[0][i];
		u_star[(Mp2-2)*Mp2+i] += mu*oneoverh2s*u[0][(Mp2-2)*Mp2+i];
	}
	// Solve the second M systems
	for(uint_32 i=1; i<Mp2-1; ++i) {
		solveB(u[0], u_star, i, mu*oneoverh2s, M);
	}

	// Solver nominal loop

	double factorizer_factor = mu*mu*oneoverh1s*oneoverh2s;
	for(uint_32 n=1; n<N; ++n) {
		double t = (1+n)*tau;
		double ts = t*t;

		//u_star boundaries
		for(uint_32 j=1; j<Mp2-1; ++j) {
			u_star_oj[j] = mid_y_fac*u[n][j*Mp2]+off_y_fac*u[n][(j-1)*Mp2]+off_y_fac*u[n][(j+1)*Mp2];
			u_star_mj[j] = mid_y_fac*u[n][j*Mp2+(Mp2-1)]+off_y_fac*u[n][(j-1)*Mp2+(Mp2-1)]+off_y_fac*u[n][(j+1)*Mp2+(Mp2-1)];
		}
	
		for(uint_32 j=1; j<Mp2-1; ++j) {
			uint_32 indy = j*Mp2;
			for(uint_32 i=1;i<Mp2-1; ++i) {
				// Source term
				b[indy+i] = mu*f[n][indy+i];
				// Derivative terms
				for(uint_32 k=0; k<n; ++k) {
					double term = (ak[n-k-1]-ak[n-k])*u[k][indy+i];
					b[indy+i] += term;
				}
				// Extra Divergence terms
				uint_32 jm1 = j-1;
				uint_32 jp1 = j+1;
				uint_32 im1 = i-1;
				uint_32 ip1 = i+1;
				double term = factorizer_factor*(
					(u[n-1][jm1*(Mp2)+im1]+u[n-1][jp1*(Mp2)+im1]+
					 u[n-1][jm1*(Mp2)+ip1]+u[n-1][jp1*(Mp2)+ip1])
					-2*(u[n-1][jm1*(Mp2)+i]+u[n-1][jp1*(Mp2)+i]+
					    u[n-1][j*(Mp2)+im1]+u[n-1][j*(Mp2)+ip1])
					+4*u[n-1][j*Mp2+i]);
				b[indy+i] += term;
			}
			//Boundary conditions for u_star
			b[indy+1] += mu*oneoverh1s*u_star_oj[j];
			b[indy+Mp2-2] += mu*oneoverh1s*u_star_mj[j];
		}

		// Solve the first M systems
		for(uint_32 j=1; j<Mp2-1; ++j) {
			solveA(u_star, b, j, mu*oneoverh1s, M);
		}

		// Handle intermediate boundary conditions
		for(uint_32 i=1; i<Mp2-1; ++i) {
			u_star[i] += mu*oneoverh2s*u[n][i];
			u_star[(Mp2-2)*Mp2+i] += mu*oneoverh2s*u[n][(Mp2-2)*Mp2+i];
		}
		// Solve the second M systems
		for(uint_32 i=1; i<Mp2-1; ++i) {
			solveB(u[n], u_star, i, mu*oneoverh2s, M);
		}
	}

	// Solver finished
	t2=clock();

	double norm1, error, *tempsolver=new double[Ms], max=0, max1;
	norm1=inf_norm_real(exactsol, Ms);
	norm1*=(T*T);

	for (uint_32 n=0; n<N; n++) {
		t=(n+1)*tau;
		ts=t*t;

		for (uint_32 j=0; j<Mp2; ++j) {
			uint_32 indy = j*Mp2;
			for(uint_32 i=0; i<Mp2; ++i) {
				tempsolver[indy+i]=u[n][indy+i]-ts*exactsol[indy+i];
			}
		}

		max1=inf_norm_real(tempsolver, Ms);

		if(max<max1) {
			max=max1;
		}
	}

	error=max/norm1;
	infnorm_error = error;
	s=(double)(t2-t1)/CLOCKS_PER_SEC;
	cpu_time = s;

	delete[]exactsol;
	delete[]tempsolver;
	delete[] ak;

	for (uint_32 i=0; i<N; i++) {
		delete [] u[i];
		delete [] f[i];
	}
	delete [] u;
	delete [] f;
	delete [] u_star;
	delete [] b;
	delete [] u_star_oj;
	delete [] u_star_mj;
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
double inf_norm_real_grid(double* vec, uint_32 veclen) {
	double max=0.;


	for (uint_32 i=1; i<veclen; i++) {
		if (max<fabs(vec[i])) {
			max=fabs(vec[i]);
		}
	}

	return max;
}
inline double postv_func(double x, double y) {
	return 1.;
}
