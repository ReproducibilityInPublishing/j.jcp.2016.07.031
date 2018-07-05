#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include "argument_parser.h"
#include "defs.h"
#include "example2_BFSMGM.h"
#include "reporter.h"

int main(int argc, char** argv) {
	uint_32 level=2;
	uint_32 N=15000;
	double X_L=-1;
	double X_R=1;
	double Y_Low=-1;
	double Y_Upp=1;
	double T=0.5;
	double alpha=0.01;
	double eps=0.5e-8;
	double tol=1.0e-8;

	int s;
	if ((s = argument_parser(level, N, X_L, X_R, Y_Low, Y_Upp, T, alpha, eps, tol, argc, argv)) != 0) {
		return s;
	}

	double average_iter = 0.;
	double cpu_time = 0.;
	double infnorm_error = 0.;
	example2_BFSMGM(average_iter, cpu_time, infnorm_error, level, N, X_L, X_R, Y_Low, Y_Upp, T, alpha, eps, tol);

	std::cerr << "Driver parameter summary" << std::endl;
	std::cerr << "level: " << level << std::endl;
	std::cerr << "N: " << N << std::endl;
	std::cerr << "alpha: " << alpha << std::endl;
	std::cerr << "Spatial Domain: (X_L,Y_Low,X_R,Y_Upp) = (" << X_L << "," << Y_Low << "," << X_R << "," << Y_Upp << ")" << std::endl;
	std::cerr << "Time Domain: [0,T] = [0," << T << "]" << std::endl;
	std::cerr << "Approximation parameter epsilon: " << eps << std::endl;
	std::cerr << "Iteration Tolerance: " << tol << std::endl;

	std::cerr << std::endl;
	std::cerr << "Results:" << std::endl;
	std::cerr<<"averageiter: "<<average_iter<<std::endl;
	std::cerr.setf(std::ios::scientific);
	std::cerr<<"the relative error under infinite norm is: "<<infnorm_error<<std::endl;
	std::cerr.unsetf(std::ios::scientific);
	std::cerr.setf(std::ios::fixed);
	std::cerr<<std::setprecision(7)<<"the running time is : "<<cpu_time<<std::endl;

	reporter(level, N, X_L, X_R, Y_Low, Y_Upp, T, alpha, eps, tol, average_iter, cpu_time, infnorm_error);
}
