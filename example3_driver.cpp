#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include "argument_parser.h"
#include "defs.h"
#include "example3.h"

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

	double average_iter;
	double cpu_time;
	double infnorm_error;
	example3(average_iter, cpu_time, infnorm_error, level, N, X_L, X_R, Y_Low, Y_Upp, T, alpha, eps, tol);

	std::cout << "Driver parameter summary" << std::endl;
	std::cout << "level: " << level << std::endl;
	std::cout << "N: " << N << std::endl;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "Spatial Domain: (X_L,Y_Low,X_R,Y_Upp) = (" << X_L << "," << Y_Low << "," << X_R << "," << Y_Upp << ")" << std::endl;
	std::cout << "Time Domain: [0,T] = [0," << T << "]" << std::endl;
	std::cout << "Approximation parameter epsilon: " << eps << std::endl;
	std::cout << "Iteration Tolerance: " << tol << std::endl;

	std::cout << std::endl;
	std::cout << "Results:" << std::endl;
	std::cout<<"averageiter: "<<average_iter<<std::endl;
	std::cout.setf(std::ios::scientific);
	std::cout<<"the relative error under infinite norm is: "<<infnorm_error<<std::endl;
	std::cout.unsetf(std::ios::scientific);
	std::cout.setf(std::ios::fixed);
	std::cout<<std::setprecision(7)<<"the running time is : "<<cpu_time<<std::endl;
}
