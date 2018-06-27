#include <iostream>
#include <iomanip>
#include "reporter.h"
#include <cmath>

void reporter(const uint_32 level, const uint_32 N, const double X_L, const double X_R, const double Y_Low, const double Y_Upp, const double T, const double alpha, const double eps, const double tol, const double average_iter, const double cpu_time, const double infnorm_error) {
	std::cout.unsetf(std::ios_base::floatfield);
	std::cout << X_L << "," << Y_Low << "," << X_R << "," << Y_Upp << "," << T << "," << eps << "," << tol << "," << pow(2,level) << "," << alpha << "," << N << "," << average_iter << "," << cpu_time << ",";
	std::cout << std::scientific << infnorm_error << std::endl;
	std::cout.unsetf(std::ios_base::floatfield);
}
