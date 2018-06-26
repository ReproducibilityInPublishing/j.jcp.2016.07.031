#include "ArgParseStandalone.h"
#include <iostream>
#include <cstdio>
typedef unsigned int uint_32;

int main(int argc, char** argv) {
	uint_32 level=2;
	uint_32 N=8192;
	double X_L=0;
	double X_R=M_PI;
	double Y_Low=0;
	double Y_Upp=M_PI;
	double T=0.5;
	double alpha=0.01;
	double eps=0.5e-8;
	double tol=1.0e-8;

	ArgParse::ArgParser Parser("Driver program for executing examples for 10.1016/j.jcp.2016.07.031");
	Parser.AddArgument("-l/--level", "Specify spatial grid resolution: M+1 here: M+1 = 2^level", &level, ArgParse::Argument::Optional);
	Parser.AddArgument("-N", "Specify time resolution here. N time steps.", &N, ArgParse::Argument::Optional);
	Parser.AddArgument("-a/--alpha", "Specify the value for alpha", &alpha);
	Parser.AddArgument("-XL", "Specify the left x bound", &X_L);
	Parser.AddArgument("-XR", "Specify the right x bound", &X_R);
	Parser.AddArgument("-YL", "Specify the lower y bound", &Y_Low);
	Parser.AddArgument("-YU", "Specify the upper y bound", &Y_Upp);
	Parser.AddArgument("-T", "Specify the upper time bound", &T);
	Parser.AddArgument("-e/--epsilon", "Specify the value for epsilon", &eps);
	Parser.AddArgument("-t/--tolerance", "Specify the tolerance value", &tol);

	if(Parser.ParseArgs(argc,argv) < 0) {
		fprintf(stderr, "There was a problem parsing the arguments!\n");
		return -1;
	}

	if (Parser.HelpPrinted()) {
		return 0;
	}

	std::cout << "Driver parameter summary" << std::endl;
	std::cout << "level: " << level << std::endl;
	std::cout << "N: " << N << std::endl;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "Spatial Domain: (X_L,Y_Low,X_R,Y_Upp) = (" << X_L << "," << Y_Low << "," << X_R << "," << Y_Upp << ")" << std::endl;
	std::cout << "Time Domain: [0,T] = [0," << T << "]" << std::endl;
	std::cout << "Approximation parameter epsilon: " << eps << std::endl;
	std::cout << "Iteration Tolerance: " << tol << std::endl;
}
