#include "argument_parser.h"
#include "ArgParseStandalone.h"

int argument_parser(uint_32& level, uint_32& N, double& X_L, double& X_R, double& Y_Low, double& Y_Upp, double& T, double& alpha, double& eps, double& tol, int argc, char** argv) {
	ArgParse::ArgParser Parser("An Example inside article 10.1016/j.jcp.2016.07.031");
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
		return 1;
	}
	return 0;
}
