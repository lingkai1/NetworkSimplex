
#include <fstream>
#include <sstream>
#include <getopt.h>

#include "nmfs.hpp"

using namespace std;


struct Main {

	char* title;

	int main(int argCount, char** argVars) {

		// parse arguments
		title = argVars[0];
		struct option options[] = {
				{"help", no_argument, 0, 0},
				{"inputfile", required_argument, 0, 0},
				{"outputfile", required_argument, 0, 0},
				{0,0,0,0}
		};
		int optionIndex;
		while (getopt_long_only(argCount, argVars, "", options, &optionIndex) != -1) {
			switch (optionIndex) {
			case 0:
				printHelp();
				return mainExit(2);
			case 1:
				//cout << "i:" << optarg << endl;
				//optarg = argVars[optind]; optind++;
				//cout << "i:" << optarg << endl;
				if (freopen(optarg, "r", stdin) == NULL)
					mainExit(1);
				break;
			case 2:
				if (freopen(optarg, "w", stdout) == NULL)
					mainExit(1);
				break;
			default:
				return mainExit(1);
			}
		}

		// build network simplex solver
		cout << "Network Max Flow Simplex" << endl;
		NetworkMaxFlowSimplex nmfs(stdin, NetworkMaxFlowSimplex::FORMAT_DIMACS, 1);

		cout << endl;
		cout << "Solving..." << endl;
		// call solver

		nmfs.solve();

		// output results
		cout << endl;
		cout << "Max Flow: " << nmfs.flow << endl;
		cout << endl;
#ifdef USE_STATS_TIME
		nmfs.printTimeStats();
#endif
		cout << endl;
#ifdef USE_STATS_COUNT
		nmfs.printCountStats();
#endif

		return mainExit(0);
	}

	int mainExit(int code) {
		// error message
		switch (code) {
		case 1: cerr << "Error: Bad argument" << endl; break;
		default:
			break;
		}
		// cleanup
		switch (code) {
		case 0:
			break;
		default:
			break;
		}
		return code;
	}

	void printHelp() {
		cerr << "Usage: " << endl;
		cerr << title << " [--help]" << endl;
	}


};

int main(int argCount, char** argVars) {
	return Main().main(argCount, argVars);
}
