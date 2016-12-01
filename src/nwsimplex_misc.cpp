#include "nwsimplex.hpp"

using namespace std;

// set the object in a state which makes destructor defined under any circumstance
// this is not needed if we don't use pointers
void NetworkMaxFlowSimplex::initialize() {
	//nodes = nullptr;
	//arcs = nullptr;
}

// constructors
NetworkMaxFlowSimplex::NetworkMaxFlowSimplex(istream& is, int format, int verbose) {
	initialize();
	this->verbose = verbose;
#ifdef USE_STATS_TIME
	t_parse = timer();
#endif
	switch (format) {
	case FORMAT_DIMACS: {
		constructorDimacs(is);
	}	break;
	default:
		break;
	}
#ifdef USE_STATS_TIME
	t_parse = timer() - t_parse;
#endif
}

// destructor
NetworkMaxFlowSimplex::~NetworkMaxFlowSimplex() {
	//if (nodes != nullptr) delete[] nodes;
	//if (arcs != nullptr) delete[] arcs;
}

#ifdef USE_STATS_TIME
void NetworkMaxFlowSimplex::printTimeStats() {
	cout << "Parsing time: " << t_parse << endl;
	cout << "Build initial basis time: " << t_buildInitialBasis << endl;
	cout << "Solving time: " << t_solve << endl;
}
#endif







