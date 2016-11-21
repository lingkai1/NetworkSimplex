#include "nwsimplex.hpp"

using namespace std;

// set the object in a state which makes destructor defined under any circumstance
// this is not needed if we don't use pointers
void NWS::initialize() {
	//nodes = nullptr;
	//arcs = nullptr;
}

// constructors
NWS::NWS(istream& is, int format) {
	initialize();
	switch (format) {
	case FORMAT_DIMACS: {
		constructorDimacs(is);
	}	break;
	default:
		break;
	}
}

// destructor
NWS::~NWS() {
	//if (nodes != nullptr) delete[] nodes;
	//if (arcs != nullptr) delete[] arcs;
}


// miscellaneous graph functions

// unused stuff.
ArcID NWS::findArc(NodeID u, NodeID v, bool binarySearch) {
	typename vector<Arc>::iterator it;
	if (!binarySearch) {
		it = find_if(arcs.begin()+nodes[u].first, arcs.begin()+nodes[u+1].first, [&v](Arc arc){ return v == arc.head; });
		if (it != arcs.begin()+nodes[u+1].first)
			return distance(arcs.begin(), it);
		else
			return UNDEF_ARC;
	}
	else {
		Arc refArc;
		refArc.head = v;
		auto pair = equal_range(arcs.begin()+nodes[u].first, arcs.begin()+nodes[u+1].first, refArc, [&](Arc arc1, Arc arc2){ return arc1.head < arc2.head; });
		it = pair.first;
		if (pair.second > pair.first)
			return distance(arcs.begin(), it);
		else
			return UNDEF_ARC;
	}
}




