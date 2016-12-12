#include "nwsimplex.hpp"

using namespace std;

// set the object in a state which makes destructor defined under any circumstance
// this is not needed if we don't use pointers
void NetworkMaxFlowSimplex::initialize() {
	//nodes = nullptr;
	//arcs = nullptr;
}

// constructors
NetworkMaxFlowSimplex::NetworkMaxFlowSimplex(istream& is, int format, int verbose) : listPivots(*this), listRelabel(*this)
{
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
	cout << "Parsing time: " << t_parse << " s" << endl;
	cout << "Build initial basis time: " << t_buildInitialBasis << " s" << endl;
	cout << "Solving time: " << t_solve << " s" << endl;
	cout << "Total solving time: " << t_total << " s" << endl;
}
#endif

#ifdef USE_STATS_COUNT
void NetworkMaxFlowSimplex::printCountStats() {
	cout << "Iterations: " << c_basisIter << endl;
	cout << "S growths: " << c_basisGrowS << endl;
	cout << "T growths: " << c_basisGrowT << endl;
	cout << "Basis unchanged: " << c_basisNoChange << endl;
	cout << "Size of trees moved from S to T: " << c_StoTMoves << endl;
	cout << "Size of trees moved from T to S: " << c_TtoSMoves << endl;
	cout << "Total length of augmenting paths: " << c_augPathTotalLen << endl;
	cout << "Pivots inserted: " << c_pivotsInserted << endl;
	cout << "Pivots deleted: " << c_pivotsDeleted << endl;
	cout << "Make current calls: " << c_makeCur << endl;
	cout << "Relabel calls: " << c_relabel << endl;
	cout << "Relabel arc scans: " << c_relabelArcScans << endl;
	cout << "Global updates: " << c_globalUpdate << endl;
	cout << "Global updates arc scans: " << c_guArcScans << endl;
	cout << "Gaps: " << c_gap << endl;
}
#endif


// print a subtree with a custom printNode function
template <typename PrintNode>
void NetworkMaxFlowSimplex::printSubTree(NodeID root, const PrintNode& print) {
	forAllSubTree(root, [&](NodeID u){
		Node& nu = nodes[u];
		if (u == root) {
			print(root);
		}
		else {
			cout << ",";
			print(u);
		}
	});
}
void NetworkMaxFlowSimplex::printSubTree(NodeID root) {
	return printSubTree(root, [&](NodeID u){ cout << u; });
}

void NetworkMaxFlowSimplex::printS() {
	cout << "S tree rooted at " << source << endl;
	printSubTree(source, [&](NodeID u){ cout << u << " (d=" << nodes[u].d << ")"; });
	cout << endl;
	cout << "S arcs:" << endl;
	forAllSubTree(source, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (u != source && r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << p << "->" << u << " rc: " << ai.resCap << endl;
		}
	});
}
void NetworkMaxFlowSimplex::printT() {
	cout << "T tree rooted at " << sink << endl;
	printSubTree(sink, [&](NodeID u){ cout << u << " (d=" << nodes[u].d << ")"; });
	cout << endl;
	cout << "T arcs:" << endl;
	forAllSubTree(sink, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (u != sink && r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << u << "->" << p << " rc: " << ai.resCap << endl;
		}
	});
}




