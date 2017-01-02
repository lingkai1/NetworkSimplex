#include "nmfs.hpp"
#include "util.hpp"

using namespace std;


// constructors
NetworkMaxFlowSimplex::NetworkMaxFlowSimplex(FILE* f, int format, int verbose) : qr(*this)
, qt(*this)
{
	this->verbose = verbose;
#ifdef USE_STATS_TIME
	t_parse = timer();
#endif
	switch (format) {
	case FORMAT_DIMACS: {
		constructorDimacs(f);
	}	break;
	default:
		break;
	}
#ifdef USE_STATS_TIME
	t_parse = timer() - t_parse;
#endif
	// allocate
	qt.initialize();
	qr.initialize();
	d.resize(nSentinel+1); d.shrink_to_fit();
	fill(d.begin(), d.end(), 0);
	color.resize(nMax+1); color.shrink_to_fit();
	dfsiv.resize(nMax+1); dfsiv.shrink_to_fit();
}

// destructor
NetworkMaxFlowSimplex::~NetworkMaxFlowSimplex() {

}

void NetworkMaxFlowSimplex::addReverseArcs(vector<NodeID>& tails) {

	// add reverse arcs
	ArcID vu = mMax+1;
	for (ArcID uv = mMin; uv <= mMax; uv++) {
		Arc& auv = arcs[uv];
		NodeID v = auv.head; Node& nv = nodes[v];
		NodeID u = tails[uv]; Node& nu = nodes[u];
		Arc avu;
		auv.rev = vu;
		avu.head = u;
		avu.resCap = 0;
		avu.rev = uv;
		tails.push_back(v);
		arcs.push_back(avu);
		vu++;
	}

	m = (ArcID)arcs.size();
	mMin = 0;
	mMax = m - 1;
}


void NetworkMaxFlowSimplex::sortArcs(vector<NodeID>& tails) {

	// sort arcs
	vector<ArcID>& p = *new vector<ArcID>(arcs.size());
	sortPermutation(p, [&](ArcID i, ArcID j){ Arc& ai = arcs[i], & aj = arcs[j];
	return (tails[i] < tails[j] || (tails[i] == tails[j] && ai.head < aj.head));
	} );
	applyPermutation(arcs, p);
	applyPermutation(tails, p);

	// determine reverse arcs
	// compute p inverse
	vector<ArcID>& q = *new vector<ArcID>(p.size());
	for (ArcID i = mMin; i <= mMax; i++)
		q[p[i]] = i;
	delete &p;
	// fix reverse arcs
	for (ArcID i = mMin; i <= mMax; i++)
		if (arcs[i].rev != UNDEF_ARC)
			arcs[i].rev = q[arcs[i].rev];
	delete &q;
}

void NetworkMaxFlowSimplex::mergeAddParallelArcs(vector<NodeID>& tails) {
	// merge parallel arcs and fix reverse arcs again
	ArcID i;
	ArcID ni;
	for (i = 0, ni = 0; i <= mMax; ) {
		Cap totCap = 0;

		ArcID revMin = mMax+1;
		do {
			if (arcs[i].rev < revMin)
				revMin = arcs[i].rev;
			totCap += arcs[i].resCap;
			i++;
		} while (i <= mMax && tails[i] == tails[ni] && arcs[i].head == arcs[ni].head);

		if (tails[ni] < arcs[ni].head) {
			if (revMin <= mMax) {
				arcs[revMin].rev = ni;
				//arcs[ni].rev = revMin;
			}
		}
		else {
			arcs[arcs[ni].rev].rev = ni;
		}

		arcs[ni].resCap = totCap;
		ni++;
		if (i <= mMax) {
			arcs[ni] = arcs[i];
			tails[ni] = tails[i];
		}
	}
	mMax = ni - 1;
	m = mMax - mMin + 1;
	arcs.resize(mMax+1+1);
}

void NetworkMaxFlowSimplex::setFirstArcs(vector<NodeID>& tails) {

	// determine first arcs
	for (ArcID i = mMin; i <= mMax; i++) {
		NodeID u = tails[i];
		if (nodes[u].first == UNDEF_ARC)
			nodes[u].first = i;
	}
	// sentinel node
	Node& node = nodes[nMax+1];
	node.first = mMax+1;
	for (NodeID u = nMax; u >= 0; u--) {
		if (nodes[u].first == UNDEF_ARC)
			nodes[u].first = nodes[u+1].first;
	}
	node.d = nMax+1;

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




