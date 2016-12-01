#include "misc.hpp"
#include <sstream>
#include "nwsimplex.hpp"

using namespace std;

void NetworkMaxFlowSimplex::constructorDimacs(istream& is) {
	stringstream ss;
	string line;
	int lineNum;
	string s;
	char c;
	NodeID u;
	ArcID i;
	// not sure this is the most efficient way of parsing text files in C++. May be improved.

	if (verbose >= 1) cout << "Parsing dimacs format" << endl;

	lineNum = 0;
	do {
		if (!getline(is, line)) return constructorDimacsFail(lineNum, 1);
		ss.str(""); ss.clear(); ss << line; lineNum++;
		if (!(ss >> c)) return constructorDimacsFail(lineNum, 2);
	} while (c != 'p');
	if (!(ss >> s)) return constructorDimacsFail(lineNum, 2);
	if (s.compare("max")) return constructorDimacsFail(lineNum, 3);
	if (verbose >= 1) cout << "Max-flow problem" << endl;
	if (!(ss >> nMax)) return constructorDimacsFail(lineNum, 2);
	if (verbose >= 1) cout << "Node max id: " << nMax << endl;
	if (!(ss >> m)) return constructorDimacsFail(lineNum, 2);
	if (verbose >= 1) cout << "Arcs in file: " << m << endl;

	n = nMax+1;
	arcs.reserve(2*m);
	nodes.resize(n+1);

	forAllNodes(u) {
		Node& node = nodes[u];
		node.first = UNDEF_ARC;
	}


	vector<NodeID>& tails = *new vector<NodeID>();
	tails.reserve(2*m);
	Arc arc;
	i = 0;
	bool sourceSet = false, sinkSet = false;
	ArcID nArcsProcessed = 0;
	while (!sourceSet || !sinkSet || nArcsProcessed < m) {
		if (!getline(is, line)) return constructorDimacsFail(lineNum, 1);
		ss.str(""); ss.clear(); ss << line; lineNum++;
		if (ss >> c)
			switch (c) {
			case 'c': {
				// comment
			} break;
			case 'n': {
				NodeID id;
				if (!(ss >> id)) return constructorDimacsFail(lineNum, 2);
				if (!(ss >> c)) return constructorDimacsFail(lineNum, 2);
				if (c == 's') {
					source = id; sourceSet = true;
				}
				else if (c == 't') {
					sink = id; sinkSet = true;
				}
				else return constructorDimacsFail(lineNum, 2);
			} break;
			case 'a': {
				if (!(ss >> u)) return constructorDimacsFail(lineNum, 2);
				if (!(ss >> arc.head)) return constructorDimacsFail(lineNum, 2);
				if (!(ss >> arc.resCap)) return constructorDimacsFail(lineNum, 2);
				arc.rev = i + 1;
				tails.push_back(u);
				arcs.push_back(arc);
				i++;

				// add reverse arc
				swap(u, arc.head);
				arc.resCap = 0;
				arc.rev = i - 1;
				tails.push_back(u);
				arcs.push_back(arc);
				i++;

				nArcsProcessed++;
			} break;
			default:
				break;
			}
	}

	assert(arcs.size() == tails.size());
	m = (ArcID)arcs.size();

	// sort arcs
	vector<ArcID>& p = *new vector<ArcID>(m);
	sortPermutation(p, [&](ArcID i, ArcID j){
		Arc& arc1 = arcs[i], & arc2 = arcs[j];
		return ((tails[i] < tails[j]) || (tails[i] == tails[j] && arc1.head < arc2.head));
	} );
	applyPermutation(arcs, p);
	applyPermutation(tails, p);

	// determine reverse arcs
	// compute p inverse
	vector<ArcID>& q = *new vector<ArcID>(m);
	for (i = 0; i < m; i++) q[p[i]] = i;
	// fix reverse arcs
	for (i = 0; i < m; i++) {
		arcs[i].rev = q[arcs[i].rev];
	}
	delete &q;
	delete &p;

	// eliminate parallel arcs and fix reverse arcs again
	ArcID na;
	for (i = 0, na = 0; i < m; ) {
		ArcID revMin = arcs[i].rev;
		Cap totCap = 0;
		do {
			if (arcs[i].rev < revMin)
				revMin = arcs[i].rev;
			totCap += arcs[i].resCap;
			i++;
		} while (i < m && tails[i] == tails[na] && arcs[i].head == arcs[na].head);
		if (tails[na] < arcs[na].head)
			arcs[revMin].rev = na;
		else
			arcs[arcs[na].rev].rev = na;
		arcs[na].resCap = totCap;
		na++;
		if (i < m) {
			arcs[na] = arcs[i];
			tails[na] = tails[i];
		}
	}
	m = na;

	// determine first arcs
	for (i = 0; i < m; i++) {
		u = tails[i];
		if (nodes[u].first == UNDEF_ARC)
			nodes[u].first = i;
	}
	// sentinel node
	nSentinel = nMax+1;
	Node& node = nodes[nSentinel];
	node.first = m;
	for (u = nSentinel; u >= 0; u--) {
		if (nodes[u].first == UNDEF_ARC)
			nodes[u].first = nodes[u+1].first;
	}

	delete &tails;

	arcs.resize(m);
	arcs.shrink_to_fit();

	bToRelabel.resize(n);
	bPivots.resize(n);


	if (verbose >= 1) cout << "n: " << n << " m: " << m << endl;

	if (verbose >= 2) {
		forAllArcs(u, i, is) {
			Arc& arc = arcs[i];
			cout << "Arc: " << i << " (" << u << ", " << arc.head << ") c: " << arc.resCap << " rev: " << arc.rev << endl;
		}
	}

}

void NetworkMaxFlowSimplex::constructorDimacsFail(int lineNum, int code) {
	cerr << "Error: ";
	switch (code) {
	case 1: cerr << "Cannot read line " << lineNum << endl; break;
	case 2: cerr << "Cannot parse line " << lineNum << endl; break;
	case 3: cerr << "Not a max flow problem" << endl; break;
	default: cerr << "Unknown failure" << endl; break;
	}
}
