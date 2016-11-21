#include "misc.hpp"
#include <sstream>
#include "nwsimplex.hpp"

using namespace std;

void NWS::constructorDimacs(istream& is) {
	stringstream ss;
	string line;
	int lineNum;
	string s;
	char c;
	NodeID u;
	ArcID a;
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
	a = 0;
	bool sourceSet = false, sinkSet = false;
	ArcID nArcsProcessed = 0;
	while (!sourceSet || !sinkSet || nArcsProcessed < m) {
		if (!getline(is, line)) return constructorDimacsFail(lineNum, 1);
		ss.str(""); ss.clear(); ss << line; lineNum++;
		if (ss >> c)
			switch (c) {
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
				arc.rev = a + 1;
				tails.push_back(u);
				arcs.push_back(arc);
				a++;

				// add reverse arc
				swap(u, arc.head);
				arc.resCap = 0;
				arc.rev = a - 1;
				tails.push_back(u);
				arcs.push_back(arc);
				a++;

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
	for (a = 0; a < m; a++) q[p[a]] = a;
	// fix reverse arcs
	for (a = 0; a < m; a++) {
		arcs[a].rev = q[arcs[a].rev];
	}
	delete &q;
	delete &p;

	// eliminate parallel arcs and fix reverse arcs again
	ArcID na;
	for (a = 0, na = 0; a < m; ) {
		ArcID revMin = arcs[a].rev;
		Cap totCap = 0;
		do {
			if (arcs[a].rev < revMin)
				revMin = arcs[a].rev;
			totCap += arcs[a].resCap;
			a++;
		} while (a < m && tails[a] == tails[na] && arcs[a].head == arcs[na].head);
		if (tails[na] < arcs[na].head)
			arcs[revMin].rev = na;
		else
			arcs[arcs[na].rev].rev = na;
		arcs[na].resCap = totCap;
		na++;
		if (a < m) {
			arcs[na] = arcs[a];
			tails[na] = tails[a];
		}
	}
	m = na;

	// determine first arcs
	for (a = 0; a < m; a++) {
		u = tails[a];
		if (nodes[u].first == UNDEF_ARC)
			nodes[u].first = a;
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


	bToRelabelFirst.resize(n);
	bPivotArcFirst.resize(n);


	if (verbose >= 1) cout << "n: " << n << " m: " << m << endl;

	if (verbose >= 2) {
		forAllArcs(u, a, aStop) {
			Arc& arc = arcs[a];
			cout << "Arc: " << a << " (" << u << ", " << arc.head << ") c: " << arc.resCap << " rev_id: " << arc.rev << endl;
		}
	}

}

void NWS::constructorDimacsFail(int lineNum, int code) {
	cerr << "Error: ";
	switch (code) {
	case 1: cerr << "Cannot read line " << lineNum << endl; break;
	case 2: cerr << "Cannot parse line " << lineNum << endl; break;
	case 3: cerr << "Not a max flow problem" << endl; break;
	default: cerr << "Unknown failure" << endl; break;
	}
}
