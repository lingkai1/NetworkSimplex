#include <sstream>

#include "nmfs.hpp"
#include "parser.hpp"

using namespace std;



void NetworkMaxFlowSimplex::constructorDimacs(FILE* f) {
	string line;
	int lineNum;
	string s;
	char c;
	NodeID u;
	ArcID i;

	if (verbose >= 1) cout << "Parsing dimacs format" << endl;

	lineNum = 0;
	Parser& pr = *new Parser(1024, f);

	while (true) {
		lineNum++;
		if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
		c = pr.s[0];
		if (c == 'p') break;
		if (!pr.gotoNextLine()) return constructorDimacsFail(lineNum, 2);
	}

	if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
	//if (strcmp(s, "max")) return constructorDimacsFail(lineNum, 3);
	if (verbose >= 1) cout << "Max flow format" << endl;
	if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
	nMax = strtoll(pr.s, nullptr, 10);
	if (verbose >= 1) cout << "Node max id: " << nMax << endl;
	if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
	m = strtoll(pr.s, nullptr, 10);
	if (verbose >= 1) cout << "Arcs in file: " << m << endl;

	nMin = 0;
	n = nMax - nMin + 1;
	nSentinel = nMax+1;

	mMin = 0;
	mMax = m - 1;

	arcs.reserve(2*(mMax+1+1));
	nodes.resize(nMax+1+1);

	forAllNodes(u) {
		Node& node = nodes[u];
		node.first = UNDEF_ARC;
	}


	vector<NodeID>& tails = *new vector<NodeID>();
	tails.reserve(2*(mMax+1+1));
	Arc a;
	i = 0;
	bool sourceSet = false, sinkSet = false;
	ArcID mRead = 0;
	while (!sourceSet || !sinkSet || mRead < m) {
		lineNum++;
		if (!pr.gotoNextLine()) break;
		if (!pr.nextString()) break;
		c = *pr.s;
		switch (c) {
		case 'c': {
			// comment
		} break;
		case 'n': {
			NodeID id;
			if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
			id = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
			c = pr.s[0];
			if (c == 's') {
				source = id; sourceSet = true;
			}
			else if (c == 't') {
				sink = id; sinkSet = true;
			}
			else return constructorDimacsFail(lineNum, 2);
		} break;
		case 'a': {
			if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
			u = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
			a.head = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorDimacsFail(lineNum, 2);
			a.resCap = strtoll(pr.s, nullptr, 10);

			a.rev = UNDEF_ARC;
			tails.push_back(u);
			arcs.push_back(a);
			i++;

			mRead++;
		} break;
		default:
			break;
		}
	}
	m = (ArcID)arcs.size();
	mMin = 0;
	mMax = m - 1;

	delete &pr;

	addReverseArcs(tails);

	sortArcs(tails);

	mergeAddParallelArcs(tails);

	setFirstArcs(tails);

	delete &tails;

	nodes.shrink_to_fit();
	arcs.shrink_to_fit();

	if (verbose >= 1) cout << "n: " << n << " m: " << m << endl;

	if (verbose >= 4) {
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
