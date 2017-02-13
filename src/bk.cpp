#include <sstream>

#include "nmfs.hpp"
#include "parser.hpp"

using namespace std;



void NetworkMaxFlowSimplex::constructorBk(FILE* f) {
	string line;
	int lineNum;
	char c;
	NodeID u;

	if (verbose >= 1) cout << "Parsing dimacs format" << endl;

	lineNum = 0;
	Parser& pr = *new Parser(1024, f);

	while (true) {
		lineNum++;
		if (!pr.nextString()) return constructorBkFail(lineNum, 2);
		c = pr.s[0];
		if (c == 'p') break;
		if (!pr.gotoNextLine()) return constructorBkFail(lineNum, 2);
	}

	//if (!pr.nextString()) return constructorBkFail(lineNum, 2);
	//if (strcmp(s, "max")) return constructorBkFail(lineNum, 3);
	if (verbose >= 1) cout << "Max flow format" << endl;
	if (!pr.nextString()) return constructorBkFail(lineNum, 2);
	nMax = strtoll(pr.s, nullptr, 10);
	source = nMax++; // add source at the end
	sink = nMax++; // add sink at the end
	if (verbose >= 1) cout << "Node max id: " << nMax << endl;
	if (!pr.nextString()) return constructorBkFail(lineNum, 2);
	m = strtoll(pr.s, nullptr, 10);
	if (verbose >= 1) cout << "Arcs in file: " << m << endl;

	nMin = 0;
	n = nMax - nMin + 1;
	nSentinel = nMax+1;

	mMin = 0;
	mMax = m - 1;
	mMax += n;

	arcs.reserve(2*(mMax+1+1));
	nodes.resize(nMax+1+1);

	forAllNodes(u) {
		Node& node = nodes[u];
		node.first = UNDEF_ARC;
	}


	vector<NodeID>& tails = *new vector<NodeID>();
	tails.reserve(2*(mMax+1+1));
	Arc a;
	ArcID mRead = 0;
	while (mRead < m) {
		lineNum++;
		if (!pr.gotoNextLine()) break;
		if (!pr.nextString()) break;
		c = *pr.s;
		switch (c) {
		case 'c': {

		} break;
		case 'n': {
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			NodeID id = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			int capFromSource = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			int capToSink = strtoll(pr.s, nullptr, 10);

			//cout << "id:" << id << " csrc:" << capFromSource << " csink:" << capToSink << endl;

			if (capFromSource > 0) {
				a.head = id;
				a.resCap = capFromSource;
				a.rev = UNDEF_ARC;
				tails.push_back(source);
				arcs.push_back(a);
			}
			if (capToSink > 0) {
				a.head = sink;
				a.resCap = capToSink;
				a.rev = UNDEF_ARC;
				tails.push_back(id);
				arcs.push_back(a);
			}
		} break;
		case 'a': {
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			u = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			a.head = strtoll(pr.s, nullptr, 10);
			if (!pr.nextString()) return constructorBkFail(lineNum, 2);
			a.resCap = strtoll(pr.s, nullptr, 10);

			a.rev = UNDEF_ARC;
			tails.push_back(u);
			arcs.push_back(a);

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

void NetworkMaxFlowSimplex::constructorBkFail(int lineNum, int code) {
	cerr << "Error: ";
	switch (code) {
	case 1: cerr << "Cannot read line " << lineNum << endl; break;
	case 2: cerr << "Cannot parse line " << lineNum << endl; break;
	case 3: cerr << "Not a max flow problem" << endl; break;
	default: cerr << "Unknown failure" << endl; break;
	}
}
