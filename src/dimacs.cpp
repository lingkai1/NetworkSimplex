#include <sstream>

#include "nmfs.hpp"
#include "util.hpp"
#include "parser.hpp"

using namespace std;

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

	/*// eliminate parallel arcs and fix reverse arcs again
	ArcID i, ni;
	for (i = 0, ni = 0; i < m; ) {
		ArcID revMin = arcs[i].rev;
		Cap totCap = 0;
		do {
			if (arcs[i].rev < revMin)
				revMin = arcs[i].rev;
			totCap += arcs[i].resCap;
			i++;
		} while (i < m && tails[i] == tails[ni] && arcs[i].head == arcs[ni].head);
		if (tails[ni] < arcs[ni].head)
			arcs[revMin].rev = ni;
		else
			arcs[arcs[ni].rev].rev = ni;
		arcs[ni].resCap = totCap;
		ni++;
		if (i < m) {
			arcs[ni] = arcs[i];
			tails[ni] = tails[i];
		}
	}
	mMax = ni - 1;
	m = mMax - mMin + 1;
	arcs.resize(mMax+1+1);*/
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

void NetworkMaxFlowSimplex::constructorDimacs(istream& is) {
	string line;
	int lineNum;
	string s;
	char c;
	NodeID u;
	ArcID i;
	// not sure this is the most efficient way of parsing text files in C++. May be improved.

	if (verbose >= 1) cout << "Parsing dimacs format" << endl;

	lineNum = 0;
	Parser& pr = *new Parser(1024, stdin);

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

	delete &pr;

	m = (ArcID)arcs.size();
	mMin = 0;
	mMax = m - 1;

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
