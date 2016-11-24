#ifndef NWS_H_
#define NWS_H_

#define __cplusplus 201103L

#include "types.hpp"
//#define NDEBUG
#include "assert.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <stack>
#include <queue>
#include <deque>
#include <set>

#define IN_NONE 0
#define IN_S 1
#define IN_T 2

#define forAllNodes(u) for (NodeID u = 0; u <= nMax; u++)
#define forAllOutArcs(u, i, iStop) if (nodes[u].first != UNDEF_ARC)\
		for (NodeID i = nodes[u].first, iStop = nodes[u+1].first; i != iStop; i++)
#define forAllArcs(u, i, iStop) forAllNodes(u)\
		forAllOutArcs(u, i, iStop)


class NWS {
public:
	const static int FORMAT_DIMACS = 1;

	int verbose;
	NodeID n; // number of nodes
	ArcID m; // number of arcs
	NodeID nMin; // smallest node id
	NodeID nMax; // highest node id
	std::vector<Node> nodes; // array of nodes
	std::vector<Arc> arcs; // array of arcs
	NodeID source; // source node pointer
	NodeID sink; // sink node pointer
	Flow flow; // flow value
	NodeID nSentinel; // end of the node list marker

	// buckets
	// todo implement
	std::vector<BucketToRelabel> bToRelabel; // for each level d first node to relabel
	std::vector<BucketPivot> bPivots; // for each level d first pivot arc
	std::stack<NodeID> toRelabel; // todo remove this and change implementation

	NWS(std::istream& is, int format = FORMAT_DIMACS);
	//NWS(int n, int m);
	~NWS();

	void solve();
	void buildInitialBasisFromBFS();
	void buildInitialBasisFromDFS();
	void prepareNextPivot();

	bool makeCur(NodeID v);
	void relabel(NodeID v);



	// todo: change implementation (use all buckets) now only bucket[d=0] is used as a FIFO queue
	void pivotInsert(ArcID i);
	bool pivotDelete(ArcID i);
	bool pivotExtractMin(ArcID& i);

	template <typename Op> inline void forAllSubTree(NodeID root, const Op& op) {
		NodeID u = root;
		for (int k = 0; k < nodes[root].stSize; k++) {
			op(u);
			u = nodes[u].next;
		}
	}

	void printSubTree(NodeID root);
	void printS();
	void printT();

	ArcID findArc(NodeID u, NodeID v, bool binarySearch = true);

#include "dfsbfs.hpp"

private:
	void initialize();
	void constructorDimacs(std::istream& is);
	void constructorDimacsFail(int lineNum, int code);

};

#include "dfsbfs.cpp"
#endif
