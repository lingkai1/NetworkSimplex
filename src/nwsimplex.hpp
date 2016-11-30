#ifndef NWS_H_
#define NWS_H_

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
#include <list>

#define USE_STATS_TIME
#define USE_STATS_COUNT

#ifdef USE_STATS_TIME
double timer();
#endif

#define IN_NONE 0
#define IN_S 1
#define IN_T 2

#define PIVOTS_QUEUE
//#define PIVOTS_BUCKETS // todo: implement the pivot list with buckets



#define forAllNodes(u) for (NodeID u = 0; u <= nMax; u++)
#define forAllOutArcs(u, i, iStop) if (nodes[u].first != UNDEF_ARC)\
		for (NodeID i = nodes[u].first, iStop = nodes[u+1].first; i != iStop; i++)
#define forAllArcs(u, i, iStop) forAllNodes(u)\
		forAllOutArcs(u, i, iStop)


class NetworkMaxFlowSimplex {
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


	// statistics
	// times
#ifdef USE_STATS_TIME
	double t_parse;
	double t_buildInitialBasis;
	double t_solve;
#endif
	// counters


	NetworkMaxFlowSimplex(std::istream& is, int format = FORMAT_DIMACS, int verbose = 0);
	//NWS(int n, int m);
	~NetworkMaxFlowSimplex();

	void solve();
	void buildInitialBasis();
	void prepareNextPivot();

	bool makeCur(NodeID v);
	void relabel(NodeID v);



	// todo: change implementation (use all buckets)
#ifdef PIVOTS_QUEUE
	std::list<ArcID> pivots;
#endif
	void pivotInsert(ArcID i);
	bool pivotDelete(ArcID i);
	bool pivotExtractMin(ArcID& i);



	void printSubTree(NodeID root);
	void printS();
	void printT();

	ArcID findArc(NodeID u, NodeID v, bool binarySearch = true);

#ifdef USE_STATS_TIME
	void printTimeStats();
#endif


#include "dfsbfs.hpp"
	void testBFS();
	void testDFS();

private:
	void initialize();
	void constructorDimacs(std::istream& is);
	void constructorDimacsFail(int lineNum, int code);

	template <typename Op> inline void forAllSubTree(NodeID root, const Op& op) {
		NodeID u = root;
		for (NodeID k = 0; k < nodes[root].stSize; k++) {
			op(u);
			u = nodes[u].next;
		}
	}

	void getSubtreeLastNode(NodeID u, NodeID& uLast);
	void deleteSubtree(NodeID q, NodeID u, NodeID uLast);
	void addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u);
	void changeRoot(NodeID q, NodeID r, NodeID& rLast);
};

#include "dfsbfs.cpp"
#endif
