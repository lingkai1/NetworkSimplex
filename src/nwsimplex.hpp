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

//#define PLAIN_SIMPLEX
#define LAZY_SIMPLEX // todo: implement the pivot list with buckets



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
	//std::stack<NodeID> toRelabel; // todo remove this and change implementation

	NetworkMaxFlowSimplex(std::istream& is, int format = FORMAT_DIMACS, int verbose = 0);
	//NWS(int n, int m);
	~NetworkMaxFlowSimplex();

	void buildInitialBasis();
	void solve();

	// todo: change implementation (use all buckets)
#if defined(PLAIN_SIMPLEX)
	std::list<ArcID> pivots;
	void pivotsInsert(ArcID i);
	ArcID pivotsExtractMin();
	bool pivotsDelete(ArcID i);
#elif defined(LAZY_SIMPLEX)
	bool pivotsContains(NodeID v);
	void pivotsInsert(NodeID v, Dist d);
	ArcID pivotsExtractMin();
	NodeID pivotsDelete(NodeID v);
#endif

#if defined(LAZY_SIMPLEX)
	bool toRelabelContains(NodeID v);
	void toRelabelInsert(NodeID v, Dist d);
	NodeID toRelabelDelete(NodeID v);
	NodeID toRelabelExtractMin();
	bool toRelabelEmpty();
#endif

#if defined(LAZY_SIMPLEX)
	bool makeCur(NodeID v);
	void relabel(NodeID v);
#endif

	template <typename PrintNode> void printSubTree(NodeID root, const PrintNode& print);
	void printSubTree(NodeID root);

	void printS();
	void printT();

#ifdef USE_STATS_TIME
	void printTimeStats();
#endif
#ifdef USE_STATS_COUNT
	void printCountStats();
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

	// op can return false to break loop
	template <typename Op> inline void forAllOutPivots(NodeID u, const Op& op) {
		Node& nu = nodes[u];
		assert(nu.tree == IN_S);
		forAllOutArcs(u, uv, is) {
			Arc& auv = arcs[uv];
			NodeID v = auv.head; Node& nv = nodes[v];
			if (nv.tree == IN_T) {
				if (auv.resCap > 0) {
					if (!op(uv)) break;
				}
			}
		}
	}
	// op can return false to break loop
	template <typename Op> inline void forAllInPivots(NodeID u, const Op& op) {
		Node& nu = nodes[u];
		assert(nu.tree == IN_T);
		forAllOutArcs(u, uv, is) {
			Arc& auv = arcs[uv];
			NodeID v = auv.head; Node& nv = nodes[v];
			ArcID vu = auv.rev; Arc& avu = arcs[vu];
			if (nv.tree == IN_S) {
				if (avu.resCap > 0) {
					if (!op(vu)) break;
				}
			}
		}
	}

	void getSubtreeLastNode(NodeID u, NodeID& uLast);
	void deleteSubtree(NodeID q, NodeID u, NodeID uLast);
	void addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u);
	void changeRoot(NodeID q, NodeID r, NodeID& rLast);

	bool hasOutPivots(NodeID v);
	bool hasInPivots(NodeID v);

	// statistics
	// times
#ifdef USE_STATS_TIME
	double t_parse;
	double t_buildInitialBasis;
	double t_solve;
#endif
	// counters
#ifdef USE_STATS_COUNT
	long long c_basisIter;
	long long c_basisGrowS;
	long long c_basisGrowT;
	long long c_basisNoChange;
	long long c_pivotsInserted;
	long long c_pivotsDeleted;
	long long c_makeCur;
	long long c_relabel;
#endif

};

#include "dfsbfs.cpp"
#endif
