#ifndef NWS_H_
#define NWS_H_

#include "types.hpp"
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

#define IN_NONE -1
#define IN_S 0
#define IN_T 1

#define GGT_RELABEL // Use GGT algorithm (relabeling)
#define LAZY_RELABEL // use lazy relabeling heuristic
//#define GAP_RELABEL // use gap relabeling heuristic
//#define FORCE_STRICT_PIVOTS // force the pivots list to strictly contain pivots

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

	std::vector<BucketToRelabel> bToRelabel; // for each level d first node to relabel
	Dist bToRelabelMinD, bToRelabelMaxD;
	std::vector<BucketPivot> bPivots; // for each level d first pivot arc
	Dist bPivotsMinD, bPivotsMaxD;

	NetworkMaxFlowSimplex(std::istream& is, int format = FORMAT_DIMACS, int verbose = 0);
	~NetworkMaxFlowSimplex();

	void solve();

	void buildInitialBasis();

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

	template <typename Op> inline void forAllSubTree(NodeID root, const Op& op) {
		NodeID u = root;
		for (NodeID k = 0; k < nodes[root].stSize; k++) {
			op(u);
			u = nodes[u].next;
		}
	}

	bool hasOutPivots(NodeID u);
	bool hasInPivots(NodeID u);

	bool isTreeArc(Node& nu, Node& nv, ArcID uv, ArcID vu) { return nv.parent == vu || nu.parent == uv; }
	bool isResidual(Arc& auv) { return auv.resCap > 0; }
	bool isResidualOrTreeArc(Node& nu, Node& nv, ArcID uv, ArcID vu, Arc& auv) { return isResidual(auv) || isTreeArc(nu, nv, uv, vu); }
	bool isPivot(Node& nu, Node& nv, Arc& auv) { return isResidual(auv) && nu.tree == IN_S && nv.tree == IN_T; }
	bool isInPivot(Node& nv, Arc& avu) { return isResidual(avu) && nv.tree == IN_S; }
	bool isOutPivot(Node& nv, Arc& auv) { return isResidual(auv) && nv.tree == IN_T; }

private:
	void initialize();
	void constructorDimacs(std::istream& is);
	void constructorDimacsFail(int lineNum, int code);

	bool pivotsContains(NodeID v);
	void pivotsInsert(NodeID v, Dist d);
	ArcID pivotsExtractMin();
	NodeID pivotsDelete(NodeID v);

#if defined(GGT_RELABEL)
	NodeID globalRelabelWork;
	double globalRelabelFreq;
	NodeID globalRelabelThreshold;

	bool toRelabelContains(NodeID v);
	void toRelabelInsert(NodeID v, Dist d);
	NodeID toRelabelDelete(NodeID v);
	NodeID toRelabelExtractMin();
	bool toRelabelEmpty();
	bool toRelabelProcessed();

	bool makeCur(NodeID v);
	void relabel(NodeID v);
	void toRelabelFlush();
	void doGlobalRelabel();
#endif

	void getSubtreeLastNode(NodeID u, NodeID& uLast);
	void deleteSubtree(NodeID q, NodeID u, NodeID uLast);
	void addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u);
	void changeRoot(NodeID q, NodeID r, NodeID& rLast);

#include "dfsbfs.hpp"
	void testBFS();
	void testDFS();

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
