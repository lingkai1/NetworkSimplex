#ifndef NMFS_H_
#define NMFS_H_

#include "types.hpp"
#include "assert.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <stack>
#include <queue>

#define USE_STATS_TIME
#define USE_STATS_COUNT

#define GGT_RELABEL // Use GGT algorithm (enable relabeling)
#ifdef GGT_RELABEL
//#define LAZY_RELABEL // use lazy relabeling heuristic
#if defined(LAZY_RELABEL)
//#define GLOBAL_RELABEL // use global update heuristic
#endif
#define GAP_RELABEL // use gap relabeling heuristic
#endif


//#define USE_PIVOTS_V_QUEUE
#define USE_CUR_T_QUEUE

#if defined(USE_PIVOTS_V_QUEUE)
//#define FORCE_STRICT_PIVOTS // force the pivots list to strictly contain pivots
#endif


#ifdef USE_STATS_TIME
double timer();
#endif

class NetworkMaxFlowSimplex {
public:
	const static int FORMAT_DIMACS = 1;

	int verbose;
	NodeID n; // number of nodes
	ArcID m; // number of arcs
	NodeID nMin; // smallest node id
	NodeID nMax; // highest node id
	NodeID nSentinel; // end of the node list marker == nMax+1
	std::vector<Node> nodes; // array of nodes
	std::vector<Arc> arcs; // array of arcs
	NodeID source; // source node pointer
	NodeID sink; // sink node pointer
	Flow flow; // flow value
	std::vector<NodeID> dc; // for each level d maintains the number of nodes with label == d

	// pivot vertex queue (vertices (of S) with at least 1 outgoing pivot)
#if defined(USE_PIVOTS_V_QUEUE)
#define bpq(name) qp##name
#define BPQ_PIVOTS_V
#include "bpq.hpp"
#undef BPQ_PIVOTS_V
#undef bpq
#endif
	// current t vertices queue
#if defined(USE_CUR_T_QUEUE)
#define bpq(name) qt##name
#define BPQ_CUR_T
#include "bpq.hpp"
#undef BPQ_CUR_T
#undef bpq
#endif
	// to relabel queue
#define bpq(name) qr##name
#define BPQ_RELABEL
#include "bpq.hpp"
#undef BPQ_RELABEL
#undef bpq

	NetworkMaxFlowSimplex(std::istream& is, int format = FORMAT_DIMACS, int verbose = 0);
	~NetworkMaxFlowSimplex();

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
		for (NodeID k = 0; k < nodes[root].size; k++) {
			op(u);
			u = nodes[u].next;
		}
	}

	void solve();
	void buildInitialBasis();

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


	ArcID extractMinPivot();
	void getSubtreeLastNode(NodeID u, NodeID& uLast);
	void deleteSubtree(NodeID q, NodeID u, NodeID uLast);
	void addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u);
	void changeRoot(NodeID q, NodeID r, NodeID& rLast);

	void setLabel(NodeID v, Dist k);
	bool isCurrent(NodeID v);

#if defined(GGT_RELABEL)
#if defined(GLOBAL_RELABEL)
	NodeID globalRelabelWork;
	double globalRelabelFreq;
	NodeID globalRelabelThreshold;
#endif
	void doGlobalRelabel();

	bool makeCur(NodeID v);
	void relabel(NodeID v);
#if defined(GAP_RELABEL)
	void gap(Dist k);
#endif
#endif


	bool checkValidCurArc(NodeID v);
	bool checkCurrent(NodeID v);


#include "dfsbfs.hpp"

	// statistics
	// times
#ifdef USE_STATS_TIME
	double t_parse;
	double t_buildInitialBasis;
	double t_solve;
	double t_total;
#endif
	// counters
#ifdef USE_STATS_COUNT
	long long c_basisIter;
	long long c_basisGrowS;
	long long c_basisGrowT;
	long long c_basisNoChange;
	long long c_StoTMoves;
	long long c_TtoSMoves;
	long long c_augPathTotalLen;
	long long c_pivotsInserted;
	long long c_pivotsDeleted;
	long long c_makeCur;
	long long c_relabel;
	long long c_relabelArcScans;
	long long c_globalUpdate;
	long long c_guArcScans;
	long long c_gap;
#endif
};

#include "dfsbfs.cpp"
#endif
