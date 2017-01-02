#ifndef NMFS_H_
#define NMFS_H_

//#define NDEBUG
#include "assert.h"

#include "types.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <stack>
#include <queue>

#define USE_STATS_TIME
#define USE_STATS_COUNT

#define GGT_RELABEL // Use GGT algorithm (enable relabeling)
#ifdef GGT_RELABEL
#define GAP_RELABEL // use gap relabeling heuristic
#define LAZY_RELABEL // use lazy relabeling heuristic
#define GLOBAL_RELABEL // use global update heuristic
#define GU_FREQ 0.2
#define RELABEL_OPT // relabel optimization
// choose only 1 of the 4
//#define LEAVING_ARC_NEAR_S
//#define LEAVING_ARC_NEAR_T
#define LEAVING_ARC_NEAR_PIVOT // may be faster
//#define LEAVING_ARC_FARFROM_PIVOT
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
	ArcID mMin;
	ArcID mMax;
	NodeID nMin; // smallest node id
	NodeID nMax; // highest node id
	NodeID nSentinel; // end of the node list marker == nMax+1
	std::vector<Node> nodes; // array of nodes
	std::vector<Arc> arcs; // array of arcs
	NodeID source; // source node pointer
	NodeID sink; // sink node pointer
	Flow flow; // flow value
	std::vector<NodeID> d; // for each level k maintains the number of nodes with label == k
	std::vector<int> color;
	std::queue<NodeID> bfsq;
	std::stack<NodeID> dfss;
	std::vector<ArcID> dfsiv;

	// pivot vertex queue (vertices (of S) with at least 1 outgoing pivot)
	// current t vertices queue
#define bpq(name) qt##name
#define BPQ_QT
#include "bpq.hpp"
#undef BPQ_QT
#undef bpq
	// to relabel queue
#define bpq(name) qr##name
#define BPQ_QR
#include "bpq.hpp"
#undef BPQ_QR
#undef bpq

	NetworkMaxFlowSimplex(FILE* f, int format = FORMAT_DIMACS, int verbose = 0);
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

	bool isTreeArc(Node& nu, Node& nv, ArcID uv, ArcID vu) { return nv.parent == vu || nu.parent == uv; }
	bool isResidual(Arc& auv) { return auv.resCap > 0; }
	bool isPseudoResidual(Node& nu, Node& nv, ArcID uv, ArcID vu, Arc& auv) { return isResidual(auv) || isTreeArc(nu, nv, uv, vu); }
	bool isPivot(Node& nu, Node& nv, Arc& auv) { return isResidual(auv) && nu.tree == IN_S && nv.tree == IN_T; }
	bool isInPivot(Node& nv, Arc& avu) { return isResidual(avu) && nv.tree == IN_S; }
	bool isOutPivot(Node& nv, Arc& auv) { return isResidual(auv) && nv.tree == IN_T; }

private:
	void addReverseArcs(std::vector<NodeID>& tails);
	void sortArcs(std::vector<NodeID>& tails);
	void mergeAddParallelArcs(std::vector<NodeID>& tails);
	void setFirstArcs(std::vector<NodeID>& tails);
	void constructorDimacs(FILE* f);
	void constructorDimacsFail(int lineNum, int code);


	ArcID extractMinPivot();
	void getSubtreeLastNode(NodeID u, NodeID& uLast);
	void deleteSubtree(NodeID q, NodeID u, NodeID uLast);
	void addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u);
	void changeRoot(NodeID q, NodeID r, NodeID& rLast);

	void setLabel(NodeID v, Dist k);
	bool isCurrent(NodeID v);

#if defined(GLOBAL_RELABEL)
	NodeID globalRelabelWork;
	double globalRelabelFreq;
	NodeID globalRelabelThreshold;
#endif
	void globalUpdate();

	void relabelPhase();
	bool makeCur(NodeID v);
	void relabel(NodeID v);
#if defined(GAP_RELABEL)
	void gap(Dist k);

#endif


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
	long long c_makeCur;
	long long c_relabel;
	long long c_relabelArcScans;
	long long c_globalUpdate;
	long long c_guArcScans;
	long long c_gap;
#endif


	// some debug or unused functions
	// returns true if u has outgoing pivots
	// u is assumed to be in S (it doesn't have to be)
	bool hasOutPivots(NodeID u) {
		Node& nu = nodes[u];
		forAllOutArcs(u, uv, is) {
			Arc& auv = arcs[uv];
			NodeID v = auv.head; Node& nv = nodes[v];
			if (isOutPivot(nv, auv))
				return true;
		}
		return false;
	}
	// returns true if u has incoming pivots
	// u is assumed to be in T (it doesn't have to be)
	bool hasInPivots(NodeID u) {
		Node& nu = nodes[u];
		forAllOutArcs(u, uv, is) {
			Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
			ArcID vu = auv.rev; Arc& avu = arcs[vu];
			if (isInPivot(nv, avu))
				return true;
		}
		return false;
	}

	bool checkValidCurArc(NodeID v) {
		bool r = true;
		assert(v != source);
		if (nodes[v].d == n) return true;
		forAllOutArcs(v, vu, is) {
			Node& nv = nodes[v];
			Arc& avu = arcs[vu]; NodeID u = avu.head; Node& nu = nodes[u];
			ArcID uv = avu.rev; Arc& auv = arcs[uv];
			Dist luv = isPseudoResidual(nu, nv, uv, vu, auv) ? 1 : INF_DIST;
			//		if (u != sink) {
			if (nv.d > nu.d + luv) {
				r = false;
				std::cerr<<"invalid labeling check for (u,v)  vu="<<vu<<std::endl;
				std::cerr<<"d(v="<<v<<")="<<nv.d<<" > "<<"d(u="<<u<<")="<<nu.d<<" + "<<luv<<std::endl;
			}
			if (nv.d == nu.d + luv && vu < nv.cur) {
				r = false;
				std::cerr<<"invalid current arc nv.cur="<<nv.cur<<" > vu="<<vu<<" but"<<std::endl;
				std::cerr<<"d(v="<<v<<")="<<nv.d<<" == "<<"d(u="<<u<<")="<<nu.d<<" + "<<luv<<std::endl;
			}
			//		}
		}
		if (!r) fflush(stderr);
		return r;
	}

	bool checkCurrent(NodeID v) {
		assert(v != source);
		if (nodes[v].first >= nodes[v+1].first) return true;
		if (nodes[v].d == n) return true;
		Node& nv = nodes[v];
		assert(nv.first <= nv.cur && nv.cur < nodes[v+1].first);
		ArcID vu = nv.cur; Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		if (isPseudoResidual(nu, nv, uv, vu, auv)) {
			if (nv.d != nu.d + 1) {
				std::cerr<<"v:"<<v<<", u:"<<u<<" ";
				std::cerr<<"dv="<<nv.d;
				std::cerr<<" != 1+du="<<(nu.d+1);
				std::cerr<<std::endl;
			}
			return checkValidCurArc(v) && nv.d == nu.d + 1;
		}
		return checkValidCurArc(v);
	}
};

#include "dfsbfs.cpp"
#endif
