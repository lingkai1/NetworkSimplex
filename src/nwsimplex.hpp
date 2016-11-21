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


#define COLOR_WHITE 0
#define COLOR_GREY  1
#define COLOR_BLACK 2
#define COLOR_BLUE 3


#define forAllNodes(u) for (NodeID u = 0; u <= nMax; u++)
#define forAllOutArcs(u, a, aStop) if (nodes[u].first != UNDEF_ARC)\
		for (NodeID a = nodes[u].first, aStop = nodes[u+1].first; a != aStop; a++)
#define forAllArcs(u, a, aStop) forAllNodes(u)\
		forAllOutArcs(u, a, aStop)


#define nNode( u ) ( u + nMin )
#define nArc( a )  ( ( a == NULL )? -1 : a )

class NWS {
public:
	const static int FORMAT_DIMACS = 1;

	int verbose = 1;
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
	std::vector<NodeID> bToRelabelFirst; // for each level d first node to relabel
	std::vector<ArcID> bPivotArcFirst; // for each level d first pivot arc
	std::stack<NodeID> toRelabel; // todo remove this and change implementation

	NWS(std::istream& is, int format = FORMAT_DIMACS);
	//NWS(int n, int m);
	~NWS();

	void solve();
	void buildInitialBasis();
	void prepareNextPivot();

	bool makeCur(NodeID v);
	void relabel(NodeID v);


	// only level0 of buckets is used now for plain simplex.
	// todo: change later
	void insertPivotArc(NodeID u, ArcID a);
	bool extractMinPivotArc(NodeID& u, ArcID& a);



	ArcID findArc(NodeID u, NodeID v, bool binarySearch = true);

#include "dfsbfs.hpp"

private:
	void initialize();
	void constructorDimacs(std::istream& is);
	void constructorDimacsFail(int lineNum, int code);

};

#include "dfsbfs.cpp"
#endif
