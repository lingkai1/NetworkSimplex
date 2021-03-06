#ifndef NMFS_TYPES_H_
#define NMFS_TYPES_H_

#ifdef CAPACITY64
typedef unsigned long long Cap;
#else
typedef unsigned long Cap;
#endif

typedef Cap Flow;

typedef int NodeID;
typedef long long ArcID;
#define UNDEF_NODE (NodeID)-1
#define UNDEF_ARC (ArcID)-1

#define IN_NONE -1
#define IN_S 0
#define IN_T 1

typedef NodeID Dist;
#define INF_DIST (Dist)1000000000

struct Arc;
struct Node;


struct Arc {
	Cap resCap; // residual capacity
	NodeID head; // arc head
	ArcID rev; // reverse arc
};

struct Node {
	ArcID first; // first outgoing arc
	ArcID cur; // current outgoing arc. if this node is v we are interested in (u,v) the reverse arc of cur
	Dist d; // distance label

	// data to represent the basis tree (S or T) this node belongs to
	ArcID parent; // arc pointing to the parent of this node
	// (Doubly linked list in DFS order)
	NodeID next; // next node in the tree
	NodeID prev; // prev node in the tree
	NodeID size; // size of the subtree rooted at this node
	int tree; // current tree (simplex multiplier)

	// doubly linked to allow simple deletion
	NodeID qtPrev, qtNext;
	NodeID qrPrev, qrNext;
};

#define forAllNodes(u) for (NodeID u = 0; u < nSentinel; u++)
#define forAllOutArcs(u, i, iStop) if (nodes[u].first != UNDEF_ARC)\
		for (NodeID i = nodes[u].first, iStop = nodes[u+1].first; i != iStop; i++)
#define forAllArcs(u, i, iStop) forAllNodes(u)\
		forAllOutArcs(u, i, iStop)


#endif
