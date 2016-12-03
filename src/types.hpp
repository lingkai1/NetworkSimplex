#ifndef NWS_TYPES_H_
#define NWS_TYPES_H_

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

typedef int Dist;
#define INF_DIST (Dist)1000000000;

struct Arc;
struct Node;


struct Arc {
	union { // use anonymous union for alternate names -> different meanings for field reuse for space optimization
		Cap resCap; // residual capasity
		Cap cap;
	};
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
	NodeID stSize; // size of the subtree rooted at this node

	int tree; // current tree (simplex multiplier)


	// doubly linked to allow deletion inside
	NodeID bToRelabelNext, bToRelabelPrev;

	NodeID bPivotsNext, bPivotsPrev;

};

// push on one side pop on the other for FIFO order
struct BucketToRelabel {
	NodeID first, last;
};

struct BucketPivot {
	NodeID first, last;
};

#endif
