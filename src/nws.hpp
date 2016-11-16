#ifndef NWS_H_
#define NWS_H_

#include "types.hpp"

struct NWS {
	const static int FORMAT_DIMACS = 1;

	NodeID n; // number of nodes
	ArcID m; // number of arcs
	NodeID nMin; // smallest node id
	Node* nodes; // array of nodes
	Arc* arcs; // array of arcs
	Cap* cap; // array of capacities
	Node* source; // source node pointer
	Node* sink; // sink node pointer
	Flow flow; // flow value
	Node* sentinelNode; // end of the node list marker
	Arc* stopA; // used in forAllArcs

	//NWS(std::istream is, int format = FORMAT_DIMACS);
	//NWS(int n, int m);

};

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode( i ) ( (i) - nodes + nMin )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )



#endif
