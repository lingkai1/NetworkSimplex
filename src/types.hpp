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

struct Arc;
struct Node;

struct Arc
{
   union {
	   Cap resCap; // residual capasity
	   Cap cap;
   };
   Node& head; // arc head
   Arc& rev; // reverse arc
};

struct Node
{
   Arc* first; // first outgoing arc
   Arc* current; // current outgoing arc
   long d; // distance label
};


#endif
