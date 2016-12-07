#ifndef bpq
#define bpq(a) listPivots
#define __cplusplus 201103L
#include "types.hpp"
#include "assert.h"
#include <vector>
#include <list>
#define Prev
#define Next
#endif

struct  bpq() { // bucket based priority queue
	NetworkMaxFlowSimplex& p;
	NodeID& nSentinel;
	std::vector<Node>& nodes;
	std::vector<Bucket> b;
	NodeID dmin; // min d in the queue
	NodeID dmax; // max d in the queue

	bpq()(NetworkMaxFlowSimplex& p) : p(p), nSentinel(p.nSentinel), nodes(p.nodes) {}

	void initialize() {
		b.resize(nSentinel+1); b.shrink_to_fit();
		forAllNodes(v) {
			Node& nv = nodes[v];
			nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;
		}
		for (Dist d = 0; d <= nSentinel; d++) {
			Bucket& bd = b[d];
			bd.first = nSentinel;
		}
		Node& nv = nodes[nSentinel];
		nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;
		Bucket& bv = b[nSentinel];
		bv.first = nSentinel;
		dmin = nSentinel;
		dmax = 0;
	}

	bool contains(NodeID v) {
		Node& nv = nodes[v];
		assert((nv.bpq(Prev) == UNDEF_NODE) == (nv.bpq(Next) == UNDEF_NODE));
		return nv.bpq(Prev) != UNDEF_NODE;
	}

	void insert(NodeID v, Dist k) {
		assert(!contains(v));
		assert(0 <= k); assert(k <= nSentinel);
		Node& nv = nodes[v]; Bucket& bv = b[k];
		NodeID w = bv.first; Node& nw = nodes[w];

		nv.bpq(Next) = w;
		nv.bpq(Prev) = nSentinel;
		nw.bpq(Prev) = v;
		bv.first = v;

		if (k < dmin)
			dmin = k;
		if (k > dmax)
			dmax = k;
	}

	void remove(NodeID v) {
		assert(contains(v));
		Node& nv = nodes[v]; Bucket& bv = b[nv.d];
		NodeID u = nv.bpq(Prev); Node& nu = nodes[u];
		NodeID w = nv.bpq(Next); Node& nw = nodes[w];

		if (bv.first == v) {
			bv.first = w;
			nw.bpq(Prev) = nSentinel;
		}
		else {
			nu.bpq(Next) = w;
			nw.bpq(Prev) = u;
		}

		nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;
	}

	void update(NodeID v, int k) {
		if (contains(v)) {
			remove(v);
			insert(v, k);
		}
	}

	void offer(NodeID v, int k) {
		if (!contains(v))
			insert(v, k);
	}

	NodeID extractMin() {
		while (b[dmin].first == nSentinel && dmin <= dmax)
			dmin++;
		if (dmin <= dmax) {
			NodeID v = b[dmin].first;
			remove(v);
			assert(b[dmin].first != v);
			return v;
		}
		else {
			return UNDEF_NODE;
		}
	}

	NodeID peekMin() {
		while (b[dmin].first == nSentinel && dmin <= dmax)
			dmin++;
		if (dmin <= dmax) {
			return b[dmin].first;
		}
		else {
			return UNDEF_NODE;
		}
	}

	bool empty() {
		return peekMin() == UNDEF_NODE;
	}

} bpq();

