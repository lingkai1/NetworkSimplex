
struct  bpq() { // bucket based priority queue
	NetworkMaxFlowSimplex& p;
	NodeID& nSentinel;
	std::vector<Node>& nodes;
	std::vector<NodeID> first;
	NodeID dmin; // min d in the queue
	NodeID dmax; // max d in the queue

	bpq()(NetworkMaxFlowSimplex& p) : p(p), nSentinel(p.nSentinel), nodes(p.nodes) {}

	void initialize() {
		first.resize(nSentinel+1); first.shrink_to_fit();
		forAllNodes(v) {
			Node& nv = nodes[v];
			nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;
		}
		for (Dist d = 0; d <= nSentinel; d++) {
			first[d] = nSentinel;
		}
		Node& nv = nodes[nSentinel];
		nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;
		first[nSentinel] = nSentinel;
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
#if defined(BPQ_CUR_T)
		//assert(p.isCurrent(v));
#endif
		Node& nv = nodes[v];
		NodeID w = first[k]; Node& nw = nodes[w];

		nv.bpq(Next) = w;
		nv.bpq(Prev) = nSentinel;
		nw.bpq(Prev) = v;
		first[k] = v;

		if (k < dmin)
			dmin = k;
		if (k > dmax)
			dmax = k;

	}

	void remove(NodeID v) {
		assert(contains(v));
		Node& nv = nodes[v]; Dist d = nv.d;
		NodeID u = nv.bpq(Prev); Node& nu = nodes[u];
		NodeID w = nv.bpq(Next); Node& nw = nodes[w];

		if (first[d] == v) {
			first[d] = w;
			nw.bpq(Prev) = nSentinel;
		}
		else {
			nu.bpq(Next) = w;
			nw.bpq(Prev) = u;
		}

		nv.bpq(Prev) = nv.bpq(Next) = UNDEF_NODE;

		while (first[dmin] == nSentinel && dmin <= dmax)
			dmin++;

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
		while (first[dmin] == nSentinel && dmin <= dmax)
			dmin++;
		if (dmin <= dmax) {
			NodeID v = first[dmin];
			assert(v != nSentinel);
			remove(v);
			return v;
		}
		else {
			return UNDEF_NODE;
		}
	}

	NodeID peekMin() {
		while (first[dmin] == nSentinel && dmin <= dmax)
			dmin++;
		if (dmin <= dmax) {
			return first[dmin];
		}
		else {
			return UNDEF_NODE;
		}
	}

	bool empty() {
		return peekMin() == UNDEF_NODE;
	}

	void flush() {
		/*forAllNodes(u) {
			if (contains(u))
				remove(u);
		}*/
		for (Dist d=dmin; d <= dmax; d++)
			while (first[d] != nSentinel) {
				remove(first[d]);
			}
		dmin = nSentinel;
		dmax = 0;
	}

#if defined(BPQ_RELABEL)
	// returns true if qr is empty (or does not need to be processed further in the lazy relabeling case)
	bool processed() {
		while (first[dmin] == nSentinel && dmin <= dmax) {
#if defined(LAZY_RELABEL)
			if (dmin >= p.qt.dmin) {
				//dmax = dmin;
				flush();
				return true;
			}
#endif
			dmin++;
		}
		if (dmin <= dmax) {
			return false;
		}
		else {
			return true;
		}
	}
#endif

} bpq();

