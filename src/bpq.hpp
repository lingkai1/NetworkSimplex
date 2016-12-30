
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
		return nv.bpq(Prev) != UNDEF_NODE;
	}

	void insert(NodeID v, Dist k) {
		assert(!contains(v));
		Node& nv = nodes[v];
#if defined(RELABEL_OPT) && defined(BPQ_QR)
		k += 1;
		p.d[nv.d]--;
		nv.d = k;
		p.d[nv.d]++;
#endif
		assert(0 <= k); assert(k <= nSentinel);
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
		if (!contains(v)) return;
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

		nv.bpq(Prev) = UNDEF_NODE;

		while (first[dmin] == nSentinel && dmin <= dmax)
			dmin++;
		if (dmin > dmax) {
			dmin = nSentinel;
			dmax = 0;
		}
	}

	void flush() {
		while (dmin <= dmax) {
			while (first[dmin] != nSentinel)
				remove(first[dmin]);
			dmin++;
		}
		dmin = nSentinel;
		dmax = 0;
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
			dmin = nSentinel;
			dmax = 0;
			return UNDEF_NODE;
		}
	}

} bpq();

