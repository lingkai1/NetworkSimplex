#include "nwsimplex.hpp"

using namespace std;


void NetworkMaxFlowSimplex::buildInitialBasis() {
	bToRelabel.resize(n+1);
	bPivots.resize(n+1);

	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.d = u != sink ? 0 : 1; // every node starts with d = 0 except sink (d = 1)
		nu.cur = nu.first;
		nu.tree = IN_NONE;
		nu.bToRelabelNext = nu.bToRelabelPrev = UNDEF_NODE;
		nu.bPivotsNext = nu.bPivotsPrev = UNDEF_NODE;
	}
	forAllArcs(u, i, is) {
		Arc& ai = arcs[i];
	}
	for (NodeID d = 0; d <= n; d++) {
		bToRelabel[d].first = bToRelabel[d].last = UNDEF_NODE;
		bPivots[d].first = bPivots[d].last = UNDEF_NODE;
	}

	vector<int> color(n); // todo: to optimize (reuse another field maybe)
	fill(color.begin(), color.end(), COLOR_WHITE);
	if (verbose >= 1) cout << "Initial BFS" << endl;
	bfs(source,
			[&](NodeID u){
		// pre node
		Node& nu = nodes[u];
		nu.tree = IN_S;
	},
	[&](NodeID u, bool leaf){ // post node
	},
	[&](NodeID u, ArcID i, int headColor){
		// process BFS arc
		Arc& ai = arcs[i];
		if (ai.resCap > 0 && headColor == COLOR_WHITE) { // forward residual tree arc
			if (ai.head != sink) {
				if (verbose >= 3) cout << "tree arc: " << u << "->"  << ai.head << " rescap: " << ai.resCap << endl;
				// set parent
				Node& z = nodes[ai.head];
				z.parent = ai.rev;
				return true; // accept this choice as a tree arc
			}
			else {
				// this arc is an initial pivot arc.
				if (verbose >= 3) cout << "pivot arc: " << u << "->"  << ai.head << " rescap: " << ai.resCap << endl;
				// don't do pivot insertion now.
				// do it after initial relabeling
				return false; // reject so that we do not search further than sink
			}
		}
		else if (headColor == COLOR_BLACK) { // backward tree arc
			return false;
		}
		else if (ai.resCap > 0) { // residual non-tree arc
			return false;
		}
		else {
			return false; // reject other types of arcs
		}
	},
	color);

	//do DFS to initialize subtree fields of each vertex
	NodeID prevu = source;
	nodes[source].prev = UNDEF_NODE;
	fill(color.begin(), color.end(), COLOR_WHITE);
	if (verbose >= 1) cout << "Initial DFS" << endl;
	dfs_i(source,
			[&](NodeID u){
		// pre node
		// link nodes in DFS order as i doubly linked list
		Node& prevnu = nodes[prevu];
		prevnu.next = u;
		Node& nu = nodes[u];
		nu.prev = prevu;
		nu.stSize = 1; // initialize subtree size
		prevu = u;
	},
	[&](NodeID u, bool leaf){
		// post node
		if (u == source) {
			nodes[prevu].next = UNDEF_NODE;
		}
	},
	[&](NodeID u, ArcID i, int headColor) {
		// arc choice (pre arc)
		Arc& ai = arcs[i];
		ArcID r = ai.rev;
		NodeID v = ai.head; Node& nv = nodes[v];
		if (headColor == COLOR_WHITE && v != sink && nv.parent == r) // follow arcs (u,v) of S (such that u=parent(v))
			return true;
		else
			return false;
	},
	[&](NodeID u, ArcID i) {
		// post arc
		Arc& ai = arcs[i];
		Node& nu = nodes[u];
		Node& nh = nodes[ai.head];
		nu.stSize += nh.stSize; // compute subtree sizes in post visit order
	},
	color);


	Node& nt = nodes[sink];
	nt.parent = UNDEF_ARC;
	nt.next = UNDEF_NODE;
	nt.stSize = 1;
	nt.tree = IN_T;

	Node& ns = nodes[source];
	ns.parent = UNDEF_ARC;

}



// main network simplex function
void NetworkMaxFlowSimplex::solve() {
#ifdef USE_STATS_TIME
	t_solve = timer();
#endif
#ifdef USE_STATS_COUNT
	c_basisIter = 0;
	c_basisGrowS = 0;
	c_basisGrowT = 0;
	c_basisNoChange = 0;
	c_pivotsInserted = 0;
	c_pivotsDeleted = 0;
	c_makeCur = 0;
	c_relabel = 0;
#endif

	flow = 0;

#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer();
#endif
	buildInitialBasis();
#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer() - t_buildInitialBasis;
#endif


#if defined(GGH_SIMPLEX)
	globalRelabelFreq = 0.2;
	globalRelabelThreshold = n;
	globalRelabelWork = 0;
	doGlobalRelabel();
#endif

	// search initial pivots
	forAllSubTree(source, [&](NodeID u){
		Node& nu = nodes[u];
		if (hasOutPivots(u)) {
			assert(!pivotsContains(u));
			pivotsInsert(u, nu.d);
		}
	});

	double t_refresh=timer();
	while (true) {
		assert(nodes[source].parent == UNDEF_ARC);
		assert(nodes[sink].parent == UNDEF_ARC);
		if (timer() - t_refresh > 1.0) { // print flow every sec
			cout << "makecurs: " << c_makeCur << endl;
			cout << "relabels: " << c_relabel << endl;
			cout << "Flow: " << flow << endl;
			fflush(stdout);
			t_refresh = timer();
		}


		if (verbose >= 3) {
			printS();
			printT();
		}

		if (verbose >= 1) cout << "Selecting Pivot..." << endl;
		NodeID v, w;
		ArcID vw, wv;
		vw = pivotsExtractMin();
		if (vw == UNDEF_ARC) {
			// no pivot available.
			// Simplex method ended
			if (verbose >= 1) cout << "No pivots left. Simplex method ended." << endl;

#ifdef USE_STATS_TIME
			t_solve = timer() - t_solve;
#endif
			return;
		}

		Arc& avw = arcs[vw];
		wv = avw.rev;
		Arc& awv = arcs[wv];
		v = awv.head;
		w = avw.head;

		if (verbose >= 1) cout << "Pivot (entering) vw arc: " << v << "->" << w << endl;

		assert(nodes[v].tree == IN_S);
		assert(nodes[w].tree == IN_T);

		if (verbose >= 1) cout << "Searching bottleneck..." << endl;
		NodeID u;
		Cap delta;

		NodeID x, y;
		ArcID xy, yx;

		int xyTree = IN_NONE;
		// initialize at pivot arc
		xy = vw;
		delta = avw.resCap;
		// path s->v
		u = v;
		while (u != source) {
			Node& nu = nodes[u];
			Arc& ar = arcs[nu.parent]; // (u,p)
			ArcID i = ar.rev; Arc& ai = arcs[i]; // (p,u)
			NodeID p = ar.head;
			if (ai.resCap <= delta) { // <= so we save the closest to the source
				delta = ai.resCap;
				xy = i;
				xyTree = IN_S;
			}
			u = p;
		}
		// path w->t
		u = w;
		while (u != sink) {
			Node& nu = nodes[u];
			ArcID i = nu.parent; Arc& ai = arcs[nu.parent]; // (u,p)
			NodeID p = ai.head;
			if (ai.resCap < delta) { // need strict <
				delta = ai.resCap;
				xy = i;
				xyTree = IN_T;
			}
			u = p;
		}


		Arc& axy = arcs[xy];
		yx = axy.rev;
		Arc& ayx = arcs[yx];
		x = ayx.head;
		y = axy.head;
		if (verbose >= 1) cout << "Bottleneck (leaving) xy arc: " << x << "->" << y << " delta: " << delta << endl;

		if (verbose >= 1) cout << "Augmenting flow..." << endl;
		// pivot (v,w)
		avw.resCap -= delta;
		awv.resCap += delta;
		// path s->v
		u = v;
		while (u != source) {
			Node& nu = nodes[u];
			Arc& ar = arcs[nu.parent]; // (u,p)
			ArcID i = ar.rev;
			Arc& ai = arcs[i]; // (p,u)
			NodeID p = ar.head;

			ai.resCap -= delta;
			ar.resCap += delta;

			u = p;
		}
		// path w->t
		u = w;
		while (u != sink) {
			Node& nu = nodes[u];
			ArcID i = nu.parent; Arc& ai = arcs[nu.parent]; // (u,p)
			NodeID p = ai.head;
			Arc& ar = arcs[ai.rev];

			ai.resCap -= delta;
			ar.resCap += delta;

			u = p;
		}
		flow += delta;
		assert(axy.resCap == 0);


		if (verbose >= 1) cout << "Updating basis..." << endl;
		if (xyTree == IN_S) {
			if (verbose >= 1) cout << "xy in S => T grows" << endl;

			// Q the subtree rooted at y in S moves from S to T
			// R will represent Q rooted at v in T
			NodeID nLast;
			getSubtreeLastNode(y, nLast);
			deleteSubtree(source, y, nLast); // disconnect Q from S
			changeRoot(y, v, nLast); // this transforms Q into R (updates parents and the doubly linked list)
			nodes[v].parent = vw;
			addSubtreeAsChild(sink, v, nLast, w); // connect R to T


			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_T;
			});

			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				if (pivotsContains(u))
					pivotsDelete(u);

				// scan for new pivots to u
				forAllOutArcs(u, uv, is) {
					Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					if (isInPivot(nv, avu)) {
						if (!pivotsContains(v))
							pivotsInsert(v, nv.d);
					}
				}
			});

#if defined(GGH_SIMPLEX)
			makeCur(y);
			//forAllNodes(v) makeCur(v);
			//globalRelabelWork = 0;
			while (!toRelabelEmpty())
				relabel(toRelabelExtractMin());
#endif
#ifdef USE_STATS_COUNT
			c_basisGrowT++;
#endif
		}
		else if (xyTree == IN_T) {
			if (verbose >= 1) cout << "xy in T => S grows" << endl;

			// Q the subtree rooted at x in T moves from T to S
			// R will represent Q rooted at w in S
			NodeID nLast;
			getSubtreeLastNode(x, nLast);
			deleteSubtree(sink, x, nLast); // disconnect Q from T
			changeRoot(x, w, nLast); // this transforms Q into R (updates parents and the doubly linked list)
			nodes[w].parent = wv;
			addSubtreeAsChild(source, w, nLast, v); // connect R to S

			forAllSubTree(w, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_S;
			});

			forAllSubTree(w, [&](NodeID u){
				Node& nu = nodes[u];

				/*nu.tree = IN_T;
				forAllOutArcs(u, uv, is) {
					Node& nu = nodes[u]; Arc& auv = arcs[uv];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					NodeID v = auv.head; Node& nv = nodes[v];
					nu.tree = IN_S;
					if (pivotsContains(v))
						if (!hasOutPivots(v))
							pivotsDelete(v);
				}
				nu.tree = IN_S;*/

				if (hasOutPivots(u))
					pivotsInsert(u, nu.d);
			});

#if defined(GGH_SIMPLEX)
			makeCur(y);
			//forAllNodes(v) makeCur(v);
			//globalRelabelWork = 0;
			while (!toRelabelEmpty())
				relabel(toRelabelExtractMin());
#endif
#ifdef USE_STATS_COUNT
			c_basisGrowS++;
#endif
		}
		else {
			if (verbose >= 1) cout << "xy == vw => Unchanged basis" << endl;
			// pivot is the leaving arc. keep the same basis
			if (!hasOutPivots(v))
				pivotsDelete(v);

#if defined(GGH_SIMPLEX)
			makeCur(y);
			//forAllNodes(v) makeCur(v);
			//globalRelabelWork = 0;
			while (!toRelabelEmpty())
				relabel(toRelabelExtractMin());
#endif

#ifdef USE_STATS_COUNT
			c_basisNoChange++;
#endif
		}

#ifdef USE_STATS_COUNT
		c_basisIter++;
#endif
	}

}

bool NetworkMaxFlowSimplex::hasOutPivots(NodeID u) {
	Node& nu = nodes[u];
	forAllOutArcs(u, uv, is) {
		Arc& auv = arcs[uv];
		NodeID v = auv.head; Node& nv = nodes[v];
		if (isOutPivot(nv, auv))
			return true;
	}
	return false;
}
bool NetworkMaxFlowSimplex::hasInPivots(NodeID u) {
	Node& nu = nodes[u];
	forAllOutArcs(u, uv, is) {
		Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
		ArcID vu = auv.rev; Arc& avu = arcs[vu];
		if (isInPivot(nv, avu))
			return true;
	}
	return false;
}

void NetworkMaxFlowSimplex::getSubtreeLastNode(NodeID u, NodeID& uLast) {
	uLast = u;
	for (NodeID k = 0; k < nodes[u].stSize - 1; k++) {
		assert(uLast != UNDEF_NODE);
		uLast = nodes[uLast].next;
	}
}

// delete from Q a subtree rooted at u
// output the subtree last element in uLast so we can reuse it without retraversing the tree again
void NetworkMaxFlowSimplex::deleteSubtree(NodeID q, NodeID u, NodeID uLast) {
	//getSubtreeLastNode(u, uLast);// should be done outside so we can save uLast for reuse when adding subtrees

	Node& nFirst = nodes[u], & nLast = nodes[uLast];

	//start at u, follow next pointers for the SIZE steps to find the last vertex, remove the sublist corresponding to the subtree
	if (nFirst.prev != UNDEF_NODE)
		nodes[nFirst.prev].next = nLast.next;
	if (nLast.next != UNDEF_NODE)
		nodes[nLast.next].prev = nFirst.prev;
	nFirst.prev = nLast.next = UNDEF_NODE;

	//Then go from u up to the root of Q updating subtree sizes.
	//At this point of the algorithm, the parents on q,u path are old, so it is still possible to do it.
	//Updating the subtree sizes now is probably faster if graphs are large and dense.
	//Another solution would be to wait and update only roots of subtrees who moved later. But it would require exploring a whole subtree unlike this.
	NodeID uSize = nodes[u].stSize;
	NodeID z = u;
	while (z != q) {
		z = arcs[nodes[z].parent].head;
		nodes[z].stSize -= uSize;
	}

	if (verbose >= 4) {cout<<"SUBTREE deleted: "; printSubTree(u); cout << " from: "; printSubTree(q); cout<<endl;}
}

// add to R a subtree rooted at p as a child of u
void NetworkMaxFlowSimplex::addSubtreeAsChild(NodeID r, NodeID p, NodeID pLast, NodeID u) {

	//insert the subtree list in the tree list immediately after u
	nodes[p].prev = u;
	nodes[pLast].next = nodes[u].next;
	if (nodes[u].next != UNDEF_NODE)
		nodes[nodes[u].next].prev = pLast;
	nodes[u].next = p;

	// go up to R's root updating subtree sizes.
	//At this point of the algorithm, the parents on r,u path are new (going to r)
	int pSize = nodes[p].stSize;
	NodeID z = u;
	while (z != r) {
		nodes[z].stSize += pSize;
		z = arcs[nodes[z].parent].head;
	}
	nodes[r].stSize += pSize;

	if (verbose >= 4) {cout<<"SUBTREE added to: "; printSubTree(r); cout<<endl;}
}

// transforms Q into R, outputs the last node of R
void NetworkMaxFlowSimplex::changeRoot(NodeID q, NodeID r, NodeID& rLast) {
	ArcID ip, itmp;
	NodeID p, pLast;
	NodeID z = r;
	ip = nodes[z].parent;
	getSubtreeLastNode(z, pLast);
	deleteSubtree(q, z, pLast);
	while (z != q) {
		p = arcs[ip].head;
		getSubtreeLastNode(p, pLast);
		deleteSubtree(q, p, pLast);
		//make z the parent of p
		assert(arcs[arcs[nodes[z].parent].rev].head == z);
		itmp = nodes[p].parent;
		nodes[p].parent = arcs[ip].rev;
		ip = itmp;
		assert(arcs[nodes[p].parent].head == z);
		assert(arcs[arcs[nodes[p].parent].rev].head == p);
		addSubtreeAsChild(r, p, pLast, z);
		z = p;
	}

	getSubtreeLastNode(r, rLast);
}

template <typename PrintNode>
void NetworkMaxFlowSimplex::printSubTree(NodeID root, const PrintNode& print) {
	forAllSubTree(root, [&](NodeID u){
		Node& nu = nodes[u];
		if (u == root) {
			print(root);
		}
		else {
			cout << ",";
			print(u);
		}
	});
}
void NetworkMaxFlowSimplex::printSubTree(NodeID root) {
	return printSubTree(root, [&](NodeID u){ cout << u; });
}

void NetworkMaxFlowSimplex::printS() {
	cout << "S tree rooted at " << source << endl;
	printSubTree(source, [&](NodeID u){ cout << u << " (d=" << nodes[u].d << ")"; });
	cout << endl;
	cout << "S arcs:" << endl;
	forAllSubTree(source, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (u != source && r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << p << "->" << u << " rc: " << ai.resCap << endl;
		}
	});
}
void NetworkMaxFlowSimplex::printT() {
	cout << "T tree rooted at " << sink << endl;
	printSubTree(sink, [&](NodeID u){ cout << u << " (d=" << nodes[u].d << ")"; });
	cout << endl;
	cout << "T arcs:" << endl;
	forAllSubTree(sink, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (u != sink && r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << u << "->" << p << " rc: " << ai.resCap << endl;
		}
	});
}


// todo: at first put all vertices in the bucket, then be smarter
bool NetworkMaxFlowSimplex::pivotsContains(NodeID v) {
	Node& nv = nodes[v];
	return !(nv.bPivotsPrev == UNDEF_NODE && nv.bPivotsNext == UNDEF_NODE && bPivots[nv.d].first != v);
}
void NetworkMaxFlowSimplex::pivotsInsert(NodeID v, Dist d) {
	assert(!pivotsContains(v));
	//cout<<"PIVOT INSERT "<<v<<" at d="<<d<<endl;
	Node& nv = nodes[v];
	assert(nv.tree == IN_S);
	BucketPivot& ll = bPivots[d];
	if (ll.last != UNDEF_NODE) {
		nodes[ll.last].bPivotsNext = v;
		nv.bPivotsPrev = ll.last;
		nv.bPivotsNext = UNDEF_NODE;
		ll.last = v;
	}
	else {
		nv.bPivotsNext = nv.bPivotsPrev = UNDEF_NODE;
		ll.first = ll.last = v;
	}
#ifdef USE_STATS_COUNT
	c_pivotsInserted++;
#endif
}
NodeID NetworkMaxFlowSimplex::pivotsDelete(NodeID v) {
	Node& nv = nodes[v]; Dist d = nv.d;
	//cout<<"PIVOT DELETE "<<v<<" at d="<<d<<endl;
	BucketPivot& ll = bPivots[d];
	if (nv.bPivotsPrev != UNDEF_NODE)
		nodes[nv.bPivotsPrev].bPivotsNext = nv.bPivotsNext;
	else
		ll.first = nv.bPivotsNext;
	if (nv.bPivotsNext != UNDEF_NODE)
		nodes[nv.bPivotsNext].bPivotsPrev = nv.bPivotsPrev;
	else
		ll.last = nv.bPivotsPrev;
	nv.bPivotsNext = nv.bPivotsPrev = UNDEF_NODE;
#ifdef USE_STATS_COUNT
	c_pivotsDeleted++;
#endif
	return v;
}
// todo: add pointer to the first non empty bucket level
ArcID NetworkMaxFlowSimplex::pivotsExtractMin() {
	for (Dist d = 0; d <= n; d++) {
		BucketPivot& ll = bPivots[d];
		while (ll.first != UNDEF_NODE) {
			NodeID v = ll.first; Node& nv = nodes[v];
			assert(nv.tree == IN_S);
			//scan the arc list to find an outgoing pivot arc.
			forAllOutArcs(v, vu, is) {
				Arc& avu = arcs[vu];
				NodeID u = avu.head; Node& nu = nodes[u];
				ArcID uv = avu.rev; Arc& auv = arcs[uv];
				if (isOutPivot(nu, avu))
					return vu;
			}
			// could not find at least one pivot arc for this pivot vertex
			//assert(false); // there can be such vertices
			if (pivotsDelete(ll.first) == ll.first) {
				assert(false);
			}
		}
	}
	// no pivots found
	return UNDEF_ARC;
}


#if defined(GGH_SIMPLEX)
bool NetworkMaxFlowSimplex::toRelabelContains(NodeID v) {
	Node& nv = nodes[v];
	return !(nv.bToRelabelPrev == UNDEF_NODE && nv.bToRelabelNext == UNDEF_NODE && bToRelabel[nv.d].first != v);
}
void NetworkMaxFlowSimplex::toRelabelInsert(NodeID v, Dist d) {
	Node& nv = nodes[v];
	BucketToRelabel& ll = bToRelabel[d];
	if (ll.last != UNDEF_NODE) {
		nodes[ll.last].bToRelabelNext = v;
		nv.bToRelabelPrev = ll.last;
		nv.bToRelabelNext = UNDEF_NODE;
		ll.last = v;
	}
	else {
		nv.bToRelabelNext = nv.bToRelabelPrev = UNDEF_NODE;
		ll.first = ll.last = v;
	}
}
NodeID NetworkMaxFlowSimplex::toRelabelDelete(NodeID v) {
	Node& nv = nodes[v]; Dist d = nv.d;
	BucketToRelabel& ll = bToRelabel[d];
	if (nv.bToRelabelPrev != UNDEF_NODE)
		nodes[nv.bToRelabelPrev].bToRelabelNext = nv.bToRelabelNext;
	else
		ll.first = nv.bToRelabelNext;
	if (nv.bToRelabelNext != UNDEF_NODE)
		nodes[nv.bToRelabelNext].bToRelabelPrev = nv.bToRelabelPrev;
	else
		ll.last = nv.bToRelabelPrev;
	nv.bToRelabelNext = nv.bToRelabelPrev = UNDEF_NODE;
	return v;
}
// todo: add pointer to the first non empty bucket level
NodeID NetworkMaxFlowSimplex::toRelabelExtractMin() {
	for (Dist d = 0; d <= n; d++) {
		BucketToRelabel& ll = bToRelabel[d];
		if (ll.first != UNDEF_NODE) {
			if (nodes[ll.first].d != d) {
				cout << "ll.first.d: " << nodes[ll.first].d << " != d:" << d << endl;
			}
			assert(nodes[ll.first].d == d);
			return toRelabelDelete(ll.first);
		}
	}
	assert(false);
	return UNDEF_NODE;
}
bool NetworkMaxFlowSimplex::toRelabelEmpty() {
	for (Dist d = 0; d <= n; d++) {
		BucketToRelabel& ll = bToRelabel[d];
		if (ll.first != UNDEF_NODE)
			return false;
	}
	return true;
}

bool NetworkMaxFlowSimplex::makeCur(NodeID v) {
#ifdef USE_STATS_COUNT
	c_makeCur++;
#endif
	if (v == source) return true;
	Node& nv = nodes[v];
	if (nv.first >= nodes[v+1].first) // no arc incident arc to v
		return true;
	while (nv.cur < nodes[v+1].first) {
		ArcID vu = nv.cur; Arc& avu = arcs[vu];
		Node& nu = nodes[avu.head];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		Dist luv = isResidualOrTreeArc(nu, nv, uv, vu, auv) ? 1 : INF_DIST;
		if (nv.d == nu.d + luv) {
			if (verbose >= 3) cout << v << " (d=" << nv.d << ") made current to " << avu.head << endl;
			return true;
		}
		nv.cur++;
	}
	// add v to to_relabel list
	if (nv.d < n)
		toRelabelInsert(v, nv.d);
	return false;
}

void NetworkMaxFlowSimplex::relabel(NodeID v) {
	globalRelabelWork++;
	if (globalRelabelWork * globalRelabelFreq > globalRelabelThreshold) {
		doGlobalRelabel();
		globalRelabelWork = 0;
		assert(toRelabelEmpty());
		return;
	}
#ifdef USE_STATS_COUNT
	c_relabel++;
#endif
	Node& nv = nodes[v];
	assert(nv.d >= 0);
	assert(nv.d < n);
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist newD = INF_DIST;
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		Dist luv = isResidualOrTreeArc(nu, nv, uv, vu, auv) ? 1 : INF_DIST;
		if (nu.d + luv < newD) {
			newD = nu.d + luv;
			nv.cur = vu;
		}
	}

	if (newD <= nv.d) cerr << "For v:"<<v<<" newD:"<<newD<<" <= d(v):"<<nv.d<<endl;
	assert(newD > nv.d);

	if (newD >= n) {
		nv.d = newD;
		if (verbose >= 3) cout<<"d("<<v<<") >= n so delete it"<<endl;
		return;
	}

	if (pivotsContains(v)) {
		pivotsDelete(v);
		pivotsInsert(v, newD);
	}
	nv.d = newD; // increase key of v

	if (verbose >= 3) cout<<v<<" (d="<<nodes[v].d<<") relabeled now current to "<<arcs[nv.cur].head<<endl;
	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev;
		if (nu.cur == uv) {
			if (verbose >= 3) cout << u << " (d=" << nu.d << ") made non-current by relabeling of " << v << endl;
			makeCur(u);
		}
	}

}

void NetworkMaxFlowSimplex::doGlobalRelabel() {
	static int k = 0; k++; cout << "Global relabeling: "<<k<<endl;

	// empty toRelabel list
	forAllNodes(u) {
		if (toRelabelContains(u))
			toRelabelDelete(u);
	}

	// do BFS from source to compute exact labels
	nodes[source].d = 0;
	vector<int> color(n);
	fill(color.begin(), color.end(), COLOR_WHITE);
	bfs(source,
			[&](NodeID u){ // pre node
	},
	[&](NodeID u, bool leaf){ // post node
	},
	[&](NodeID u, ArcID uv, int headColor){ // process arc
		//cout<<"BFS uv="<<uv<<endl;
		Node& nu = nodes[u];
		Arc& auv = arcs[uv]; ArcID vu = auv.rev;
		NodeID v = auv.head; Node& nv = nodes[v];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
			if (nv.d > nu.d + 1 || headColor == COLOR_WHITE) {
				if (pivotsContains(v)) {
					pivotsDelete(v);
					pivotsInsert(v, nu.d + 1);
				}
				nv.d = nu.d + 1;
				nv.cur = vu;
			}
			else if (nv.d == nu.d + 1) {
				if (nv.cur > vu) {
					nv.cur = vu;
				}
			}
			return headColor == COLOR_WHITE && u != sink; // respect BFS order and do not search further than sink
		}
		else
			return false;

	},
	color);

	forAllNodes(v) {
		forAllOutArcs(v, vu, is) {
			Node& nv = nodes[v];
			Arc& avu = arcs[vu]; NodeID u = avu.head; Node& nu = nodes[u];
			ArcID uv = avu.rev; Arc& auv = arcs[uv];
			if (nv.d == nu.d + 1 && isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
				if (nv.cur != vu) {
					cout<<"v="<<v<<endl;
					cout<<"nv.cur="<<nv.cur<<endl;
					cout<<"vu="<<vu<<endl;
				}
				assert(nv.cur == vu);
				break;
			}
		}
	}
}


#endif

//static double t_=timer(); if (timer()-t_>1.0) { cout<<"event"<<endl; t_=timer(); }


void NetworkMaxFlowSimplex::testBFS() {
	vector<int> color(n, COLOR_WHITE);

	cout << "BFS" << endl;
	bfs(source,
			[&](NodeID u){
		cout << "pre:  " << u << endl;
	},
	[&](NodeID u, bool leaf){
		cout << "post: " << u << endl;
	},
	[&](NodeID u, ArcID i, int headColor){
		switch (headColor) {
		case COLOR_WHITE:
			cout << "forward arc: " << u << "->"  << arcs[i].head << endl;
			break;
		case COLOR_GREY:
			cout << "cross arc: " << u << "->"  << arcs[i].head << endl;
			break;
		case COLOR_BLACK:
			cout << "backward arc: " << u << "->"  << arcs[i].head << endl;
			break;
		default:
			break;
		}
		return true;
	}, color);

}

void NetworkMaxFlowSimplex::testDFS() {
	vector<int> color(n, COLOR_WHITE);

	cout << "DFS" << endl;
	dfs_i(source,
			[&](NodeID u){
		cout << "pre:  " << u << endl;
	},
	[&](NodeID u, bool leaf){ cout << "post: " << u << endl; },
	[&](NodeID u, ArcID i, int headColor){
		switch (headColor) {
		case COLOR_WHITE:
			cout << "forward pre arc: " << u << "->"  << arcs[i].head << endl;
			break;
		case COLOR_BLACK:
			cout << "cross pre arc: " << u << "->"  << arcs[i].head << endl;
			break;
		case COLOR_GREY:
			cout << "backward pre arc: " << u << "->"  << arcs[i].head << endl;
			break;
		default:
			break;
		}
		return true;
	},
	[&](NodeID u, ArcID a){
		cout << "post arc: " << u << "->"  << arcs[a].head << endl;
	}, color);
}


