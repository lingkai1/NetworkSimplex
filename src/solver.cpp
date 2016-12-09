#include "nwsimplex.hpp"

using namespace std;

// initialize the algorithm
// build S = V - {t} and T={t} in a bfs
// do dfs to construct an initial dfs order of nodes
// do a global relabeling to initialize d labels. (this could be optimized and done in the first bfs)
// search for initial pivots
void NetworkMaxFlowSimplex::buildInitialBasis() {
	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.d = u != sink ? 0 : 1; // every node starts with d = 0 except sink (d = 1)
		nu.cur = nu.first;
		nu.tree = IN_NONE;
	}
	forAllArcs(u, i, is) {
		Arc& ai = arcs[i];
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

	nodes[source].parent = UNDEF_ARC;

	listPivots.initialize();
	listRelabel.initialize();

	dc.resize(nSentinel+1);
	fill(dc.begin(), dc.end(), 0);
	forAllNodes(u) {
		dc[nodes[u].d]++;
	}

#if defined(GLOBAL_RELABEL)
	globalRelabelFreq = 0.2;
	globalRelabelThreshold = nSentinel;
	globalRelabelWork = 0;
#endif

	doGlobalRelabel();

	// search initial pivots
	forAllSubTree(source, [&](NodeID u){
		Node& nu = nodes[u];
		if (hasOutPivots(u)) {
			listPivots.insert(u, nu.d);
		}
	});

	forAllNodes(u) assert(checkValidCurArc(u));
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
	c_globalupdate = 0;
	c_gap = 0;
#endif

	flow = 0;

#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer();
#endif
	buildInitialBasis();
#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer() - t_buildInitialBasis;
#endif


	// main loop
	while (true) {
		// print some debug info every sec
		static double t_loopInfo = timer(); if (timer() - t_loopInfo > 1.0) {
			//cout << "makeCur(): " << c_makeCur << endl;
			//cout << "relabel(): " << c_relabel << endl;
			cout << "Flow: " << flow << endl;
			fflush(stdout); t_loopInfo = timer();
		}


		if (verbose >= 3) {
			printS();
			printT();
		}

		if (verbose >= 1) cout << "Selecting Pivot..." << endl;
		NodeID v, w;
		ArcID vw, wv;
		vw = extractMinPivot(); // extract a pivot with min label
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

		assert(nodes[v].d < nSentinel);

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

		/*if (vw == 455138) {
			cerr<<"vw=="<<455138<<endl;
			cerr<<"v:"<<v<<endl;
			cerr<<"w:"<<w<<endl;
			cerr<<"x:"<<x<<endl;
			cerr<<"y:"<<y<<endl;
			cerr<<"xy:"<<xy<<endl;
			cerr<<"tree(xy):"<<(xyTree==IN_S?"S":xyTree==IN_T?"T":"SxT")<<endl;
		}*/


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

			// move subtree from S to T
			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_T;
			});

			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				if (listPivots.contains(u))
					listPivots.remove(u);

				// scan for new pivots to u
				forAllOutArcs(u, uv, is) {
					Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					if (isInPivot(nv, avu)) {
						listPivots.offer(v, nv.d); // costly
					}
				}
			});

#if defined(GGT_RELABEL)

			nodes[v].cur = nodes[v].first; makeCur(v);
			makeCur(y); // only y can be made non current (if its current arc was xy (yx))
			//globalRelabelWork = 0;
			while (!listRelabel.processed())
				relabel(listRelabel.extractMin());
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

#if defined(FORCE_STRICT_PIVOTS)
				// if we strictly force the pivot list to contain vertices with at least 1 pivot.
				// we need to check neighbors, but slower
				forAllOutArcs(u, uv, is) {
					Node& nu = nodes[u]; Arc& auv = arcs[uv];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					NodeID v = auv.head; Node& nv = nodes[v];
					if (listPivots.contains(v))
						if (!hasOutPivots(v))
							listPivots.remove(v);
				}
#else
				// else we don't check. neighbors might not have pivots anymore
#endif

				if (hasOutPivots(u))
					listPivots.insert(u, nu.d);  // costly
			});

#if defined(GGT_RELABEL)

			nodes[v].cur = nodes[v].first; makeCur(v);
			makeCur(y);
			//globalRelabelWork = 0;
			while (!listRelabel.processed())
				relabel(listRelabel.extractMin());
#endif
#ifdef USE_STATS_COUNT
			c_basisGrowS++;
#endif
		}
		else {
			if (verbose >= 1) cout << "xy == vw => Unchanged basis" << endl;
			// pivot is the leaving arc. keep the same basis
			if (!hasOutPivots(v))
				listPivots.remove(v);

#if defined(GGT_RELABEL)

			nodes[v].cur = nodes[v].first; makeCur(v);
			makeCur(y);
			//globalRelabelWork = 0;
			while (!listRelabel.processed())
				relabel(listRelabel.extractMin());
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

// returns true if u has outgoing pivots
// u is assumed to be in S (it doesn't have to be)
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
// returns true if u has incoming pivots
// u is assumed to be in T (it doesn't have to be)
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


// returns a v from the pivots list with minimal d and deletes it from the list
ArcID NetworkMaxFlowSimplex::extractMinPivot() {
	do {
		NodeID v = listPivots.peekMin();
		if (v == UNDEF_NODE) // no pivots found
			return UNDEF_ARC;
		Node& nv = nodes[v];
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
#if defined(FORCE_STRICT_PIVOTS)
		assert(false); // cannot happen
#else
		listPivots.remove(v);
#endif
	} while (true);
}


// go to the last node of a subtree in DFS order
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




#if defined(GGT_RELABEL)

bool NetworkMaxFlowSimplex::checkValidCurArc(NodeID v) {
	bool r = true;
	forAllOutArcs(v, vu, is) {
		Node& nv = nodes[v];
		Arc& avu = arcs[vu]; NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		Dist luv = isResidualOrTreeArc(nu, nv, uv, vu, auv) ? 1 : INF_DIST;
		if (nv.d > nu.d + luv) {
			r = false;
			cerr<<"invalid labeling check for (u,v)  vu="<<vu<<endl;
			cerr<<"d(v="<<v<<")="<<nv.d<<" > "<<"d(u="<<u<<")="<<nu.d<<" + "<<luv<<endl;
		}
		if (nv.d == nu.d + luv && vu < nv.cur) {
			r = false;
			cerr<<"invalid current arc nv.cur="<<nv.cur<<" > vu="<<vu<<" but"<<endl;
			cerr<<"d(v="<<v<<")="<<nv.d<<" == "<<"d(u="<<u<<")="<<nu.d<<" + "<<luv<<endl;
		}
	}
	if (!r) fflush(stderr);
	return r;
}


// try to make v current, find and set its new current arc
// if v is already current, does nothing
// if v cannot be made current it is added to the toRelabel list
bool NetworkMaxFlowSimplex::makeCur(NodeID v) {
#ifdef USE_STATS_COUNT
	c_makeCur++;
#endif
	if (v == source) return true;
	Node& nv = nodes[v];
	if (nv.first >= nodes[v+1].first) // no arc incident arc to v
		return true;

	checkValidCurArc(v);

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
	if (nv.d < nSentinel) {
		listRelabel.offer(v, nv.d);
	}
	return false;
}

// relabel v
// this must result in a strict increase of v's d label based on the assumption that v's current arc was of minimal index
void NetworkMaxFlowSimplex::relabel(NodeID v) {
#if defined(GLOBAL_RELABEL)
	// try global update if a lot of work has been done
	globalRelabelWork++;
	if (globalRelabelWork * globalRelabelFreq > globalRelabelThreshold) {
		doGlobalRelabel();
		globalRelabelWork = 0;
		return;
	}
#endif
#ifdef USE_STATS_COUNT
	c_relabel++;
#endif
	Node& nv = nodes[v];
	assert(nv.d >= 0);
	assert(nv.d < nSentinel);
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist oldD = nv.d;
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

	if (newD >= nSentinel) {
		listPivots.update(v, nSentinel);
		dc[nv.d]--; dc[nSentinel]++;
		nv.d = nSentinel;
		return;
	}

	listPivots.update(v, newD);
	dc[nv.d]--; dc[newD]++;
	nv.d = newD; // increase key of v

	//if (v == 81801) { cerr<<"relabel(81801) from d="<<oldD<<" to d="<<newD<<endl;assert(checkValidCurArc(v)); }

	if (verbose >= 3) cout<<v<<" (d="<<nodes[v].d<<") relabeled now current to "<<arcs[nv.cur].head<<endl;
	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev;
		if (nu.cur == uv) {
			if (verbose >= 3) cout << u << " (d=" << nu.d << ") made non-current by relabeling of " << v << endl;
			makeCur(u); // if fail u is added in toRelabel list in makeCur()
		}
	}

#if defined(GAP_RELABEL)
	if (dc[oldD] == 0) {
		gap(oldD);
		return;
	}
#endif
}

#if defined(GAP_RELABEL)
void NetworkMaxFlowSimplex::gap(Dist k) {
	//cout<<"gap detected at d="<<k<<endl;
#ifdef USE_STATS_COUNT
	c_gap++;
#endif
	assert(dc[k] == 0);
	for (Dist d = k+1; d <= listRelabel.dmax; d++) {
		NodeID v = listRelabel.first[d];
		while (v != nSentinel) {
			assert(v != UNDEF_NODE);
			Node& nv = nodes[v];
			NodeID next = nv.listRelabelNext;
			listPivots.update(v, nSentinel);
			dc[nv.d]--; dc[nSentinel]++;
			nv.d = nSentinel;
			v = next;
		}
		listRelabel.first[d] = nSentinel;
	}
	listRelabel.dmax = k-1;
}
#endif


// do a global update
// do a BFS to compute exact d labels (BFS = dijkstra in undirected graphs)
void NetworkMaxFlowSimplex::doGlobalRelabel() {
#ifdef USE_STATS_COUNT
	c_globalupdate++;
	//cout << "Global relabeling: " << c_globalupdate << endl;
#endif

	listRelabel.flush();

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
		Node& nu = nodes[u];
		Arc& auv = arcs[uv]; ArcID vu = auv.rev;
		NodeID v = auv.head; Node& nv = nodes[v];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
			if (headColor == COLOR_WHITE && nv.d > nu.d + 1) {
				listPivots.update(v, nu.d + 1);
				dc[nv.d]--; dc[nu.d + 1]++;
				nv.d = nu.d + 1;
				nv.cur = vu;
			}
			else if (nv.d == nu.d + 1) {
				if (nv.cur > vu) {
					nv.cur = vu;
				}
			}
			return u != sink; // do not search further than sink
		}
		else
			return false;
	},
	color);


	forAllOutArcs(sink, uv, is) {
		Node& nu = nodes[sink];
		Arc& auv = arcs[uv]; ArcID vu = auv.rev;
		NodeID v = auv.head; Node& nv = nodes[v];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
			if (nv.d > nu.d + 1) {
				listPivots.update(v, nu.d + 1);
				dc[nv.d]--; dc[nu.d + 1]++;
				nv.d = nu.d + 1;
				nv.cur = vu;
			}
			else if (nv.d == nu.d + 1) {
				if (nv.cur > vu) {
					nv.cur = vu;
				}
			}
		}
	}

	// check current arcs are correct after bfs
	forAllNodes(v) assert(checkValidCurArc(v));
}


#endif

