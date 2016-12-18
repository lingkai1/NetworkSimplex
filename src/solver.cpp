#include "nmfs.hpp"
#include <stdio.h>

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

	fill(color.begin(), color.end(), COLOR_WHITE);
	if (verbose >= 2) cout << "Initial BFS" << endl;
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
	[&](){ return bfsq.empty(); },
	[&](NodeID u){ bfsq.push(u); },
	[&](){ NodeID u = bfsq.front(); bfsq.pop(); return u; },
	[&](NodeID u, int c){ color[u] = c; },
	[&](NodeID u){ return color[u]; }
	);

	forAllNodes(v) if (v!= sink && color[v] == COLOR_WHITE) {
		// delete unreachable nodes beyond sink
		Node& nv = nodes[v];
		//cout<<"Unreachable node: "<<v<<endl;
		forAllOutArcs(v, vu, is) {
			Arc& avu = arcs[vu];
			ArcID uv = avu.rev; Arc& auv = arcs[uv];
			NodeID u = avu.head; Node& nu = nodes[u];
			auv.resCap = avu.resCap = 0;
		}
		nv.cur = nv.first = nodes[v+1].first; // delete arcs
		nv.d = n;
	}

	//do DFS to initialize subtree fields of each vertex
	NodeID prevu = source;
	nodes[source].prev = UNDEF_NODE;
	fill(color.begin(), color.end(), COLOR_WHITE);
	if (verbose >= 2) cout << "Initial DFS" << endl;
	dfs_i(source,
			[&](NodeID u){
		// pre node
		// link nodes in DFS order as i doubly linked list
		Node& prevnu = nodes[prevu];
		prevnu.next = u;
		Node& nu = nodes[u];
		nu.prev = prevu;
		nu.size = 1; // initialize subtree size
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
		nu.size += nh.size; // compute subtree sizes in post visit order
	},
	[&](){ return dfss.empty(); },
	[&](NodeID pu, ArcID i){ dfss.push(i); },
	[&](NodeID& pu, ArcID& i){ i = dfss.top(); pu = arcs[arcs[i].rev].head; dfss.pop(); return arcs[i].head; },
	[&](NodeID u, int c){ color[u] = c; },
	[&](NodeID u){ return color[u]; }
	);

	Node& nt = nodes[sink];
	nt.parent = UNDEF_ARC;
	nt.next = UNDEF_NODE;
	nt.size = 1;
	nt.tree = IN_T;

	nodes[source].parent = UNDEF_ARC;

	d[0] = n - 1;
	d[1] = 1;


#if defined(GLOBAL_RELABEL)
	globalRelabelFreq = 0.2;
	globalRelabelThreshold = n;
	globalRelabelWork = 0;
#endif

	globalUpdate();
	//qt.insert(sink, nodes[sink].d); // set in globalUpdate()
}



// main network simplex function
void NetworkMaxFlowSimplex::solve() {
#ifdef USE_STATS_TIME
	t_total = timer();
#endif
#ifdef USE_STATS_COUNT
	c_basisIter = 0;
	c_basisGrowS = 0;
	c_basisGrowT = 0;
	c_basisNoChange = 0;
	c_StoTMoves = 0;
	c_TtoSMoves = 0;
	c_augPathTotalLen = 0;
	c_makeCur = 0;
	c_relabel = 0;
	c_relabelArcScans = 0;
	c_globalUpdate = 0;
	c_guArcScans = 0;
	c_gap = 0;
#endif

	flow = 0;

	if (verbose >= 1) cout << "Building initial basis..." << endl;

#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer();
#endif
	buildInitialBasis();
#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer() - t_buildInitialBasis;
#endif

#ifdef USE_STATS_TIME
	t_solve = timer();
#endif

	if (verbose >= 1) cout << "Simplex method..." << endl;

	// main loop
	while (true) {
		// print some debug info every sec
		static double t_loopInfo = timer(); if (timer() - t_loopInfo > 1.0) {
			//cout << "makeCur(): " << c_makeCur << endl;
			//cout << "relabel(): " << c_relabel << endl;
			cout << "Flow: " << flow << endl;
			fflush(stdout); t_loopInfo = timer();
		}

		//if (c_gap == 0) forAllNodes(u) if (nodes[u].tree == IN_T) if (!qt.contains(u) && !qr.contains(u)) { cout<<"problem with u: " << u<<endl; assert(false); }

		if (verbose >= 3) {
			printS();
			printT();
		}

		if (verbose >= 2) cout << "Selecting Pivot..." << endl;
		NodeID v, w;
		ArcID vw, wv;
		vw = extractMinPivot(); // extract a pivot with min label
		if (vw == UNDEF_ARC) {
			// no pivot available.
			// Simplex method ended
			if (verbose >= 1) cout << "No pivots left. Simplex method ended." << endl;

#ifdef USE_STATS_TIME
			t_solve = timer() - t_solve;
			t_total= timer() - t_total;
#endif
			return;
		}

		Arc& avw = arcs[vw];
		wv = avw.rev;
		Arc& awv = arcs[wv];
		v = awv.head; Node& nv = nodes[v];
		w = avw.head; Node& nw = nodes[w];

		assert(nv.d < n);
		//assert(nw.d == nv.d + 1);

		if (verbose >= 2) cout << "Pivot (entering) vw arc: " << v << "->" << w << endl;

		assert(nv.tree == IN_S);
		assert(nw.tree == IN_T);

		if (verbose >= 2) cout << "Searching bottleneck..." << endl;
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
#ifdef USE_STATS_COUNT
			c_augPathTotalLen++;
#endif
		}
#ifdef USE_STATS_COUNT
		c_augPathTotalLen++;
#endif
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
#ifdef USE_STATS_COUNT
			c_augPathTotalLen++;
#endif
		}


		Arc& axy = arcs[xy];
		yx = axy.rev;
		Arc& ayx = arcs[yx];
		x = ayx.head; Node& nx = nodes[x];
		y = axy.head; Node& ny = nodes[y];
		if (verbose >= 2) cout << "Bottleneck (leaving) xy arc: " << x << "->" << y << " delta: " << delta << endl;

		if (verbose >= 2) cout << "Augmenting flow..." << endl;
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



		if (verbose >= 2) cout << "Updating basis..." << endl;
		if (xyTree == IN_S) {
			if (verbose >= 2) cout << "xy in S => T grows" << endl;

			// Q the subtree rooted at y in S moves from S to T
			// R will represent Q rooted at v in T
			NodeID nLast;
			getSubtreeLastNode(y, nLast);
			deleteSubtree(source, y, nLast); // disconnect Q from S
			changeRoot(y, v, nLast); // this transforms Q into R (updates parents and the doubly linked list)
			nv.parent = vw;
			addSubtreeAsChild(sink, v, nLast, w); // connect R to T

			// move subtree from S to T
			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_T;
				if (u != y) {
#if defined(LAZY_RELABEL)
					if (isCurrent(u))
#else
						//assert(isCurrent(u));
#endif
						qt.insert(u, nu.d);
				}
			});

			// reinsert w because it was removed when pivot was extracted
			// w is already current
			//assert(isCurrent(w));
			qt.insert(w, nw.d);

			//assert(ny.tree == IN_T);
			if (makeCur(y)) {
				qt.insert(y, ny.d);
			} else {
				qt.remove(y);
#if defined(LAZY_RELABEL)
				if (!qr.contains(y)) // y can possibly be in Qr whether it was in S or T
#else
					assert(!qr.contains(y)); // y was in S so it cannot be in Qr
#endif
				qr.insert(y, ny.d);
			}
			processQr();

#ifdef USE_STATS_COUNT
			c_StoTMoves += nv.size;
			c_basisGrowT++;
#endif
		}
		else if (xyTree == IN_T) {
			if (verbose >= 2) cout << "xy in T => S grows" << endl;

			// Q the subtree rooted at x in T moves from T to S
			// R will represent Q rooted at w in S
			NodeID nLast;
			getSubtreeLastNode(x, nLast);
			deleteSubtree(sink, x, nLast); // disconnect Q from T
			changeRoot(x, w, nLast); // this transforms Q into R (updates parents and the doubly linked list)
			nw.parent = wv;
			addSubtreeAsChild(source, w, nLast, v); // connect R to S

			forAllSubTree(w, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_S;
				qt.remove(u);
			});

			qt.remove(y);

			//assert(ny.tree == IN_T);
			if (makeCur(y)) {
				qt.insert(y, ny.d);
			} else {
				qt.remove(y);
#if defined(LAZY_RELABEL) || defined(LAZY_RELABEL2)
				if (!qr.contains(y))
#else
					//assert(!qr.contains(y)); Qr is empty between pivots
#endif
					qr.insert(y, nodes[y].d);
			}
			processQr();

#ifdef USE_STATS_COUNT
			c_TtoSMoves += nw.size;
			c_basisGrowS++;
#endif
		}
		else {
			if (verbose >= 2) cout << "xy == vw => Unchanged basis" << endl;

			//qt.insert(w, nw.d); // not needed because y == w

			//assert(ny.tree == IN_T);
			if (makeCur(y)) {
				qt.insert(y, ny.d);
			} else {
				qt.remove(y);
#if defined(LAZY_RELABEL) || defined(LAZY_RELABEL2)
				if (!qr.contains(y))
#else
					//assert(!qr.contains(y)); Qr is empty between pivots
#endif
					qr.insert(y, ny.d);
			}
			processQr();

#ifdef USE_STATS_COUNT
			c_basisNoChange++;
#endif
		}

		//forAllNodes(v) if (v != source) assert(checkCurrent(v));

#ifdef USE_STATS_COUNT
		c_basisIter++;
#endif
	}

}

void NetworkMaxFlowSimplex::processQr() {
	while (!qr.processed())
		relabel(qr.extractMin());

#if defined(LAZY_RELABEL2)
	for (Dist k=qr.dmin; k<=qr.dmax; k++)
		for (NodeID v = qr.first[k]; v != nSentinel; ) {
			NodeID w = nodes[v].qrNext;
			if (nodes[v].tree == IN_S) {
				qr.remove(v);
				relabel(v);
			}
			v = w;
		}
#endif
}


// returns a v from the pivots list with minimal d and deletes it from the list
ArcID NetworkMaxFlowSimplex::extractMinPivot() {
	NodeID w = qt.extractMin();
	if (w == UNDEF_NODE) // no pivots found
		return UNDEF_ARC;
	Node& nw = nodes[w];
	assert(nw.tree == IN_T);
	ArcID wv = nw.cur; Arc& awv = arcs[wv];
	ArcID vw = awv.rev; Arc& avw = arcs[vw];
	NodeID v = awv.head; Node& nv = nodes[v];
	//assert(isCurrent(w));
	assert(nw.d == nv.d + 1);
	if (nv.tree != IN_S) {
		cerr<<"w:"<<w<<"->v:"<<v<<" not in S"<<endl;
	}
	assert(nv.tree == IN_S);
	return vw;
}


// go to the last node of a subtree in DFS order
void NetworkMaxFlowSimplex::getSubtreeLastNode(NodeID u, NodeID& uLast) {
	uLast = u;
	for (NodeID k = 0; k < nodes[u].size - 1; k++) {
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
	NodeID uSize = nodes[u].size;
	NodeID z = u;
	while (z != q) {
		z = arcs[nodes[z].parent].head;
		nodes[z].size -= uSize;
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
	int pSize = nodes[p].size;
	NodeID z = u;
	while (z != r) {
		nodes[z].size += pSize;
		z = arcs[nodes[z].parent].head;
	}
	nodes[r].size += pSize;

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

void NetworkMaxFlowSimplex::setLabel(NodeID v, Dist k) {
	Node& nv = nodes[v];
	if (qt.contains(v)) {
		qt.remove(v);
		if (k < n)
			qt.insert(v, k);
	}
	d[nv.d]--; d[k]++;
	nv.d = k;
}

bool NetworkMaxFlowSimplex::isCurrent(NodeID v) {
	return !qr.contains(v);
	/*Node& nv = nodes[v];
	forAllOutArcs(v, vu, is) {
		Node& nv = nodes[v];
		Arc& avu = arcs[vu]; NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
			assert(nv.d <= nu.d + 1);
			if (nv.d == nu.d + 1) {
				return nv.cur == vu;
			}
		}
	}
	return false;*/
}


// try to make v current, find and set its new current arc
// if v is already current, does nothing
// if v cannot be made current it is added to the toRelabel list
bool NetworkMaxFlowSimplex::makeCur(NodeID v) {
#ifdef USE_STATS_COUNT
	c_makeCur++;
#endif
	if (v == source) return true;
	assert(v != source);
	Node& nv = nodes[v];
	assert(nv.d < n); //if (nv.d >= n) return true;
	if (nv.first >= nodes[v+1].first) return false;
	//checkValidCurArc(v);
	assert(nv.first <= nv.cur && nv.cur <= nodes[v+1].first);
	while (nv.cur < nodes[v+1].first) {
		ArcID vu = nv.cur; Arc& avu = arcs[vu];
		Node& nu = nodes[avu.head];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv) && nv.d == nu.d + 1) {
			if (verbose >= 3) cout << v << " (d=" << nv.d << ") made current to " << avu.head << endl;
			return true;
		}
		nv.cur++;
	}
	return false;
}
// in the general case after makeCur we want to do:
/*
			if (makeCur(u)) {
				if (nu.tree == IN_T) qt.insert(u, nu.d);
			} else {
				if (nu.tree == IN_T) qt.remove(u);
				if (!qr.contains(u)) qr.insert(u, nu.d);
			}
 */

// relabel v
// this must result in a strict increase of v's d label based on the assumption that v's current arc was of minimal index
void NetworkMaxFlowSimplex::relabel(NodeID v) {
#if defined(GLOBAL_RELABEL)
	// try global update if a lot of work has been done
	globalRelabelWork++;
	if (globalRelabelWork * globalRelabelFreq > globalRelabelThreshold) {
		globalUpdate();
		globalRelabelWork = 0;
		return;
	}
#endif
#ifdef USE_STATS_COUNT
	c_relabel++;
#endif
	Node& nv = nodes[v];
	assert(nv.d >= 0);
	assert(nv.d < n);
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist oldD = nv.d;
	Dist newD = INF_DIST;
	nv.cur = nv.first;
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		if (isResidualOrTreeArc(nu, nv, uv, vu, auv)) {
			if (nu.d + 1 < newD) {
				newD = nu.d + 1;
				nv.cur = vu;
			}
		}
#ifdef USE_STATS_COUNT
		c_relabelArcScans++;
#endif
	}

	if (newD <= nv.d) cerr << "For v:"<<v<<" newD:"<<newD<<" <= d(v):"<<nv.d<<endl;
	assert(newD > nv.d);
	if (newD > n)
		newD = n;

	setLabel(v, newD);
	if (verbose >= 3) cout<<v<<" (d="<<nodes[v].d<<") relabeled now current to "<<arcs[nv.cur].head<<endl;

	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev;
		if (nu.cur == uv && nu.d < n) {
			if (verbose >= 3) cout << u << " (d=" << nu.d << ") made non-current by relabeling of " << v << endl;
			//if (nu.tree == IN_T) assert(qt.contains(u)); //qt.remove(u);
			if (makeCur(u)) {
				//if (nu.tree == IN_T) qt.insert(u, nu.d);
			} else {
				if (nu.tree == IN_T) qt.remove(u);
				//assert(!qr.contains(u));
				qr.insert(u, nu.d);
			}
		}
	}

	if (newD == n) {
		if (verbose >= 3) cout<<v<<" (d="<<nodes[v].d<<") relabeled now ignored"<<endl;
		return;
	}
	else {
		if (nv.tree == IN_T) {
			qt.insert(v, nv.d);
		}
	}

#if defined(GAP_RELABEL)
	if (d[oldD] == 0) {
		gap(oldD);
		return;
	}
#endif
}

#if defined(GAP_RELABEL)

void NetworkMaxFlowSimplex::gap(Dist k) {
#ifdef USE_STATS_COUNT
	c_gap++;
#endif
	if (verbose >= 1) cout<<"Gap detected at d="<<k<<endl;
	qt.flush();
	qr.flush();
	return;
	/*
	assert(dc[k] == 0);
	for (Dist d = k+1; d <= qr.dmax; d++) {
		NodeID v = qr.first[d];
		while (v != nSentinel) {
			assert(v != UNDEF_NODE);
			Node& nv = nodes[v];
			NodeID next = nv.qrNext;
			setLabel(v, n);
			v = next;
		}
		qr.first[d] = nSentinel;
	}
	qr.dmax = k-1;
	 */
}
#endif

// do a global update
// do a BFS to compute exact d labels (BFS = dijkstra in unweighted graphs)
void NetworkMaxFlowSimplex::globalUpdate() {
#ifdef USE_STATS_COUNT
	c_globalUpdate++;
	if (verbose >= 1) cout << "Global relabeling: " << c_globalUpdate << endl;
#endif

	qr.flush();
	qt.flush(); // we need to rebuild qt

	// do BFS from source to compute exact labels
	nodes[source].d = 0;

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
			if (headColor == COLOR_WHITE) {
				setLabel(v, nu.d + 1);
				nv.cur = vu;
				if (nv.tree == IN_T)
					qt.insert(v, nv.d);
			}
			else if (nv.d == nu.d + 1) {
				if (nv.cur > vu) {
					nv.cur = vu;
				}
			}
#ifdef USE_STATS_COUNT
			c_guArcScans++;
#endif
			return true; // BFS in pseudo-residual graph
		}
		else
			return false;
	},
	[&](){ return bfsq.empty(); },
	[&](NodeID u){ bfsq.push(u); },
	[&](){ NodeID u = bfsq.front(); bfsq.pop(); return u; },
	[&](NodeID u, int c){ color[u] = c; },
	[&](NodeID u){ return color[u]; }
	);

	forAllNodes(v) if (color[v] == COLOR_WHITE) {
		// unreachable nodes.
		setLabel(v, n);
	}

	// check current arcs are correct after bfs
	// forAllNodes(v) if (v != source) assert(checkCurrent(v));
}




