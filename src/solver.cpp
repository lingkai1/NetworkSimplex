#include "nwsimplex.hpp"

using namespace std;


void NetworkMaxFlowSimplex::buildInitialBasis() {
	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.d = u != sink ? 0 : 1; // every node starts with d = 0 except sink (d = 1)
		nu.cur = nu.first;
		nu.tree = IN_NONE;
		nu.bToRelabelNext = nu.bToRelablePrev = UNDEF_NODE;
	}
	forAllArcs(u, i, is) {
		Arc& ai = arcs[i];
	}
	for (NodeID d = 0; d < n; d++) {
		bToRelabel[d].first = UNDEF_NODE;
		bToRelabel[d].last = UNDEF_NODE;
		bPivots[d].first = UNDEF_ARC;
		bPivots[d].last = UNDEF_ARC;
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
				// this arc is an initial pivot arc. add it to the priority queue
				if (verbose >= 3) cout << "pivot arc: " << u << "->"  << ai.head << " rescap: " << ai.resCap << endl;
				pivotInsert(i);
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

	forAllNodes(v) {
		if (v != sink)
			makeCur(v);
	}

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
#endif

	flow = 0;

#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer();
#endif
	buildInitialBasis();
#ifdef USE_STATS_TIME
	t_buildInitialBasis = timer() - t_buildInitialBasis;
#endif

	double t_refresh=timer();
	while (true) {
		if (timer() - t_refresh > 1.0) { // print flow every sec
			cout << "Flow: " << flow << endl;
			fflush(stdout);
			t_refresh = timer();
		}


		if (verbose >= 3) {
			cout<<"TREE S: "; printSubTree(source); cout<<endl;
			cout<<"TREE T: "; printSubTree(sink); cout<<endl;
			printS();
			printT();
		}

		if (verbose >= 1) cout << "Selecting Pivot..." << endl;
		NodeID v, w;
		ArcID vw, wv;
		if (!pivotExtractMin(vw)) {
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


			// update pivots list
			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_T;
			});
			// examine neighbors, this is probably costly
			forAllSubTree(v, [&](NodeID u){
				Node& nu = nodes[u];
				forAllOutArcs(u, uv, is) {
					Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					if (nv.tree == IN_T) {
						if (auv.resCap > 0) {
							bool ok = pivotDelete(uv); //assert(ok);
						}
					}
					else if (nv.tree == IN_S) {
						if (avu.resCap > 0) {
							pivotInsert(vu);
						}
					}
				}
			});
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

			// update pivots list
			forAllSubTree(w, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_S;
			});
			// examine neighbors
			forAllSubTree(w, [&](NodeID u){
				Node& nu = nodes[u];
				forAllOutArcs(u, uv, is) {
					Arc& auv = arcs[uv]; NodeID v = auv.head; Node& nv = nodes[v];
					ArcID vu = auv.rev; Arc& avu = arcs[vu];
					if (nv.tree == IN_S) {
						if (avu.resCap > 0) {
							bool ok = pivotDelete(vu); //assert(ok);
						}
					}
					else if (nv.tree == IN_T) {
						if (auv.resCap > 0) {
							pivotInsert(uv);
						}
					}
				}
			});
#ifdef USE_STATS_COUNT
			c_basisGrowS++;
#endif
		}
		else {
			if (verbose >= 1) cout << "xy == vw => Unchanged basis" << endl;
			// pivot is the leaving arc. keep the same basis
#ifdef USE_STATS_COUNT
			c_basisNoChange++;
#endif
		}

#ifdef USE_STATS_COUNT
		c_basisIter++;
#endif
	}

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

	if (verbose >= 3) {cout<<"SUBTREE deleted: "; printSubTree(u); cout << " from: "; printSubTree(q); cout<<endl;}
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

	if (verbose >= 3) {cout<<"SUBTREE added to: "; printSubTree(r); cout<<endl;}
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

void NetworkMaxFlowSimplex::printSubTree(NodeID root) {
	forAllSubTree(root, [&](NodeID u){
		Node& nu = nodes[u];
		if (u == root)
			cout << root;
		else
			cout << "," << u;
	});
}

void NetworkMaxFlowSimplex::printS() {
	cout << "S tree rooted at " << source << endl;
	forAllSubTree(source, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << u << "->" << p << " rc: " << ai.resCap << endl;
		}
	});
}
void NetworkMaxFlowSimplex::printT() {
	cout << "T tree rooted at " << sink << endl;
	forAllSubTree(sink, [&](NodeID u){
		ArcID r = nodes[u].parent;
		if (r != UNDEF_ARC) {
			Arc& ar = arcs[r]; ArcID i = arcs[r].rev; Arc& ai = arcs[i];
			NodeID p = ar.head;
			cout << u << "->" << p << " rc: " << ai.resCap << endl;
		}
	});
}


// todo: change implementation to use buckets
#ifdef PIVOTS_QUEUE
void NetworkMaxFlowSimplex::pivotInsert(ArcID i) {
	pivots.push_back(i);
#ifdef USE_STATS_COUNT
	c_pivotsInserted++;
#endif
}
bool NetworkMaxFlowSimplex::pivotDelete(ArcID i) {
	auto it = find(pivots.begin(), pivots.end(), i);
	if (it == pivots.end())
		return false;
	else {
		pivots.erase(it);
#ifdef USE_STATS_COUNT
		c_pivotsDeleted++;
#endif
		return true;
	}
}
bool NetworkMaxFlowSimplex::pivotExtractMin(ArcID& i) {
	if (pivots.empty())
		return false;
	else {
		//i = pivots.front();
		//pivots.pop_front();
		auto it = pivots.begin();
		advance(it, rand()%(int)pivots.size());
		i = *it;
		pivots.erase(it);
		return true;
	}
}
#else
#endif






// not yet used
void NetworkMaxFlowSimplex::prepareNextPivot() {
	// to call after an augmentation to make vertices current, relabel to prepare for next pivot

	forAllNodes(v) {
		if (v != sink)
			makeCur(v);
	}

	while (!toRelabel.empty()) {
		NodeID v = toRelabel.top();
		toRelabel.pop();
		relabel(v);
		//cout<<v<<" relabeled at "<<nodes[v].d<<endl;
	}


}


// todo: implement current arcs and relabeling in algorithm
// not yet used
bool NetworkMaxFlowSimplex::makeCur(NodeID v) {
	Node& nv = nodes[v];
	if (nv.first >= nodes[v+1].first) // no arc incident arc to v
		return true;
	while (nv.cur < nodes[v+1].first) {
		ArcID vu = nv.cur; Arc& avu = arcs[vu];
		Node& nu = nodes[avu.head];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		Dist l = (auv.resCap > 0 || nu.parent == uv || nv.parent == vu) ? 1 : INF_DIST;
		if (nv.d == nu.d + l) {
			return true;
		}
		nv.cur++;
	}
	// add v to to_relabel list
	toRelabel.push(v);
	return false;
}
// todo: implement relabeling
// not yet used
void NetworkMaxFlowSimplex::relabel(NodeID v) {
	Node& nv = nodes[v];
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist oldD = nv.d;
	nv.d = INF_DIST;
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev; Arc& auv = arcs[uv];
		Dist l = (auv.resCap > 0 || nu.parent == uv || nv.parent == vu) ? 1 : INF_DIST;
		Dist newDist = nu.d + l;
		if (newDist < nv.d) {
			nv.d = newDist;
			nv.cur = vu;
		}
	}
	assert(nv.d > oldD);
	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, vu, is) {
		Arc& avu = arcs[vu];
		NodeID u = avu.head; Node& nu = nodes[u];
		ArcID uv = avu.rev;
		if (nu.cur == uv) {
			makeCur(u);
		}
	}
}





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


