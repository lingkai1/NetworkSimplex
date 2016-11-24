#include "nwsimplex.hpp"

using namespace std;



// main network simplex function
void NWS::solve() {

	flow = 0;

	buildInitialBasis();


	while (true) {

		if (verbose >= 1)  cout << "Flow: " << flow << endl;

		if (verbose >= 1) {
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


			return;
		}

		Arc& avw = arcs[vw];
		wv = avw.rev;
		Arc& awv = arcs[wv];
		v = awv.head;
		w = avw.head;

		if (verbose >= 1) cout << "Pivot (entering) vw arc: " << v << "->" << w << endl;


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
			Arc& ar = arcs[nu.parent]; // (u,pu)
			ArcID i = ar.rev; Arc& ai = arcs[i]; // (pu,u)
			NodeID pu = ar.head;
			if (ai.resCap <= delta) { // <= so we save the closest to the source
				delta = ai.resCap;
				xy = i;
				xyTree = IN_S;
			}
			u = pu;
		}
		// path w->t
		u = w;
		while (u != sink) {
			Node& nu = nodes[u];
			ArcID i = nu.parent; Arc& ai = arcs[nu.parent]; // (u,pu)
			NodeID pu = ai.head;
			if (ai.resCap < delta) { // need strict <
				delta = ai.resCap;
				xy = i;
				xyTree = IN_T;
			}
			u = pu;
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
			Arc& ar = arcs[nu.parent]; // (u,pu)
			ArcID i = ar.rev;
			Arc& ai = arcs[i]; // (pu,u)
			NodeID pu = ar.head;

			ai.resCap -= delta;
			ar.resCap += delta;

			u = pu;
		}
		// path w->t
		u = w;
		while (u != sink) {
			Node& nu = nodes[u];
			ArcID i = nu.parent; Arc& ai = arcs[nu.parent]; // (u,pu)
			NodeID pu = ai.head;
			Arc& ar = arcs[ai.rev];

			ai.resCap -= delta;
			ar.resCap += delta;

			u = pu;
		}
		flow += delta;
		assert(axy.resCap == 0);


		if (verbose >= 1) cout << "Updating basis..." << endl;
		if (xyTree == IN_S) {
			if (verbose >= 1) cout << "xy in S => T grows" << endl;
			// reverse parent pointers on the y->v path and make v->w
			u = v;
			ArcID pi = vw;
			while (u != y) {
				Node& nu = nodes[u];
				ArcID r = nu.parent;
				nu.parent = pi;
				Arc& ar = arcs[r];
				pi = ar.rev;
				u = ar.head;
			}
			nodes[u].parent = pi;
			// update pivots list
			// subtree(y) moves from S to T
			forAllSubTree(y, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_T;
			});
			// examine neighbors
			forAllSubTree(y, [&](NodeID u){
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
		}
		else if (xyTree == IN_T) {
			if (verbose >= 1) cout << "xy in T => S grows" << endl;
			// reverse parent pointers on the w->x path and make w->v
			u = w;
			ArcID pi = wv;
			while (u != x) {
				Node& nu = nodes[u];
				ArcID i = nu.parent;
				nu.parent = pi;
				Arc& ai = arcs[i];
				pi = ai.rev;
				u = ai.head;
			}
			nodes[u].parent = pi;
			// update pivots list
			// subtree(x) moves from T to S
			forAllSubTree(x, [&](NodeID u){
				Node& nu = nodes[u];
				nu.tree = IN_S;
			});
			// examine neighbors
			forAllSubTree(x, [&](NodeID u){
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
		}
		else {
			if (verbose >= 1) cout << "xy == vw => Unchanged basis" << endl;
			// pivot is the leaving arc. keep the same basis
		}

	}

}

void NWS::printSubTree(NodeID root) {
	forAllSubTree(root, [&](NodeID u){
		Node& nu = nodes[u];
		if (u != root) {
			cout << arcs[nu.parent].head << "->" << u << " rc: " << arcs[arcs[nu.parent].rev].resCap << endl;
		}
	});
}

void NWS::printS() {
	cout << "S tree rooted at " << source << endl;
	vector<int> color(n);
	fill(color.begin(), color.end(), COLOR_WHITE);
	dfs_i(source,
			[&](NodeID u){},
			[&](NodeID u, bool leaf){},
			[&](NodeID u, ArcID i, int headColor) {
				// arc choice (pre arc)
				Arc& ai = arcs[i];
				ArcID r = ai.rev;
				NodeID v = ai.head; Node& nv = nodes[v];
				if (headColor == COLOR_WHITE && v != sink && nv.parent == r) {
					cout << v << "->" << u << " rc: " << ai.resCap << endl;
					return true;
				}
				else
					return false;
			},
			[&](NodeID u, ArcID i){},
			color);
}
void NWS::printT() {
	cout << "T tree rooted at " << sink << endl;
	vector<int> color(n);
	fill(color.begin(), color.end(), COLOR_WHITE);
	dfs_i(sink,
			[&](NodeID v){},
			[&](NodeID v, bool leaf){},
			[&](NodeID v, ArcID r, int headColor) {
				// arc choice (pre arc)
				Arc& ar = arcs[r];
				ArcID i = ar.rev; Arc& ai = arcs[i];
				NodeID u = ar.head; Node& nu = nodes[u];
				if (headColor == COLOR_WHITE && u != source && nu.parent == i) {
					cout << u << "->" << v << " rc: " << ai.resCap << endl;
					return true;
				}
				else
					return false;
			},
			[&](NodeID v, ArcID i){},
			color);
}

void NWS::buildInitialBasis() {
	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.d = u != sink ? 0 : 1; // every node starts with d = 0 except sink (d = 1)
		nu.cur = nu.first;
		nu.tree = IN_NONE;
		nu.bToRelabelNext = nu.bToRelablePrev = UNDEF_NODE;
	}
	forAllArcs(u, i, is) {
		Arc& ai = arcs[i];
		ai.bPivotNext = ai.bPivotPrev = UNDEF_ARC;
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
				if (verbose >= 2) cout << "tree arc: " << u << "->"  << ai.head << " rescap: " << ai.resCap << endl;
				ai.l = 1; //0;
				// set parent
				Node& z = nodes[ai.head];
				z.parent = ai.rev;
				return true; // accept this choice as i tree arc
			}
			else {
				// this arc is an initial pivot arc. add it to the priority queue
				if (verbose >= 2) cout << "pivot arc: " << u << "->"  << ai.head << " rescap: " << ai.resCap << endl;
				pivotInsert(i);
				ai.l = INF_DIST;
				return false; // reject so that we do not search further than sink
			}
		}
		else if (headColor == COLOR_BLACK) { // backward tree arc
			ai.l = 1; //0;
			return false;
		}
		else if (ai.resCap > 0) {
			// residual non-tree arc
			ai.l = 1;
			return false;
		}
		else {
			ai.l = INF_DIST;
			return false; // reject other types of arcs
		}
	},
	color);

	//do DFS to initialize subtree fields of each vertex
	NodeID prevu = source;
	fill(color.begin(), color.end(), COLOR_WHITE);
	if (verbose >= 1) cout << "Initial DFS" << endl;
	dfs_i(source,
			[&](NodeID u){
		// pre node
		// link nodes in DFS order as i singly linked list
		Node& prevnu = nodes[prevu];
		prevnu.next = u;
		prevu = u;
		Node& nu = nodes[u];
		nu.stSize = 1; // initialize subtree size
	},
	[&](NodeID u, bool leaf){
		// post node
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



// not yet used
void NWS::prepareNextPivot() {
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
bool NWS::makeCur(NodeID v) {
	Node& nv = nodes[v];
	if (nv.first >= nodes[v+1].first) // no arc incident arc to v
		return true;
	while (nv.cur < nodes[v+1].first) {
		Arc& ai = arcs[nv.cur];
		Node& nu = nodes[ai.head];
		if (nv.d == nu.d + arcs[ai.rev].l) {
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
void NWS::relabel(NodeID v) {
	Node& nv = nodes[v];
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist oldD = nv.d;
	nv.d = INF_DIST;
	forAllOutArcs(v, i, is) {
		Arc& ai = arcs[i];
		NodeID u = ai.head; Node& nu = nodes[u];
		ArcID r = ai.rev; Arc& ar = arcs[r];
		Dist newDist = nu.d + ar.l;
		if (newDist < nv.d) {
			nv.d = newDist;
			nv.cur = i;
		}
	}
	assert(nv.d > oldD);
	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, i, is) {
		Arc& ai = arcs[i];
		NodeID u = ai.head; Node& nu = nodes[u];
		if (nu.cur == ai.rev) {
			makeCur(u);
		}
	}
}


// now only bucket[d=0] is used
void NWS::pivotInsert(ArcID i) {
	Arc& ai = arcs[i];
	ArcID j = bPivots[0].last;
	if (j != UNDEF_ARC)
		arcs[j].bPivotNext = i;
	else
		bPivots[0].first = i;
	ai.bPivotPrev = j;
	ai.bPivotNext = UNDEF_ARC;
	bPivots[0].last = i;
	if (verbose >= 2) cout <<"Pivot inserted: ("<<arcs[arcs[i].rev].head<<","<<arcs[i].head<<")"<<endl;
}
// O(n) deletion since only level d=0 is used now
bool NWS::pivotDelete(ArcID i) {
	ArcID j = bPivots[0].first;
	while (j != UNDEF_ARC) {
		Arc& aj = arcs[j];
		if (j == i) {
			ArcID x = aj.bPivotPrev, y = aj.bPivotNext;
			if (x != UNDEF_ARC) {
				Arc& ax = arcs[x];
				ax.bPivotNext = y;
			}
			else
				bPivots[0].first = y;
			if (y != UNDEF_ARC) {
				Arc& ay = arcs[y];
				ay.bPivotPrev = x;
			}
			else
				bPivots[0].last = x;
			aj.bPivotNext = aj.bPivotPrev = UNDEF_ARC;
			if (verbose >= 2) cout <<"Pivot removed: ("<<arcs[arcs[i].rev].head<<","<<arcs[i].head<<")"<<endl;
			return true;
		}
		j = aj.bPivotNext;
	}
	return false;
}
bool NWS::pivotExtractMin(ArcID& i) {
	i = bPivots[0].first;
	if (i != UNDEF_ARC) {
		Arc& ai = arcs[i];
		ArcID j = ai.bPivotNext;
		if (j != UNDEF_NODE)
			arcs[j].bPivotPrev = UNDEF_ARC;
		else
			bPivots[0].last = UNDEF_NODE;
		bPivots[0].first = j;
		ai.bPivotNext = ai.bPivotPrev = UNDEF_ARC;
		return true;
	}
	else
		return false;
}





void NWS::testBFS() {
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

void NWS::testDFS() {
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


