#include "nwsimplex.hpp"

using namespace std;


/*
 //print subtree
NodeID u0 = 3;
cout<<"subtree("<<u0<<")"<<endl;
NodeID u = u0;
for (int i = 0; i < nodes[u0].stSize; i++) {
	u = nodes[u].next;
	cout<<u<<endl;
}
 */


// main network simplex function
void NWS::solve() {

	buildInitialBasis();



	if (verbose >= 1) cout << "Selecting Pivot" << endl;
	NodeID v;
	ArcID avw;
	if (!extractMinPivotArc(v, avw)) {
		// no pivot available.
		// Simplex method ended
		if (verbose >= 1) cout << "Simplex method ended. No more pivot arcs" << endl;

		return;
	}
	Arc& vw = arcs[avw];
	NodeID w = vw.head;

	cout << "Pivot: " << v << "->" << w << endl;



}

void NWS::buildInitialBasis() {
	// construct a basis S (BFS) and T

	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.d = u != sink ? 0 : 1; // every node starts with d = 0 except sink (d = 1)
		nu.cur = nu.first;
	}


	forAllNodes(u) {
		Node& nu = nodes[u];
		nu.bToRelabelNext = nu.bToRelablePrev = UNDEF_NODE;
	}

	forAllArcs(u, a, aStop) {
		Arc& aa = arcs[a];
		aa.bPivotArcNext = aa.bPivotArcPrev = UNDEF_ARC;
	}

	for (NodeID d = 0; d < n; d++) {
		bToRelabelFirst[d] = UNDEF_NODE;
		bPivotArcFirst[d] = UNDEF_ARC;
	}


	//do BFS to build tree


	//do DFS to initialize subtree fields of each vertex
	NodeID prevu = source;
	vector<int> color(n, COLOR_WHITE); // todo: to optimize (reuse another field maybe)
	if (verbose >= 1) cout << "Initial DFS" << endl;
	dfs_i(source,
			[&](NodeID u){
		// pre node
		if (verbose>=2) cout<<"pre node: "<<u<<endl;
		// link nodes in DFS order as a singly linked list
		Node& prevnu = nodes[prevu];
		prevnu.next = u;
		prevu = u;
	},
	[&](NodeID u, int color){
		// post node
		if (color == COLOR_BLACK) { // if leaf node
			Node& nu = nodes[u];
			nu.stSize = 0;
		}
		if (verbose>=2) cout<<"post node: "<<u<<" subtree size: "<<nodes[u].stSize<<endl;
	},
	[&](NodeID u, ArcID a, int headColor) {
		// arc choice (pre arc)
		Arc& aa = arcs[a];
		if (aa.resCap > 0 && headColor == COLOR_WHITE) { // forward residual tree arc
			// process forward DFS arc
			if (aa.head != sink) {
				if (verbose >= 2) cout << "tree arc: " << u << "->"  << aa.head << " rescap: " << aa.resCap << endl;
				aa.l = 1; //0;
				// set parent
				Node& z = nodes[aa.head];
				z.parent = u;
				return true; // accept this choice as a tree arc
			}
			else {
				// this arc is an initial pivot arc. add it to the priority queue
				if (verbose >= 2) cout << "pivot arc: " << u << "->"  << aa.head << " rescap: " << aa.resCap << endl;
				insertPivotArc(u, a);
				aa.l = INF_DIST;
				return false; // reject so that we do not search further than sink
			}
		}
		else if (headColor == COLOR_GREY) { // backward tree arc
			aa.l = 1; //0;
			return false;
		}
		else if (aa.resCap > 0) {
			// residual non-tree arc
			aa.l = 1;
			return false;
		}
		else {
			aa.l = INF_DIST;
			return false; // reject other types of arcs
		}
	},
	[&](NodeID u, ArcID a) {
		// post arc
		Arc& aa = arcs[a];
		//if (verbose >= 2) cout << "post arc: " << u << "->"  << arc.head << endl;
		Node& nu = nodes[u];
		Node& nh = nodes[aa.head];
		nu.stSize += nh.stSize + 1;
	},
	color);


	forAllNodes(v) {
		if (v != sink)
			makeCur(v);
	}

}

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

bool NWS::makeCur(NodeID v) {
	Node& nv = nodes[v];
	if (nv.first >= nodes[v+1].first) // no arc incident arc to v
		return true;
	while (nv.cur < nodes[v+1].first) {
		Arc& aa = arcs[nv.cur];
		Node& nu = nodes[aa.head];
		if (nv.d == nu.d + arcs[aa.rev].l) {
			return true;
		}
		nv.cur++;
	}
	// add v to to_relabel list
	toRelabel.push(v);
	return false;
}

void NWS::relabel(NodeID v) {
	Node& nv = nodes[v];
	assert(nv.first < nodes[v+1].first); // there is at least one incident arc to v
	Dist oldD = nv.d;
	nv.d = INF_DIST;
	forAllOutArcs(v, a, aStop) {
		Arc& aa = arcs[a];
		Node& nu = nodes[aa.head];
		Dist newDist = nu.d + arcs[aa.rev].l;
		if (newDist < nv.d) {
			nv.d = newDist;
			nv.cur = a;
		}
	}
	assert(nv.d > oldD);
	//After relabel(v) increases d(v), we check all arcs (v,u). If u was current before the relabeling and cur(u) = (v,u), we make u non-current.
	//During this processing, if relabel(v) makes u non-current, we apply make_cur(u) and if it fails, we add u to the relabel list.
	forAllOutArcs(v, a, aStop) {
		Arc& aa = arcs[a];
		NodeID u = aa.head;
		Node& nu = nodes[u];
		if (nu.cur == aa.rev) {
			makeCur(u);
		}
	}
}



void NWS::insertPivotArc(NodeID u, ArcID a) {
	if (bPivotArcFirst[0] == UNDEF_ARC) {
		bPivotArcFirst[0] = a;
	}
	else {
		ArcID b = bPivotArcFirst[0];
		Arc& na = arcs[a];
		Arc& nb = arcs[b];
		bPivotArcFirst[0] = a;
		na.bPivotArcPrev = UNDEF_ARC;
		na.bPivotArcNext = b;
		nb.bPivotArcPrev = a;
	}
}
bool NWS::extractMinPivotArc(NodeID& u, ArcID& a) {
	if (bPivotArcFirst[0] == UNDEF_ARC) {
		return false;
	}
	else {
		a = bPivotArcFirst[0];
		Arc& na = arcs[a];
		u = arcs[na.rev].head;
		NodeID b = na.bPivotArcNext;
		bPivotArcFirst[0] = b;
		if (b != UNDEF_ARC) {
			Arc& nb = arcs[b];
			nb.bPivotArcPrev = UNDEF_ARC;
		}
		return true;
	}
}
