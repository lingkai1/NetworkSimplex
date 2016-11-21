#ifndef DFSBFS_CPP_
#define DFSBFS_CPP_
// only template funcs in this file

#include "nwsimplex.hpp"


// Iterative implementation of BFS with pre and post processing of nodes and arcs defined
// based on a node queue
template <
typename PreProcessNode,
typename PostProcessNode,
typename ProcessArc,
typename QueueIsEmpty,
typename QueuePush,
typename QueuePop,
typename ColorSet,
typename ColorGet
>
inline void NWS::bfs(NodeID s,
		const PreProcessNode& preProcessNode,
		const PostProcessNode& postProcessNode,
		const ProcessArc& processArc,
		const QueueIsEmpty& queueIsEmpty,
		const QueuePush& queuePush,
		const QueuePop& queuePop,
		const ColorSet& colorSet,
		const ColorGet& colorGet
) {
	// queue: queue of next nodes to visit
	// white: unvisited, untouched node
	// grey:  unvisited node in queue
	// black: visited node
	queuePush(s);
	colorSet(s, COLOR_GREY);
	while (!queueIsEmpty()) {
		NodeID u = queuePop();
		preProcessNode(u);
		colorSet(u, COLOR_BLACK);
		forAllOutArcs(u, a, aStop) {
			NodeID v = arcs[a].head;
			if (processArc(u, a, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				queuePush(v);
				colorSet(v, COLOR_GREY);
			}
		}
		postProcessNode(u);
	}
}
template <
typename PreProcessNode,
typename PostProcessNode,
typename ProcessArc>
inline void NWS::bfs(NodeID s,
		const PreProcessNode& preProcessNode,
		const PostProcessNode& postProcessNode,
		const ProcessArc& processArc,
		std::vector<int>& color
) {
	std::queue<NodeID> q;
	return bfs(s,
			preProcessNode,
			postProcessNode,
			processArc,
			[&](){ return q.empty(); },
			[&](NodeID u){ q.push(u); },
			[&](){ NodeID u = q.front(); q.pop(); return u; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
}








// Recursive implementation of DFS with pre and post processing of nodes and arcs defined
// advantage: no stack needed (call stack implicitly used), disadvantage: may need very large call stack
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc,
typename ColorSet,
typename ColorGet
>
inline void NWS::dfs_r(NodeID u,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		const ColorSet& colorSet,
		const ColorGet& colorGet
) {
	// white: unvisited, untouched node
	// grey: ancestor visited node
	// black: visited node in another branch or descendant.
	preOrderNode(u);
	colorSet(u, COLOR_GREY);
	forAllOutArcs(u, a, aStop) {
		NodeID v = arcs[a].head;
		if (preOrderArc(u, a, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
			dfs_r(v, preOrderNode, postOrderNode, preOrderArc, postOrderArc, colorSet, colorGet);
			colorSet(v, COLOR_BLACK);
			postOrderArc(u, a);
		}
	}
	postOrderNode(u);
}
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc
>
inline void NWS::dfs_r(NodeID u,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		std::vector<int>& color
) {
	return dfs_r(u,
			preOrderNode,
			postOrderNode,
			preOrderArc,
			postOrderArc,
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
}

// simple iterative node stack based implementation of DFS with only pre order processing of node and arc defined
template <
typename PreOrderNode,
typename PreOrderArc,
typename PreNodeStackIsEmpty,
typename PreNodeStackPush,
typename PreNodeStackPop,
typename ColorSet,
typename ColorGet
>
inline void NWS::dfs_ip(NodeID s,
		const PreOrderNode& preOrderNode,
		const PreOrderArc& preOrderArc,
		const PreNodeStackIsEmpty& preStackIsEmpty,
		const PreNodeStackPush& preStackPush,
		const PreNodeStackPop& preStackPop,
		const ColorSet& colorSet,
		const ColorGet& colorGet
) {
	// prestack: stack of next nodes to visit
	// white: unvisited, untouched node
	// grey: ancestor visited node
	// black: visited node in another branch or descendant.
	preStackPush(s);
	colorSet(s, COLOR_GREY);
	while (!preStackIsEmpty()) {
		NodeID u = preStackPop();
		preOrderNode(u);
		colorSet(u, COLOR_BLACK);
		forAllOutArcs(u, a, aStop) {
			NodeID v = arcs[a].head;
			if (preOrderArc(u, a, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				colorSet(u, COLOR_GREY);
				preStackPush(v);
				colorSet(v, COLOR_GREY);
			}
		}
	}
}
template <
typename PreOrderNode,
typename PreOrderArc
>
inline void NWS::dfs_ip(NodeID s,
		const PreOrderNode& preOrderNode,
		const PreOrderArc& preOrderArc,
		std::vector<int>& color
) {
	std::stack<NodeID> pre;
	return dfs_ip(s,
			preOrderNode,
			preOrderArc,
			[&](){ return pre.empty(); },
			[&](NodeID u){ pre.push(u); },
			[&](){ NodeID u = pre.top(); pre.pop(); return u; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
}

// Iterative implementation of DFS with pre and post processing of nodes and arcs defined
// based on two stacks of arcs
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc,
typename PreStackIsEmpty,
typename PreStackPush,
typename PreStackPop,
typename PostStackIsEmpty,
typename PostStackPush,
typename PostStackPop,
typename ColorSet,
typename ColorGet
>
inline void NWS::dfs_i(NodeID s,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		const PreStackIsEmpty& preStackIsEmpty,
		const PreStackPush& preStackPush,
		const PreStackPop& preStackPop,
		const PostStackIsEmpty& postStackIsEmpty,
		const PostStackPush& postStackPush,
		const PostStackPop& postStackPop,
		const ColorSet& colorSet,
		const ColorGet& colorGet
) {
	// prestack: stack of next arcs to visit
	// poststack: stack of visited arcs
	// white: unvisited, untouched node
	// blue: unvisited node, saved in the prestack
	// grey: ancestor visited node
	// black: visited node in another branch or descendant.
	NodeID pu, u; ArcID a;
	u = s;
	do {
		preOrderNode(u);
		colorSet(u, COLOR_BLACK); // regard u as a leaf by default
		forAllOutArcs(u, a, aStop) {
			NodeID v = arcs[a].head;
			if (preOrderArc(u, a, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				colorSet(u, COLOR_GREY); // u is not a leaf
				preStackPush(u, a);
				colorSet(v, COLOR_BLUE); // marked as in the pre stack
			}
			//else if (colorGet(v) != COLOR_BLUE) preOrderArc(u, a, colorGet(v));
		}

		if (colorGet(u) == COLOR_BLACK) { // if u is a leaf go backwards until next branch
			if (!postStackIsEmpty()) {
				pu = u;

				NodeID pv;
				if (!preStackIsEmpty()) {
					preStackPop(pv, a); preStackPush(pv, a); // pv := next arc tail
				}
				else pv = s; // dfs done, go back to root

				while (pu != pv) {
					assert(!postStackIsEmpty());
					assert(pu != s);
					postOrderNode(pu, colorGet(pu));
					colorSet(pu, COLOR_BLACK);
					u = postStackPop(pu, a);
					postOrderArc(pu, a);
				}
			}
			else assert(u == s);
		}


		if (!preStackIsEmpty()) {
			u = preStackPop(pu, a); // where we go next
			postStackPush(pu, a); // save where we come from
			//preOrderArc(pu, a, COLOR_WHITE);
		} else {
			postOrderNode(s, colorGet(s));
			break;
		}
	} while (true);
}
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc
>
inline void NWS::dfs_i(NodeID s,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		std::vector<int>& color
) {
	// implementation with double ended stack
	std::deque<ArcID> q;
	int pre = 0, post = 0;
	return dfs_i(s,
			preOrderNode,
			postOrderNode,
			preOrderArc,
			postOrderArc,
			[&](){ return pre == 0; },
			[&](NodeID pu, ArcID a){ pre++; q.push_back(a); },
			[&](NodeID& pu, ArcID& a){ pre--; a = q.back(); pu = arcs[arcs[a].rev].head; q.pop_back(); return arcs[a].head; },
			[&](){ return post == 0; },
			[&](NodeID pu, ArcID a){ post++; q.push_front(a); },
			[&](NodeID& pu, ArcID& a){ post--; a = q.front(); pu = arcs[arcs[a].rev].head; q.pop_front(); return arcs[a].head; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
	/*
	// implementation with two stacks
	std::stack<std::pair<NodeID, ArcID>> pre;
	std::stack<std::pair<NodeID, ArcID>> post;
	return dfs_i(s,
			preOrderNode,
			postOrderNode,
			preOrderArc,
			postOrderArc,
			[&](){ return pre.empty(); },
			[&](NodeID pu, ArcID a){ pre.push(std::pair<NodeID, ArcID>(pu, a)); },
			[&](NodeID& pu, ArcID& a){ pu = pre.top().first; a = pre.top().second; pre.pop(); return arcs[a].head; },
			[&](){ return post.empty(); },
			[&](NodeID pu, ArcID a){ post.push(std::pair<NodeID, ArcID>(pu, a)); },
			[&](NodeID& pu, ArcID& a){ pu = post.top().first; a = post.top().second; post.pop(); return arcs[a].head; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
	 */
}



#endif
