#ifndef DFSBFS_CPP_
#define DFSBFS_CPP_
// only template funcs in this file

#include "nwsimplex.hpp"

#define COLOR_WHITE 0
#define COLOR_GREY  1
#define COLOR_BLACK 2
#define COLOR_BLUE 3


// Iterative implementation of BFS with pre and post processing of nodes and arcs defined
// based on i node queue
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
		colorSet(u, COLOR_BLACK); // regard u as leaf by default (u's color unused for other purposes in the next lines, is used to hold this information now)
		forAllOutArcs(u, i, is) {
			NodeID v = arcs[i].head;
			if (processArc(u, i, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				colorSet(u, COLOR_GREY); // u is not a leaf
				queuePush(v);
				colorSet(v, COLOR_GREY);
			}
		}
		postProcessNode(u, colorGet(u)); // color to determine u is a leaf or not (GREY=not leaf, BLACK=leaf)
		colorSet(u, COLOR_BLACK); // now u's black color only mean it is processed (does not mean it is necessarily a leaf)
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
	colorSet(u, COLOR_BLACK);
	forAllOutArcs(u, i, is) {
		NodeID v = arcs[i].head;
		if (preOrderArc(u, i, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
			colorSet(u, COLOR_GREY);
			dfs_r(v, preOrderNode, postOrderNode, preOrderArc, postOrderArc, colorSet, colorGet);
			postOrderArc(u, i);
		}
	}
	postOrderNode(u, colorGet(u));
	colorSet(u, COLOR_BLACK);
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
	NodeID pu, u; ArcID i;
	u = s;
	do {
		preOrderNode(u);
		colorSet(u, COLOR_BLACK); // regard u as i leaf by default
		forAllOutArcs(u, i, is) {
			NodeID v = arcs[i].head;
			if (preOrderArc(u, i, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				colorSet(u, COLOR_GREY); // u is not i leaf
				preStackPush(u, i);
				colorSet(v, COLOR_BLUE); // marked as in the pre stack
			}
		}

		if (colorGet(u) == COLOR_BLACK) { // if u is i leaf go backwards until next branch
			if (!postStackIsEmpty()) {
				pu = u;

				NodeID pv;
				if (!preStackIsEmpty()) {
					preStackPop(pv, i); preStackPush(pv, i); // pv := next arc tail
				}
				else pv = s; // dfs done, go back to root

				while (pu != pv) {
					assert(!postStackIsEmpty());
					assert(pu != s);
					postOrderNode(pu, colorGet(pu));
					colorSet(pu, COLOR_BLACK);
					u = postStackPop(pu, i);
					postOrderArc(pu, i);
				}
			}
			else assert(u == s);
		}


		if (!preStackIsEmpty()) {
			u = preStackPop(pu, i); // where we go next
			postStackPush(pu, i); // save where we come from
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
	// implementation with double ended queue (two stacks)
	std::deque<ArcID> q;
	int pre = 0, post = 0;
	return dfs_i(s,
			preOrderNode,
			postOrderNode,
			preOrderArc,
			postOrderArc,
			[&](){ return pre == 0; },
			[&](NodeID pu, ArcID i){ pre++; q.push_back(i); },
			[&](NodeID& pu, ArcID& i){ pre--; i = q.back(); pu = arcs[arcs[i].rev].head; q.pop_back(); return arcs[i].head; },
			[&](){ return post == 0; },
			[&](NodeID pu, ArcID i){ post++; q.push_front(i); },
			[&](NodeID& pu, ArcID& i){ post--; i = q.front(); pu = arcs[arcs[i].rev].head; q.pop_front(); return arcs[i].head; },
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
			[&](NodeID pu, ArcID i){ pre.push(std::pair<NodeID, ArcID>(pu, i)); },
			[&](NodeID& pu, ArcID& i){ pu = pre.top().first; i = pre.top().second; pre.pop(); return arcs[i].head; },
			[&](){ return post.empty(); },
			[&](NodeID pu, ArcID i){ post.push(std::pair<NodeID, ArcID>(pu, i)); },
			[&](NodeID& pu, ArcID& i){ pu = post.top().first; i = post.top().second; post.pop(); return arcs[i].head; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
	 */
}



#endif
