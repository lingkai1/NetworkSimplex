#ifndef DFSBFS_CPP_
#define DFSBFS_CPP_
// only template funcs in this file

#include "nmfs.hpp"

#define COLOR_WHITE 0
#define COLOR_GREY  1
#define COLOR_BLACK 2



// Iterative implementation of BFS with pre and post processing of nodes and arcs user-defined
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
inline void NetworkMaxFlowSimplex::bfs(NodeID s,
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
	bool leaf;
	queuePush(s);
	colorSet(s, COLOR_GREY);
	while (!queueIsEmpty()) {
		NodeID u = queuePop();
		preProcessNode(u);
		leaf = true; // regard u as leaf by default
		forAllOutArcs(u, i, is) {
			NodeID v = arcs[i].head;
			if (processArc(u, i, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				leaf = false; // u is not a leaf
				queuePush(v);
				colorSet(v, COLOR_GREY);
			}
		}
		postProcessNode(u, leaf);
		colorSet(u, COLOR_BLACK);
	}
}
template <
typename PreProcessNode,
typename PostProcessNode,
typename ProcessArc>
inline void NetworkMaxFlowSimplex::bfs(NodeID s,
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







// Iterative implementation of DFS with pre and post processing of nodes and arcs user-defined
// based on 1 stack of arcs
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc,
typename PostStackIsEmpty,
typename PostStackPush,
typename PostStackPop,
typename ColorSet,
typename ColorGet
>
inline void NetworkMaxFlowSimplex::dfs_i(NodeID s,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		const PostStackIsEmpty& postStackIsEmpty,
		const PostStackPush& postStackPush,
		const PostStackPop& postStackPop,
		const ColorSet& colorSet,
		const ColorGet& colorGet
) {
	// poststack: stack of visited arcs
	// white: unvisited, untouched node
	// grey: ancestor visited node
	// black: visited node in another branch or descendant.
	fill(dfsiv.begin(), dfsiv.end(), 0); // save arc index position to avoid going back to the beginning every time

	NodeID pu, u; ArcID i;
	bool goingDown = true;
	u = s;
	do {
		if (colorGet(u) == COLOR_WHITE)
			preOrderNode(u);
		colorSet(u, COLOR_BLACK); // no more children to explore
		for (NodeID i = nodes[u].first + dfsiv[u], is = nodes[u+1].first; i != is; i++, dfsiv[u]++) {
			NodeID v = arcs[i].head;
			if (preOrderArc(u, i, colorGet(v)) && colorGet(v) == COLOR_WHITE) {
				colorSet(u, COLOR_GREY); // u has more children to explore
				postStackPush(u, i); // save where we come from
				u = v;
				goingDown = true;
				break;
			}
		}

		if (colorGet(u) == COLOR_BLACK) { // if u has no more children to explore
			if (!postStackIsEmpty()) {
				postOrderNode(u, goingDown);
				postStackPop(pu, i);
				postOrderArc(pu, i);
				u = pu;
				goingDown = false;
			}
			else break;
		}
	} while (true);
	postOrderNode(s, colorGet(s) == COLOR_BLACK);
}
template <
typename PreOrderNode,
typename PostOrderNode,
typename PreOrderArc,
typename PostOrderArc
>
inline void NetworkMaxFlowSimplex::dfs_i(NodeID s,
		const PreOrderNode& preOrderNode,
		const PostOrderNode& postOrderNode,
		const PreOrderArc& preOrderArc,
		const PostOrderArc& postOrderArc,
		std::vector<int>& color
) {
	std::stack<ArcID> post;
	return dfs_i(s,
			preOrderNode,
			postOrderNode,
			preOrderArc,
			postOrderArc,
			[&](){ return post.empty(); },
			[&](NodeID pu, ArcID i){ post.push(i); },
			[&](NodeID& pu, ArcID& i){ i = post.top(); pu = arcs[arcs[i].rev].head; post.pop(); return arcs[i].head; },
			[&](NodeID u, int c){ color[u] = c; },
			[&](NodeID u){ return color[u]; }
	);
}


#endif
