
// BFS
template <typename PreProcessNode, typename PostProcessNode, typename ProcessArc, typename QueueIsEmpty, typename QueuePush, typename QueuePop, typename ColorSet, typename ColorGet>
inline void bfs(NodeID s, const PreProcessNode& preProcessNode, const PostProcessNode& postProcessNode, const ProcessArc& processArc, const QueueIsEmpty& queueIsEmpty, const QueuePush& queuePush, const QueuePop& queuePop, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreProcessNode, typename PostProcessNode, typename ProcessArc>
inline void bfs(NodeID s, const PreProcessNode& preProcessNode, const PostProcessNode& postProcessNode, const ProcessArc& processArc, std::vector<int>& color);

// DFS
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc, typename ColorSet, typename ColorGet>
inline void dfs_r(NodeID u, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc>
inline void dfs_r(NodeID u, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, std::vector<int>& color);
template <typename PreOrderNode, typename PreOrderArc, typename PreNodeStackIsEmpty, typename PreNodeStackPush, typename PreNodeStackPop, typename ColorSet, typename ColorGet>
inline void dfs_ip(NodeID s, const PreOrderNode& preOrderNode, const PreOrderArc& preOrderArc, const PreNodeStackIsEmpty& preStackIsEmpty, const PreNodeStackPush& preStackPush, const PreNodeStackPop& preStackPop, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreOrderNode, typename PreOrderArc>
inline void dfs_ip(NodeID s, const PreOrderNode& preOrderNode, const PreOrderArc& preOrderArc, std::vector<int>& color);
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc, typename PreStackIsEmpty, typename PreStackPush, typename PreStackPop, typename PostStackIsEmpty, typename PostStackPush, typename PostStackPop, typename ColorSet, typename ColorGet>
inline void dfs_i(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, const PreStackIsEmpty& preStackIsEmpty, const PreStackPush& preStackPush, const PreStackPop& preStackPop, const PostStackIsEmpty& postStackIsEmpty, const PostStackPush& postStackPush, const PostStackPop& postStackPop, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc>
inline void dfs_i(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, std::vector<int>& color);

/*
void testBFS() {
	std::vector<int> color(n, COLOR_WHITE);

	std::cout << "BFS" << std::endl;
	bfs(source,
			[&](NodeID u){
		std::cout << "pre:  " << u << std::endl;
	},
	[&](NodeID u){
		std::cout << "post: " << u << std::endl;
	},
	[&](NodeID u, ArcID a, int headColor){
		switch (headColor) {
		case COLOR_WHITE:
			std::cout << "forward arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		case COLOR_GREY:
			std::cout << "cross arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		case COLOR_BLACK:
			std::cout << "backward arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		default:
			break;
		}
		return true;
	}, color);

}
*/
/*
void testDFS() {
	std::vector<int> color(n, COLOR_WHITE);

	std::cout << "DFS" << std::endl;
	dfs_i(source,
			[&](NodeID u){
		std::cout << "pre:  " << u << std::endl;
	},
	[&](NodeID u){ std::cout << "post: " << u << std::endl; },
	[&](NodeID u, ArcID a, int headColor){
		switch (headColor) {
		case COLOR_WHITE:
			std::cout << "forward pre arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		case COLOR_BLACK:
			std::cout << "cross pre arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		case COLOR_GREY:
			std::cout << "backward pre arc: " << u << "->"  << arcs[a].head << std::endl;
			break;
		default:
			break;
		}
		return true;
	},
	[&](NodeID u, ArcID a){
		std::cout << "post arc: " << u << "->"  << arcs[a].head << std::endl;
	}, color);
}
*/
