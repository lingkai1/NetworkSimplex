
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
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc, typename PostStackIsEmpty, typename PostStackPush, typename PostStackPop, typename ColorSet, typename ColorGet>
inline void dfs_i(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, const PostStackIsEmpty& postStackIsEmpty, const PostStackPush& postStackPush, const PostStackPop& postStackPop, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc>
inline void dfs_i(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, std::vector<int>& color);


/*
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
*/
