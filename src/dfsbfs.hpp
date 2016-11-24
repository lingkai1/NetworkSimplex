
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
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc, typename PreStackIsEmpty, typename PreStackPush, typename PreStackPop, typename PostStackIsEmpty, typename PostStackPush, typename PostStackPop, typename ColorSet, typename ColorGet>
inline void dfs_i2(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, const PreStackIsEmpty& preStackIsEmpty, const PreStackPush& preStackPush, const PreStackPop& preStackPop, const PostStackIsEmpty& postStackIsEmpty, const PostStackPush& postStackPush, const PostStackPop& postStackPop, const ColorSet& colorSet, const ColorGet& colorGet);
template <typename PreOrderNode, typename PostOrderNode, typename PreOrderArc, typename PostOrderArc>
inline void dfs_i2(NodeID s, const PreOrderNode& preOrderNode, const PostOrderNode& postOrderNode, const PreOrderArc& preOrderArc, const PostOrderArc& postOrderArc, std::vector<int>& color);

