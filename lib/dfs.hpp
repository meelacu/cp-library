// @name Depth-First Search (DFS)
// @prefix dfs
// @description Depth-First Search (Recursive)

// DFS
 vb seen(G.size());
auto dfs = [&](this auto&& self, int v) -> void {
    seen[v] = true;
    for (int nxt : G[v]) {
        if (seen[nxt]) continue;
        self(nxt);
    }
};
dfs(0);$0
