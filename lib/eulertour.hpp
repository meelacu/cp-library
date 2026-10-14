// @name EulerTour
// @prefix eulertour
// @description 木（森）の行きがけ順 tin/tout/order（部分木=[tin,tout)・heavy-first・非再帰・森対応・edge_child・is_anc）

// EulerTour（行きがけ順・heavy-first・非再帰・森対応）
//   部分木 v = [tin[v], tout[v])、order[t] = 時刻 t の頂点、is_anc(u, v): u が v の祖先（u = v を含む）
//   森は root の後、未訪問で番号最小の頂点を根にする。edge_child[i] = 辺 i の子側の頂点
// 使い方メモ
//   初期配列: init[t] = a[order[t]]
//   S(v) = 根→v の和（点更新）: x の更新を [tin[x], tout[x]) への区間加算にし、tin[v] を点取得（差分 BIT はサイズ n+1）
//   u-v パス: S(u) + S(v) - 2S(w) + a[w]（w = LCA）。辺重みは edge_child に載せ、+a[w] を外す
//   木 DP: order を逆順に回すと子→親の順（par[v] に集約）
//   2N 時刻（木上 Mo など）: 入 = 2*tin[v] - dep[v]、出 = 2*tout[v] - dep[v] - 1
//   LCA などの上位層は par / dep を受け取って作る
struct EulerTour {
    int n;
    vector<int> tin, tout, order, par, dep, edge_child;
    vector<pair<int, int>> edges;

    explicit EulerTour(int n_) : n(n_) {}

    int add_edge(int u, int v) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        edges.emplace_back(u, v);
        return int(edges.size()) - 1;
    }

    void build(int root = 0) {
        assert(0 <= root && root < n);
        int m = int(edges.size());
        vector<int> st(n + 1, 0);  // CSR: v の隣接は adj[st[v], st[v+1])、(相手, 辺 id)
        vector<pair<int, int>> adj(2 * m);
        for (auto [u, v] : edges) st[u]++, st[v]++;
        for (int v = 0; v < n; v++) st[v + 1] += st[v];
        for (int i = m - 1; i >= 0; i--) {
            auto [u, v] = edges[i];
            adj[--st[u]] = {v, i}, adj[--st[v]] = {u, i};
        }
        par.assign(n, -1), dep.assign(n, -1), edge_child.assign(m, -1);  // dep = -1 は未訪問
        vector<int> bfs;
        bfs.reserve(n);
        for (int k = -1; k < n; k++) {
            int r = k < 0 ? root : k;
            if (dep[r] >= 0) continue;
            dep[r] = 0, bfs.push_back(r);
            for (int i = int(bfs.size()) - 1; i < int(bfs.size()); i++) {
                int v = bfs[i];
                for (int j = st[v]; j < st[v + 1]; j++) {
                    auto [c, e] = adj[j];
                    if (dep[c] >= 0) continue;
                    par[c] = v, dep[c] = dep[v] + 1, edge_child[e] = c;
                    bfs.push_back(c);
                }
            }
        }
        vector<int> sz(n, 1), heavy(n, -1);
        int roots = 0;
        for (int i = n - 1; i >= 0; i--) {
            int v = bfs[i], p = par[v];
            if (p < 0) { roots++; continue; }
            sz[p] += sz[v];
            if (heavy[p] < 0 || sz[heavy[p]] < sz[v]) heavy[p] = v;
        }
        assert(m == n - roots);  // 閉路・多重辺・自己ループがない（森である）
        tin.assign(n, 0), tout.assign(n, 0), order.assign(n, 0);
        vector<int> nxt(n);  // nxt[p] = p の次の軽い子に渡す時刻
        int offset = 0;
        for (int v : bfs) {  // 親の区間を heavy 子 → 軽い子（隣接順）の順に切り分ける
            int p = par[v];
            if (p < 0) tin[v] = offset, offset += sz[v];
            else if (v == heavy[p]) tin[v] = tin[p] + 1;
            else tin[v] = nxt[p], nxt[p] += sz[v];
            nxt[v] = tin[v] + 1 + (heavy[v] < 0 ? 0 : sz[heavy[v]]);
            tout[v] = tin[v] + sz[v], order[tin[v]] = v;
            // HLD に拡張するときはここで head[v] = (p >= 0 && v == heavy[p]) ? head[p] : v;
        }
    }

    pair<int, int> subtree(int v) const { return {tin[v], tout[v]}; }
    int size(int v) const { return tout[v] - tin[v]; }
    bool is_anc(int u, int v) const { return tin[u] <= tin[v] && tin[v] < tout[u]; }
};
