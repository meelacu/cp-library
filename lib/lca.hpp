// @name LCA (doubling)
// @prefix lca
// @description LCA(ダブリング) テンプレ: lca/la/dist/dist_w/jump/is_ancestor/on_path/meet, 森対応, 非再帰

struct LCA {
    int n, LOG;
    vector<vector<pair<int, long long>>> g;
    vector<vector<int>> up;   // up[k][v] : v の 2^k 個上の祖先 (無ければ -1)
    vector<int> dep;          // 根からの辺数
    vector<long long> dw;     // 根からの重み付き距離

    LCA(int n_ = 0) { init(n_); }
    void init(int n_) {
        n = n_;
        LOG = 1;
        while ((1 << LOG) < n) ++LOG;
        g.assign(n, {});
    }
    void add_edge(int u, int v, long long w = 1) {
        g[u].push_back({v, w});
        g[v].push_back({u, w});
    }
    // 根を指定して前処理 (非再帰BFSなのでスタックオーバーフローしない)
    void build(int root = 0) {
        up.assign(LOG + 1, vector<int>(n, -1));
        dep.assign(n, -1);
        dw.assign(n, 0);
        vector<int> q;
        q.reserve(n);
        auto bfs = [&](int s) {
            dep[s] = 0;
            dw[s] = 0;
            int head = q.size();
            q.push_back(s);
            while (head < (int)q.size()) {
                int v = q[head++];
                for (auto [to, w] : g[v])
                    if (dep[to] < 0) {
                        dep[to] = dep[v] + 1;
                        dw[to] = dw[v] + w;
                        up[0][to] = v;
                        q.push_back(to);
                    }
            }
        };
        if (0 <= root && root < n) bfs(root);
        for (int v = 0; v < n; ++v)
            if (dep[v] < 0) bfs(v);   // 森でも動く (連結成分ごとに根を取る)
        for (int k = 0; k < LOG; ++k)
            for (int v = 0; v < n; ++v)
                up[k + 1][v] = (up[k][v] < 0 ? -1 : up[k][up[k][v]]);
    }
    // v の k 個上の祖先 (存在しなければ -1)
    int la(int v, int k) const {
        if (k < 0 || k > dep[v]) return -1;
        for (int i = 0; k; k >>= 1, ++i)
            if (k & 1) v = up[i][v];
        return v;
    }
    // 最小共通祖先 (別成分なら -1)
    int lca(int u, int v) const {
        if (dep[u] > dep[v]) swap(u, v);
        v = la(v, dep[v] - dep[u]);
        if (u == v) return u;
        for (int k = LOG; k >= 0; --k)
            if (up[k][u] != up[k][v]) { u = up[k][u]; v = up[k][v]; }
        return up[0][u];
    }
    int dist(int u, int v) const {          // 辺数
        int l = lca(u, v);
        return l < 0 ? -1 : dep[u] + dep[v] - 2 * dep[l];
    }
    long long dist_w(int u, int v) const {  // 重み付き距離
        int l = lca(u, v);
        return l < 0 ? -1 : dw[u] + dw[v] - 2 * dw[l];
    }
    // u -> v のパス上で u から k 本進んだ頂点 (範囲外なら -1)
    int jump(int u, int v, int k) const {
        int l = lca(u, v);
        if (l < 0) return -1;
        int du = dep[u] - dep[l], dv = dep[v] - dep[l];
        if (k < 0 || k > du + dv) return -1;
        return k <= du ? la(u, k) : la(v, du + dv - k);
    }
    bool is_ancestor(int a, int v) const {  // a は v の祖先か (a == v も true)
        return dep[a] <= dep[v] && la(v, dep[v] - dep[a]) == a;
    }
    bool on_path(int u, int v, int x) const {  // x が u-v パス上にあるか
        return dist(u, x) + dist(x, v) == dist(u, v);
    }
    int meet(int a, int b, int c) const {  // 3頂点のパスが集まる合流点
        int x = lca(a, b), y = lca(b, c), z = lca(c, a);
        if (dep[x] < dep[y]) x = y;
        if (dep[x] < dep[z]) x = z;
        return x;
    }
};
