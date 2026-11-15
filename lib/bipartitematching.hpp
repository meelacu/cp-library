// @name BipartiteMatching
// @prefix bipartitematching
// @description 二部マッチング（Hopcroft-Karp, 非再帰）と König の定理: add_edge/solve/ml/mr/min_vertex_cover/max_independent_set

// BipartiteMatching（二部グラフの最大マッチング）: Hopcroft-Karp O(E sqrt(V))、非再帰、これ単体で動く
//   BipartiteMatching g(nl, nr);
//   add_edge(l, r)           : 左 l と右 r を結ぶ辺
//   solve()                  : 最大マッチングの大きさ。辺を足してから呼び直すと前回のマッチングの続きから計算する
//   ml[l] / mr[r]            : マッチ相手（いなければ -1）
//   min_vertex_cover()       : 最小点被覆 {左の頂点, 右の頂点}（König の定理。大きさ = 最大マッチング）
//   max_independent_set()    : 最大独立集合 {左の頂点, 右の頂点}（最小点被覆の補集合。大きさ = nl + nr - 最大マッチング）
//   最小辺被覆の大きさ = nl + nr - 最大マッチング（孤立点がないとき）
//   DAG の最小パス被覆 = 頂点数 - 最大マッチング（各頂点を「出る側」「入る側」に分けた二部グラフで）
//   重み付きのマッチングは MinCostFlow を使う
struct BipartiteMatching {
    BipartiteMatching(int nl, int nr) : ml(nl, -1), mr(nr, -1), L(nl), R(nr) {}
    void reserve(int m) { eu.reserve(m); ev.reserve(m); }
    void add_edge(int l, int r) {
        assert(0 <= l && l < L && 0 <= r && r < R);
        eu.push_back(l); ev.push_back(r);
        built = false;
    }

    // 最大マッチングの大きさ。辺を足してから呼び直すと、前回のマッチングの続きから計算する
    int solve() {
        build();
        for (; seen < int(eu.size()); seen++)   // 新しい辺で貪欲に初期マッチングを取る
            if (ml[eu[seen]] < 0 && mr[ev[seen]] < 0) { ml[eu[seen]] = ev[seen]; mr[ev[seen]] = eu[seen]; size++; }
        while (bfs()) {
            for (int v = 0; v < L; v++) it[v] = st[v];
            for (int v = 0; v < L; v++) if (ml[v] < 0 && dfs(v)) size++;
        }
        return size;
    }
    std::vector<int> ml, mr;   // ml[l] = 左 l のマッチ相手（いなければ -1）。mr[r] も同様

    // König の定理: 最小点被覆 = (Z に入らない左頂点) + (Z に入る右頂点)。Z = 未マッチの左頂点から交互路で届く頂点
    std::pair<std::vector<int>, std::vector<int>> min_vertex_cover() {
        auto [zl, zr] = reachable();
        std::vector<int> cl, cr;
        for (int v = 0; v < L; v++) if (!zl[v]) cl.push_back(v);
        for (int u = 0; u < R; u++) if (zr[u]) cr.push_back(u);
        return {cl, cr};
    }
    std::pair<std::vector<int>, std::vector<int>> max_independent_set() {
        auto [zl, zr] = reachable();
        std::vector<int> il, ir;
        for (int v = 0; v < L; v++) if (zl[v]) il.push_back(v);
        for (int u = 0; u < R; u++) if (!zr[u]) ir.push_back(u);
        return {il, ir};
    }

  private:
    int L, R, seen = 0, size = 0;
    std::vector<int> eu, ev, st, adj, lev, it, stk, que;
    bool built = false;

    void build() {
        if (built) return;
        built = true;
        st.assign(L + 1, 0);
        for (int l : eu) st[l + 1]++;
        for (int v = 0; v < L; v++) st[v + 1] += st[v];
        std::vector<int> cur(st.begin(), st.end() - 1);
        adj.resize(eu.size());
        for (int i = 0; i < int(eu.size()); i++) adj[cur[eu[i]]++] = ev[i];
        lev.assign(L, -1); it.assign(L, 0); que.resize(L);
    }

    bool bfs() {
        int qt = 0;
        for (int v = 0; v < L; v++) { lev[v] = ml[v] < 0 ? 0 : -1; if (lev[v] == 0) que[qt++] = v; }
        bool found = false;
        for (int qh = 0; qh < qt; qh++) {
            int v = que[qh];
            for (int i = st[v]; i < st[v + 1]; i++) {
                int w = mr[adj[i]];
                if (w < 0) found = true;
                else if (lev[w] < 0) { lev[w] = lev[v] + 1; que[qt++] = w; }
            }
        }
        return found;
    }

    // 非再帰: 交互路が長くてもスタックオーバーフローしない
    bool dfs(int root) {
        stk.assign(1, root);
        while (!stk.empty()) {
            int v = stk.back();
            bool advanced = false;
            for (int& i = it[v]; i < st[v + 1]; i++) {
                int u = adj[i], w = mr[u];
                if (w < 0) {                       // スタックを逆にたどって増加させる
                    for (int k = int(stk.size()) - 1; k >= 0; k--) {
                        int x = stk[k], y = (k + 1 == int(stk.size())) ? u : adj[it[x]];
                        ml[x] = y; mr[y] = x; lev[x] = -1;
                    }
                    return true;
                }
                if (lev[w] == lev[v] + 1) { stk.push_back(w); advanced = true; break; }
            }
            if (!advanced) { lev[v] = -1; stk.pop_back(); }
        }
        return false;
    }

    std::pair<std::vector<char>, std::vector<char>> reachable() {
        assert(built && seen == int(eu.size()) && "call solve() first");
        std::vector<char> zl(L, 0), zr(R, 0);
        std::vector<int> q;
        for (int v = 0; v < L; v++) if (ml[v] < 0) { zl[v] = 1; q.push_back(v); }
        for (int qi = 0; qi < int(q.size()); qi++) {
            int v = q[qi];
            for (int i = st[v]; i < st[v + 1]; i++) {
                int u = adj[i];
                if (zr[u]) continue;
                zr[u] = 1;
                int w = mr[u];
                if (w >= 0 && !zl[w]) { zl[w] = 1; q.push_back(w); }
            }
        }
        return {zl, zr};
    }
};
$0
