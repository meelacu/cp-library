// @name SCC
// @prefix scc
// @description 強連結成分分解（非再帰Tarjan, CSR）: add_edge/build/comp/num/groups/dag

// SCC（強連結成分分解）: 非再帰 Tarjan + CSR
//   add_edge(u, v) をすべて呼んだあと build()
//   comp[v] / scc[v] : v の属する成分番号（トポロジカル順: 辺 u->v で comp[u]!=comp[v] なら comp[u] < comp[v]）
//   num              : 成分数
//   groups()         : 成分ごとの頂点リスト（成分番号順）
//   dag(dedup)       : 縮約 DAG の隣接リスト（dedup=true で多重辺除去）
struct SCC {
    int n, num = 0;
    vector<pair<int, int>> edges;
    vector<int> start, elist, comp;
    SCC(int n) : n(n), comp(n, -1) {}
    void add_edge(int u, int v) { edges.emplace_back(u, v); }
    int operator[](int v) const { return comp[v]; }
    int build() {
        // CSR 構築
        start.assign(n + 1, 0);
        elist.resize(edges.size());
        for (auto [u, v] : edges) start[u + 1]++;
        for (int i = 0; i < n; i++) start[i + 1] += start[i];
        vector<int> it(start.begin(), start.end() - 1);
        for (auto [u, v] : edges) elist[it[u]++] = v;
        for (int i = 0; i < n; i++) it[i] = start[i];
        // 非再帰 Tarjan
        vector<int> ord(n, -1), low(n), st, call;
        st.reserve(n), call.reserve(n);
        int now = 0;
        num = 0;
        for (int s = 0; s < n; s++) {
            if (ord[s] != -1) continue;
            ord[s] = low[s] = now++;
            st.push_back(s), call.push_back(s);
            while (!call.empty()) {
                int v = call.back();
                if (it[v] < start[v + 1]) {
                    int u = elist[it[v]++];
                    if (ord[u] == -1) {
                        ord[u] = low[u] = now++;
                        st.push_back(u), call.push_back(u);
                    } else if (comp[u] == -1) {
                        low[v] = min(low[v], ord[u]);
                    }
                } else {
                    call.pop_back();
                    if (!call.empty()) low[call.back()] = min(low[call.back()], low[v]);
                    if (low[v] == ord[v]) {
                        while (true) {
                            int u = st.back();
                            st.pop_back();
                            comp[u] = num;
                            if (u == v) break;
                        }
                        num++;
                    }
                }
            }
        }
        for (int v = 0; v < n; v++) comp[v] = num - 1 - comp[v];
        return num;
    }
    vector<vector<int>> groups() const {
        vector<vector<int>> res(num);
        for (int v = 0; v < n; v++) res[comp[v]].push_back(v);
        return res;
    }
    vector<vector<int>> dag(bool dedup = true) const {
        vector<vector<int>> res(num);
        for (int v = 0; v < n; v++)
            for (int i = start[v]; i < start[v + 1]; i++)
                if (comp[v] != comp[elist[i]]) res[comp[v]].push_back(comp[elist[i]]);
        if (dedup)
            for (auto& e : res) {
                sort(e.begin(), e.end());
                e.erase(unique(e.begin(), e.end()), e.end());
            }
        return res;
    }
};
$0
