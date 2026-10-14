// @name Rerooting (全方位木DP)
// @prefix rerooting
// @description 全方位木DP O(N)。本体のみ。使い方は rerooting-ex

// 全方位木DP (Rerooting)  O(N)  —— 使い方は rerooting-ex
template <class S, class Merge, class AddEdge, class AddVert>
struct Rerooting {
    int n;
    std::vector<std::vector<std::pair<int, int>>> g;
    S id;
    Merge mg;      // S(S, S)
    AddEdge ae;    // S(S, int eid)
    AddVert av;    // S(S, int v)
    bool vertex_in_answer=false;

    Rerooting(int n_, S id_, Merge mg_, AddEdge ae_, AddVert av_)
        : n(n_), g(n_), id(id_), mg(mg_), ae(ae_), av(av_) {}

    void add(int u, int v, int eid = -1) {
        g[u].emplace_back(v, eid);
        g[v].emplace_back(u, eid);
    }

    std::vector<S> solve(int root = 0) {
        std::vector<int> par(n, -1), ord;
        ord.reserve(n);
        {   // 非再帰DFSで行きがけ順を作る
            std::vector<char> seen(n, 0);
            std::vector<int> st{root};
            seen[root] = 1;
            while (!st.empty()) {
                int v = st.back(); st.pop_back();
                ord.push_back(v);
                for (auto [u, eid] : g[v]) if (!seen[u]) {
                    seen[u] = 1; par[u] = v; st.push_back(u);
                }
            }
        }

        // 1回目: 葉 -> 根
        std::vector<S> dp(n, id);
        for (int i = n - 1; i >= 0; --i) {
            int v = ord[i];
            S acc = id;
            for (auto [u, eid] : g[v])
                if (u != par[v]) acc = mg(acc, ae(dp[u], eid));
            dp[v] = av(acc, v);
        }

        // 2回目: 根 -> 葉
        std::vector<S> up(n, id), ans(n, id);
        for (int v : ord) {
            int d = (int)g[v].size();
            std::vector<S> val(d, id), pre(d + 1, id), suf(d + 1, id);
            for (int i = 0; i < d; ++i) {
                auto [u, eid] = g[v][i];
                val[i] = (u == par[v]) ? up[v] : ae(dp[u], eid);
            }
            for (int i = 0; i < d; ++i)      pre[i + 1] = mg(pre[i], val[i]);
            for (int i = d - 1; i >= 0; --i) suf[i]     = mg(val[i], suf[i + 1]);

            ans[v] = vertex_in_answer ? av(pre[d], v) : pre[d];

            for (int i = 0; i < d; ++i) {
                auto [u, eid] = g[v][i];
                if (u == par[v]) continue;
                up[u] = ae(av(mg(pre[i], suf[i + 1]), v), eid);
            }
        }
        return ans;
    }
};
$0
