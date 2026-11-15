// @name ProjectSelection
// @prefix projectselection
// @description 燃やす埋める・最大重み閉包（要 maxflow スニペット）: add_profit/add_cost/add_penalty/add_implication/force/add_profit_if_all/add_cost_if_any/solve/feasible/assignment

// ProjectSelection（燃やす埋める / 最大重み閉包）（MaxFlow スニペットが必要）
//   各頂点を group 0（s 側）か group 1（t 側）に振り分けて、利得の合計を最大化する
//   ProjectSelection<T> g(n);
//   add_profit(v, group, w)             : v が group なら +w（w < 0 なら自動で費用として扱う）
//   add_cost(v, group, c)               : v が group なら -c
//   add_penalty(u, v, c)                : u が 0 かつ v が 1 なら -c（c >= 0）。逆向きは引数を入れ替える
//   add_implication(u, v)               : u が 0 なら v も 0（「v が 1 なら u も 1」と同じ）
//   force(v, group)                     : v を group に固定
//   add_profit_if_all(us, group, w)     : us が全て group なら +w（w >= 0、補助頂点を 1 個使う）
//   add_profit_if_both(u, v, group, w)  : その 2 頂点版
//   add_cost_if_any(us, group, c)       : us のどれかが group なら -c（c >= 0）
//   solve()                             : 最大利得（= 利得の総和 - 最小カット）。項を足してからもう一度呼んでもよい
//   feasible()                          : force / add_implication が矛盾していたら false（そのとき solve() の値は無意味）
//   assignment()                        : 最適解での各頂点の group（0 / 1）
//   mf                                  : 内部の MaxFlow
//   最大重み閉包（プロジェクト選択）は「選ぶ = group 0」として add_profit(v, 0, w[v]) と add_implication(v, 前提)
//   注意: 「u が 0 かつ v が 1 なら利得」のような条件（優モジュラ）は最小カットでは表せない
//   注意: 有限の重みの総和は T の最大値の 1/4 未満にする（T は long long 推奨）
template <class T = long long> struct ProjectSelection {
    explicit ProjectSelection(int n) : mf(n + 2), _n(n), S(n), TT(n + 1) {}

    void add_profit(int v, int group, T w) {          // v が group なら +w
        if (w < 0) { add_cost(v, group, -w); return; }
        base += w; add_cost(v, 1 - group, w);
    }
    void add_cost(int v, int group, T c) {             // v が group なら -c
        if (c < 0) { add_profit(v, group, -c); return; }
        if (group == 1) mf.add_edge(S, v, c);
        else mf.add_edge(v, TT, c);
    }
    void add_penalty(int u, int v, T c) {              // u が 0 かつ v が 1 なら -c
        assert(c >= 0);
        mf.add_edge(u, v, c);
    }
    void add_implication(int u, int v) {               // u が 0 なら v も 0
        mf.add_edge(u, v, INF);
    }
    void force(int v, int group) { add_cost(v, 1 - group, INF); }   // v を group に固定
    // us が全て group なら +w（w >= 0）。補助頂点を 1 個使う
    void add_profit_if_all(const std::vector<int>& us, int group, T w) {
        assert(w >= 0);
        if (us.size() == 1) { add_profit(us[0], group, w); return; }
        base += w;
        if (us.empty() || w == 0) return;
        int x = mf.add_vertex();
        if (group == 0) { mf.add_edge(S, x, w); for (int v : us) mf.add_edge(x, v, INF); }
        else { mf.add_edge(x, TT, w); for (int v : us) mf.add_edge(v, x, INF); }
    }
    void add_profit_if_both(int u, int v, int group, T w) { add_profit_if_all({u, v}, group, w); }
    // us のどれかが group なら -c（c >= 0）
    void add_cost_if_any(const std::vector<int>& us, int group, T c) {
        assert(c >= 0);
        base -= c;
        add_profit_if_all(us, 1 - group, c);
    }

    // 最大利得。項を足してからもう一度呼んでもよい
    T solve() {
        if (cut < INF) cut += mf.flow(S, TT, INF - cut);
        solved = true;
        return base - cut;
    }
    // force() / add_implication() が矛盾していたら false（そのとき solve() の値は無意味）
    bool feasible() const { return cut < INF; }
    std::vector<int> assignment() {                    // 0 = s 側、1 = t 側
        assert(solved);
        auto vis = mf.min_cut(S);
        std::vector<int> res(_n);
        for (int v = 0; v < _n; v++) res[v] = vis[v] ? 0 : 1;
        return res;
    }

    MaxFlow<T> mf;   // 内部の MaxFlow（min_cut() などを使いたいとき用に公開）

  private:
    static constexpr T INF = std::numeric_limits<T>::max() / 4;
    int _n, S, TT;
    T base = 0, cut = 0;
    bool solved = false;
};
$0
