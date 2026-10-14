// @name TwoSat
// @prefix twosat
// @description 2-SAT（要 scc スニペット）: add_clause/implies/nand/set_val/eq/ne/at_most_one/satisfiable

// 2-SAT（SCC スニペットが必要）
//   変数 x_0..x_{n-1}（bool）。add_var() で後から変数を追加できる
//   ノード 2*i+f が「x_i == f」を表す
//   add_clause(i,f,j,g)    : (x_i==f) || (x_j==g)
//   add_implies(i,f,j,g)   : (x_i==f) -> (x_j==g)
//   add_nand(i,f,j,g)      : !((x_i==f) && (x_j==g))
//   set_val(i,f)           : x_i == f
//   add_eq(i,j) / add_ne(i,j) : x_i == x_j / x_i != x_j
//   at_most_one({{i,f},...}) : 指定リテラルのうち true は高々1つ（補助変数 k-1 個、辺 O(k)）
//   satisfiable()          : 充足可能なら true を返し ans[i] に x_i の値を格納（補助変数分も含む）
struct TwoSat {
    int n;
    vector<pair<int, int>> es;
    vector<bool> ans;
    TwoSat(int n = 0) : n(n) {}
    int add_var() { return n++; }
    void add_clause(int i, bool f, int j, bool g) {
        int a = 2 * i + f, b = 2 * j + g;
        es.emplace_back(a ^ 1, b), es.emplace_back(b ^ 1, a);
    }
    void add_implies(int i, bool f, int j, bool g) { add_clause(i, !f, j, g); }
    void add_nand(int i, bool f, int j, bool g) { add_clause(i, !f, j, !g); }
    void set_val(int i, bool f) { add_clause(i, f, i, f); }
    void add_eq(int i, int j) { add_implies(i, true, j, true), add_implies(i, false, j, false); }
    void add_ne(int i, int j) { add_implies(i, true, j, false), add_implies(i, false, j, true); }
    void at_most_one(const vector<pair<int, bool>>& lits) {
        int k = lits.size();
        if (k == 2) { add_nand(lits[0].first, lits[0].second, lits[1].first, lits[1].second); return; }
        // prev : lits[0..t-1] のどれかが true
        int prev = -1;
        for (int t = 0; t < k; t++) {
            auto [i, f] = lits[t];
            if (prev != -1) add_implies(prev, true, i, !f);
            if (t + 1 < k) {
                int s = add_var();
                add_implies(i, f, s, true);
                if (prev != -1) add_implies(prev, true, s, true);
                prev = s;
            }
        }
    }
    bool satisfiable() {
        SCC scc(2 * n);
        for (auto [a, b] : es) scc.add_edge(a, b);
        scc.build();
        ans.assign(n, false);
        for (int i = 0; i < n; i++) {
            if (scc[2 * i] == scc[2 * i + 1]) return false;
            ans[i] = scc[2 * i] < scc[2 * i + 1];
        }
        return true;
    }
};
$0
