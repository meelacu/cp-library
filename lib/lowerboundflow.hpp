// @name LowerBoundFlow
// @prefix lowerboundflow
// @description 下限付き流量（要 maxflow スニペット）: add_edge(lower, upper)/feasible/max_flow/min_flow/get_flow

// LowerBoundFlow（下限付き流量）（MaxFlow スニペットが必要）
//   LowerBoundFlow<Cap> g(n);
//   add_edge(u, v, lower, upper)  : 流量が [lower, upper] の辺。辺 id を返す（最初の問い合わせより前に全部足す）
//   feasible()                    : 全辺の下限・上限を満たす循環流があるか
//   max_flow(s, t)                : {実行可能か, 値}。下限・上限を満たす s-t フローの流量の最大値
//   min_flow(s, t)                : 同じく最小値
//   get_flow(i)                   : 最後の問い合わせの解で辺 i に流れている量
//   mf                            : 内部の MaxFlow（min_cut() などを使いたいとき）
//   max_flow と min_flow は両方呼んでよい（同じ s, t で）。feasible() も使うなら最初に呼ぶ
//   値は t -> s 向きに流せる構造だと負になりうる（それが本当の最小値）。非負に限るなら max(0, 値) とする
//   費用も付くなら MinCostFlow の add_edge(u, v, lower, upper, cost) と b_flow() を使う
template <class Cap = long long> struct LowerBoundFlow {
    explicit LowerBoundFlow(int n) : mf(n + 2), _n(n), ex(n, 0) {}

    int add_edge(int from, int to, Cap lower, Cap upper) {
        assert(!started && "add every edge before the first query");
        assert(0 <= lower && lower <= upper);
        low.push_back(lower);
        ex[from] -= lower; ex[to] += lower;
        return mf.add_edge(from, to, upper - lower);
    }
    // 全辺の [lower, upper] を満たす循環流があるか（max_flow / min_flow より前に呼ぶ）
    bool feasible() {
        assert(_s < 0 && "feasible() must come before max_flow() / min_flow()");
        return saturate();
    }
    // {実行可能か, 値}: 下限・上限を満たす s-t フローの流量の最大値 / 最小値。
    // 両方呼んでよい（同じ s, t で）。get_flow() は最後に呼んだほうの解を返す
    std::pair<bool, Cap> max_flow(int s, int t) {
        if (!open(s, t)) return {false, Cap(0)};
        cur += mf.flow(s, t);
        return {true, cur};
    }
    std::pair<bool, Cap> min_flow(int s, int t) {
        if (!open(s, t)) return {false, Cap(0)};
        cur -= mf.flow(t, s);
        return {true, cur};
    }
    Cap get_flow(int i) { return low[i] + mf.get_edge(i).flow; }
    MaxFlow<Cap> mf;   // 内部の MaxFlow（min_cut() などを使いたいとき用に公開）

  private:
    int _n, _s = -1, _t = -1;
    std::vector<Cap> ex, low;
    Cap need = 0, got = 0, cur = 0;
    bool started = false, ok = false;

    // 超頂点 S から流せるだけ流す。全部流れれば実行可能
    bool saturate() {
        if (!started) {
            started = true;
            for (int v = 0; v < _n; v++) {
                if (ex[v] > 0) { mf.add_edge(_n, v, ex[v]); need += ex[v]; }
                else if (ex[v] < 0) mf.add_edge(v, _n + 1, -ex[v]);
            }
        }
        got += mf.flow(_n, _n + 1);
        return got == need;
    }
    bool open(int s, int t) {
        if (_s < 0) {
            assert(0 <= s && s < _n && 0 <= t && t < _n && s != t);
            _s = s; _t = t;
            // t-s 間は両向きに流せる辺にする（下限のせいで s-t 流量が負になることがある）
            int e = mf.add_undirected_edge(t, s, std::numeric_limits<Cap>::max() / 2);
            ok = saturate();
            cur = mf.get_edge(e).flow;
            mf.change_edge(e, 0, 0);
        }
        assert(_s == s && _t == t && "max_flow / min_flow need the same s and t");
        return ok;
    }
};
$0
