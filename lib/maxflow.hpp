// @name MaxFlow
// @prefix maxflow
// @description 最大流・最小カット（双方向Dinic, CSR, 非再帰）: add_edge/add_undirected_edge/add_vertex/flow/min_cut/min_cut_edges/get_edge/change_edge/reset

// MaxFlow（最大流 / 最小カット）: Dinic + CSR、非再帰。O(V^2 E)、単位容量や二部グラフは O(E sqrt(V))
//   MaxFlow<Cap> g(n);               Cap の既定は long long。double も可（eps 未満の残余容量は 0 扱い、既定 1e-9）
//   add_edge(u, v, cap)              : 有向辺。辺 id を返す
//   add_undirected_edge(u, v, cap)   : 無向辺（どちら向きにも cap まで流せる）
//   add_vertex() / size()            : 頂点を 1 つ追加して番号を返す / 現在の頂点数
//   flow(s, t) / flow(s, t, limit)   : 最大流（limit まで）。もう一度呼ぶと残余グラフの続きから流し、増えた分を返す
//   min_cut(s)                       : flow() の直後に呼ぶ。最小カットの s 側なら true の vector<bool>
//   min_cut_edges(s)                 : その最小カットを跨ぐ辺 id（容量の和 = 最大流）
//   get_edge(i) / edges()            : {from, to, cap, flow}。無向辺の flow は from->to 向きが正、逆向きが負
//   change_edge(i, cap, flow)        : 辺 i の容量と流量を直接書き換える
//   reset()                          : 全辺の流量を 0 に戻す（二分探索などで同じグラフを使い回すとき）
//   scaling = true                   : 容量スケーリング O(V E log U)。容量が巨大で遅いときだけ試す
//   辺・頂点は flow() の後から追加してもよい（流れは保たれる）
//   ブロッキングフローは「t からの後退 DFS」と「s からの前進 DFS」を予算倍々で交互に進めるので、
//   多くの増加路が共有する長い一本道が s 側・t 側のどちらにあっても O(V + E) で済む
//   注意: limit で打ち切った直後の min_cut() は最小カットではない
//   注意: 無向辺を使うときは Cap を符号付きの型にし、容量は最大値の半分以下にする（残余容量が 2*cap まで増える）
template <class Cap = long long> struct MaxFlow {
    struct Edge { int from, to; Cap cap, flow; };

    MaxFlow() : _n(0) {}
    explicit MaxFlow(int n) : _n(n) {}
    void reserve(int m) { _from.reserve(m); _to.reserve(m); _cap.reserve(m); }
    int size() const { return _n; }
    int add_vertex() { built = false; return _n++; }

    int add_edge(int from, int to, Cap cap) { return add(from, to, cap, false); }
    int add_undirected_edge(int u, int v, Cap cap) { return add(u, v, cap, true); }

    // 無向辺の flow は from -> to 向きが正、to -> from 向きが負
    Edge get_edge(int i) {
        build();
        return Edge{_from[i], _to[i], _cap[i], _cap[i] - el[epos[i]].cap};
    }
    std::vector<Edge> edges() {
        build();
        std::vector<Edge> res(_from.size());
        for (int i = 0; i < int(_from.size()); i++) res[i] = get_edge(i);
        return res;
    }
    void change_edge(int i, Cap new_cap, Cap new_flow) {
        build();
        _cap[i] = new_cap;
        el[epos[i]].cap = new_cap - new_flow;
        el[el[epos[i]].rev].cap = _und[i] ? new_cap + new_flow : new_flow;
    }
    // 全辺の流量を 0 に戻す（同じグラフを使い回すとき）
    void reset() {
        build();
        for (int i = 0; i < int(_from.size()); i++) {
            el[epos[i]].cap = _cap[i];
            el[el[epos[i]].rev].cap = _und[i] ? _cap[i] : Cap(0);
        }
    }

    bool scaling = false;  // 容量スケーリング O(V E log U)。容量が巨大なとき用
    Cap eps = std::is_floating_point<Cap>::value ? Cap(1e-9) : Cap(1);  // eps 未満の残余容量は 0 とみなす

    Cap flow(int s, int t) { return flow(s, t, std::numeric_limits<Cap>::max()); }
    Cap flow(int s, int t, Cap flow_limit) {
        assert(0 <= s && s < _n && 0 <= t && t < _n && s != t);
        build();
        level.assign(_n, -1); que.resize(_n);
        fv.resize(_n + 1); frem.resize(_n + 1); fgot.resize(_n + 1);
        Cap flow = 0, lim = eps;
        if (scaling) { Cap mx = 0; for (auto& e : el) mx = std::max(mx, e.cap); while (lim <= mx / 2) lim *= 2; }
        while (true) {
            while (flow_limit - flow >= eps && bfs(s, t, lim)) flow += blocking(s, t, flow_limit - flow, lim);
            if (!scaling || lim / 2 < eps) break;
            lim /= 2;
        }
        return flow;
    }

    // 最小カットの s 側の頂点が true（flow() の直後に呼ぶ）
    std::vector<bool> min_cut(int s) {
        build();
        std::vector<bool> vis(_n, false);
        std::vector<int> stk{s};
        vis[s] = true;
        while (!stk.empty()) {
            int v = stk.back(); stk.pop_back();
            for (int i = start[v]; i < start[v + 1]; i++)
                if (el[i].cap >= eps && !vis[el[i].to]) { vis[el[i].to] = true; stk.push_back(el[i].to); }
        }
        return vis;
    }
    // その最小カットを跨ぐ辺の id（容量の和 = 最大流）
    std::vector<int> min_cut_edges(int s) {
        auto vis = min_cut(s);
        std::vector<int> res;
        for (int i = 0; i < int(_from.size()); i++)
            if (vis[_from[i]] != vis[_to[i]] && (vis[_from[i]] || _und[i])) res.push_back(i);
        return res;
    }

  private:
    struct _edge { int to, rev; Cap cap; };
    int _n;
    bool built = false;
    std::vector<int> _from, _to;
    std::vector<Cap> _cap;
    std::vector<bool> _und;
    std::vector<int> start, epos, level, iter[2], que, fv;
    std::vector<Cap> frem, fgot;
    std::vector<_edge> el;   // 始点ごとにまとめた弧（CSR）。el[el[a].rev] が逆辺

    int add(int from, int to, Cap cap, bool und) {
        assert(0 <= from && from < _n && 0 <= to && to < _n && 0 <= cap);
        _from.push_back(from); _to.push_back(to); _cap.push_back(cap); _und.push_back(und);
        built = false;
        return int(_from.size()) - 1;
    }

    void build() {
        if (built) return;
        int m = int(_from.size()), old = int(epos.size());
        std::vector<Cap> fw(m), bw(m);
        for (int i = 0; i < old; i++) { fw[i] = el[epos[i]].cap; bw[i] = el[el[epos[i]].rev].cap; }
        for (int i = old; i < m; i++) { fw[i] = _cap[i]; bw[i] = _und[i] ? _cap[i] : Cap(0); }
        start.assign(_n + 1, 0);
        for (int i = 0; i < m; i++) { start[_from[i] + 1]++; start[_to[i] + 1]++; }
        for (int v = 0; v < _n; v++) start[v + 1] += start[v];
        std::vector<int> cur(start.begin(), start.end() - 1);
        el.assign(2 * m, _edge{0, 0, Cap(0)});
        epos.assign(m, 0);
        for (int i = 0; i < m; i++) {
            int a = cur[_from[i]]++, b = cur[_to[i]]++;
            el[a] = _edge{_to[i], b, fw[i]};
            el[b] = _edge{_from[i], a, bw[i]};
            epos[i] = a;
        }
        built = true;
    }

    bool bfs(int s, int t, Cap lim) {
        std::fill(level.begin(), level.end(), -1);
        level[s] = 0;
        int qh = 0, qt = 0;
        que[qt++] = s;
        while (qh < qt) {
            int v = que[qh++];
            for (int i = start[v]; i < start[v + 1]; i++) {
                int w = el[i].to;
                if (el[i].cap < lim || level[w] >= 0) continue;
                level[w] = level[v] + 1;
                if (w == t) return true;
                que[qt++] = w;
            }
        }
        return false;
    }

    // 現在のレベルグラフのブロッキングフロー。非再帰の多重 push DFS を 2 方向用意し、手数の予算を
    // 倍々にしながら交互に進める: t からの後退（多くの増加路が共有する末尾は 1 回しか歩かない）と
    // s からの前進（共有する先頭は 1 回しか歩かない）。グラフに合う向きがフェーズを終わらせるので、
    // 長い一本道が s 側・t 側のどちらにあっても O(V E) ではなく O(V + E) で済む
    Cap blocking(int s, int t, Cap up, Cap lim) {
        iter[1].assign(start.begin(), start.end() - 1);
        bool fw = false, used_fwd = false;
        Cap res = 0;
        for (long long budget = _n + (long long)el.size();; budget *= 2, fw = !fw) {
            if (fw && !used_fwd) { iter[0].assign(start.begin(), start.end() - 1); used_fwd = true; }
            res += push(fw, fw ? s : t, fw ? t : s, up - res, lim, budget);
            if (level[fw ? s : t] == -2 || up - res < eps) return res;
        }
    }

    // src から dst へ向かう DFS を 1 回走らせる。budget 手を超えたら（流した分は保ったまま）中断する
    Cap push(bool fw, int src, int dst, Cap up, Cap lim, long long budget) {
        std::vector<int>& it = iter[fw ? 0 : 1];
        const int step = fw ? 1 : -1;
        bool quit = false;
        int k = 0;
        fv[0] = src; frem[0] = up; fgot[0] = 0;
        while (true) {
            int v = fv[k];
            Cap d;
            if (v == dst) d = frem[k];
            else {
                int& i = it[v];
                const int e = start[v + 1], i0 = i;
                while (i < e && (level[el[i].to] != level[v] + step || (fw ? el[i].cap : el[el[i].rev].cap) < lim)) i++;
                budget -= i - i0 + 1;
                if (i < e && budget >= 0) {
                    fv[k + 1] = el[i].to; frem[k + 1] = std::min(frem[k], fw ? el[i].cap : el[el[i].rev].cap); fgot[k + 1] = 0;
                    k++;
                    continue;
                }
                if (i < e) quit = true;   // 予算切れ: 何も dead にせず巻き戻す
                else level[v] = -2;       // 使える弧がもうない: このフェーズ中は dead
                d = fgot[k];
            }
            while (true) {                // フレーム k が流せた量 d を親に返す
                if (k == 0) return d;
                k--;
                _edge& a = fw ? el[it[fv[k]]] : el[el[it[fv[k]]].rev];
                a.cap -= d; el[a.rev].cap += d;
                frem[k] -= d; fgot[k] += d;
                if (frem[k] >= eps && !quit) break;
                d = fgot[k];
            }
        }
    }
};
$0
