// @name MinCostFlow
// @prefix mincostflow
// @description 最小費用流（primal-dual, CSR, 負辺・b-flow・下限付き辺・slope・双対解）: add_edge/flow/slope/min_cost_any/b_flow/potential

// MinCostFlow（最小費用流）: primal-dual 法。1 ラウンド = Dijkstra + タイトな辺だけの小さなグラフ上の最大流
//   MinCostFlow<Cap, Cost> g(n);          既定はどちらも long long。Cap は符号付き整数型、Cost は double でもよい
//   add_edge(u, v, cap, cost)             : 有向辺。辺 id を返す。cost は負でもよい（負閉路がなければ）
//   add_edge(u, v, lower, upper, cost)    : 流量が [lower, upper] の辺（b_flow() で解く）
//   add_vertex() / size()                 : 頂点を 1 つ追加して番号を返す / 現在の頂点数
//   flow(s, t) / flow(s, t, limit)        : {流量, 費用}。最大流（limit まで）を最小費用で流す
//   slope(s, t) / slope(s, t, limit)      : (流量, 最小費用) の折れ線の端点。先頭は (0, 0) で下に凸
//   min_cost_any(s, t)                    : {流量, 費用}。費用が最小になる流量で止める（流すと損になる手前まで）
//   add_supply(v, b) / add_demand(v, b)   : 頂点 v の供給量 / 需要量
//   b_flow()                              : {実行可能か, 費用}。供給・需要・下限をすべて満たす最小費用流。負閉路があってもよい。1 回だけ呼べる
//   get_edge(i) / edges()                 : {from, to, cap, flow, cost}
//   potential()                           : 双対解 p。余裕のある辺で cost + p[from] - p[to] >= 0、流れている辺で <= 0
//   flow / slope / min_cost_any は何度呼んでもよい（s, t を変えてもよい）。現在の流れの続きから流し、増えた分を返す
//   辺・頂点は流した後から追加してもよい。b_flow() の後に flow() で追加で流すこともできる
//   ラウンド数は「s-t 最短路長の種類数」以下。割当問題のように同じ長さの増加路が多いほど速い
//   負辺があるときの最初のポテンシャル計算は、残余グラフが DAG なら O(V + E)、そうでなければ SPFA（最悪 O(V E)）
//   注意: flow / slope / min_cost_any は負閉路があると assert で止まる。流した後に辺を足したせいで負閉路が
//         できる（= 今の流れが最適でなくなる）場合も同じ。負閉路があるなら b_flow() を使う
//   注意: 供給・需要・下限付き辺を入れたら、b_flow() を呼ぶまで flow() 系は使えない（assert）
template <class Cap = long long, class Cost = long long> struct MinCostFlow {
    struct Edge { int from, to; Cap cap, flow; Cost cost; };

    MinCostFlow() : _n(0) {}
    explicit MinCostFlow(int n) : _n(n), excess(n, 0), node(n) {}
    void reserve(int m) { _from.reserve(m); _to.reserve(m); _cap.reserve(m); _cost.reserve(m); }
    int size() const { return _n; }
    int add_vertex() { built = false; excess.push_back(Cap(0)); node.emplace_back(); return _n++; }

    int add_edge(int from, int to, Cap cap, Cost cost) {
        assert(0 <= from && from < _n && 0 <= to && to < _n && 0 <= cap);
        _from.push_back(from); _to.push_back(to); _cap.push_back(cap); _cost.push_back(cost);
        built = false;
        // 新しい弧が現在のポテンシャルで負の被約費用になるなら、次に流す前にポテンシャルを直す
        if (cap > 0 && cost + node[from].p - node[to].p < 0) dirty = true;
        return int(_from.size()) - 1;
    }
    // 流量が [lower, upper] でなければならない辺（b_flow() で解く）
    int add_edge(int from, int to, Cap lower, Cap upper, Cost cost) {
        assert(lower <= upper);
        int m = add_edge(from, to, upper - lower, cost);
        low.resize(_from.size(), Cap(0));
        low[m] = lower;
        if (lower != Cap(0)) {
            excess[from] -= lower; excess[to] += lower;
            base_cost += lower * cost;
            has_excess = true;
        }
        return m;
    }
    void add_supply(int v, Cap b) { excess[v] += b; has_excess = true; }
    void add_demand(int v, Cap b) { excess[v] -= b; has_excess = true; }

    Edge get_edge(int i) {
        build();
        Cap fw = el[epos[i]].cap, bw = el[el[epos[i]].rev].cap;
        Cap l = i < int(low.size()) ? low[i] : Cap(0);
        return Edge{_from[i], _to[i], fw + bw + l, bw + l, _cost[i]};
    }
    std::vector<Edge> edges() {
        build();
        std::vector<Edge> res(_from.size());
        for (int i = 0; i < int(_from.size()); i++) res[i] = get_edge(i);
        return res;
    }

    // どの呼び出しも現在の流れの続きから流し、今回増えた分だけを返す
    std::pair<Cap, Cost> flow(int s, int t) { return flow(s, t, std::numeric_limits<Cap>::max()); }
    std::pair<Cap, Cost> flow(int s, int t, Cap flow_limit) { return st(s, t, flow_limit, false, false).back(); }
    std::vector<std::pair<Cap, Cost>> slope(int s, int t) { return slope(s, t, std::numeric_limits<Cap>::max()); }
    std::vector<std::pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) { return st(s, t, flow_limit, true, false); }
    // それ以上流すと費用が増えるところで止める（流量を問わない最小費用）
    std::pair<Cap, Cost> min_cost_any(int s, int t) { return st(s, t, std::numeric_limits<Cap>::max(), false, true).back(); }

    // 供給・需要・下限をすべて満たす最小費用流。負閉路があってもよい
    std::pair<bool, Cost> b_flow() {
        assert(!b_done && "b_flow() may be called only once");
        b_done = true;
        Cap total = 0;
        for (int v = 0; v < _n; v++) total += excess[v];
        if (total != Cap(0)) return {false, Cost(0)};   // 供給と需要の総量が合わない
        build();
        Cost base = base_cost;
        std::vector<Cap> ex = excess;
        // 負の弧を全部流し切る。以降、残余グラフの弧はすべて cost >= 0
        for (int v = 0; v < _n; v++)
            for (int i = start[v]; i < start[v + 1]; i++) {
                _edge& e = el[i];
                if (e.cost < 0 && e.cap > 0) {
                    base += e.cap * e.cost;
                    ex[v] -= e.cap; ex[e.to] += e.cap;
                    el[e.rev].cap += e.cap; e.cap = 0;
                }
            }
        // 超頂点 S, T を足したコピーの上で過不足を流す
        const int m = int(_from.size()), S = _n, T = _n + 1;
        MinCostFlow aux(_n + 2);
        for (int i = 0; i < m; i++) aux.add_edge(_from[i], _to[i], el[epos[i]].cap + el[el[epos[i]].rev].cap, _cost[i]);
        Cap need = 0;
        for (int v = 0; v < _n; v++) {
            if (ex[v] > 0) { aux.add_edge(S, v, ex[v], Cost(0)); need += ex[v]; }
            else if (ex[v] < 0) aux.add_edge(v, T, -ex[v], Cost(0));
        }
        aux.build();
        for (int i = 0; i < m; i++) {
            _edge& a = aux.el[aux.epos[i]];
            a.cap = el[epos[i]].cap; aux.el[a.rev].cap = el[el[epos[i]].rev].cap;
        }
        aux.dirty = false;   // 残余費用はすべて >= 0 なので、ポテンシャルは 0 のままで正しい
        auto res = aux.run(S, T, std::numeric_limits<Cap>::max(), false, false).back();
        for (int i = 0; i < m; i++) {
            const _edge& a = aux.el[aux.epos[i]];
            el[epos[i]].cap = a.cap; el[el[epos[i]].rev].cap = aux.el[a.rev].cap;
        }
        for (int v = 0; v < _n; v++) node[v].p = aux.node[v].p;
        dirty = false; has_excess = false;
        return {res.first == need, base + res.second};
    }

    // 現在の流れに対する双対解 p: 余裕のある辺では cost + p[from] - p[to] >= 0、
    // 流れている辺では <= 0
    std::vector<Cost> potential() {
        prepare();
        std::vector<Cost> res(_n);
        for (int v = 0; v < _n; v++) res[v] = node[v].p;
        return res;
    }

  private:
    struct _edge { int to, rev; Cap cap; Cost cost; };
    struct _node { Cost p = 0, dist = 0; bool done = false; };
    int _n;
    bool built = false, dirty = false, has_excess = false, b_done = false;
    std::vector<int> _from, _to;
    std::vector<Cap> _cap, excess, low, frem, fgot;
    std::vector<Cost> _cost;
    Cost base_cost = 0;
    struct _cand { int arc, from, to; Cost nd; };   // 弧 arc を通ると to に距離 nd で届く
    std::vector<_cand> cand;
    std::vector<int> start, epos, zero, tstart, tadj, level, iter[2], que, fv;
    std::vector<_edge> el;      // 始点ごとにまとめた弧（CSR）。el[el[a].rev] が逆辺
    std::vector<_node> node;    // ポテンシャルと Dijkstra の状態（キャッシュ効率のため同じ構造体に置く）
    std::vector<std::pair<Cost, int>> heap;

    void build() {
        if (built) return;
        int m = int(_from.size()), old = int(epos.size());
        std::vector<Cap> fw(m), bw(m, Cap(0));
        for (int i = 0; i < old; i++) { fw[i] = el[epos[i]].cap; bw[i] = el[el[epos[i]].rev].cap; }
        for (int i = old; i < m; i++) fw[i] = _cap[i];
        start.assign(_n + 1, 0);
        for (int i = 0; i < m; i++) { start[_from[i] + 1]++; start[_to[i] + 1]++; }
        for (int v = 0; v < _n; v++) start[v + 1] += start[v];
        std::vector<int> cur(start.begin(), start.end() - 1);
        el.assign(2 * m, _edge{0, 0, Cap(0), Cost(0)});
        epos.assign(m, 0);
        for (int i = 0; i < m; i++) {
            int a = cur[_from[i]]++, b = cur[_to[i]]++;
            el[a] = _edge{_to[i], b, fw[i], _cost[i]};
            el[b] = _edge{_from[i], a, bw[i], -_cost[i]};
            epos[i] = a;
        }
        built = true;
    }

    std::vector<std::pair<Cap, Cost>> st(int s, int t, Cap flow_limit, bool want_slope, bool stop_when_positive) {
        assert(0 <= s && s < _n && 0 <= t && t < _n && s != t);
        assert(!has_excess && "supplies or lower bounds present: use b_flow()");
        prepare();
        return run(s, t, flow_limit, want_slope, stop_when_positive);
    }

    // ポテンシャルを正しい状態に戻す: 残余グラフの全弧で cost + p[u] - p[v] >= 0
    void prepare() {
        build();
        if (!dirty) return;
        dirty = false;
        bool ok = fix_potential();
        assert(ok && "negative cycle: use b_flow()");
        (void)ok;
    }
    bool fix_potential() {
        // 残余グラフが DAG（流す前はたいていそう）なら、トポロジカル順に 1 回なめるだけ
        std::vector<int> deg(_n, 0), ord;
        ord.reserve(_n);
        for (auto& e : el) if (e.cap > 0) deg[e.to]++;
        for (int v = 0; v < _n; v++) if (deg[v] == 0) ord.push_back(v);
        for (int qi = 0; qi < int(ord.size()); qi++)
            for (int i = start[ord[qi]]; i < start[ord[qi] + 1]; i++)
                if (el[i].cap > 0 && --deg[el[i].to] == 0) ord.push_back(el[i].to);
        auto relax = [&](int v, int i) {
            if (el[i].cap <= 0 || !(node[v].p + el[i].cost < node[el[i].to].p)) return false;
            node[el[i].to].p = node[v].p + el[i].cost;
            return true;
        };
        if (int(ord.size()) == _n) {
            for (int v : ord) for (int i = start[v]; i < start[v + 1]; i++) relax(v, i);
            return true;
        }
        // そうでなければキュー版 Bellman-Ford（最悪 O(V E)）
        std::vector<int> cnt(_n, 0);
        std::vector<char> inq(_n, 1);
        ord.resize(_n);
        for (int v = 0; v < _n; v++) ord[v] = v;
        for (int head = 0, size = _n; size > 0;) {
            int v = ord[head];
            head = head + 1 == _n ? 0 : head + 1; size--; inq[v] = 0;
            for (int i = start[v]; i < start[v + 1]; i++) {
                if (!relax(v, i) || inq[el[i].to]) continue;
                if (++cnt[el[i].to] > _n) return false;   // 負閉路
                int tail = head + size; if (tail >= _n) tail -= _n;
                ord[tail] = el[i].to; size++; inq[el[i].to] = 1;
            }
        }
        return true;
    }

    // primal-dual 法。1 ラウンド = 被約費用での Dijkstra + タイトな弧（最短路上の弧）だけを使った最大流。
    // これで s-t 距離はラウンドごとに必ず増える
    std::vector<std::pair<Cap, Cost>> run(int s, int t, Cap flow_limit, bool want_slope, bool stop_when_positive) {
        const Cost inf = std::numeric_limits<Cost>::max();
        level.resize(_n); que.resize(_n); tstart.resize(_n + 1);
        fv.resize(_n + 1); frem.resize(_n + 1); fgot.resize(_n + 1);
        Cap flow = 0;
        Cost cost = 0, prev_d = 0;
        bool first = true;
        std::vector<std::pair<Cap, Cost>> result{{Cap(0), Cost(0)}};
        while (flow < flow_limit) {
            for (auto& x : node) { x.dist = inf; x.done = false; }
            heap.clear(); zero.clear(); cand.clear();
            node[s].dist = 0;
            zero.push_back(s);
            bool hit = false;   // t が確定した。以降は t と同じ距離の頂点だけ見ればよい
            while (true) {
                int v;
                if (!zero.empty()) { v = zero.back(); zero.pop_back(); }   // 直前に取り出した頂点と同じ距離
                else if (!heap.empty() && !(hit && heap.front().first > node[t].dist)) {
                    std::pop_heap(heap.begin(), heap.end(), std::greater<>());
                    v = heap.back().second; heap.pop_back();
                } else break;
                if (node[v].done) continue;
                node[v].done = true;
                if (v == t) { hit = true; continue; }
                const Cost dv = node[v].dist, pv = node[v].p;
                for (int i = start[v]; i < start[v + 1]; i++) {
                    const _edge& e = el[i];
                    if (e.cap <= 0) continue;
                    _node& w = node[e.to];
                    Cost nd = dv + (e.cost + pv - w.p);
                    if (nd > w.dist) continue;                    // ほとんどの弧はここで終わる
                    if (w.done) {                                 // 確定済みの頂点は書き換えない
                        if (nd == w.dist) cand.push_back({i, v, e.to, nd});
                        continue;
                    }
                    if (hit && nd > node[t].dist) continue;
                    cand.push_back({i, v, e.to, nd});
                    if (nd < w.dist) {
                        w.dist = nd;
                        if (nd == dv) zero.push_back(e.to);
                        else { heap.push_back({nd, e.to}); std::push_heap(heap.begin(), heap.end(), std::greater<>()); }
                    }
                }
            }
            if (!hit) break;
            // タイトな弧 = 終点がちょうどその距離で確定した候補。保存した値どうしを比べるので、
            // 実数コストでも Dijkstra の木の辺は必ず残る（毎ラウンド必ず 1 本は流せる）
            std::fill(tstart.begin(), tstart.end(), 0);
            int nt = 0;
            for (const _cand& c : cand) {
                if (!node[c.to].done || c.nd != node[c.to].dist) continue;
                cand[nt++] = c;
                tstart[c.from + 1]++; tstart[c.to + 1]++;
            }
            cand.resize(nt);
            for (int v = 0; v < _n; v++) tstart[v + 1] += tstart[v];
            tadj.resize(2 * nt);
            iter[0].assign(tstart.begin(), tstart.end() - 1);
            for (const _cand& c : cand) { tadj[iter[0][c.from]++] = c.arc; tadj[iter[0][c.to]++] = el[c.arc].rev; }

            const Cost dt = node[t].dist;
            for (auto& x : node) if (x.done) x.p -= dt - x.dist;
            Cost d = node[t].p - node[s].p;   // タイトな s-t パスの 1 単位あたりの本当の費用
            if (stop_when_positive && d > 0) break;

            Cap c = 0;
            while (c < flow_limit - flow && bfs(s, t)) c += blocking(s, t, flow_limit - flow - c);
            flow += c; cost += c * d;
            if (want_slope && !first && prev_d == d) result.pop_back();
            result.push_back({flow, cost});
            prev_d = d; first = false;
        }
        return result;
    }

    // ---- タイトな弧の上の Dinic: tadj[tstart[v] .. tstart[v+1]) が v から出る弧の番号 ----
    bool bfs(int s, int t) {
        std::fill(level.begin(), level.end(), -1);
        level[s] = 0;
        int qh = 0, qt = 0;
        que[qt++] = s;
        while (qh < qt) {
            int v = que[qh++];
            for (int i = tstart[v]; i < tstart[v + 1]; i++) {
                const _edge& e = el[tadj[i]];
                if (e.cap <= 0 || level[e.to] >= 0) continue;
                level[e.to] = level[v] + 1;
                if (e.to == t) return true;
                que[qt++] = e.to;
            }
        }
        return false;
    }
    // t からの後退 DFS と s からの前進 DFS を、予算を倍々にしながら交互に進める（MaxFlow と同じ）
    Cap blocking(int s, int t, Cap up) {
        iter[1].assign(tstart.begin(), tstart.end() - 1);
        bool fw = false, used_fwd = false;
        Cap res = 0;
        for (long long budget = _n + (long long)tadj.size();; budget *= 2, fw = !fw) {
            if (fw && !used_fwd) { iter[0].assign(tstart.begin(), tstart.end() - 1); used_fwd = true; }
            res += push(fw, fw ? s : t, fw ? t : s, up - res, budget);
            if (level[fw ? s : t] == -2 || up - res <= 0) return res;
        }
    }
    Cap push(bool fw, int src, int dst, Cap up, long long budget) {
        std::vector<int>& it = iter[fw ? 0 : 1];
        const int step = fw ? 1 : -1;
        bool quit = false;
        int k = 0;
        fv[0] = src; frem[0] = up; fgot[0] = 0;
        auto arc = [&](int i) -> _edge& { _edge& x = el[tadj[i]]; return fw ? x : el[x.rev]; };
        while (true) {
            int v = fv[k];
            Cap d;
            if (v == dst) d = frem[k];
            else {
                int& i = it[v];
                const int e = tstart[v + 1], i0 = i;
                while (i < e && (level[el[tadj[i]].to] != level[v] + step || arc(i).cap <= 0)) i++;
                budget -= i - i0 + 1;
                if (i < e && budget >= 0) {
                    fv[k + 1] = el[tadj[i]].to; frem[k + 1] = std::min(frem[k], arc(i).cap); fgot[k + 1] = 0;
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
                _edge& a = arc(it[fv[k]]);
                a.cap -= d; el[a.rev].cap += d;
                frem[k] -= d; fgot[k] += d;
                if (frem[k] > 0 && !quit) break;
                d = fgot[k];
            }
        }
    }
};
$0
