// @name Tree Diameter
// @prefix treediameter
// @description 木の直径(経路つき)・各頂点の最遠距離と最遠点・中心・半径、木の2点間距離とパス。重みなし/重み付き両対応

// 木の直径まわり。main の外に置く
//   G[v] = 隣接リスト。重みなしなら vvi、重み付きなら vector<vector<pair<int, W>>>（{行き先, 重み}。重みは 0 以上）
//          .to（と .cost）を持つ struct の辺でもよい。閉路・多重辺・自己ループがあると assert で止まる
//   距離の型 D は、重みなし → int、重みが整数 → long long、重みが実数 → その型
//
//   auto td = tree_diameter(G);       O(n)。森でもよい（下の「森のとき」を参照）
//     td.len               直径の長さ（重みなしなら辺数、重み付きなら重みの和）
//     td.path              直径を端から端まで通る頂点列
//     td.u, td.v           直径の両端（= td.path.front(), td.path.back()）
//     td.ecc(x)            x から最も遠い頂点までの距離          td.far(x)   その最も遠い頂点（x の木の直径の端のどちらか）
//     td.center()          ecc が最小の頂点（木の中心）           td.centers()  直径上の中心すべて（u 側から順。重みが正なら 1 個か 2 個で、これで全部）
//     td.radius()          中心の ecc
//     td.du[x], td.dv[x]   直径の端 u, v から x までの距離
//     td.pu[x], td.pv[x]   直径の端 u, v を根にしたときの x の親（根は -1）
//     森のとき: len / u / v / path / center / radius は直径が最大の木のもの。du / dv / pu / pv / ecc / far は各頂点が属する木について正しい
//               td.trees = 木ごとの {len, u, v}（最小の頂点番号の昇順）、td.comp[x] = x が属する木の番号 → x の木は td.trees[td.comp[x]]
//               k 番目の木の中心・半径・直径の経路は td.center(k) / td.centers(k) / td.radius(k) / td.path_of(k)
//   vector<D> d = tree_dist(G, s);    s から各頂点までの距離（別の木の頂点は -1）。O(n)
//   vi p = tree_path(G, a, b);        a から b までのパスの頂点列（別の木なら空）。O(n)
// 例: x から距離ちょうど k の頂点を 1 つ（重みなし）: k <= td.ecc(x) のときだけ存在し、x から td.far(x) へ k 歩進んだ頂点
//     = td.du[x] >= td.dv[x] なら pu を、そうでなければ pv を親とみて k 個上（doubling スニペット: Doubling<>(td.pu, n).jump(x, k) で O(log n)）
namespace tree_diameter_internal {
template <class E>
int to(const E& e) {
    if constexpr (is_integral_v<E>) return e;
    else if constexpr (requires { e.to; }) return e.to;
    else return get<0>(e);
}
// 重み付きの辺か: pair・tuple（2 番目が重み）か、.cost を持つ struct
template <class E>
constexpr bool has_weight = !is_integral_v<E> && (requires(const E& e) { e.cost; } || !requires(const E& e) { e.to; });
template <class E>
auto weight(const E& e) {
    if constexpr (requires { e.cost; }) return e.cost;
    else return get<1>(e);
}
template <class E>
auto dist_zero() {
    if constexpr (!has_weight<E>) return 0;
    else if constexpr (is_integral_v<decay_t<decltype(weight(declval<E>()))>>) return 0LL;
    else return decay_t<decltype(weight(declval<E>()))>(0);
}
template <class Graph> using edge_t = typename Graph::value_type::value_type;
template <class Graph> using dist_t = decltype(dist_zero<edge_t<Graph>>());

// s から届く範囲をたどる。{距離（届かない頂点は -1）, 親（s と届かない頂点は -1）}
template <class Graph>
pair<vector<dist_t<Graph>>, vector<int>> walk(const Graph& G, int s) {
    using E = edge_t<Graph>;
    using D = dist_t<Graph>;
    int n = G.size();
    assert(0 <= s && s < n);
    vector<D> d(n, D(-1));
    vector<int> par(n, -1), q{s};
    d[s] = 0;
    for (int i = 0; i < (int)q.size(); i++) {
        int v = q[i];
        for (auto& e : G[v]) {
            int u = to(e);
            if (d[u] >= 0) continue;
            if constexpr (has_weight<E>) {
                assert(weight(e) >= 0);
                d[u] = d[v] + (D)weight(e);
            } else d[u] = d[v] + 1;
            par[u] = v;
            q.push_back(u);
        }
    }
    return {d, par};
}
}  // namespace tree_diameter_internal

template <class D>
struct TreeDiameter {
    struct Tree {
        D len;
        int u, v;
    };
    D len = 0;
    int u = -1, v = -1;
    vector<int> path;
    vector<D> du, dv;
    vector<int> pu, pv, comp;
    vector<Tree> trees;
    D ecc(int x) const { return max(du[x], dv[x]); }
    int far(int x) const { return du[x] >= dv[x] ? trees[comp[x]].u : trees[comp[x]].v; }
    // 以下は k 番目の木について（k を省くと直径が最大の木）
    vector<int> path_of(int k) const {
        vector<int> res;
        for (int x = trees[k].v; x != -1; x = pu[x]) res.push_back(x);
        reverse(res.begin(), res.end());
        return res;
    }
    vector<int> centers(int k) const {  // ecc が最小の頂点は必ず直径の上にあるので、直径の上だけを見る
        vector<int> res;
        D best = ecc(trees[k].v);
        for (int x = trees[k].v; x != -1; x = pu[x]) best = min(best, ecc(x));
        for (int x = trees[k].v; x != -1; x = pu[x])
            if (ecc(x) == best) res.push_back(x);
        reverse(res.begin(), res.end());
        return res;
    }
    vector<int> centers() const { return u < 0 ? vector<int>() : centers(comp[u]); }
    int center(int k) const { return centers(k)[0]; }
    int center() const {
        assert(u >= 0);
        return center(comp[u]);
    }
    D radius(int k) const { return ecc(center(k)); }
    D radius() const { return ecc(center()); }
};

template <class Graph>
auto tree_dist(const Graph& G, int s) {
    return tree_diameter_internal::walk(G, s).first;
}

template <class Graph>
vector<int> tree_path(const Graph& G, int a, int b) {
    auto [d, par] = tree_diameter_internal::walk(G, a);
    vector<int> res;
    if (d[b] < 0) return res;
    for (int x = b; x != -1; x = par[x]) res.push_back(x);
    reverse(res.begin(), res.end());
    return res;
}

template <class Graph>
auto tree_diameter(const Graph& G) {
    using tree_diameter_internal::has_weight, tree_diameter_internal::to, tree_diameter_internal::weight;
    using E = tree_diameter_internal::edge_t<Graph>;
    using D = tree_diameter_internal::dist_t<Graph>;
    int n = G.size();
    TreeDiameter<D> t;
    t.du.assign(n, D(-1)), t.dv.assign(n, D(-1)), t.pu.assign(n, -1), t.pv.assign(n, -1), t.comp.assign(n, -1);
    vector<int> q;
    // s の木を歩いて d（-1 の頂点だけ進む）と par を埋め、s から最も遠い頂点を返す。歩いた頂点は q に残る
    auto sweep = [&](int s, vector<D>& d, vector<int>& par) {
        q.assign(1, s);
        d[s] = 0, par[s] = -1;
        int far = s;
        for (int i = 0; i < (int)q.size(); i++) {
            int x = q[i];
            if (d[x] > d[far]) far = x;
            for (auto& e : G[x]) {
                int y = to(e);
                if (d[y] >= 0) continue;
                if constexpr (has_weight<E>) {
                    assert(weight(e) >= 0);  // 負の重みがあると「最も遠い頂点は直径の端」が成り立たない
                    d[y] = d[x] + (D)weight(e);
                } else d[y] = d[x] + 1;
                par[y] = x;
                q.push_back(y);
            }
        }
        return far;
    };
    for (int s = 0; s < n; s++) {
        if (t.comp[s] >= 0) continue;
        int a = sweep(s, t.dv, t.pv);  // s から最も遠い頂点 a は直径の端（dv, pv は仮置きなので消す）
        for (int x : q) t.dv[x] = -1;
        int b = sweep(a, t.du, t.pu);  // a から最も遠い頂点 b がもう一方の端
        size_t deg = 0;
        for (int x : q) t.comp[x] = t.trees.size(), deg += G[x].size();
        assert(deg == 2 * (q.size() - 1));  // 閉路・多重辺・自己ループがない（木である）
        sweep(b, t.dv, t.pv);
        t.trees.push_back({t.du[b], a, b});
        if (t.u < 0 || t.du[b] > t.len) t.len = t.du[b], t.u = a, t.v = b;
    }
    for (int x = t.v; x != -1; x = t.pu[x]) t.path.push_back(x);
    reverse(t.path.begin(), t.path.end());
    return t;
}
$0
