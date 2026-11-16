// @name Centroid Decomposition
// @prefix centroid
// @description 木の重心・重心分解: 重心木(par/dep/order)、各重心の担当範囲の列挙(部分木ごと・距離つき)、重心祖先と距離。重みなし/重み付き・森も可

// 木の重心と重心分解。main の外に置く
//   G[v] = 隣接リスト。重みなしなら vvi、重み付きなら vector<vector<pair<int, W>>>（{行き先, 重み}）。森でもよい
//          .to（と .cost）を持つ struct の辺でもよい。閉路・多重辺・自己ループがあると assert で止まる
//   距離の型 D は、重みなし → int、重みが整数 → long long、重みが実数 → その型
//
//   vi cs = tree_centroids(G);        木の重心（1 個か 2 個、昇順）。森なら tree_centroids(G, s) で s を含む木の重心。O(n)
//
//   CentroidDecomposition cd(G);      O(n log n)。型は CentroidDecomposition<D>。作ったあとは G を使わない（G を書き換えても壊れない）
//   cd.par[v]     重心木での v の親（根 = 各木の重心は -1）   cd.dep[v]   重心木での v の深さ（根は 0、最大でも log2(n)）
//   cd.sz[v]      v の担当範囲（v が重心に選ばれた時点の連結成分 = 重心木での v の部分木）の頂点数
//   cd.order      全頂点を重心木の浅い順（dep の昇順）に並べたもの（逆順になめると重心木の葉から根への DP になる）
//   ─ 使い方 1: 重心ごとに「その重心を通るパス」をまとめて処理する（パスの数え上げなど）
//   auto& t = cd.component(c);        重心 c の担当範囲を、c を根にした木として取り出す。O(大きさ)、全頂点分で合計 O(n log n)
//       t.vs[i]     頂点（t.vs[0] = c。親が先に並ぶ）      t.par[i]    親の位置（t.vs の添字。t.par[0] = -1）
//       t.dist[i]   c からの距離                             t.dep[i]    c からの辺数
//       t.sub       c の子の部分木ごとの区間 {l, r}（t.vs[l..r) がひとつの部分木で、t.vs[l] が c の隣の頂点）   t.size()  頂点数
//       戻り値は内部バッファへの参照（次に component を呼ぶと書き換わる）。辺ごとの重みは t.dist[i] - t.dist[t.par[i]]
//   ─ 使い方 2: 頂点ごとに O(log n) 個の重心祖先をたどる（点更新 + 質問など）
//   cd.ancestors(v, f)                v の重心祖先 a（v 自身から根へ順に）ごとに f(a, d) を呼ぶ。d = v と a の距離
//                                     f(a, d, ch) と 3 引数で受けると、ch = a の子（重心木）のうち v がいる側（a = v のときは -1）
//   D cd.dist(v, a)                   v とその重心祖先 a の距離。O(1)（a が v の重心祖先でないときの値は無意味）
//   どの 2 頂点 x, y（同じ木）についても、x と y の重心木での最小共通祖先 a は x-y パス上にあり、dist(x, y) = dist(x, a) + dist(y, a)
//
// 例 1: 距離がちょうど K の頂点対の個数（重みなし）
//   ll ans = 0;
//   for (int c : cd.order) {
//       auto& t = cd.component(c);
//       vi cnt(t.size());                       // cnt[d] = c と、これまでに見た部分木のうち、c からの距離が d の頂点数
//       cnt[0] = 1;
//       for (auto [l, r] : t.sub) {
//           for (int i = l; i < r; i++) if (0 <= K - t.dist[i] && K - t.dist[i] < t.size()) ans += cnt[K - t.dist[i]];
//           for (int i = l; i < r; i++) cnt[t.dist[i]]++;
//       }
//   }
// 例 2: 印をつけた頂点のうち最も近いものまでの距離（印をつける・質問する、どちらも O(log n)）
//   vll best(n, LINF);
//   印をつける: cd.ancestors(v, [&](int a, ll d) { chmin(best[a], d); });
//   質問:       ll res = LINF; cd.ancestors(v, [&](int a, ll d) { chmin(res, best[a] + d); });
// 例 3: v から距離 r 以内の頂点数。a の担当範囲で数えてから、同じ子の側（パスが a を通らない相手）の分を引く
//   vector<vll> own(n), up(n);        // own[a] = a の担当範囲の各頂点と a の距離、up[c] = c の担当範囲の各頂点と par[c] の距離
//   rep(x, n) cd.ancestors(x, [&](int a, ll d, int ch) { own[a].push_back(d); if (ch != -1) up[ch].push_back(d); });
//   rep(a, n) sort(all(own[a])), sort(all(up[a]));
//   質問:       ll res = 0; cd.ancestors(v, [&](int a, ll d, int ch) { res += ub(own[a], r - d); if (ch != -1) res -= ub(up[ch], r - d); });
//   （重心ごとに「距離を添字にした BIT」などを持つなら、大きさは cd.sz[a] + 1 で足りる。重みなしのとき own[a] の値は sz[a] 未満、up[c] の値は sz[c] 以下）
namespace centroid_internal {
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
}  // namespace centroid_internal

template <class Graph>
vector<int> tree_centroids(const Graph& G, int s = 0) {
    using centroid_internal::to;
    int n = G.size();
    assert(0 <= s && s < n);
    vector<int> par(n, -1), sz(n, 1), mx(n, 0), q{s}, res;
    vector<char> seen(n);
    seen[s] = 1;
    for (int i = 0; i < (int)q.size(); i++)
        for (auto& e : G[q[i]]) {
            int u = to(e);
            if (!seen[u]) seen[u] = 1, par[u] = q[i], q.push_back(u);
        }
    int tot = q.size();
    for (int i = tot - 1; i >= 1; i--) {
        int v = q[i], p = par[v];
        sz[p] += sz[v];
        mx[p] = max(mx[p], sz[v]);
    }
    for (int v : q)
        if (2 * max(mx[v], tot - sz[v]) <= tot) res.push_back(v);  // v を取り除いた残りの成分がどれも半分以下
    sort(res.begin(), res.end());
    return res;
}

template <class D = int>
struct CentroidDecomposition {
    struct Component {
        vector<int> vs, par, dep;
        vector<D> dist;
        vector<pair<int, int>> sub;
        int size() const { return vs.size(); }
    };
    int n = 0;
    vector<int> par, dep, sz, order;

    CentroidDecomposition() = default;
    template <class Graph>
    explicit CentroidDecomposition(const Graph& G) { build(G); }

    D dist(int v, int a) const { return dt[(size_t)dep[a] * n + lab[v]]; }
    template <class F>
    void ancestors(int v, F f) const {
        const size_t x = lab[v];
        for (int a = v, ch = -1, k = dep[v]; a != -1; ch = a, a = par[a], k--) {  // k = dep[a]
            const D d = dt[(size_t)k * n + x];
            if constexpr (invocable<F&, int, D, int>) f(a, d, ch);
            else f(a, d);
        }
    }

    const Component& component(int c) const {
        Component& t = buf;
        const int xc = lab[c], d = ld[xc];
        t.vs.assign(1, c), t.par.assign(1, -1), t.dep.assign(1, 0), t.dist.assign(1, D(0)), t.sub.clear();
        bl.assign(1, xc);
        auto push = [&](int y, int pi, D w) {  // ラベル y の頂点を、位置 pi の頂点の子として足す
            int nd = t.dep[pi] + 1;
            D ndist = t.dist[pi] + w;
            bl.push_back(y), t.vs.push_back(orig[y]), t.par.push_back(pi), t.dep.push_back(nd), t.dist.push_back(ndist);
        };
        auto grow = [&](int l) {  // 位置 l の頂点から先へ広げて、c の子の部分木を 1 つ取り終える
            for (int i = l; i < (int)bl.size(); i++) {
                int x = bl[i], px = bl[t.par[i]];
                for (int y = cl[x]; y < cr[x]; y++)
                    if (ld[y] > d && y != px) push(y, i, wof(y));
                if (int p = P[x]; p >= 0 && p != px && ld[p] > d) push(p, i, wof(x));
            }
            t.sub.push_back({l, (int)bl.size()});
        };
        if (int p = P[xc]; p >= 0 && ld[p] > d) push(p, 0, wof(xc)), grow(1);
        for (int y = cl[xc]; y < cr[xc]; y++)
            if (ld[y] > d) {
                int l = bl.size();
                push(y, 0, wof(y));
                grow(l);
            }
        return t;
    }

   private:
    // 内部では頂点を「木ごとの BFS 順」に振り直した番号（ラベル）で持つ: 親のラベル < 子のラベル、子どうしのラベルは連続
    vector<int> lab, orig;      // lab[頂点] = ラベル、orig[ラベル] = 頂点
    vector<int> P, cl, cr, ld;  // ラベル x の親（根は -1）、子のラベルの範囲 [cl, cr)、重心木での深さ
    vector<D> wt;               // 親との辺の重み（重み付きのときだけ）
    vector<D> dt;               // dt[k * n + x] = ラベル x の頂点と、その重心祖先のうち深さ k のものとの距離
    bool weighted = false;
    mutable Component buf;
    mutable vector<int> bl;     // component 用: buf.vs と同じ並びのラベル
    D wof(int x) const { return weighted ? wt[x] : D(1); }

    template <class Graph>
    void build(const Graph& G) {
        using centroid_internal::has_weight, centroid_internal::to, centroid_internal::weight;
        using E = centroid_internal::edge_t<Graph>;
        n = G.size();
        weighted = has_weight<E>;
        lab.assign(n, -1), orig.reserve(n), sz.assign(n, 0);
        P.assign(n, -1), cl.assign(n, 0), cr.assign(n, 0), ld.assign(n, -1);
        if (weighted) wt.assign(n, D(0));
        vector<array<int, 4>> tasks;  // {ord 上の区間 L, R, 重心木での親（ラベル）, 重心木での深さ}。先頭から順に処理する
        // 1. 木ごとに BFS 順のラベルを振る。ついでに木であること（閉路・多重辺・自己ループがない）を確かめる
        for (int s = 0; s < n; s++) {
            if (lab[s] >= 0) continue;
            int start = orig.size();
            lab[s] = start, orig.push_back(s);
            size_t deg = 0;
            for (int x = start; x < (int)orig.size(); x++) {
                int v = orig[x];
                deg += G[v].size();
                cl[x] = orig.size();
                for (auto& e : G[v]) {
                    int u = to(e);
                    if (lab[u] >= 0) continue;
                    int y = orig.size();
                    lab[u] = y, orig.push_back(u), P[y] = x;
                    if constexpr (has_weight<E>) wt[y] = (D)weight(e);
                }
                cr[x] = orig.size();
            }
            assert(deg == 2 * (orig.size() - start - 1));
            tasks.push_back({start, (int)orig.size(), -1, 0});
        }
        // 2. 分解。成分はいつも「ラベルの昇順に並んだ ord の区間」で持つ（先頭 = 成分内で一番上の頂点、他の頂点の親も成分内）
        //    → 部分木の大きさも重心からの距離も、隣接リストなしで区間をなめるだけで求まる
        vector<int> ord(n), tmp(n), ssz(n), mx(n), cid(n), cnt, cpar(n, -1);
        vector<char> anc(n);
        iota(ord.begin(), ord.end(), 0);
        order.reserve(n);
        dt.reserve((size_t)bit_width((unsigned)n) * n);  // 重心木の深さは log2(n) 以下。先に確保しておき、伸ばすたびの再確保（メモリの山）を避ける
        for (size_t ti = 0; ti < tasks.size(); ti++) {
            auto [L, R, cp, d] = tasks[ti];
            const int tot = R - L, m = ord[L];
            for (int i = L; i < R; i++) ssz[ord[i]] = 1, mx[ord[i]] = 0;
            for (int i = R - 1; i > L; i--) {
                int x = ord[i], p = P[x];
                ssz[p] += ssz[x];
                mx[p] = max(mx[p], ssz[x]);
            }
            int c = -1;
            for (int i = L; i < R && c < 0; i++)
                if (int x = ord[i]; 2 * max(mx[x], tot - ssz[x]) <= tot) c = x;
            ld[c] = d, cpar[c] = cp, order.push_back(c), sz[orig[c]] = tot;
            // 重心 c からの距離（深さ d の段）と、c を取り除いたあとの成分番号 cid を決める
            if (dt.size() < (size_t)(d + 1) * n) dt.resize((size_t)(d + 1) * n);
            D* row = dt.data() + (size_t)d * n;
            row[c] = 0;
            for (int x = c; x != m; x = P[x]) row[P[x]] = row[x] + wof(x), anc[P[x]] = 1;  // c の祖先たち
            int k = (c != m);  // 成分の個数。c より上の側があればそれが 0 番
            for (int i = L; i < R; i++) {
                int x = ord[i];
                if (x == c) continue;
                if (anc[x]) {
                    cid[x] = 0;
                    continue;
                }
                int p = P[x];
                row[x] = row[p] + wof(x);
                cid[x] = (p == c) ? k++ : cid[p];
            }
            for (int x = c; x != m; x = P[x]) anc[P[x]] = 0;
            // ord[L..R) を「c、成分 0、成分 1、…」の順に並べ替える（成分の中はラベルの昇順のまま）
            cnt.assign(k + 1, 0);
            for (int i = L; i < R; i++)
                if (ord[i] != c) cnt[cid[ord[i]] + 1]++;
            for (int g = 0; g < k; g++) cnt[g + 1] += cnt[g];
            for (int i = L; i < R; i++)
                if (int x = ord[i]; x != c) tmp[L + 1 + cnt[cid[x]]++] = x;
            ord[L] = c;
            copy(tmp.begin() + L + 1, tmp.begin() + R, ord.begin() + L + 1);
            for (int g = 0, s = L + 1; g < k; g++) {
                int e = L + 1 + cnt[g];
                tasks.push_back({s, e, c, d + 1});
                s = e;
            }
        }
        // 3. 頂点番号に戻す
        par.assign(n, -1), dep.assign(n, 0);
        for (int x = 0; x < n; x++) {
            dep[orig[x]] = ld[x];
            if (cpar[x] >= 0) par[orig[x]] = orig[cpar[x]];
        }
        for (int& x : order) x = orig[x];
    }
};
template <class Graph>
CentroidDecomposition(const Graph&) -> CentroidDecomposition<centroid_internal::dist_t<Graph>>;
$0
