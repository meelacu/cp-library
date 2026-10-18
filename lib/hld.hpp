// @name HLD
// @prefix hld
// @description 重軽分解（要 eulertour スニペット）: lca/dist/la/jump/on_path/meet/path（経路順・非可換可）/subtree_rooted/auxiliary_tree

// HLD（重軽分解）: EulerTour を継承（eulertour スニペットが必要）。森対応・非再帰
//   EulerTour の tin は heavy-first の行きがけ順なので、heavy path は tin 上の連続区間になる
//   build(root)      : EulerTour::build + head の計算 O(n)。以下のクエリは O(log n)
//   head[v]          : v を含む heavy path の最上点
//   lca(u, v)        : 別の木なら -1
//   dist(u, v)       : 辺数（別の木なら -1）
//   la(v, k)         : v の k 個上の祖先（無ければ -1）
//   jump(u, v, k)    : u → v のパスで u から k 本進んだ頂点（範囲外・別の木なら -1）
//   on_path(u, v, x) : x が u-v パス上にあるか
//   meet(a, b, c)    : 3 頂点を結ぶパスの合流点 = c を根にしたときの lca(a, b)（別の木を含めば -1）
//   path(u, v, edge) : u → v のパスを tin の半開区間 {l, r, rev} に分解し経路順に返す（O(log n) 個、u, v は同じ木）
//                      rev = true の区間は経路上で r-1 → l の向きに通る。edge = true なら LCA の位置を除く
//   subtree(v, edge) : 部分木の区間 [l, r)（edge = true なら v 自身の位置を除く）
//   subtree_rooted(v, r, edge) : 根を r にしたときの v の部分木（空でない区間を高々 2 個）
//   edge_tin(e)      : 辺 e の位置 = tin[edge_child[e]]（辺 (x, par[x]) を位置 tin[x] に置く）
//   auxiliary_tree(xs) : 補助木 {ws, p}。ws = xs（重複可）とその LCA を tin 順に並べたもの（親が先）、
//                        p[i] = ws[i] の親の添字（根は -1）。O(k log n)。辺の長さは dep[ws[i]] - dep[ws[p[i]]]
// 使い方メモ（seg は位置の昇順に積を取るセグ木、seg_rev は op の引数を逆にしたセグ木）
//   初期配列: 頂点 init[tin[v]] = a[v]、辺 init[edge_tin(e)] = w[e]（根の位置は単位元のまま）
//   可換  : for (auto [l, r, rev] : hld.path(u, v)) res = op(res, seg.prod(l, r));
//   非可換: 上の seg.prod(l, r) を rev ? seg_rev.prod(l, r) : seg.prod(l, r) にする
//   二分探索: acc（ここまでの積）を持って区間を順に見て、pred(op(acc, 区間積)) が偽になる区間で
//             f(x) = pred(op(acc, x)) とすると、初めて偽になる位置は
//             !rev なら seg.max_right(l, f)、rev なら seg_rev.min_left(r, f) - 1
struct HLD : EulerTour {
    vector<int> head;

    using EulerTour::EulerTour;

    void build(int root = 0) {
        EulerTour::build(root);
        head.resize(n);
        for (int v : order) {  // heavy 子は tin が親 + 1
            int p = par[v];
            head[v] = (p >= 0 && tin[v] == tin[p] + 1) ? head[p] : v;
        }
    }

    int lca(int u, int v) const {
        while (head[u] != head[v]) {  // tin が大きい側の head は LCA より真に下
            if (tin[u] > tin[v]) swap(u, v);
            v = par[head[v]];
            if (v < 0) return -1;
        }
        return tin[u] < tin[v] ? u : v;
    }

    int dist(int u, int v) const {
        int w = lca(u, v);
        return w < 0 ? -1 : dep[u] + dep[v] - 2 * dep[w];
    }

    int la(int v, int k) const {
        if (k < 0 || dep[v] < k) return -1;
        while (tin[v] - tin[head[v]] < k) k -= tin[v] - tin[head[v]] + 1, v = par[head[v]];
        return order[tin[v] - k];
    }

    int jump(int u, int v, int k) const {
        int w = lca(u, v);
        if (w < 0) return -1;
        int du = dep[u] - dep[w], dv = dep[v] - dep[w];
        if (k < 0 || du + dv < k) return -1;
        return k <= du ? la(u, k) : la(v, du + dv - k);
    }

    bool on_path(int u, int v, int x) const {
        int w = lca(u, v);
        return w >= 0 && is_anc(w, x) && (is_anc(x, u) || is_anc(x, v));
    }

    int meet(int a, int b, int c) const {
        int x = lca(a, b), y = lca(b, c), z = lca(c, a);
        return x < 0 || y < 0 ? -1 : x ^ y ^ z;  // 3 つのうち 2 つは一致し、残りが答え
    }

    struct Seg { int l, r; bool rev; };

    vector<Seg> path(int u, int v, bool edge = false) const {
        Seg buf[64];  // light edge は片側 log2(n) 本以下なので区間は 61 個以下
        int k = 0, j = 64;  // u 側（上り）は前から、v 側（下り）は後ろから詰める
        while (head[u] != head[v]) {
            if (tin[u] > tin[v]) buf[k++] = {tin[head[u]], tin[u] + 1, true}, u = par[head[u]];
            else buf[--j] = {tin[head[v]], tin[v] + 1, false}, v = par[head[v]];
            assert(u >= 0 && v >= 0);  // u, v は同じ木
        }
        if (tin[u] > tin[v]) buf[k++] = {tin[v] + edge, tin[u] + 1, true};
        else if (tin[u] + edge <= tin[v]) buf[k++] = {tin[u] + edge, tin[v] + 1, false};
        while (j < 64) buf[k++] = buf[j++];
        return vector<Seg>(buf, buf + k);
    }

    pair<int, int> subtree(int v, bool edge = false) const { return {tin[v] + edge, tout[v]}; }
    int edge_tin(int e) const { return tin[edge_child[e]]; }

    vector<pair<int, int>> subtree_rooted(int v, int r, bool edge = false) const {
        vector<pair<int, int>> res;
        auto add = [&](int a, int b) { if (a < b) res.emplace_back(a, b); };
        if (v != r && !is_anc(v, r)) {  // r が v の部分木の外（別の木を含む）なら普通の部分木
            add(tin[v] + edge, tout[v]);
        } else {
            int t = la(v, dep[v]);  // 木の根
            if (v == r) {
                add(tin[t] + edge, tout[t]);
            } else {
                int c = la(r, dep[r] - dep[v] - 1);  // v の子のうち r 側
                add(tin[t] + edge, tin[c]), add(tout[c], tout[t]);
            }
        }
        return res;
    }

    pair<vector<int>, vector<int>> auxiliary_tree(vector<int> xs) const {
        auto by_tin = [&](int a, int b) { return tin[a] < tin[b]; };
        sort(xs.begin(), xs.end(), by_tin);
        xs.erase(unique(xs.begin(), xs.end()), xs.end());
        for (int i = int(xs.size()) - 1; i > 0; i--) {  // tin 順で隣り合う 2 点の LCA を足せば閉じる
            int w = lca(xs[i - 1], xs[i]);
            if (w >= 0) xs.push_back(w);
        }
        sort(xs.begin(), xs.end(), by_tin);
        xs.erase(unique(xs.begin(), xs.end()), xs.end());
        vector<int> p(xs.size(), -1), st;
        for (int i = 0; i < int(xs.size()); i++) {
            while (!st.empty() && !is_anc(xs[st.back()], xs[i])) st.pop_back();
            if (!st.empty()) p[i] = st.back();
            st.push_back(i);
        }
        return {std::move(xs), std::move(p)};
    }
};
