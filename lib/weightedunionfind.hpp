// @name Weighted UnionFind
// @prefix weightedunionfind
// @description 重み付きUnionFind(ポテンシャル): 差分制約の追加と矛盾判定・diff。足し算(pair・arrayは成分ごと)・xor・mod・非可換群

// 重み付き UnionFind（ポテンシャル付き UnionFind）
//   各頂点に未知の値 w(v) があり、「w(y) = w(x) + d」の形の関係を足していく。値そのものは決めず、差だけを管理する
//   WeightedUnionFind<> uf(n);                 long long の足し算（既定）
//   WeightedUnionFind<T> uf(n);                T の足し算。mint・double のほか、pair・array は成分ごとに足す
//                                              （2 次元の相対位置なら WeightedUnionFind<pair<ll, ll>>、unite(x, y, {dx, dy})）
//   XorUnionFind<> uf(n);                      xor 版（w(y) = w(x) ^ d）。XorUnionFind<ll> なども可
//                                              unite(x, y, 1) =「x と y は違う」、unite(x, y, 0) =「同じ」で二部グラフ判定・2 色の整合性
//   WeightedUnionFind<T, op, e, inv> uf(n);    群 (T, op, e, inv) を自分で渡す（mod の掛け算・置換・アフィン変換など。非可換でもよい）
//                                              非可換のとき「w(y) = w(x) + d」は「w(y) = op(w(x), d)」と読む
//   矛盾判定は == で行う。double で誤差を許して判定したいときは、same(x, y) なら diff(x, y) を自分で比べ、違う成分のときだけ unite する
//
//   bool unite(x, y, d)       w(y) = w(x) + d を追加する。既に同じ成分なら矛盾しないかを確かめるだけ
//                             矛盾したら false を返し、何も変えない（新しく繋いだ・矛盾しなかった なら true）
//   bool unite(x, y, d, f)    実際に繋ぐ直前に f(root, child, g) を呼ぶ（成分ごとのデータの統合用）
//                             g = 繋いだ後の child のポテンシャル。child 側の頂点 v のポテンシャルは op(g, 今の pot(v)) に変わる
//                             f の中で size(root) / size(child) は繋ぐ前の値
//   T    diff(x, y)           w(y) - w(x)（op(inv(w(x)), w(y))）。x, y は同じ成分であること
//   T    pot(x)               w(x) - w(根)。成分ごとに根を 0 とした値の割り当てになる
//   int  find(x) / bool same(x, y) / int size(x) / int count()（成分数）/ vector<vector<int>> groups()
namespace weighted_uf_internal {
// 既定の群 = 足し算。+ や単項 - が無い pair・array は成分ごとに計算する
template <class T> T add(T a, T b) {
    if constexpr (requires { a + b; }) return a + b;
    else if constexpr (requires { a.first; }) return {add(a.first, b.first), add(a.second, b.second)};
    else {
        for (size_t i = 0; i < a.size(); i++) a[i] = add(a[i], b[i]);
        return a;
    }
}
template <class T> T zero() { return T(); }
template <class T> T neg(T a) {
    if constexpr (requires { -a; }) return -a;
    else if constexpr (requires { a.first; }) return {neg(a.first), neg(a.second)};
    else {
        for (auto& x : a) x = neg(x);
        return a;
    }
}
template <class T> T bxor(T a, T b) { return a ^ b; }
template <class T> T ident(T a) { return a; }
}  // namespace weighted_uf_internal

template <class T = long long, auto op = weighted_uf_internal::add<T>, auto e = weighted_uf_internal::zero<T>,
          auto inv = weighted_uf_internal::neg<T>>
struct WeightedUnionFind {
    int n, cnt;
    vector<int> p;  // 根なら -size、そうでなければ親
    vector<T> w;    // 親から見たポテンシャル: pot(v) = op(pot(親), w[v])。根では e()

    WeightedUnionFind(int n_ = 0) : n(n_), cnt(n_), p(n_, -1), w(n_, e()) {}

    // 根を返す。経路上の頂点を根の直下に付け替え、w を根から見た値にする（union by size なので深さ O(log n)）
    int find(int x) {
        assert(0 <= x && x < n);
        if (p[x] < 0) return x;
        int r = find(p[x]);
        w[x] = op(w[p[x]], w[x]);
        return p[x] = r;
    }
    T pot(int x) {
        find(x);
        return w[x];
    }
    bool same(int x, int y) { return find(x) == find(y); }
    T diff(int x, int y) {
        assert(same(x, y));
        return op(inv(pot(x)), pot(y));
    }

    template <class F>
    bool unite(int x, int y, const T& d, F f) {
        int rx = find(x), ry = find(y);
        T px = w[x], py = w[y];  // 根から見た x, y（x, y が根なら e()）
        if (rx == ry) return op(inv(px), py) == d;
        if (-p[rx] >= -p[ry]) {  // ry を rx の下に: pot(ry) = px + d - py
            w[ry] = op(op(px, d), inv(py));
            f(rx, ry, w[ry]);
            p[rx] += p[ry];
            p[ry] = rx;
        } else {                 // rx を ry の下に: pot(rx) = py - d - px
            w[rx] = op(op(py, inv(d)), inv(px));
            f(ry, rx, w[rx]);
            p[ry] += p[rx];
            p[rx] = ry;
        }
        cnt--;
        return true;
    }
    bool unite(int x, int y, const T& d) { return unite(x, y, d, [](int, int, const T&) {}); }

    int size(int x) { return -p[find(x)]; }
    int count() const { return cnt; }
    vector<vector<int>> groups() {
        vector<vector<int>> g(n);
        for (int i = 0; i < n; i++) g[find(i)].push_back(i);
        g.erase(remove_if(g.begin(), g.end(), [](const vector<int>& v) { return v.empty(); }), g.end());
        return g;
    }
};
template <class T = int>
using XorUnionFind = WeightedUnionFind<T, weighted_uf_internal::bxor<T>, weighted_uf_internal::zero<T>, weighted_uf_internal::ident<T>>;
$0
