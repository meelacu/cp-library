// @name Doubling
// @prefix doubling
// @description ダブリング: k歩先(k<=1e18)・通った値のモノイド積・条件を満たす間進める最大歩数、全頂点を一斉にk歩。行き止まり(-1)も可

// ダブリング（各頂点から行き先が 1 つ決まっている遷移 nxt を、何歩も進める）。main の外に置く
//   Doubling<> db(nxt, max_k);              行き先だけ。nxt[v] = v の 1 歩先。-1 は行き止まり（次がない）
//   Doubling<S, op, e> db(nxt, w, max_k);   1 歩ごとの値つき。w[v] = 「v → nxt[v]」の 1 歩の値、(S, op, e) はモノイド
//       max_k : 進める歩数の最大値（1e18 まで可）。表は log2(max_k) 段 × (n + 1) 個（n = 2e5, max_k = 1e18 で int 48MB + S の分）
//               2^(段数 - 1) 歩分の積まで作るので、和などが溢れないよう必要以上に大きくしない
//
//   int  jump(v, k)                 v から k 歩先の頂点。途中で行き止まりより先へ出るなら -1          O(log k)
//   pair<int, S> prod(v, k)         {k 歩先の頂点, 通った k 歩の値の積（進んだ順）}                     O(log k)
//                                   行き止まりより先へ出るなら {-1, 行き止まりまでの積}
//   ll   max_steps(v, f)            f が true のまま進める最大の歩数（max_k 以下）                        O(log max_k)
//   ll   max_steps(v, f, limit)     同上で limit 歩まで
//        f(to, acc) : to = 進んだ先の頂点、acc = そこまでの積（値なしなら f(to)）
//        歩数について単調（ある歩数までは true、それ以降は false）で、0 歩（f(v, e())）は true であること
//        行き止まりより先へは進まない。ずっと true なら、行き止まりまでの歩数（なければ max_k / limit）を返す
//   static vector<int> jump_all(nxt, k)              全頂点を一斉に k 歩進めた先。表を作らない（メモリ O(n)、時間 O(n log k)）
//   static pair<vector<int>, vector<S>> prod_all(nxt, w, k)    同上で、通った値の積も
//        呼び方: Doubling<>::jump_all(nxt, k) / Doubling<S, op, e>::prod_all(nxt, w, k)
// 例: 高橋君とホテル（ARC060 E）: nxt[i] = 1 日で行ける最も遠いホテル（右向き）、a < b の答え = db.max_steps(a, [&](int to) { return to < b; }) + 1
//     木の親 par（根は -1）を渡すと、jump(v, k) = k 個上の祖先（なければ -1）、max_steps(v, [](int) { return true; }) = v の深さ
namespace doubling_internal {
struct None {};
inline None none_op(None, None) { return {}; }
inline None none_e() { return {}; }
}  // namespace doubling_internal

template <class S = doubling_internal::None, auto op = doubling_internal::none_op, auto e = doubling_internal::none_e>
struct Doubling {
    static constexpr bool AGG = !is_empty_v<S>;
    int n = 0, LOG = 1;
    long long kmax = 0;  // 構築時の max_k
    vector<int> up;  // up[j * (n + 1) + v] : v から 2^j 歩先。添字 n は「行き止まりの先」を表す番兵（n の先は n）
    vector<S> val;   // val[j * (n + 1) + v] : その 2^j 歩の値の積（値つきのときだけ。番兵は e()）

    Doubling() = default;
    Doubling(const vector<int>& nxt, long long max_k) {
        static_assert(!AGG, "値つきは Doubling(nxt, w, max_k)");
        build(nxt, {}, max_k);
    }
    Doubling(const vector<int>& nxt, const vector<S>& w, long long max_k) { build(nxt, w, max_k); }

    int jump(int v, long long k) const {
        assert(0 <= v && v < n && 0 <= k && k <= cap());
        const size_t m = n + 1;
        for (int j = 0; k; j++, k >>= 1)
            if (k & 1) v = up[j * m + v];
        return v == n ? -1 : v;
    }
    pair<int, S> prod(int v, long long k) const {
        assert(0 <= v && v < n && 0 <= k && k <= cap());
        const size_t m = n + 1;
        S acc = e();
        for (int j = 0; k; j++, k >>= 1)
            if (k & 1) {
                if constexpr (AGG) acc = op(acc, val[j * m + v]);
                v = up[j * m + v];
            }
        return {v == n ? -1 : v, acc};
    }
    template <class F>
    long long max_steps(int v, F f, long long limit) const {
        assert(0 <= v && v < n && 0 <= limit && limit <= cap());
        const size_t m = n + 1;
        long long k = 0;
        S acc = e();
        for (int j = LOG - 1; j >= 0; j--) {
            if (k + (1LL << j) > limit) continue;
            size_t id = j * m + v;
            if (up[id] == n) continue;  // 行き止まりより先へ出る
            if constexpr (AGG) {
                S nacc = op(acc, val[id]);
                if (f(up[id], nacc)) k += 1LL << j, v = up[id], acc = nacc;
            } else {
                if (f(up[id])) k += 1LL << j, v = up[id];
            }
        }
        return k;
    }
    template <class F>
    long long max_steps(int v, F f) const { return max_steps(v, f, kmax); }

    static vector<int> jump_all(const vector<int>& nxt, long long k) {
        assert(k >= 0);
        int m = nxt.size();
        vector<int> r(m), b(m + 1, m), t(m + 1, m);  // 添字 m は番兵
        for (int v = 0; v < m; v++)
            if (nxt[v] >= 0) b[v] = nxt[v];
        iota(r.begin(), r.end(), 0);
        for (; k; k >>= 1) {
            if (k & 1)
                for (int v = 0; v < m; v++) r[v] = b[r[v]];
            if (k > 1) {
                for (int v = 0; v < m; v++) t[v] = b[b[v]];
                swap(b, t);
            }
        }
        for (int& x : r)
            if (x == m) x = -1;
        return r;
    }
    static pair<vector<int>, vector<S>> prod_all(const vector<int>& nxt, const vector<S>& w, long long k) {
        assert(k >= 0 && nxt.size() == w.size());
        int m = nxt.size();
        vector<int> r(m), b(m + 1, m), t(m + 1, m);
        vector<S> rv(m, e()), bv(m + 1, e()), tv(m + 1, e());
        for (int v = 0; v < m; v++)
            if (nxt[v] >= 0) b[v] = nxt[v], bv[v] = w[v];
        iota(r.begin(), r.end(), 0);
        for (; k; k >>= 1) {
            if (k & 1)
                for (int v = 0; v < m; v++) {
                    rv[v] = op(rv[v], bv[r[v]]);
                    r[v] = b[r[v]];
                }
            if (k > 1) {
                for (int v = 0; v < m; v++) {
                    t[v] = b[b[v]];
                    tv[v] = op(bv[v], bv[b[v]]);
                }
                swap(b, t);
                swap(bv, tv);
            }
        }
        for (int& x : r)
            if (x == m) x = -1;
        return {r, rv};
    }

   private:
    long long cap() const { return LOG >= 63 ? LLONG_MAX : (1LL << LOG) - 1; }  // 表で進める歩数の上限（>= max_k）
    void build(const vector<int>& nxt, const vector<S>& w, long long max_k) {
        assert(max_k >= 0);
        kmax = max_k;
        n = nxt.size();
        LOG = max(1, (int)bit_width((unsigned long long)max_k));
        const size_t m = n + 1;
        up.assign(LOG * m, n);
        for (int v = 0; v < n; v++) {
            assert(-1 <= nxt[v] && nxt[v] < n);
            if (nxt[v] >= 0) up[v] = nxt[v];
        }
        if constexpr (AGG) {
            assert((int)w.size() == n);
            val.assign(LOG * m, e());
            for (int v = 0; v < n; v++)
                if (nxt[v] >= 0) val[v] = w[v];
        }
        for (int j = 0; j + 1 < LOG; j++) {
            const int* a = up.data() + j * m;
            int* b = up.data() + (j + 1) * m;
            for (int v = 0; v < n; v++) b[v] = a[a[v]];
            if constexpr (AGG) {
                const S* x = val.data() + j * m;
                S* y = val.data() + (j + 1) * m;
                for (int v = 0; v < n; v++) y[v] = op(x[v], x[a[v]]);
            }
        }
    }
};
$0
