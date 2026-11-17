// @name Zeta / Mobius transform
// @prefix zeta
// @prefix mobius
// @prefix sos
// @description 高速ゼータ・メビウス変換（部分集合 / 上位集合 / 約数 / 倍数 / b 進数・多次元）。演算は max・xor などに差し替え可。畳み込みは zetaconv

// ==== 高速ゼータ変換・メビウス変換 ====
// どれも f を上書きする。f は vector / array / 生配列 / span など（f[i] と std::size(f) が使えるもの）
//
//   何を集めるか            ゼータ（集める）        メビウス（ゼータの逆変換）
//   部分集合 T ⊆ S          subset_zeta(f)          subset_mobius(f)          g[S] = Σ_{T ⊆ S} f[T]
//   上位集合 T ⊇ S          superset_zeta(f)        superset_mobius(f)        g[S] = Σ_{T ⊇ S} f[T]
//   約数     d | n          divisor_zeta(f)         divisor_mobius(f)         g[n] = Σ_{d | n} f[d]
//   倍数     n | m          multiple_zeta(f)        multiple_mobius(f)        g[n] = Σ_{n | m} f[m]
//   全桁が以下 y ≤ x        prefix_zeta(f, b)       prefix_mobius(f, b)       添字を b 進数で見て、全部の桁で y_k ≤ x_k
//   全桁が以上 y ≥ x        suffix_zeta(f, b)       suffix_mobius(f, b)       （b = 2 なら subset / superset と同じ）
//
// ・長さは 2 冪（b 冪）でなくてよい。長さ n なら添字 n 未満だけを対象にした和になる
// ・divisor / multiple は f[i] が整数 i に対応する（N まで扱うなら長さ N + 1）。f[0] は使わない・触らない
// ・prefix / suffix は b の代わりに {d0, d1, ...}（vector<int> でも可）を渡すと桁ごとに長さを変えられる（多次元の累積和）
//     添字 i = x0 + d0 * (x1 + d1 * (x2 + ...))、0 ≤ x_k < d_k。d の積が長さと一致すること
//     H × W の表を a[i * W + j] で持っているなら {W, H}（下の桁が先）
// ・演算の差し替え: 最後の引数に op を渡す。op は可換・結合的であること
//     a = op(a, b) になる関数（セグ木用の op_max などもそのまま渡せる）か、第 1 引数を T& / auto& で受けて直接書き換える void の関数
//     メビウスには逆演算（引き算にあたるもの）を渡す。max のように逆がない演算はゼータだけ
// ・配列の一部（連続する区間）だけにかけるなら span を渡す: subset_zeta(span(a, n))、subset_zeta(span(f).subspan(l, len))
//
// 例
//   subset_zeta(f, [](ll a, ll b) { return max(a, b); });                  // 部分集合の中の最大値
//   subset_zeta(f, [](auto& a, const auto& b) { ... });                    // a を直接更新（上位 2 個の pair の併合など）
//   subset_zeta(f, bit_xor<>()); subset_mobius(f, bit_xor<>());            // xor は自分自身が逆演算
//   prefix_zeta(f, 3);                                                     // 3 進数の各桁が以下のものの和
//   vector<int> mu(n + 1); mu[1] = 1; divisor_mobius(mu);                  // メビウス関数 μ(1..n)
//   vector<ll> phi(n + 1); iota(all(phi), 0); divisor_mobius(phi);         // オイラーの φ(1..n)
//   multiple_zeta(c);  rep1(d, n) c[d] = c[d] * (c[d] - 1) / 2;  multiple_mobius(c);
//                                  // c[x] = x の個数 → d の倍数の個数 → 両方 d の倍数の組 → gcd がちょうど d の組の数
//   N の約数だけの上でやる（N が大きくても約数の個数ぶんで済む）:
//     vector<int> dims; vector<ll> ds{1};                                  // ds[i] = 添字 i が表す約数
//     for (auto [p, e] : factorize(N)) { dims.push_back(e + 1); int m = ds.size(); rep(i, e * m) ds.push_back(ds[i] * p); }
//     prefix_zeta(f, dims);                                                // f[i] = Σ_{ds[j] | ds[i]} f[j]（suffix なら倍数側）
//
// 計算量: subset / superset O(n log n)、prefix / suffix O(n × 桁数)、divisor / multiple O(n log log n)
// 目安:   subset_zeta は n = 2^20 で 5〜7ms・2^24 で 0.15〜0.2 秒、divisor_zeta は N = 10^7 で 0.1〜0.15 秒（初回だけ素数の篩で +0.05 秒）
namespace zeta_internal {
struct Add { template <class T> void operator()(T& a, const T& b) const { a += b; } };
struct Sub { template <class T> void operator()(T& a, const T& b) const { a -= b; } };
// a に b を取り込む。op が値を返すなら a = op(a, b)、void なら op(a, b) が a を直接書き換える
template <class Op, class T>
inline void upd(Op& op, T& a, const T& b) {
    using R = std::invoke_result_t<Op&, T&, const T&>;
    if constexpr (std::is_void_v<R>) {
        static_assert(!std::is_invocable_v<Op&, T, const T&>, "void を返す op は第 1 引数を T& か auto& で受けて書き換えること（auto&& は不可）。値を返すつもりなら return の書き忘れ");
        op(a, b);
    } else {
        static_assert(std::is_same_v<T, bool> || !std::is_same_v<std::decay_t<R>, bool>, "op が bool を返している（chmin / chmax の戻り値を return していないか）");
        a = op(a, b);
    }
}
// f[0, n) を長さ blk ごとに区切り、各区間の中で j と j - s を結ぶ（s = この軸の歩幅、blk = 1 つ上の軸の歩幅）
//   UP  : true なら f[j] ← f[j - s]（下から上へ集める）、false なら f[j - s] ← f[j]
//   ASC : j を昇順に回すか。ゼータは集める向きに、メビウスはその逆に回す（軸の長さが 3 以上のときに効く）
template <bool UP, bool ASC, class It, class Op>
void axis(It f, int n, int s, int blk, Op& op) {
    for (int b = 0; b < n; b += blk) {
        const int lo = b + s, hi = std::min(b + blk, n);
        if (ASC) {
            for (int j = lo; j < hi; j++) {
                if (UP) upd(op, f[j], f[j - s]);
                else upd(op, f[j - s], f[j]);
            }
        } else {
            for (int j = hi - 1; j >= lo; j--) {
                if (UP) upd(op, f[j], f[j - s]);
                else upd(op, f[j - s], f[j]);
            }
        }
    }
}
// f[0, n) に軸 0 .. k-1 の変換をかける。st[a] = 軸 a の歩幅（st[0] = 1 < st[1] < ...、最後は n 以上）
// 軸の順番は自由なので、大きい区間は一番上の軸で切り分けて下の軸を先に済ませる（下の軸の処理がキャッシュに収まる）
// 長さが LEAF 以下になったら、軸ごとに全体をなめる普通のループで処理する
constexpr int LEAF = 1 << 10;
template <bool UP, bool ASC, class It, class Op>
void rec(It f, int n, const std::vector<int>& st, int k, Op& op) {
    while (k > 0 && st[k - 1] >= n) k--;
    if (n <= LEAF) {
        for (int a = 0; a < k; a++) axis<UP, ASC>(f, n, st[a], st[a + 1], op);
        return;
    }
    const int s = st[k - 1], step = s < LEAF ? LEAF / s * s : s;  // 切り分けた 1 個が小さいときは LEAF 以下でまとめる
    for (int i = 0; i < n; i += step) rec<UP, ASC>(f + i, std::min(step, n - i), st, k - 1, op);
    axis<UP, ASC>(f, n, s, st[k], op);
}
// 歩幅の列を作る。dims が整数ならその進数（全桁同じ長さ）、列なら桁ごとの長さ
template <class D>
std::vector<int> strides(int n, const D& dims) {
    std::vector<int> st{1};
    auto push = [&](long long d) { st.push_back((int)std::min<long long>(st.back() * d, 1 << 30)); };
    if constexpr (std::is_integral_v<D>) {
        assert(dims >= 2);
        while (st.back() < n) push(dims);
    } else {
        [[maybe_unused]] long long prod = 1;
        for (long long d : dims) {
            assert(d >= 0);
            d = std::min(d, 1LL << 31);
            prod = std::min(prod * d, 1LL << 31);
            if (d >= 2 && st.back() < n) push(d);  // 長さ 1 の桁は変換に関係ない
        }
        assert(prod == n);  // dims の積が配列の長さと合っていない
    }
    return st;
}
template <bool UP, bool ASC, class V, class D, class Op>
void grid(V&& f, const D& dims, Op& op) {
    assert(std::size(f) <= (size_t(1) << 30));
    const int n = (int)std::size(f);
    const std::vector<int> st = strides(n, dims);
    rec<UP, ASC>(std::begin(f), n, st, (int)st.size() - 1, op);
}
// 素数の昇順（n 以下は必ず全部入っている）。一度作ったら使い回し、足りなくなったら作り直す
inline const std::vector<int>& primes(int n) {
    static std::vector<int> ps;
    static int lim = 1;
    if (n > lim) {
        lim = n;
        ps.assign(1, 2);
        std::vector<char> comp(lim / 2 + 1);  // comp[i] : 2i+1 が合成数
        for (int i = 1; 2 * i + 1 <= lim; i++) {
            if (comp[i]) continue;
            const int p = 2 * i + 1;
            ps.push_back(p);
            for (long long j = (long long)p * p / 2; j <= lim / 2; j += p) comp[j] = 1;
        }
    }
    return ps;
}
// 素数 p ごとに k と kp を結ぶ（素因数ごとの指数を軸とみた多次元の変換。UP / ASC の意味は axis と同じ）
template <bool UP, bool ASC, class V, class Op>
void divs(V&& f, Op& op) {
    const int n = (int)std::size(f) - 1;
    for (int p : primes(n)) {
        if (p > n) break;
        if (ASC) {
            for (int k = 1, j = p; j <= n; k++, j += p) {
                if (UP) upd(op, f[j], f[k]);
                else upd(op, f[k], f[j]);
            }
        } else {
            for (int k = n / p, j = k * p; k >= 1; k--, j -= p) {
                if (UP) upd(op, f[j], f[k]);
                else upd(op, f[k], f[j]);
            }
        }
    }
}
}  // namespace zeta_internal

// grid / divs の <UP, ASC> は、下から集めるゼータ = <true, true>、その逆変換 = <true, false>、
// 上から集めるゼータ = <false, false>、その逆変換 = <false, true>（長さ 2 の軸は回す向きを問わないので全部昇順にしている）
template <class V, class Op = zeta_internal::Add>
void subset_zeta(V&& f, Op op = Op()) { zeta_internal::grid<true, true>(f, 2, op); }
template <class V, class Op = zeta_internal::Sub>
void subset_mobius(V&& f, Op op = Op()) { zeta_internal::grid<true, true>(f, 2, op); }
template <class V, class Op = zeta_internal::Add>
void superset_zeta(V&& f, Op op = Op()) { zeta_internal::grid<false, true>(f, 2, op); }
template <class V, class Op = zeta_internal::Sub>
void superset_mobius(V&& f, Op op = Op()) { zeta_internal::grid<false, true>(f, 2, op); }

template <class V, class Op = zeta_internal::Add>
void divisor_zeta(V&& f, Op op = Op()) { zeta_internal::divs<true, true>(f, op); }
template <class V, class Op = zeta_internal::Sub>
void divisor_mobius(V&& f, Op op = Op()) { zeta_internal::divs<true, false>(f, op); }
template <class V, class Op = zeta_internal::Add>
void multiple_zeta(V&& f, Op op = Op()) { zeta_internal::divs<false, false>(f, op); }
template <class V, class Op = zeta_internal::Sub>
void multiple_mobius(V&& f, Op op = Op()) { zeta_internal::divs<false, true>(f, op); }

template <class V, class D, class Op = zeta_internal::Add>
void prefix_zeta(V&& f, const D& dims, Op op = Op()) { zeta_internal::grid<true, true>(f, dims, op); }
template <class V, class D, class Op = zeta_internal::Sub>
void prefix_mobius(V&& f, const D& dims, Op op = Op()) { zeta_internal::grid<true, false>(f, dims, op); }
template <class V, class D, class Op = zeta_internal::Add>
void suffix_zeta(V&& f, const D& dims, Op op = Op()) { zeta_internal::grid<false, false>(f, dims, op); }
template <class V, class D, class Op = zeta_internal::Sub>
void suffix_mobius(V&& f, const D& dims, Op op = Op()) { zeta_internal::grid<false, true>(f, dims, op); }
// {d0, d1, ...} と直接書いたとき用
template <class V, class Op = zeta_internal::Add>
void prefix_zeta(V&& f, std::initializer_list<int> dims, Op op = Op()) { zeta_internal::grid<true, true>(f, dims, op); }
template <class V, class Op = zeta_internal::Sub>
void prefix_mobius(V&& f, std::initializer_list<int> dims, Op op = Op()) { zeta_internal::grid<true, false>(f, dims, op); }
template <class V, class Op = zeta_internal::Add>
void suffix_zeta(V&& f, std::initializer_list<int> dims, Op op = Op()) { zeta_internal::grid<false, false>(f, dims, op); }
template <class V, class Op = zeta_internal::Sub>
void suffix_mobius(V&& f, std::initializer_list<int> dims, Op op = Op()) { zeta_internal::grid<false, true>(f, dims, op); }
