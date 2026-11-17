// @name Zeta convolution (OR / AND / GCD / LCM / subset)
// @prefix zetaconv
// @prefix orconv
// @prefix andconv
// @prefix gcdconv
// @prefix lcmconv
// @prefix subsetconv
// @description OR・AND・GCD・LCM・subset 畳み込み（ゼータ変換 → 各点の積 → メビウス変換）。要 zeta スニペット

// ==== ゼータ変換でできる畳み込み ====（zeta スニペットを先に貼ること）
//   or_convolution(a, b)       c[k] = Σ_{i | j = k} a[i] * b[j]
//   and_convolution(a, b)      c[k] = Σ_{i & j = k} a[i] * b[j]
//   gcd_convolution(a, b)      c[k] = Σ_{gcd(i, j) = k} a[i] * b[j]     添字 0 は使わない（c[0] = 0）
//   lcm_convolution(a, b)      c[k] = Σ_{lcm(i, j) = k} a[i] * b[j]     添字 0 は使わない。lcm が長さ以上になる組は捨てる
//   subset_convolution(a, b)   c[S] = Σ_{T ⊆ S} a[T] * b[S \ T]         i & j = 0 の組だけを集める OR 畳み込み
// ・a, b は vector<T>、新しい vector<T> を返す。T は + - * のある型（ll, mint, double など）で、T() が 0 であること
// ・a と b の長さは違ってもよい。結果は長い方の長さ（or / subset は 2 冪に切り上げ）
// ・同じ配列どうし（or_convolution(a, a) など）は変換を 1 回省く
// 計算量: or / and O(n log n)、gcd / lcm O(n log log n)、subset O(n log^2 n)
//         subset_convolution は作業用に T を n (log n + 2) 個使う（n = 2^20 の mint で 90MB）
// 目安:   n = 2^20 で or / and 15〜70ms、subset 0.3〜0.7 秒。N = 10^6 で gcd / lcm 15〜30ms
namespace zeta_internal {
using zeta_snippet_required = Add;  // ここでエラーになるなら zeta スニペットが貼られていない（先に貼る）
// ゼータ → 各点の積 → メビウス。a, b は長さ n にそろえる
template <class T, class Z, class M>
std::vector<T> conv(const std::vector<T>& a, const std::vector<T>& b, size_t n, Z zeta, M mobius) {
    std::vector<T> fa(a);
    fa.resize(n), zeta(fa);
    if (&a == &b) {
        for (T& x : fa) x *= T(x);
    } else {
        std::vector<T> fb(b);
        fb.resize(n), zeta(fb);
        for (size_t i = 0; i < n; i++) fa[i] *= fb[i];
    }
    mobius(fa);
    return fa;
}
// ランク付きゼータ（INV = false）/ メビウス（INV = true）。f = 行 off の先頭、行 off .. off + n - 1 を処理する
// 行 s は長さ |s| + 1 で、行を詰めて並べる。ランク r = 「|T| = r の T ⊆ s についての和」
//   ゼータ中 : 行 s の x[r] がランク r（0 ≤ r ≤ |s|。それより上は 0）
//   積のあと : 行 s の x[t] がランク |s| + t（答えに要るのは |s| 以上だけで、2|s| より上は 0）
template <bool INV, class T>
void ranked(T* f, int off, int n, int lg) {
    if (n == 1) return;
    const int h = n / 2;
    T* g = f + (size_t)h * (__builtin_popcount(off) + 1) + (size_t)(h / 2) * __builtin_ctz(h);  // 行 off + h の先頭
    ranked<INV>(f, off, h, lg), ranked<INV>(g, off + h, h, lg);
    T *lo = f, *hi = g;  // 行 off + j と、その 1 ビット上の行 off + h + j
    for (int j = 0; j < h; j++) {
        const int pc = __builtin_popcount(off + j);
        if (!INV) for (int r = 0; r <= pc; r++) hi[r] += lo[r];
        else for (int t = 0, e = std::min(pc, lg - pc); t < e; t++) hi[t] -= lo[t + 1];
        lo += pc + 1, hi += pc + 2;
    }
}
}  // namespace zeta_internal

template <class T>
std::vector<T> or_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    const size_t m = std::max(a.size(), b.size());
    size_t n = m ? 1 : 0;
    while (n < m) n <<= 1;
    return zeta_internal::conv(a, b, n, [](auto& f) { subset_zeta(f); }, [](auto& f) { subset_mobius(f); });
}
template <class T>
std::vector<T> and_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    return zeta_internal::conv(a, b, std::max(a.size(), b.size()), [](auto& f) { superset_zeta(f); }, [](auto& f) { superset_mobius(f); });
}
template <class T>
std::vector<T> gcd_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    auto c = zeta_internal::conv(a, b, std::max(a.size(), b.size()), [](auto& f) { multiple_zeta(f); }, [](auto& f) { multiple_mobius(f); });
    if (!c.empty()) c[0] = T();
    return c;
}
template <class T>
std::vector<T> lcm_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    auto c = zeta_internal::conv(a, b, std::max(a.size(), b.size()), [](auto& f) { divisor_zeta(f); }, [](auto& f) { divisor_mobius(f); });
    if (!c.empty()) c[0] = T();
    return c;
}
template <class T>
std::vector<T> subset_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    const size_t m = std::max(a.size(), b.size());
    if (m == 0) return {};
    int lg = 0;
    while ((size_t(1) << lg) < m) lg++;
    const int n = 1 << lg;
    auto lift = [&](const std::vector<T>& v) {  // ランクごとに分けてゼータ変換
        std::vector<T> r((size_t)n + (size_t)(n / 2) * lg);
        T* x = r.data();
        for (int s = 0; s < n; s++) {
            const int pc = __builtin_popcount(s);
            if (s < (int)v.size()) x[pc] = v[s];
            x += pc + 1;
        }
        zeta_internal::ranked<false>(r.data(), 0, n, lg);
        return r;
    };
    const bool sq = &a == &b;
    std::vector<T> fa = lift(a), fb;
    if (!sq) fb = lift(b);
    T* x = fa.data();
    const T* y = sq ? x : fb.data();
    for (int s = 0; s < n; s++) {  // 行ごとに多項式として掛ける。x[r - pc] に書いても、あとで読むのはそれより後ろだけ
        const int pc = __builtin_popcount(s);
        for (int r = pc, e = std::min(2 * pc, lg); r <= e; r++) {
            T v = T();
            for (int p = r - pc; p <= pc; p++) v += x[p] * y[r - p];
            x[r - pc] = v;
        }
        x += pc + 1, y += pc + 1;
    }
    zeta_internal::ranked<true>(fa.data(), 0, n, lg);
    std::vector<T> c(n);
    x = fa.data();
    for (int s = 0; s < n; s++) c[s] = *x, x += __builtin_popcount(s) + 1;
    return c;
}
