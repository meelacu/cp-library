// @name CumSum (n-dimensional prefix sum / imos)
// @prefix cumsum
// @description n 次元累積和 CumSum<T, D>：直方体の和 O(2^D)、1 点加算といもす（直方体への加算）を build なしで混ぜられる。XOR などの可換群にも対応

// ==== n 次元累積和・いもす ====
// CumSum<T, D = 1, Add = std::plus<>, Sub = std::minus<>> : D 次元。添字は 0-indexed、区間は半開
// 作る
//   CumSum S(a)                         : 入れ子の vector から（D は入れ子の深さ。int などの整数は long long で持つ）
//   CumSum<T, D> S(n0, ..., n_{D-1})    : 全部 0
//   CumSum<T, D> S(n0, ..., n_{D-1}, f) : a[i0][i1].. = f(i0, i1, ..) で作る。例: [&](int i, int j) { return s[i][j] == '#'; }
// 更新（何度でも・どの順でも。build は不要）
//   S.add(i0, .., v)                     : a[i] += v
//   S.imos(l0, .., l_{D-1}, r0, .., r_{D-1}, v) : 直方体 [l, r) の全要素に += v（いもす）
// 参照
//   S.sum(l0, .., l_{D-1}, r0, .., r_{D-1}) : 直方体 [l, r)（各軸 k で l_k ≤ i_k < r_k）の和。S.sum() で全体の和
//   S.get(i0, ..)                        : a[i]
//   例 2 次元: S.sum(x1, y1, x2, y2) は [x1, x2) × [y1, y2) の和。添字は配列でも渡せる: S.sum({x1, y1}, {x2, y2})
// 計算量: 構築 O(N)（N は全要素数）、各操作 O(2^D)
//   ただし更新のあと最初の参照と、参照のあと最初の更新で O(N D) の作り直しが入る
//   → 更新をまとめてから参照する使い方向け。更新と参照が交互に大量に来るなら BIT / セグ木を使う
// 演算は可換群なら差し替えられる（T() が単位元、Sub(a, b) = a - b）。XOR: CumSum<ull, D, std::bit_xor<>, std::bit_xor<>>
template <class T, int D = 1, class Add = std::plus<>, class Sub = std::minus<>>
struct CumSum {
    static_assert(D >= 1);
    using Idx = std::array<int, D>;
    Idx n{};
    std::array<size_t, D> st{};  // 各軸に 0 番の余白を足した配列での歩幅
    mutable std::vector<T> d;    // a[i] は d[(i0 + 1, i1 + 1, ..)] に置く
    mutable int lv = 0;          // d は a を lv 回累積したもの（-1 は差分）
    mutable bool z = true;       // d は全部 0（どの lv ともみなせる）

    explicit CumSum(const Idx& sz) { init(sz); }
    template <class... A>
    explicit CumSum(const A&... a) {
        static_assert(sizeof...(A) == D || sizeof...(A) == D + 1);
        auto t = std::forward_as_tuple(a...);
        init(take<0, D>(t));
        if constexpr (sizeof...(A) == D + 1) {
            const auto& f = std::get<D>(t);
            size_t cells = 1;
            for (int k = 0; k < D; k++) cells *= n[k];
            Idx i{};
            for (size_t c = 0; c < cells; c++) {
                d[at1(i)] = T(std::apply(f, i));
                for (int k = D - 1; k >= 0 && ++i[k] == n[k]; k--) i[k] = 0;
            }
            z = false;
        }
    }
    template <class V>
    explicit CumSum(const std::vector<V>& a) {
        Idx sz{};
        dims<0>(a, sz);
        init(sz);
        load<0>(a, 0);
        z = false;
    }

    void add(const Idx& i, const T& v) {
        in_range(i);
        if (lv < 0) return imos(i, plus1(i), v);
        to(0);
        T& x = d[at1(i)];
        x = Add{}(x, v);
    }
    void imos(const Idx& l, const Idx& r, const T& v) {
        in_range(l, r);
        to(-1);
        for (int m = 0; m < (1 << D); m++) {  // m の k ビット目が 1 なら軸 k は r 側（符号が反転）
            size_t p = 0;
            bool ok = true;
            for (int k = 0; k < D && ok; k++) {
                int x = (m >> k & 1) ? r[k] : l[k];
                ok = x < n[k];  // 右端を越える角は配列の外なので捨てる
                p += (x + 1) * st[k];
            }
            if (ok) d[p] = (__builtin_popcount(m) & 1) ? Sub{}(d[p], v) : Add{}(d[p], v);
        }
    }
    T get(const Idx& i) const {
        in_range(i);
        if (lv < 0) to(0);
        return lv == 0 ? d[at1(i)] : sum(i, plus1(i));
    }
    T sum(const Idx& l, const Idx& r) const {
        in_range(l, r);
        to(1);
        T pos = T(), neg = T();
        for (int m = 0; m < (1 << D); m++) {  // m の k ビット目が 1 なら軸 k は r 側。l 側の個数が奇数なら引く
            size_t p = 0;
            for (int k = 0; k < D; k++) p += ((m >> k & 1) ? r[k] : l[k]) * st[k];
            if ((D - __builtin_popcount(m)) & 1) neg = Add{}(neg, d[p]);
            else pos = Add{}(pos, d[p]);
        }
        return Sub{}(pos, neg);
    }

    template <class... A> void add(const A&... a) {
        static_assert(sizeof...(A) == D + 1);
        auto t = std::forward_as_tuple(a...);
        add(take<0, D>(t), T(std::get<D>(t)));
    }
    template <class... A> void imos(const A&... a) {
        static_assert(sizeof...(A) == 2 * D + 1);
        auto t = std::forward_as_tuple(a...);
        imos(take<0, D>(t), take<D, D>(t), T(std::get<2 * D>(t)));
    }
    template <class... I> T get(I... i) const {
        static_assert(sizeof...(I) == D);
        return get(Idx{int(i)...});
    }
    template <class... I> T sum(I... lr) const {
        static_assert(sizeof...(I) == 0 || sizeof...(I) == 2 * D);
        if constexpr (sizeof...(I) == 0) return sum(Idx{}, n);
        else {
            auto t = std::forward_as_tuple(lr...);
            return sum(take<0, D>(t), take<D, D>(t));
        }
    }

  private:
    void init(const Idx& sz) {
        n = sz;
        size_t tot = 1;
        for (int k = D - 1; k >= 0; k--) {
            assert(n[k] >= 0);
            st[k] = tot;
            tot *= n[k] + 1;
        }
        d.assign(tot, T());
    }
    template <int K, class V> void dims(const V& a, Idx& sz) {
        if constexpr (K < D) {
            sz[K] = (int)a.size();
            if (!a.empty()) dims<K + 1>(a[0], sz);
        }
    }
    template <int K, class V> void load(const V& a, size_t p) {
        if constexpr (K == D) d[p] = T(a);
        else {
            assert((int)a.size() == n[K]);  // 長方形であること
            for (int i = 0; i < n[K]; i++) load<K + 1>(a[i], p + (i + 1) * st[K]);
        }
    }
    template <size_t O, class Tup, size_t... k> static Idx take_(const Tup& t, std::index_sequence<k...>) { return Idx{int(std::get<O + k>(t))...}; }
    template <size_t O, size_t K, class Tup> static Idx take(const Tup& t) { return take_<O>(t, std::make_index_sequence<K>()); }
    static Idx plus1(Idx i) { for (auto& x : i) x++; return i; }
    size_t at1(const Idx& i) const { size_t p = 0; for (int k = 0; k < D; k++) p += (i[k] + 1) * st[k]; return p; }
    void in_range(const Idx& i) const { for (int k = 0; k < D; k++) assert(0 <= i[k] && i[k] < n[k]); }
    void in_range(const Idx& l, const Idx& r) const { for (int k = 0; k < D; k++) assert(0 <= l[k] && l[k] <= r[k] && r[k] <= n[k]); }
    void up() const {  // 各軸で累積
        for (int k = 0; k < D; k++) {
            size_t s = st[k], blk = s * (n[k] + 1);
            for (size_t b = 0; b < d.size(); b += blk)
                for (size_t i = b + s; i < b + blk; i++) d[i] = Add{}(d[i], d[i - s]);
        }
    }
    void down() const {  // up の逆（各軸で差分）
        for (int k = 0; k < D; k++) {
            size_t s = st[k], blk = s * (n[k] + 1);
            for (size_t b = 0; b < d.size(); b += blk)
                for (size_t i = b + blk; i-- > b + s;) d[i] = Sub{}(d[i], d[i - s]);
        }
    }
    void to(int k) const {
        if (z) { z = false, lv = k; return; }
        for (; lv < k; lv++) up();
        for (; lv > k; lv--) down();
    }
};
template <class V> struct cumsum_elem { using type = V; static constexpr int dim = 0; };
template <class V> struct cumsum_elem<std::vector<V>> { using type = typename cumsum_elem<V>::type; static constexpr int dim = cumsum_elem<V>::dim + 1; };
template <class V> CumSum(const std::vector<V>&) -> CumSum<
    std::conditional_t<std::is_integral_v<typename cumsum_elem<V>::type> && sizeof(typename cumsum_elem<V>::type) < 8, long long, typename cumsum_elem<V>::type>,
    cumsum_elem<std::vector<V>>::dim>;
