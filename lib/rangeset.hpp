// @name RangeSet
// @prefix rangeset
// @description 区間をsetで管理（RangeSet / 値付き区間・ODT・フック対応）

// RangeSet<T, V = std::monostate> : 区間を std::set で管理する
//   半開区間 [l, r) を互いに素・左端昇順で保持し、接する同値区間は自動で併合する（正規形）
//   V 省略 → 区間の集合（被覆の有無だけ）  V 指定 → 区間ごとに値を持つ（区間代入・ODT）。V には == が必要
//   計算量: 更新は償却 O(log n)、クエリは O(log n)（列挙系は + 列挙数）
//   長さは D 型で返る（整数なら long long 以上、実数なら T）。被覆長の合計が D に収まること
//
// 更新
//   insert(l, r) / erase(l, r) : [l, r) を被覆 / 未被覆にし、新たに覆った / 消した長さを返す
//   insert(x) / erase(x)       : 1点版（整数のみ）
//   assign(l, r, v)            : [l, r) を v で上書きし、新たに覆った長さを返す
//   apply(l, r, f)             : [l, r) 内の被覆部分の値を f(v&) か f(l, r, v&) で書き換え、同値の隣接を併合
//   split(x)                   : x を区間の境界にし、左端 >= x の最初の区間を返す
//   clear()
//   insert / erase / assign / split は末尾にフック on_add, on_del、clear は on_del を渡せる（apply は非対応）
// クエリ
//   find(x)  : x を含む区間（なければ end()）
//   lower(x) : x を含む区間、なければ x より右で最初の区間。prev(lower(x)) は r <= x で最後の区間
//   contains(x)   get(x, dflt)   covers(l, r) : [l, r) がすべて覆われているか   intersects(l, r)
//   mex(x = 0) : x 以上で覆われていない最小値   prev_uncovered(x) : x 以下で覆われていない最大値（整数のみ）
//   for_each(l, r, f)     : [l, r) 内の被覆部分を f(l, r, v) で列挙（[l, r) に切り詰めて渡す）
//   for_each_gap(l, r, f) : [l, r) 内の未被覆部分を f(l, r) で列挙
//   size() : 区間数   len : 被覆長の合計   for (auto& [l, r, v] : rs) で列挙   cout << rs でデバッグ出力
//
// フック on_add(l, r, v) / on_del(l, r, v) : 内部の区間が追加 / 削除されるたびに呼ばれる
//   分割・併合も del → add として通知。1回の操作では del がすべて先。変化しない区間には通知しない
//   (l, r) の2引数でもよい。v は値か const& で受ける。フックや f の中で rs 自身を変更しないこと
//   1点版にフックを付けるときは insert(x, x + 1, on_add, on_del) と書く
//
// 使用例
//   色ごとの個数 : auto add = [&](int a, int b, int c) { cnt[c] += b - a; };  auto del = [&](int a, int b, int c) { cnt[c] -= b - a; };
//                  rs.assign(l, r, c, add, del);  同色の連結成分ごと塗る: auto it = rs.find(x); rs.assign(it->l, it->r, c, add, del);
//   最長の連続区間: multiset<ll> ms;  rs.insert(l, r, [&](ll a, ll b) { ms.insert(b - a); }, [&](ll a, ll b) { ms.erase(ms.find(b - a)); });
//   未処理の部分だけ: rs.for_each_gap(l, r, [&](ll a, ll b) { /* [a, b) を処理 */ });  rs.insert(l, r);
//   ODT の区間加算: rs.apply(l, r, [&](ll& v) { v += d; });   区間和: rs.for_each(l, r, [&](ll a, ll b, ll v) { s += (b - a) * v; });
// 注意: 1点版の x + 1、prev_uncovered の l - 1 は型の端でオーバーフローする。イテレータ経由で r を書き換えない
template <class T, class V = std::monostate>
struct RangeSet {
    using D = std::conditional_t<std::is_integral_v<T>, std::common_type_t<T, long long>, T>;
    struct Node {
        T l;
        mutable T r;
        [[no_unique_address]] mutable V v;
        friend bool operator<(const Node& a, const Node& b) { return a.l < b.l; }
        friend bool operator<(const Node& a, const T& x) { return a.l < x; }
        friend bool operator<(const T& x, const Node& a) { return x < a.l; }
    };
    using Set = std::set<Node, std::less<>>;
    using iterator = typename Set::const_iterator;
    struct Nop { template <class... A> void operator()(A&&...) const {} };

    Set st;
    D len = 0;

    iterator begin() const { return st.begin(); }
    iterator end() const { return st.end(); }
    int size() const { return (int)st.size(); }
    bool empty() const { return st.empty(); }

    iterator lower(T x) const {
        auto it = st.upper_bound(x);
        if (it != st.begin() && x < std::prev(it)->r) --it;
        return it;
    }
    iterator find(T x) const {
        auto it = lower(x);
        return it != st.end() && !(x < it->l) ? it : st.end();
    }
    bool contains(T x) const { return find(x) != st.end(); }
    V get(T x, V dflt = V()) const {
        auto it = find(x);
        return it == st.end() ? dflt : it->v;
    }
    bool covers(T l, T r) const {
        if (!(l < r)) return true;
        auto it = find(l);
        if (it == st.end()) return false;
        T y = it->r;
        for (++it; y < r && it != st.end() && !(y < it->l); ++it) y = it->r;
        return !(y < r);
    }
    bool intersects(T l, T r) const {
        auto it = lower(l);
        return l < r && it != st.end() && it->l < r;
    }
    T mex(T x = T(0)) const {
        auto it = find(x);
        if (it == st.end()) return x;
        T y = it->r;
        for (++it; it != st.end() && !(y < it->l); ++it) y = it->r;
        return y;
    }
    T prev_uncovered(T x) const {
        auto it = find(x);
        if (it == st.end()) return x;
        T y = it->l;
        while (it != st.begin() && !(std::prev(it)->r < y)) y = (--it)->l;
        return y - 1;
    }
    template <class F> void for_each(T l, T r, F&& f) const {
        for (auto it = lower(l); l < r && it != st.end() && it->l < r; ++it)
            call_(f, std::max(it->l, l), std::min(it->r, r), it->v);
    }
    template <class F> void for_each_gap(T l, T r, F&& f) const {
        if (!(l < r)) return;
        T x = l;
        for (auto it = lower(l); it != st.end() && it->l < r; ++it) {
            if (x < it->l) f(x, it->l);
            x = it->r;
        }
        if (x < r) f(x, r);
    }

    template <class FA = Nop, class FD = Nop>
    D assign(T l, T r, V v, FA&& on_add = FA(), FD&& on_del = FD()) {
        if (!(l < r)) return D(0);
        auto it = st.upper_bound(l);
        if (it != st.begin()) {
            auto p = std::prev(it);  // 左端 <= l で最後の区間
            if (l < p->r) {          // p が l を含む
                if (!(p->r < r) && p->v == v) return D(0);  // 変化なし
                it = p;
                if (!(p->l < l) && it != st.begin() && !(std::prev(it)->r < l)) --it;  // l で接する左隣
            } else if (!(p->r < l)) it = p;  // p が l で接する
        }
        T nl = l, nr = r;
        D cov = 0;
        std::optional<Node> kl, kr;  // 値が違う区間の [l, r) からはみ出す部分
        while (it != st.end() && !(r < it->l)) {
            const T a = it->l, b = it->r;
            const bool same = it->v == v, ov = a < r && l < b;
            if (!ov && !same) {
                if (a < l) { ++it; continue; }  // 左で接するだけ
                break;                          // 右で接するだけ
            }
            if (ov) cov += width_(std::max(a, l), std::min(b, r));
            if (same) nl = std::min(nl, a), nr = std::max(nr, b);
            else {
                if (a < l) kl = Node{a, l, it->v};
                if (r < b) kr = Node{r, b, it->v};
            }
            it = del_(it, on_del);
        }
        if (kl) add_(it, kl->l, kl->r, std::move(kl->v), on_add);
        add_(it, nl, nr, std::move(v), on_add);
        if (kr) add_(it, kr->l, kr->r, std::move(kr->v), on_add);
        return width_(l, r) - cov;
    }
    template <class FA = Nop, class FD = Nop>
    D insert(T l, T r, FA&& on_add = FA(), FD&& on_del = FD()) { return assign(l, r, V(), on_add, on_del); }
    template <class FA = Nop, class FD = Nop>
    D erase(T l, T r, FA&& on_add = FA(), FD&& on_del = FD()) {
        if (!(l < r)) return D(0);
        auto it = lower(l);
        D rem = 0;
        std::optional<Node> kl, kr;
        while (it != st.end() && it->l < r) {
            rem += width_(std::max(it->l, l), std::min(it->r, r));
            if (it->l < l) kl = Node{it->l, l, it->v};
            if (r < it->r) kr = Node{r, it->r, it->v};
            it = del_(it, on_del);
        }
        if (kl) add_(it, kl->l, kl->r, std::move(kl->v), on_add);
        if (kr) add_(it, kr->l, kr->r, std::move(kr->v), on_add);
        return rem;
    }
    D insert(T x) { return insert(x, x + 1); }
    D erase(T x) { return erase(x, x + 1); }
    template <class F> void apply(T l, T r, F&& f) {
        if (!(l < r)) return;
        auto itr = split(r), it = split(l);
        for (auto jt = it; jt != itr; ++jt) {
            if constexpr (std::is_invocable_v<F&, T, T, V&>) f(jt->l, jt->r, jt->v);
            else f(jt->v);
        }
        auto stop = itr == st.end() ? itr : std::next(itr);  // r から始まる区間までを併合対象に
        if (it != st.begin()) --it;
        while (it != stop) {
            auto nx = std::next(it);
            if (nx == stop) break;
            if (!(it->r < nx->l) && it->v == nx->v) it->r = nx->r, st.erase(nx);
            else it = nx;
        }
    }
    template <class FA = Nop, class FD = Nop>
    iterator split(T x, FA&& on_add = FA(), FD&& on_del = FD()) {
        auto it = lower(x);
        if (it == st.end() || !(it->l < x)) return it;
        call_(on_del, it->l, it->r, it->v);
        T b = it->r;
        it->r = x;
        call_(on_add, it->l, x, it->v);
        auto jt = st.emplace_hint(std::next(it), Node{x, b, it->v});
        call_(on_add, x, b, jt->v);
        return jt;
    }
    template <class FD = Nop>
    void clear(FD&& on_del = FD()) {
        if constexpr (!std::is_same_v<std::decay_t<FD>, Nop>)
            for (auto& n : st) call_(on_del, n.l, n.r, n.v);
        st.clear(), len = 0;
    }
    friend std::ostream& operator<<(std::ostream& os, const RangeSet& s) {
        os << '{';
        bool first = true;
        for (auto& n : s.st) {
            os << (first ? "" : ", ") << '[' << n.l << ", " << n.r << ')';
            if constexpr (!std::is_same_v<V, std::monostate>) os << ':' << n.v;
            first = false;
        }
        return os << '}';
    }

    // 内部用
    static D width_(T l, T r) { return D(r) - D(l); }
    template <class F> static void call_(F& f, T l, T r, const V& v) {
        if constexpr (std::is_invocable_v<F&, T, T, const V&>) f(l, r, v);
        else f(l, r);
    }
    template <class F> iterator add_(iterator hint, T l, T r, V v, F& f) {
        auto it = st.emplace_hint(hint, Node{l, r, std::move(v)});
        len += width_(l, r);
        call_(f, l, r, it->v);
        return it;
    }
    template <class F> iterator del_(iterator it, F& f) {
        call_(f, it->l, it->r, it->v);
        len -= width_(it->l, it->r);
        return st.erase(it);
    }
};
