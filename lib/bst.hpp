// @name BST (weight-balanced tree)
// @prefix bst
// @description 重み平衡二分探索木: k番目/順位/前後探索/モノイド集約(prod・max_right)/split・merge・meld

// ============================================================================
//  BST<K, S, op, e> : 重み平衡二分探索木（Weight-Balanced Tree, Δ=3, Γ=2）
//  「常にソートされた列（多重集合）」として使う。各要素は key(K) と val(S) を持ち、
//  S がモノイド (op, e) なら添字区間・key 区間の総積が取れる（ACL segtree と同じ op / e）。
//  S を省略すると集約なしの順序付き多重集合になる。全インスタンスでノードプールを共有。
//  計算量は特記なきものは O(log n)。添字はすべて 0-indexed・半開区間。
// ----------------------------------------------------------------------------
//  [構築]  BST<K> t;  /  BST<K, S, op, e> t;
//          BST t(vector<K>) / BST t(vector<pair<K, S>>)  O(n)（未ソートなら内部でソート）
//          BST<...>::reserve(n)                           ノードを先に確保（任意）
//  [挿入・削除]
//    int  insert(key) / insert(key, val)   挿入位置を返す（同じ key の中では末尾に入る）
//         ※ insert(key) の val は S(key)（S が K から構築できるとき）
//    bool erase(key)            key を1個削除（無ければ false）
//    K    erase_at(k)           k 番目を削除してその key を返す
//    int  erase_range(lo, hi)   key ∈ [lo, hi) を全削除、個数を返す   O(log n + 削除数)
//    void clear()               O(n)
//  [探索]
//    int  size() / bool empty()                 O(1)
//    K    kth(k)  /  min_elem()  /  max_elem()
//    int  count_less(x) / count_leq(x)          x 未満 / x 以下の個数
//    int  lower_bound(x) / upper_bound(x)       上と同じ（添字として使う用）
//    int  count(x) / count_range(lo, hi) / bool contains(x)
//    optional<K> lt(x) / le(x) / gt(x) / ge(x)  x 未満の最大 / 以下の最大 / 超の最小 / 以上の最小
//  [集約]（S がモノイドのとき）
//    S    all_prod()                  O(1)
//    S    prod(l, r)                  添字 [l, r) の総積
//    S    prod_key(lo, hi)            key ∈ [lo, hi) の総積
//    int  max_right(l, f) / min_left(r, f)   ACL segtree と同じ二分探索（f(e()) = true が必要）
//    S    val_at(k) / void set_at(k, v)       k 番目の val の取得 / 書き換え
//    S    get(key) / void set(key, v) / void update(key, f)
//         map 風。get は無ければ e()。set・update は無ければ挿入。update は f(S&) で書き換え
//         （key が重複しているときはそのうち1個が対象）
//  [分割・併合]
//    BST  split_key(x)    this = x 未満,   戻り値 = x 以上
//    BST  split_at(k)     this = 先頭 k 個, 戻り値 = 残り
//    void merge(r)        this の後ろに r を連結（max(this) <= min(r) が前提）。r は空になる
//    void meld(o)         任意の2つの和集合（マージ）。o は空になる  O(m log(n/m + 1)), m は小さい方
//  [走査]  for_each(f) : 昇順に f(key, val)  O(n)  /  vector<K> keys()  O(n)
//  注意: コピー不可（ムーブのみ）。破棄時にノードをプールへ返す。
// ============================================================================
namespace bst_internal {
struct None {};
inline None none_op(None, None) { return {}; }
inline None none_e() { return {}; }
}  // namespace bst_internal
template <class K, class S = bst_internal::None, auto op = bst_internal::none_op, auto e = bst_internal::none_e>
struct BST {
   private:
    static constexpr bool AGG = !std::is_empty_v<S>;
    // ホット部（探索で読む）とコールド部（集約値）を別配列に置く。子の要約は親が持つ
    struct Hot {
        int ch[2], cnt[2];  // 子, 子の部分木サイズ
        K key;
    };
    struct Cold {
        S val, agg[2];  // 自分の値, 子の部分木の総積
    };
    static inline std::vector<Hot>* hv = nullptr;  // 破棄順問題を避けるため解放しない
    static inline std::vector<Cold>* cv = nullptr;
    static inline Hot* d = nullptr;
    static inline Cold* cd = nullptr;
    static inline Cold dummy{};
    static inline int freelist = 0;
    int root = 0;

    struct RootTag {};
    BST(int rt, RootTag) : root(rt) {}
    static void init_pool() {
        hv = new std::vector<Hot>(1);
        (*hv)[0].cnt[0] = -1;  // 番兵 0: size = cnt[0] + cnt[1] + 1 = 0
        cv = new std::vector<Cold>(AGG ? 1 : 0, Cold{e(), {e(), e()}});
        refresh();
    }
    static void refresh() {
        d = hv->data();
        if constexpr (AGG) cd = cv->data();
    }
    static Cold& C(int t) {
        if constexpr (AGG) return cd[t];
        else return dummy;
    }
    static S make_val(const K& k) {
        if constexpr (std::is_constructible_v<S, const K&>) return S(k);
        else {
            static_assert(!AGG, "insert(key) には S が K から構築できる必要があります。insert(key, val) を使ってください");
            return e();
        }
    }
    static int newnode(const K& k, const S& v) {
        if (!hv) init_pool();
        if (freelist) {
            int t = freelist;
            freelist = d[t].ch[0], d[t] = Hot{{0, 0}, {0, 0}, k};
            if constexpr (AGG) cd[t] = Cold{v, {e(), e()}};
            return t;
        }
        hv->push_back(Hot{{0, 0}, {0, 0}, k});
        if constexpr (AGG) cv->push_back(Cold{v, {e(), e()}});
        refresh();
        return (int)hv->size() - 1;
    }
    static void freenode(int t) { d[t].ch[0] = freelist, freelist = t; }
    static void free_tree(int t) {
        if (!t) return;
        free_tree(d[t].ch[0]), free_tree(d[t].ch[1]), freenode(t);
    }
    // 両方の子を先読み（方向の分岐予測が外れても次ノードが届いている）。
    // 普通の関数だと GCC が「副作用なし」とみなして呼び出しごと消すので always_inline
    __attribute__((always_inline)) static void pf(int t) {
        __builtin_prefetch(d + d[t].ch[0]), __builtin_prefetch(d + d[t].ch[1]);
    }
    __attribute__((always_inline)) static void pfa(int t) {  // 集約も読む走査用
        int a = d[t].ch[0], b = d[t].ch[1];
        __builtin_prefetch(d + a), __builtin_prefetch(d + b);
        if constexpr (AGG) __builtin_prefetch(cd + a), __builtin_prefetch(cd + b);
    }
    __attribute__((always_inline)) static void pfc(int t) {  // 後で経路を更新するノード
        if constexpr (AGG) __builtin_prefetch(cd + t, 1);
    }
    static int size_of(int t) { return d[t].cnt[0] + d[t].cnt[1] + 1; }
    static S total(int t) {
        const Cold& c = C(t);
        return op(op(c.agg[0], c.val), c.agg[1]);
    }
    static void set_child(int t, int c, int x) {
        d[t].ch[c] = x, d[t].cnt[c] = size_of(x);
        if constexpr (AGG) cd[t].agg[c] = total(x);
    }
    template <class G>
    static int build(int lo, int hi, const G& make) {
        if (lo >= hi) return 0;
        int mid = (lo + hi) >> 1, t = make(mid);
        int L = build(lo, mid, make), R = build(mid + 1, hi, make);
        set_child(t, 0, L), set_child(t, 1, R);
        return t;
    }
    // ---- 平衡（重み = サイズ + 1。Δ=3, Γ=2） ----
    static int rotate(int t, int c) {  // c 側の子を持ち上げる
        int s = d[t].ch[c], o = c ^ 1;
        Hot &T = d[t], &U = d[s];
        T.ch[c] = U.ch[o], T.cnt[c] = U.cnt[o];
        U.ch[o] = t, U.cnt[o] = T.cnt[0] + T.cnt[1] + 1;
        if constexpr (AGG) {
            Cold &TC = cd[t], &UC = cd[s];
            TC.agg[c] = UC.agg[o];
            UC.agg[o] = op(op(TC.agg[0], TC.val), TC.agg[1]);
        }
        return s;
    }
    static int lift(int p, int h) {  // h 側が重すぎる p を回転で直す（内側の孫が重ければ二重回転）
        int s = d[p].ch[h];
        if (d[s].cnt[h ^ 1] + 1 >= 2 * (d[s].cnt[h] + 1)) d[p].ch[h] = rotate(s, h ^ 1);
        return rotate(p, h);
    }
    static int balance(int t) {
        int a = d[t].cnt[0] + 1, b = d[t].cnt[1] + 1;
        return 3 * a < b ? lift(t, 1) : 3 * b < a ? lift(t, 0) : t;
    }
    static int join(int l, int m, int r) {  // l < m < r（m は単独ノード）
        int a = size_of(l) + 1, b = size_of(r) + 1;
        if (3 * a < b) {
            set_child(r, 0, join(l, m, d[r].ch[0]));
            return balance(r);
        }
        if (3 * b < a) {
            set_child(l, 1, join(d[l].ch[1], m, r));
            return balance(l);
        }
        set_child(m, 0, l), set_child(m, 1, r);
        return m;
    }
    static int pop_last(int t, int& m) {
        if (!d[t].ch[1]) return m = t, d[t].ch[0];
        set_child(t, 1, pop_last(d[t].ch[1], m));
        return balance(t);
    }
    static int merge2(int l, int r) {
        if (!l || !r) return l | r;
        int m;
        l = pop_last(l, m);
        return join(l, m, r);
    }
    static std::pair<int, int> split_lt(int t, const K& key) {  // (key 未満, key 以上)
        if (!t) return {0, 0};
        pfa(t);
        int l = d[t].ch[0], r = d[t].ch[1];
        if (d[t].key < key) {
            auto [a, b] = split_lt(r, key);
            return {join(l, t, a), b};
        }
        auto [a, b] = split_lt(l, key);
        return {a, join(b, t, r)};
    }
    static std::pair<int, int> split_k(int t, int k) {  // (先頭 k 個, 残り)
        if (!t) return {0, 0};
        pfa(t);
        int l = d[t].ch[0], r = d[t].ch[1], ls = d[t].cnt[0];
        if (k <= ls) {
            auto [a, b] = split_k(l, k);
            return {a, join(b, t, r)};
        }
        auto [a, b] = split_k(r, k - ls - 1);
        return {join(l, t, a), b};
    }
    static int unite(int a, int b) {
        if (!a || !b) return a | b;
        if (size_of(a) < size_of(b)) std::swap(a, b);
        int l = d[a].ch[0], r = d[a].ch[1];
        auto [x, y] = split_lt(b, d[a].key);
        int L = unite(l, x);
        int R = unite(r, y);
        return join(L, a, R);
    }
    // 降下で記録した経路を下から付け直す（挿入のサイズ加算は降下時に済ませてある）。
    // 挿入なら通った側 c、削除なら反対側だけが重くなりうるので片側だけ判定する
    template <bool ins>
    void fix_path(const int* path, const char* dir, int np, int child) {
        for (int i = np - 1; i >= 0; i--) {
            int p = path[i], c = dir[i], h = ins ? c : c ^ 1;
            d[p].ch[c] = child;
            if constexpr (!ins) d[p].cnt[c]--;
            if constexpr (AGG) cd[p].agg[c] = total(child);
            child = 3 * (d[p].cnt[h ^ 1] + 1) < d[p].cnt[h] + 1 ? lift(p, h) : p;
        }
        root = child;
    }
    static void fix_agg(const int* path, const char* dir, int np, int child) {
        if constexpr (AGG)
            for (int i = np - 1; i >= 0; i--) cd[path[i]].agg[(int)dir[i]] = total(child), child = path[i];
    }
    void remove_node(int t, int* path, char* dir, int np) {
        int repl;
        if (!d[t].ch[0] || !d[t].ch[1]) {
            repl = d[t].ch[0] | d[t].ch[1];
            freenode(t);
        } else {  // 後続ノードの中身を移して、後続ノードを消す
            path[np] = t, dir[np++] = 1;
            int s = d[t].ch[1];
            while (d[s].ch[0]) pfc(s), path[np] = s, dir[np++] = 0, s = d[s].ch[0];
            d[t].key = std::move(d[s].key);
            if constexpr (AGG) cd[t].val = std::move(cd[s].val);
            repl = d[s].ch[1];
            freenode(s);
        }
        fix_path<false>(path, dir, np, repl);
    }
    template <bool strict>
    static int rank_in(int t, const K& x) {  // 部分木 t の x 未満（strict）/ x 以下 の個数
        int res = 0;
        while (t) {
            pf(t);
            bool go;
            if constexpr (strict) go = d[t].key < x;
            else go = !(x < d[t].key);
            if (go) res += d[t].cnt[0] + 1, t = d[t].ch[1];
            else t = d[t].ch[0];
        }
        return res;
    }
    template <int M>  // 0: x 未満の最大, 1: x 以下の最大, 2: x 超の最小, 3: x 以上の最小
    std::optional<K> nearest(const K& x) const {
        int t = root, c = 0;
        while (t) {
            pf(t);
            const K& k = d[t].key;
            bool in = M == 0 ? k < x : M == 1 ? !(x < k) : M == 2 ? x < k : !(k < x);
            if (in) c = t;
            t = d[t].ch[(M < 2) == in];
        }
        if (!c) return std::nullopt;
        return d[c].key;
    }
    int kth_node(int k) const {
        assert(0 <= k && k < size());
        int t = root;
        while (true) {
            pf(t);
            int ls = d[t].cnt[0];
            if (k == ls) return t;
            int g = k > ls;
            k -= g ? ls + 1 : 0;
            t = d[t].ch[g];
        }
    }
    static S suffix(int t, int k) {  // 部分木 t の添字 k 以降（非空）
        S res = e();
        while (t) {
            if (k == 0) return op(total(t), res);
            pfa(t);
            int ls = d[t].cnt[0];
            if (k <= ls) {
                res = op(op(C(t).val, C(t).agg[1]), res);
                if (k == ls) return res;
                t = d[t].ch[0];
            } else k -= ls + 1, t = d[t].ch[1];
        }
        return res;
    }
    static S prefix(int t, int k) {  // 部分木 t の添字 k 未満
        S res = e();
        while (t && k > 0) {
            pfa(t);
            int ls = d[t].cnt[0];
            if (k <= ls) {
                if (k == ls) return op(res, C(t).agg[0]);
                t = d[t].ch[0];
            } else {
                res = op(op(res, C(t).agg[0]), C(t).val), k -= ls + 1;
                if (k >= d[t].cnt[1]) return op(res, C(t).agg[1]);
                t = d[t].ch[1];
            }
        }
        return res;
    }
    template <class F>
    static void dfs(int t, F& f) {
        if (!t) return;
        dfs(d[t].ch[0], f), f(d[t].key, C(t).val), dfs(d[t].ch[1], f);
    }

   public:
    BST() {}
    explicit BST(std::vector<K> ks) {
        if (!std::is_sorted(ks.begin(), ks.end())) std::sort(ks.begin(), ks.end());
        reserve((hv ? (int)hv->size() : 0) + (int)ks.size());
        root = build(0, (int)ks.size(), [&](int i) { return newnode(ks[i], make_val(ks[i])); });
    }
    explicit BST(std::vector<std::pair<K, S>> kv) {
        auto cmp = [](const auto& a, const auto& b) { return a.first < b.first; };
        if (!std::is_sorted(kv.begin(), kv.end(), cmp)) std::stable_sort(kv.begin(), kv.end(), cmp);
        reserve((hv ? (int)hv->size() : 0) + (int)kv.size());
        root = build(0, (int)kv.size(), [&](int i) { return newnode(kv[i].first, kv[i].second); });
    }
    BST(const BST&) = delete;
    BST& operator=(const BST&) = delete;
    BST(BST&& o) noexcept : root(o.root) { o.root = 0; }
    BST& operator=(BST&& o) noexcept {
        if (this != &o) clear(), root = o.root, o.root = 0;
        return *this;
    }
    ~BST() { clear(); }
    static void reserve(int n) {
        if (!hv) init_pool();
        hv->reserve(n + 1);
        if constexpr (AGG) cv->reserve(n + 1);
        refresh();
    }

    int size() const { return root ? size_of(root) : 0; }
    bool empty() const { return root == 0; }
    void clear() { free_tree(root), root = 0; }

    // ---------------- 挿入・削除 ----------------
    int insert(const K& key) { return insert(key, make_val(key)); }
    int insert(const K& key, const S& val) {
        int nw = newnode(key, val);
        if (!root) return root = nw, 0;
        int path[128], np = 0, t = root, idx = 0;
        char dir[128];
        while (t) {
            pf(t), pfc(t);
            path[np] = t;
            int g = !(key < d[t].key);
            if (g) idx += d[t].cnt[0] + 1;
            dir[np++] = (char)g, d[t].cnt[g]++;
            t = d[t].ch[g];
        }
        fix_path<true>(path, dir, np, nw);
        return idx;
    }
    bool erase(const K& key) {
        int path[128], np = 0, t = root;
        char dir[128];
        while (t) {
            pf(t), pfc(t);
            if (key < d[t].key) path[np] = t, dir[np++] = 0, t = d[t].ch[0];
            else if (d[t].key < key) path[np] = t, dir[np++] = 1, t = d[t].ch[1];
            else break;
        }
        if (!t) return false;
        remove_node(t, path, dir, np);
        return true;
    }
    K erase_at(int k) {
        assert(0 <= k && k < size());
        int path[128], np = 0, t = root;
        char dir[128];
        while (true) {
            pf(t), pfc(t);
            int ls = d[t].cnt[0];
            if (k < ls) path[np] = t, dir[np++] = 0, t = d[t].ch[0];
            else if (k > ls) k -= ls + 1, path[np] = t, dir[np++] = 1, t = d[t].ch[1];
            else break;
        }
        K res = d[t].key;
        remove_node(t, path, dir, np);
        return res;
    }
    int erase_range(const K& lo, const K& hi) {
        if (!(lo < hi)) return 0;
        auto [a, b] = split_lt(root, lo);
        auto [m, c] = split_lt(b, hi);
        int res = m ? size_of(m) : 0;
        free_tree(m);
        root = merge2(a, c);
        return res;
    }

    // ---------------- 探索 ----------------
    K kth(int k) const { return d[kth_node(k)].key; }
    K min_elem() const {
        assert(root);
        int t = root;
        while (d[t].ch[0]) t = d[t].ch[0];
        return d[t].key;
    }
    K max_elem() const {
        assert(root);
        int t = root;
        while (d[t].ch[1]) t = d[t].ch[1];
        return d[t].key;
    }
    int count_less(const K& x) const { return rank_in<true>(root, x); }
    int count_leq(const K& x) const { return rank_in<false>(root, x); }
    int lower_bound(const K& x) const { return count_less(x); }
    int upper_bound(const K& x) const { return count_leq(x); }
    int count(const K& x) const {
        int t = root;
        while (t) {
            pf(t);
            if (x < d[t].key) t = d[t].ch[0];
            else if (d[t].key < x) t = d[t].ch[1];
            else return 1 + d[t].cnt[0] - rank_in<true>(d[t].ch[0], x) + rank_in<false>(d[t].ch[1], x);
        }
        return 0;
    }
    int count_range(const K& lo, const K& hi) const {
        if (!(lo < hi)) return 0;
        int t = root;
        while (t) {
            pf(t);
            if (!(d[t].key < hi)) t = d[t].ch[0];
            else if (d[t].key < lo) t = d[t].ch[1];
            else return 1 + d[t].cnt[0] - rank_in<true>(d[t].ch[0], lo) + rank_in<true>(d[t].ch[1], hi);
        }
        return 0;
    }
    bool contains(const K& x) const {
        int t = root;
        while (t) {
            pf(t);
            if (x < d[t].key) t = d[t].ch[0];
            else if (d[t].key < x) t = d[t].ch[1];
            else return true;
        }
        return false;
    }
    std::optional<K> lt(const K& x) const { return nearest<0>(x); }
    std::optional<K> le(const K& x) const { return nearest<1>(x); }
    std::optional<K> gt(const K& x) const { return nearest<2>(x); }
    std::optional<K> ge(const K& x) const { return nearest<3>(x); }

    // ---------------- 集約 ----------------
    S all_prod() const { return root ? total(root) : e(); }
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return e();
        int t = root;
        while (true) {
            pfa(t);
            int ls = d[t].cnt[0];
            if (r <= ls) t = d[t].ch[0];
            else if (l > ls) l -= ls + 1, r -= ls + 1, t = d[t].ch[1];
            else {
                S L = l < ls ? suffix(d[t].ch[0], l) : e();
                S R = r - ls - 1 > 0 ? prefix(d[t].ch[1], r - ls - 1) : e();
                return op(op(L, C(t).val), R);
            }
        }
    }
    S prod_key(const K& lo, const K& hi) const {
        if (!(lo < hi)) return e();
        int t = root;
        while (t) {
            pfa(t);
            if (!(d[t].key < hi)) t = d[t].ch[0];
            else if (d[t].key < lo) t = d[t].ch[1];
            else break;
        }
        if (!t) return e();
        S L = e(), R = e();
        for (int u = d[t].ch[0]; u;) {
            pfa(u);
            if (!(d[u].key < lo)) L = op(op(C(u).val, C(u).agg[1]), L), u = d[u].ch[0];
            else u = d[u].ch[1];
        }
        for (int u = d[t].ch[1]; u;) {
            pfa(u);
            if (d[u].key < hi) R = op(op(R, C(u).agg[0]), C(u).val), u = d[u].ch[1];
            else u = d[u].ch[0];
        }
        return op(op(L, C(t).val), R);
    }
    template <class F>
    int max_right(int l, F f) const {
        assert(0 <= l && l <= size() && f(e()));
        int n = size();
        if (l == n) return n;
        int pv[128], pi[128], np = 0;
        int t = root, k = l, off = 0;
        while (true) {  // 位置 l への経路。l 以降の部分は「ノード + その右部分木」の列になる
            pf(t);
            int ls = d[t].cnt[0];
            if (k < ls) pfc(t), pv[np] = t, pi[np++] = off + ls, t = d[t].ch[0];
            else if (k == ls) {
                pfc(t), pv[np] = t, pi[np++] = off + ls;
                break;
            } else k -= ls + 1, off += ls + 1, t = d[t].ch[1];
        }
        S acc = e();
        while (np) {
            np--;
            int v = pv[np], idx = pi[np];
            S x = op(acc, C(v).val);
            if (!f(x)) return idx;
            acc = x, x = op(acc, C(v).agg[1]);
            if (f(x)) {
                acc = x;
                continue;
            }
            int c = d[v].ch[1];
            idx++;
            while (true) {
                pfa(c);
                x = op(acc, C(c).agg[0]);
                if (!f(x)) {
                    c = d[c].ch[0];
                    continue;
                }
                acc = x, idx += d[c].cnt[0], x = op(acc, C(c).val);
                if (!f(x)) return idx;
                acc = x, idx++, c = d[c].ch[1];
            }
        }
        return n;
    }
    template <class F>
    int min_left(int r, F f) const {
        assert(0 <= r && r <= size() && f(e()));
        if (r == 0) return 0;
        int pv[128], pi[128], np = 0;
        int t = root, k = r - 1, off = 0;
        while (true) {
            pf(t);
            int ls = d[t].cnt[0];
            if (k < ls) t = d[t].ch[0];
            else {
                pfc(t), pv[np] = t, pi[np++] = off + ls;
                if (k == ls) break;
                k -= ls + 1, off += ls + 1, t = d[t].ch[1];
            }
        }
        S acc = e();
        while (np) {
            np--;
            int v = pv[np], iv = pi[np];
            S x = op(C(v).val, acc);
            if (!f(x)) return iv + 1;
            acc = x, x = op(C(v).agg[0], acc);
            if (f(x)) {
                acc = x;
                continue;
            }
            int c = d[v].ch[0], idx = iv;
            while (true) {
                pfa(c);
                x = op(C(c).agg[1], acc);
                if (!f(x)) {
                    c = d[c].ch[1];
                    continue;
                }
                acc = x, idx -= d[c].cnt[1], x = op(C(c).val, acc);
                if (!f(x)) return idx;
                acc = x, idx--, c = d[c].ch[0];
            }
        }
        return 0;
    }
    S val_at(int k) const { return C(kth_node(k)).val; }
    void set_at(int k, const S& v) {
        assert(0 <= k && k < size());
        int path[128], np = 0, t = root;
        char dir[128];
        while (true) {
            pf(t), pfc(t);
            int ls = d[t].cnt[0];
            if (k < ls) path[np] = t, dir[np++] = 0, t = d[t].ch[0];
            else if (k > ls) k -= ls + 1, path[np] = t, dir[np++] = 1, t = d[t].ch[1];
            else break;
        }
        C(t).val = v;
        fix_agg(path, dir, np, t);
    }
    S get(const K& key) const {
        int t = root;
        while (t) {
            pf(t);
            if (key < d[t].key) t = d[t].ch[0];
            else if (d[t].key < key) t = d[t].ch[1];
            else return C(t).val;
        }
        return e();
    }
    template <class F>
    void update(const K& key, F f) {
        int path[128], np = 0, t = root;
        char dir[128];
        while (t) {
            pf(t), pfc(t);
            if (key < d[t].key) path[np] = t, dir[np++] = 0, t = d[t].ch[0];
            else if (d[t].key < key) path[np] = t, dir[np++] = 1, t = d[t].ch[1];
            else {
                f(C(t).val);
                fix_agg(path, dir, np, t);
                return;
            }
        }
        S v = e();
        f(v);
        insert(key, v);
    }
    void set(const K& key, const S& v) {
        update(key, [&](S& x) { x = v; });
    }

    // ---------------- 分割・併合 ----------------
    BST split_key(const K& x) {
        auto [a, b] = split_lt(root, x);
        root = a;
        return BST(b, RootTag{});
    }
    BST split_at(int k) {
        assert(0 <= k && k <= size());
        auto [a, b] = split_k(root, k);
        root = a;
        return BST(b, RootTag{});
    }
    void merge(BST& r) {
        if (this != &r) root = merge2(root, r.root), r.root = 0;
    }
    void meld(BST& o) {
        if (this != &o) root = unite(root, o.root), o.root = 0;
    }

    // ---------------- 走査 ----------------
    template <class F>
    void for_each(F f) const { dfs(root, f); }
    std::vector<K> keys() const {
        std::vector<K> res;
        res.reserve(size());
        for_each([&](const K& k, const S&) { res.push_back(k); });
        return res;
    }
};
