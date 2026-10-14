// @name BinaryTrie
// @prefix binarytrie
// @description Binary Trie (整数/bitset キー, 遅延 xor_all, ノード再利用, insert/erase/count/count_less/count_leq/count_range/kth/min_xor/max_xor)

// BinaryTrie<LOG, Key>
//   Key : 整数型 (0 <= x < 2^LOG, LOG <= 63) または bitset<B> (LOG <= B)
//   bit LOG-1 が最上位。全要素への XOR (xor_all) を遅延で持つ
//   要素数が減ったノードは再利用するので、ノード数は常に (要素数 * LOG + 1) 以下
template <class Key, class = void> struct TrieKeyTraits;
template <class Key>
struct TrieKeyTraits<Key, enable_if_t<is_integral_v<Key>>> {
    using U = make_unsigned_t<Key>;
    static int  get(Key x, int i) { return (int)(((U)x >> i) & 1); }
    static void set(Key& x, int i, int b) { x = (Key)(((U)x & ~((U)1 << i)) | ((U)b << i)); }
    static Key  zero() { return 0; }
};
template <size_t B>
struct TrieKeyTraits<bitset<B>, void> {
    static int  get(const bitset<B>& x, int i) { return (int)x[i]; }
    static void set(bitset<B>& x, int i, int b) { x[i] = (bool)b; }
    static bitset<B> zero() { return bitset<B>(); }
};

template <int LOG = 30, class Key = long long, class KT = TrieKeyTraits<Key>>
struct BinaryTrie {
    struct Node { int nxt[2] = {-1, -1}; int cnt = 0; };
    vector<Node> t;
    vector<int> pool;  // 再利用待ちノード
    Key lz;            // 全要素に適用中の XOR (遅延)
    BinaryTrie() : t(1), lz(KT::zero()) {}

    int  size() const { return t[0].cnt; }
    bool empty() const { return size() == 0; }
    void clear() { t.assign(1, Node()); pool.clear(); lz = KT::zero(); }
    void reserve(int n_elem) { t.reserve((size_t)n_elem * LOG + 1); }
    void xor_all(const Key& x) { lz = lz ^ x; }  // 全要素 e を e ^ x に

    void insert(Key x, int k = 1) {  // x を k 個追加
        x = x ^ lz;
        int v = 0; t[v].cnt += k;
        for (int i = LOG - 1; i >= 0; i--) {
            int b = KT::get(x, i);
            if (t[v].nxt[b] == -1) { int c = new_node(); t[v].nxt[b] = c; }
            v = t[v].nxt[b];
            t[v].cnt += k;
        }
    }
    bool erase(Key x, int k = 1) {  // x を k 個削除。足りなければ何もせず false
        x = x ^ lz;
        int path[LOG + 1], bit[LOG + 1];
        int v = 0; path[0] = 0;
        for (int d = 1; d <= LOG; d++) {
            int b = KT::get(x, LOG - d), c = t[v].nxt[b];
            if (c == -1) return false;
            bit[d] = b; v = c; path[d] = v;
        }
        if (t[v].cnt < k) return false;
        for (int d = 0; d <= LOG; d++) t[path[d]].cnt -= k;
        for (int d = LOG; d >= 1 && t[path[d]].cnt == 0; d--) {
            t[path[d - 1]].nxt[bit[d]] = -1;
            pool.push_back(path[d]);
        }
        return true;
    }
    int count(Key x) const {  // x の個数
        x = x ^ lz;
        int v = 0;
        for (int i = LOG - 1; i >= 0; i--) {
            v = t[v].nxt[KT::get(x, i)];
            if (v == -1) return 0;
        }
        return t[v].cnt;
    }
    int count_less(const Key& x) const { return count_lower(x, false); }  // e <  x の個数
    int count_leq(const Key& x) const { return count_lower(x, true); }    // e <= x の個数
    int count_range(const Key& l, const Key& r) const { return count_less(r) - count_less(l); }  // [l, r)

    Key kth(int k) const {  // 小さい方から k 番目 (0-indexed)
        assert(0 <= k && k < size());
        int v = 0; Key res = KT::zero();
        for (int i = LOG - 1; i >= 0; i--) {
            int lb = KT::get(lz, i);
            int c = t[v].nxt[lb], l = (c == -1 ? 0 : t[c].cnt);
            if (k < l) v = c;
            else { k -= l; v = t[v].nxt[lb ^ 1]; KT::set(res, i, 1); }
        }
        return res;
    }
    Key min_elem() const { return kth(0); }
    Key max_elem() const { return kth(size() - 1); }
    // x との XOR が最小 / 最大になる要素 (XOR 値自体は 戻り値 ^ x)
    Key min_xor(const Key& x) const { return xor_walk(x, 0); }
    Key max_xor(const Key& x) const { return xor_walk(x, 1); }

  private:
    int new_node() {
        if (!pool.empty()) { int v = pool.back(); pool.pop_back(); t[v] = Node(); return v; }
        t.emplace_back();
        return (int)t.size() - 1;
    }
    int count_lower(const Key& x, bool incl) const {
        int v = 0, res = 0;
        for (int i = LOG - 1; i >= 0; i--) {
            int b = KT::get(x, i), lb = KT::get(lz, i);
            if (b) { int c = t[v].nxt[lb]; if (c != -1) res += t[c].cnt; }
            v = t[v].nxt[b ^ lb];
            if (v == -1) return res;
        }
        return res + (incl ? t[v].cnt : 0);
    }
    Key xor_walk(const Key& x, int f) const {
        assert(size() > 0);
        Key y = x ^ lz;
        int v = 0; Key res = KT::zero();
        for (int i = LOG - 1; i >= 0; i--) {
            int b = KT::get(y, i) ^ f, c = t[v].nxt[b];
            if (c == -1) { b ^= 1; c = t[v].nxt[b]; }
            KT::set(res, i, b);
            v = c;
        }
        return res ^ lz;
    }
};
$0
