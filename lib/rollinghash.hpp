// @name RollingHash (mod 2^61-1)
// @prefix rollinghash
// @description RollingHash (mod 2^61-1, 基数ランダム): get / hash_of / concat / repeat / lcp / cmp / push_back・pop_back / セグ木用モノイド S, op, e, leaf。依存: bits/stdc++.h と using namespace std;

// RollingHash: mod 2^61-1, 基数ランダム（全インスタンス共通なので別々の列のハッシュ同士も比較可）
//   string / vector<int> / vector<ll> など何でも可（要素は整数型, |x| < 2^60）。0-indexed・半開区間 [l, r)
//   各要素にランダムな shift を足しているので, 0 を含む列や長さの違う列を同じ set に入れても安全
// 構築      RollingHash rh(s);  /  RollingHash rh; rh.push_back(c); rh.pop_back();   (DFS しながらの経路ハッシュ等)
// 取得      rh.get(l, r), rh.get()            表なしで全体だけ: RollingHash::hash_of(t)
// 連結      concat(h1, h2, len2)              h1 の後ろに (h2, 長さ len2) を繋ぐ。1 文字は single(c)
// 繰り返し  repeat(h, len, k)                 (h, 長さ len) を k 回。O(log k), 全長は len * k
// LCP/比較  rh.lcp(i, j), rh.lcp(l1, r1, l2, r2), lcp(a, l1, r1, b, l2, r2), cmp(同) -> -1/0/1   各 O(log n)
//           接尾辞ソート: sort(all(id), [&](int i, int j) { return rh.cmp(i, n, j, n) < 0; });   O(n log^2 n)
// 回文      RollingHash rv(逆順の s);  s[l, r) が回文 <=> rh.get(l, r) == rv.get(n - r, n - l)
// 1 点更新  segtree<RH::S, RH::op, RH::e> seg(n); seg.set(i, RH::leaf(c));
//           seg.prod(l, r).h == rh.get(l, r) と互換。.r は逆順ハッシュなので .h == .r <=> 回文
struct RollingHash {
    using u64 = unsigned long long;
    using u128 = __uint128_t;
    static constexpr u64 MOD = (1ULL << 61) - 1;

    static u64 add(u64 a, u64 b) { return (a += b) >= MOD ? a - MOD : a; }
    static u64 sub(u64 a, u64 b) { return a >= b ? a - b : a + MOD - b; }
    static u64 mul(u64 a, u64 b) {
        u128 t = (u128)a * b;
        return add((u64)(t >> 61), (u64)(t & MOD));
    }
    static u64 rnd() {
        static mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
        return rng() % (MOD - 2) + 2;
    }
    static inline const u64 base = rnd(), shift = rnd();
    static inline vector<u64> pw{1};

    // base^n。表にあれば O(1), なければ O(log n)（表は伸ばさない）
    static u64 power(long long n) {
        if (n < (long long)pw.size()) return pw[n];
        u64 r = 1, b = base;
        for (; n; n >>= 1, b = mul(b, b))
            if (n & 1) r = mul(r, b);
        return r;
    }
    static void expand(size_t n) {  // 表を base^n まで伸ばす（reserve しない: 1 個ずつ伸ばしても償却 O(1)）
        while (pw.size() <= n) pw.push_back(mul(pw.back(), base));
    }
    // 要素 -> [0, MOD) の値。長さ 1 の列のハッシュでもある
    template <class T> static u64 single(T x) {
        u64 v;
        if constexpr (is_signed_v<T>) {
            long long y = (long long)x % (long long)MOD;
            v = y < 0 ? y + MOD : y;
        } else {
            v = (u64)x % MOD;
        }
        return add(v, shift);
    }
    template <class C> static u64 hash_of(const C& s) {
        u64 r = 0;
        for (auto&& c : s) r = add(mul(r, base), single(c));
        return r;
    }
    static u64 hash_of(const char* s) { return hash_of(string(s)); }
    static u64 concat(u64 h1, u64 h2, long long len2) { return add(mul(h1, power(len2)), h2); }
    static u64 repeat(u64 h, long long len, long long k) {
        u64 p = power(len), r = 0;
        for (; k; k >>= 1, h = add(mul(h, p), h), p = mul(p, p))
            if (k & 1) r = add(mul(r, p), h);
        return r;
    }

    vector<u64> h;  // h[i] = s[0, i) のハッシュ
    RollingHash() : h(1, 0) {}
    template <class C> explicit RollingHash(const C& s) : h(1, 0) {
        h.reserve(std::size(s) + 1);
        for (auto&& c : s) push_back(c);
    }
    RollingHash(const char* s) : RollingHash(string(s)) {}

    int size() const { return (int)h.size() - 1; }
    template <class T> void push_back(T c) {
        h.push_back(add(mul(h.back(), base), single(c)));
        expand(h.size() - 1);
    }
    void pop_back() { assert(size() > 0), h.pop_back(); }

    u64 get() const { return h.back(); }
    u64 get(int l, int r) const {
        assert(0 <= l && l <= r && r <= size());
        return sub(h[r], mul(h[l], pw[r - l]));
    }
    // 元の要素 s[i] を復元（|s[i]| < 2^60 が前提）
    long long at(int i) const {
        u64 v = sub(get(i, i + 1), shift);
        return v > MOD / 2 ? (long long)v - (long long)MOD : (long long)v;
    }

    static int lcp(const RollingHash& a, int l1, int r1, const RollingHash& b, int l2, int r2) {
        int lo = 0, hi = min(r1 - l1, r2 - l2) + 1;
        while (hi - lo > 1) {
            int m = (lo + hi) / 2;
            (a.get(l1, l1 + m) == b.get(l2, l2 + m) ? lo : hi) = m;
        }
        return lo;
    }
    int lcp(int l1, int r1, int l2, int r2) const { return lcp(*this, l1, r1, *this, l2, r2); }
    int lcp(int i, int j) const { return lcp(i, size(), j, size()); }
    // 辞書順比較: a[l1, r1) < b[l2, r2) なら -1, 等しければ 0, 大きければ 1
    static int cmp(const RollingHash& a, int l1, int r1, const RollingHash& b, int l2, int r2) {
        int n1 = r1 - l1, n2 = r2 - l2, k = lcp(a, l1, r1, b, l2, r2);
        if (k == n1 || k == n2) return (n1 > n2) - (n1 < n2);
        return a.at(l1 + k) < b.at(l2 + k) ? -1 : 1;
    }
    int cmp(int l1, int r1, int l2, int r2) const { return cmp(*this, l1, r1, *this, l2, r2); }

    // セグ木用モノイド（h: 正順, r: 逆順, p: base^len）
    struct S {
        u64 h, r, p;
    };
    static S op(S a, S b) { return {add(mul(a.h, b.p), b.h), add(mul(b.r, a.p), a.r), mul(a.p, b.p)}; }
    static S e() { return {0, 0, 1}; }
    template <class T> static S leaf(T c) {
        u64 v = single(c);
        return {v, v, base};
    }
};
using RH = RollingHash;
$0
