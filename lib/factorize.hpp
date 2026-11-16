// @name Factorize / Prime (64bit)
// @prefix factorize
// @prefix isprime
// @prefix divisors
// @description 64bit整数の素数判定(Miller-Rabin)・素因数分解(Pollard rho)・約数列挙・φ・原始根、n以下の素数列挙

// 64bit 整数 1 個を扱う整数論。main の外に置く
//   bool is_prime(n)                     決定的 Miller-Rabin（long long 全域）。1 個 2µs 以下
//   vector<pair<ll, int>> factorize(n)   素因数分解 {素数, 指数}（素数の昇順）。n >= 1、factorize(1) は空
//                                        Pollard rho。1e18 以下のランダムな数で平均 10µs、最悪（31bit 素数 2 個の積）で 0.4ms 程度
//   vector<ll> prime_factors(n)          相異なる素因数（昇順）
//   vector<ll> divisors(n)               約数（昇順）
//   vector<ll> divisors(pf)              素因数分解 pf = {{素数, 指数}, ...} から約数を列挙（昇順）。分解済みの結果の使い回しや、
//                                        spf スニペットの factorize の結果、いくつかの数の積（指数を足し合わせたもの）にも使える
//   ll euler_phi(n)                      オイラーの φ(n) = 1 以上 n 以下で n と互いに素な数の個数
//   ll primitive_root(p)                 素数 p の最小の原始根
//   vector<int> primes_upto(n)           n 以下の素数の昇順（エラトステネスの篩。n = 1e8 で 0.4 秒）
// 10^7 以下の数を大量に分解・判定するなら spf スニペットや primes_upto のほうが速い
// 例: 約数の個数  ll c = 1; for (auto [p, e] : factorize(n)) c *= e + 1;
namespace factorize_internal {
using u64 = unsigned long long;
using u128 = unsigned __int128;
// 奇数 n を法とする Montgomery 乗算（R = 2^64）。値は x*R mod n の形で持つ
struct Mont {
    u64 n, ni, r2;  // ni = n^{-1} mod 2^64, r2 = R^2 mod n
    explicit Mont(u64 n_) : n(n_), ni(n_) {
        for (int i = 0; i < 5; i++) ni *= 2 - n * ni;
        u64 r1 = (u64)(-n) % n;
        r2 = (u128)r1 * r1 % n;
    }
    u64 reduce(u128 t) const {  // t * R^{-1} mod n（t < n * R）
        u64 lo = (u64)t, hi = (u64)(t >> 64);
        u64 mh = (u64)(((u128)(lo * ni) * n) >> 64);
        return hi >= mh ? hi - mh : hi - mh + n;
    }
    u64 mul(u64 a, u64 b) const { return reduce((u128)a * b); }
    u64 to(u64 a) const { return mul(a % n, r2); }
    u64 add(u64 a, u64 b) const {
        u64 c = a + b;
        return (c >= n || c < a) ? c - n : c;
    }
    u64 pow(u64 a, u64 e) const {
        u64 r = to(1);
        for (; e; e >>= 1, a = mul(a, a))
            if (e & 1) r = mul(r, a);
        return r;
    }
};
inline bool is_prime_u64(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
        if (n % p == 0) return n == p;
    if (n < 41 * 41) return true;
    Mont m(n);
    u64 d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;
    const u64 one = m.to(1), mone = m.to(n - 1);
    auto check = [&](u64 a) {
        u64 x = a % n;
        if (x == 0) return true;
        x = m.pow(m.to(x), d);
        if (x == one || x == mone) return true;
        for (int i = 1; i < s; i++) {
            x = m.mul(x, x);
            if (x == mone) return true;
        }
        return false;
    };
    if (n < 4759123141ULL) {
        for (u64 a : {2, 7, 61})
            if (!check(a)) return false;
    } else {
        for (u64 a : {2, 325, 9375, 28178, 450775, 9780504, 1795265022})
            if (!check(a)) return false;
    }
    return true;
}
// n は奇数の合成数。1 < g < n の約数 g を返す（Brent 版 Pollard rho、gcd はまとめて取る）
inline u64 pollard(u64 n) {
    Mont m(n);
    const u64 one = m.to(1);
    auto absdiff = [](u64 a, u64 b) { return a > b ? a - b : b - a; };
    for (u64 c0 = 1;; c0++) {
        const u64 c = m.to(c0);
        auto f = [&](u64 v) { return m.add(m.mul(v, v), c); };
        u64 x = 0, y = m.to(2), ys = 0, q = one, g = 1;
        const u64 M = 128;
        for (u64 r = 1; g == 1; r <<= 1) {
            x = y;
            for (u64 i = 0; i < r; i++) y = f(y);
            for (u64 k = 0; k < r && g == 1; k += M) {
                ys = y;
                for (u64 i = 0; i < M && i < r - k; i++) {
                    y = f(y);
                    q = m.mul(q, absdiff(x, y));
                }
                g = std::gcd(q, n);
            }
        }
        if (g == n) {
            do {
                ys = f(ys);
                g = std::gcd(absdiff(x, ys), n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}
inline void factor_rec(u64 n, vector<u64>& out) {
    if (n == 1) return;
    if (is_prime_u64(n)) {
        out.push_back(n);
        return;
    }
    u64 s = sqrtl((long double)n);  // 平方数（p^2 など）は rho に回すと遅いので先に割る（n < 2^63 なので (s+1)^2 は溢れない）
    while (s * s > n) s--;
    while ((s + 1) * (s + 1) <= n) s++;
    if (s * s == n) {
        factor_rec(s, out);
        factor_rec(s, out);
        return;
    }
    u64 d = pollard(n);
    factor_rec(d, out);
    factor_rec(n / d, out);
}
}  // namespace factorize_internal

inline bool is_prime(long long n) { return n >= 2 && factorize_internal::is_prime_u64(n); }

inline vector<pair<long long, int>> factorize(long long n) {
    assert(n >= 1);
    using factorize_internal::u64;
    u64 x = n;
    vector<u64> ps;
    for (u64 p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47})
        while (x % p == 0) ps.push_back(p), x /= p;
    factorize_internal::factor_rec(x, ps);
    sort(ps.begin(), ps.end());
    vector<pair<long long, int>> res;
    for (u64 p : ps) {
        if (!res.empty() && (u64)res.back().first == p) res.back().second++;
        else res.push_back({(long long)p, 1});
    }
    return res;
}

inline vector<long long> prime_factors(long long n) {
    vector<long long> res;
    for (auto [p, e] : factorize(n)) res.push_back(p);
    return res;
}

template <class P, class E>
vector<long long> divisors(const vector<pair<P, E>>& pf) {
    vector<long long> res{1};
    for (auto [p, e] : pf) {
        int sz = res.size();
        long long pw = 1;
        for (E i = 0; i < e; i++) {
            pw *= p;
            for (int j = 0; j < sz; j++) res.push_back(res[j] * pw);
        }
    }
    sort(res.begin(), res.end());
    return res;
}
inline vector<long long> divisors(long long n) { return divisors(factorize(n)); }

inline long long euler_phi(long long n) {
    long long r = n;
    for (auto [p, e] : factorize(n)) r = r / p * (p - 1);
    return r;
}

inline long long primitive_root(long long p) {
    assert(is_prime(p));
    if (p == 2) return 1;
    factorize_internal::Mont m(p);
    auto fs = factorize(p - 1);
    const unsigned long long one = m.to(1);
    for (long long g = 2;; g++) {
        unsigned long long gm = m.to(g);
        bool ok = true;
        for (auto [q, e] : fs)
            if (m.pow(gm, (p - 1) / q) == one) {
                ok = false;
                break;
            }
        if (ok) return g;
    }
}

inline vector<int> primes_upto(int n) {
    vector<int> ps;
    if (n < 2) return ps;
    ps.push_back(2);
    vector<char> comp(n / 2 + 1);  // comp[i] : 2i+1 が合成数
    for (int i = 1; 2 * i + 1 <= n; i++) {
        if (comp[i]) continue;
        int p = 2 * i + 1;
        ps.push_back(p);
        for (long long j = (long long)p * p / 2; j <= n / 2; j += p) comp[j] = 1;
    }
    return ps;
}
$0
