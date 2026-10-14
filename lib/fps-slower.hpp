// @name FPS (Formal Power Series, mod 998244353)
// @prefix fps-slower
// @description FPS: NTT畳み込み / inv / log / exp / pow / sqrt / taylor_shift / sparse mul-div / Bostan-Mori (kth_term)。依存: bits/stdc++.h と using namespace std;

// ========== Formal Power Series (mod 998244353) ==========
namespace fps {

constexpr unsigned MOD = 998244353;
struct mint {
    unsigned v = 0;
    mint() {}
    mint(long long x) { x %= (long long)MOD; if (x < 0) x += MOD; v = (unsigned)x; }
    mint& operator+=(mint r) { v += r.v; if (v >= MOD) v -= MOD; return *this; }
    mint& operator-=(mint r) { v += MOD - r.v; if (v >= MOD) v -= MOD; return *this; }
    mint& operator*=(mint r) { v = (unsigned)((unsigned long long)v * r.v % MOD); return *this; }
    mint pow(long long e) const { mint a = *this, r = 1; for (; e > 0; e >>= 1, a *= a) if (e & 1) r *= a; return r; }
    mint inv() const { return pow(MOD - 2); }
    mint& operator/=(mint r) { return *this *= r.inv(); }
    mint operator-() const { return mint() - *this; }
    friend mint operator+(mint a, mint b) { return a += b; }
    friend mint operator-(mint a, mint b) { return a -= b; }
    friend mint operator*(mint a, mint b) { return a *= b; }
    friend mint operator/(mint a, mint b) { return a /= b; }
    friend bool operator==(mint a, mint b) { return a.v == b.v; }
    friend bool operator!=(mint a, mint b) { return a.v != b.v; }
    friend ostream& operator<<(ostream& os, mint m) { return os << m.v; }
    friend istream& operator>>(istream& is, mint& m) { long long x; is >> x; m = mint(x); return is; }
};

// x^2 = a (mod MOD) の解を返す。非存在なら -1（Tonelli-Shanks）
inline long long mod_sqrt(mint a) {
    if (a.v == 0) return 0;
    if (a.pow((MOD - 1) / 2).v != 1) return -1;
    mint b = 2;
    while (b.pow((MOD - 1) / 2).v == 1) b += 1;
    unsigned long long e = 0, m = MOD - 1;
    while (m % 2 == 0) m >>= 1, e++;
    mint x = a.pow((m - 1) / 2), y = a * x * x, z = b.pow(m);
    x *= a;
    while (y.v != 1) {
        unsigned long long j = 0;
        mint t = y;
        while (t.v != 1) t *= t, j++;
        z = z.pow(1ULL << (e - j - 1));
        x *= z; z *= z; y *= z; e = j;
    }
    return x.v;
}

inline void ntt(vector<mint>& a, bool inverse) {
    int n = (int)a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (int len = 1; len < n; len <<= 1) {
        mint wlen = mint(3).pow((MOD - 1) / (2 * len));
        if (inverse) wlen = wlen.inv();
        for (int i = 0; i < n; i += 2 * len) {
            mint w = 1;
            for (int j = 0; j < len; j++) {
                mint u = a[i + j], v = a[i + j + len] * w;
                a[i + j] = u + v;
                a[i + j + len] = u - v;
                w *= wlen;
            }
        }
    }
    if (inverse) {
        mint ninv = mint(n).inv();
        for (auto& x : a) x *= ninv;
    }
}

struct FPS : vector<mint> {
    using vector<mint>::vector;

    // 先頭 n 項に切り詰め（不足分は 0 埋め）
    FPS pre(int n) const {
        FPS r(begin(), begin() + min((int)size(), n));
        r.resize(n);
        return r;
    }
    FPS operator-() const { FPS r = *this; for (auto& x : r) x = -x; return r; }
    FPS& operator+=(const FPS& g) {
        if (g.size() > size()) resize(g.size());
        for (int i = 0; i < (int)g.size(); i++) (*this)[i] += g[i];
        return *this;
    }
    FPS& operator-=(const FPS& g) {
        if (g.size() > size()) resize(g.size());
        for (int i = 0; i < (int)g.size(); i++) (*this)[i] -= g[i];
        return *this;
    }
    FPS& operator*=(mint c) { for (auto& x : *this) x *= c; return *this; }
    FPS& operator/=(mint c) { return *this *= c.inv(); }
    FPS& operator*=(const FPS& g) {
        if (empty() || g.empty()) { clear(); return *this; }
        int rs = (int)(size() + g.size()) - 1;
        if (min(size(), g.size()) <= 60) {  // 小サイズは愚直
            FPS r(rs);
            for (int i = 0; i < (int)size(); i++)
                for (int j = 0; j < (int)g.size(); j++) r[i + j] += (*this)[i] * g[j];
            return *this = r;
        }
        int n = 1;
        while (n < rs) n <<= 1;
        FPS a = pre(n), b = g.pre(n);
        ntt(a, false); ntt(b, false);
        for (int i = 0; i < n; i++) a[i] *= b[i];
        ntt(a, true);
        return *this = a.pre(rs);
    }
    // x^k 倍 / x^k で割る（切り捨て）
    FPS operator<<(int k) const { FPS r(k); r.insert(r.end(), begin(), end()); return r; }
    FPS operator>>(int k) const { return k >= (int)size() ? FPS{} : FPS(begin() + k, end()); }

    friend FPS operator+(FPS f, const FPS& g) { return f += g; }
    friend FPS operator-(FPS f, const FPS& g) { return f -= g; }
    friend FPS operator*(FPS f, const FPS& g) { return f *= g; }
    friend FPS operator*(FPS f, mint c) { return f *= c; }
    friend FPS operator*(mint c, FPS f) { return f *= c; }
    friend FPS operator/(FPS f, mint c) { return f /= c; }

    FPS diff() const {
        FPS r(max(0, (int)size() - 1));
        for (int i = 1; i < (int)size(); i++) r[i - 1] = (*this)[i] * i;
        return r;
    }
    FPS integ() const {
        int n = (int)size();
        FPS r(n + 1);
        vector<mint> iv(n + 2);
        if (n + 2 > 1) iv[1] = 1;
        for (int i = 2; i <= n + 1; i++) iv[i] = -mint(MOD / i) * iv[MOD % i];
        for (int i = 0; i < n; i++) r[i + 1] = (*this)[i] * iv[i + 1];
        return r;
    }
    // 1/f mod x^deg（f[0] != 0 が必要）
    FPS inv(int deg = -1) const {
        assert(!empty() && (*this)[0].v != 0);
        if (deg < 0) deg = (int)size();
        FPS r{(*this)[0].inv()};
        for (int d = 1; d < deg; d <<= 1)
            r = (r * (FPS{2} - (pre(2 * d) * r).pre(2 * d))).pre(2 * d);
        return r.pre(deg);
    }
    // log f mod x^deg（f[0] == 1 が必要）
    FPS log(int deg = -1) const {
        assert(!empty() && (*this)[0].v == 1);
        if (deg < 0) deg = (int)size();
        if (deg == 0) return {};
        return (diff() * inv(deg)).pre(max(0, deg - 1)).integ().pre(deg);
    }
    // exp f mod x^deg（f[0] == 0 が必要）
    FPS exp(int deg = -1) const {
        assert(empty() || (*this)[0].v == 0);
        if (deg < 0) deg = (int)size();
        if (deg == 0) return {};
        FPS r{1};
        for (int d = 1; d < deg; d <<= 1)
            r = (r * (pre(2 * d) + FPS{1} - r.log(2 * d))).pre(2 * d);
        return r.pre(deg);
    }
    // f^k mod x^deg（先頭が 0 でもよい）
    FPS pow(long long k, int deg = -1) const {
        if (deg < 0) deg = (int)size();
        if (deg == 0) return {};
        if (k == 0) { FPS r(deg); r[0] = 1; return r; }
        int p = 0;
        while (p < (int)size() && (*this)[p].v == 0) p++;
        if (p == (int)size() || p > (deg - 1) / k) return FPS(deg);  // すべて 0
        mint c = (*this)[p];
        int rest = deg - (int)(p * k);
        FPS f = (*this >> p) / c;
        FPS r = (f.log(rest) * mint(k % (long long)MOD)).exp(rest) * c.pow(k);
        return (r << (int)(p * k)).pre(deg);
    }
    // sqrt f mod x^deg。存在しなければ空を返す
    FPS sqrt(int deg = -1) const {
        if (deg < 0) deg = (int)size();
        int p = 0;
        while (p < (int)size() && (*this)[p].v == 0) p++;
        if (p == (int)size()) return FPS(deg);  // f = 0
        if (p % 2 == 1) return {};
        long long s = mod_sqrt((*this)[p]);
        if (s < 0) return {};
        int rest = deg - p / 2;
        if (rest <= 0) return FPS(deg);
        FPS f = (*this >> p) / (*this)[p];
        FPS r{1};
        mint i2 = mint(2).inv();
        for (int d = 1; d < rest; d <<= 1)
            r = (r + f.pre(2 * d) * r.inv(2 * d)).pre(2 * d) * i2;
        r = r.pre(rest) * mint(s);
        return (r << (p / 2)).pre(deg);
    }
    // f(x + c)
    FPS taylor_shift(mint c) const {
        int n = (int)size();
        if (n == 0) return {};
        vector<mint> fact(n), finv(n);
        fact[0] = 1;
        for (int i = 1; i < n; i++) fact[i] = fact[i - 1] * i;
        finv[n - 1] = fact[n - 1].inv();
        for (int i = n - 1; i > 0; i--) finv[i - 1] = finv[i] * i;
        FPS a(n), b(n);
        for (int i = 0; i < n; i++) a[n - 1 - i] = (*this)[i] * fact[i];
        mint pw = 1;
        for (int i = 0; i < n; i++) b[i] = pw * finv[i], pw *= c;
        FPS r = (a * b).pre(n);
        reverse(r.begin(), r.end());
        for (int i = 0; i < n; i++) r[i] *= finv[i];
        return r;
    }
    mint eval(mint x) const {
        mint r = 0;
        for (int i = (int)size() - 1; i >= 0; i--) r = r * x + (*this)[i];
        return r;
    }
    // 疎な g = Σ c_e x^e（(e, c) のリスト）との乗算 O(N·K)
    FPS mul_sparse(const vector<pair<int, mint>>& g, int deg = -1) const {
        if (deg < 0) deg = (int)size();
        FPS r(deg);
        for (auto [e, c] : g)
            for (int i = 0; i + e < deg && i < (int)size(); i++) r[i + e] += (*this)[i] * c;
        return r;
    }
    // 疎な g での除算 O(N·K)（g の定数項が必要）
    FPS div_sparse(vector<pair<int, mint>> g, int deg = -1) const {
        if (deg < 0) deg = (int)size();
        mint c0;
        bool has0 = false;
        vector<pair<int, mint>> h;
        for (auto [e, c] : g) {
            if (e == 0) c0 += c, has0 = true;
            else h.push_back({e, c});
        }
        assert(has0 && c0.v != 0);
        mint ic = c0.inv();
        FPS r(deg);
        for (int i = 0; i < deg; i++) {
            mint s = i < (int)size() ? (*this)[i] : mint(0);
            for (auto [e, c] : h)
                if (i - e >= 0) s -= r[i - e] * c;
            r[i] = s * ic;
        }
        return r;
    }
};

// [x^k] P(x)/Q(x) を Bostan–Mori 法で計算 O(M log M log k)
inline mint bostan_mori(FPS P, FPS Q, long long k) {
    assert(!Q.empty() && Q[0].v != 0);
    while (k > 0) {
        FPS Qm = Q;
        for (int i = 1; i < (int)Qm.size(); i += 2) Qm[i] = -Qm[i];
        FPS U = P * Qm, V = Q * Qm;
        P.clear(); Q.clear();
        for (int i = (int)(k & 1); i < (int)U.size(); i += 2) P.push_back(U[i]);
        for (int i = 0; i < (int)V.size(); i += 2) Q.push_back(V[i]);
        k >>= 1;
    }
    return P.empty() ? mint(0) : P[0] / Q[0];
}

// 線形漸化式 a[n] = Σ c[i] a[n-1-i]（初項 a[0..d-1]）の第 k 項
inline mint kth_term(const vector<mint>& a, const vector<mint>& c, long long k) {
    int d = (int)c.size();
    FPS Q(d + 1);
    Q[0] = 1;
    for (int i = 0; i < d; i++) Q[i + 1] = -c[i];
    FPS A(a.begin(), a.begin() + min(a.size(), (size_t)d));
    FPS P = (A * Q).pre(d);
    return bostan_mori(P, Q, k);
}

}  // namespace fps
using fps::FPS;
using fps::mint;
$0
