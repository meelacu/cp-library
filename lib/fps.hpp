// @name FPS (NTT / formal power series)
// @prefix fps
// @description mod 998244353 の形式的冪級数ライブラリ（NTT, inv/log/exp/pow/sqrt, taylor_shift, Bostan-Mori）

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

// ---------------- NTT（rate テーブル + DIF/DIT + radix-4）----------------
// 注意: ntt(a, false) の出力はビット反転順、ntt(a, true) はビット反転順を入力に取る。
// 畳み込み（点ごとの積を挟む用途）専用。自然順を仮定して直接呼ばないこと。
namespace ntt_internal {
constexpr int ctz_const(unsigned x) { int c = 0; while (!(x & 1u)) x >>= 1, c++; return c; }
constexpr int RANK2 = ctz_const(MOD - 1);
struct Tables {
    mint root[RANK2 + 1], iroot[RANK2 + 1];
    mint rate2[RANK2 - 1], irate2[RANK2 - 1];
    mint rate3[RANK2 - 2], irate3[RANK2 - 2];
    Tables() {
        root[RANK2] = mint(3).pow((MOD - 1) >> RANK2);
        iroot[RANK2] = root[RANK2].inv();
        for (int i = RANK2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
        { mint p = 1, ip = 1;
          for (int i = 0; i <= RANK2 - 2; i++) {
              rate2[i] = root[i + 2] * p; irate2[i] = iroot[i + 2] * ip;
              p *= iroot[i + 2]; ip *= root[i + 2];
          } }
        { mint p = 1, ip = 1;
          for (int i = 0; i <= RANK2 - 3; i++) {
              rate3[i] = root[i + 3] * p; irate3[i] = iroot[i + 3] * ip;
              p *= iroot[i + 3]; ip *= root[i + 3];
          } }
    }
};
inline const Tables& tab() { static const Tables t; return t; }
}  // namespace ntt_internal

inline void ntt(vector<mint>& a, bool inverse) {
    int n = (int)a.size();
    if (n <= 1) return;
    const auto& T = ntt_internal::tab();
    int h = __builtin_ctz((unsigned)n);
    if (!inverse) {
        for (int len = 0; len < h;) {
            if (h - len == 1) {
                int p = 1 << (h - len - 1);
                mint rot = 1;
                for (int s = 0; s < (1 << len); s++) {
                    int off = s << (h - len);
                    for (int i = 0; i < p; i++) {
                        mint l = a[i + off], r = a[i + off + p] * rot;
                        a[i + off] = l + r;
                        a[i + off + p] = l - r;
                    }
                    if (s + 1 != (1 << len)) rot *= T.rate2[__builtin_ctz(~(unsigned)s)];
                }
                len++;
            } else {
                int p = 1 << (h - len - 2);
                mint rot = 1, imag = T.root[2];
                for (int s = 0; s < (1 << len); s++) {
                    mint rot2 = rot * rot, rot3 = rot2 * rot;
                    int off = s << (h - len);
                    for (int i = 0; i < p; i++) {
                        mint a0 = a[i + off], a1 = a[i + off + p] * rot;
                        mint a2 = a[i + off + 2 * p] * rot2, a3 = a[i + off + 3 * p] * rot3;
                        mint a1na3imag = (a1 - a3) * imag;
                        mint a0a2 = a0 + a2, a0na2 = a0 - a2, a1a3 = a1 + a3;
                        a[i + off]         = a0a2 + a1a3;
                        a[i + off + p]     = a0a2 - a1a3;
                        a[i + off + 2 * p] = a0na2 + a1na3imag;
                        a[i + off + 3 * p] = a0na2 - a1na3imag;
                    }
                    if (s + 1 != (1 << len)) rot *= T.rate3[__builtin_ctz(~(unsigned)s)];
                }
                len += 2;
            }
        }
    } else {
        for (int len = h; len;) {
            if (len == 1) {
                int p = 1 << (h - len);
                mint irot = 1;
                for (int s = 0; s < (1 << (len - 1)); s++) {
                    int off = s << (h - len + 1);
                    for (int i = 0; i < p; i++) {
                        mint l = a[i + off], r = a[i + off + p];
                        a[i + off] = l + r;
                        a[i + off + p] = (l - r) * irot;
                    }
                    if (s + 1 != (1 << (len - 1))) irot *= T.irate2[__builtin_ctz(~(unsigned)s)];
                }
                len--;
            } else {
                int p = 1 << (h - len);
                mint irot = 1, iimag = T.iroot[2];
                for (int s = 0; s < (1 << (len - 2)); s++) {
                    mint irot2 = irot * irot, irot3 = irot2 * irot;
                    int off = s << (h - len + 2);
                    for (int i = 0; i < p; i++) {
                        mint a0 = a[i + off], a1 = a[i + off + p];
                        mint a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                        mint a2na3iimag = (a2 - a3) * iimag;
                        a[i + off]         = a0 + a1 + a2 + a3;
                        a[i + off + p]     = (a0 - a1 + a2na3iimag) * irot;
                        a[i + off + 2 * p] = (a0 + a1 - a2 - a3) * irot2;
                        a[i + off + 3 * p] = (a0 - a1 - a2na3iimag) * irot3;
                    }
                    if (s + 1 != (1 << (len - 2))) irot *= T.irate3[__builtin_ctz(~(unsigned)s)];
                }
                len -= 2;
            }
        }
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
        if (rs == n / 2 + 1 && n >= 8) {
            // 結果長が 2^k+1 のとき NTT 長が倍に飛ぶのを回避。
            // 長さ m=2^k の巡回畳み込みを取り、唯一回り込む x^m の係数を後で補正する。
            int m = n / 2;
            mint top = back() * g.back();
            FPS a = pre(m), b = g.pre(m);
            if ((int)size() > m) a[0] += (*this)[m];
            if ((int)g.size() > m) b[0] += g[m];
            ntt(a, false); ntt(b, false);
            for (int i = 0; i < m; i++) a[i] *= b[i];
            ntt(a, true);
            a[0] -= top;
            a.push_back(top);
            return *this = a;
        }
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

// [x^k] P(x)/Q(x) を Bostan-Mori 法で計算 O(M log M log k)
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
