// @name modint (compile-time / runtime mod)
// @prefix modint
// @description modint：ModInt<M>（コンパイル時 mod）と ModInt<0>（実行時 mod、set_mod）。mint / dmint の別名つき

// ==== modint ====
// ModInt<M>     : コンパイル時 mod（1 ≤ M < 2^63。M ≤ 2^31 なら 32bit で保持）
// ModInt<0, id> : 実行時 mod。ModInt<0, id>::set_mod(m)（1 ≤ m < 2^63、既定 998244353）
//                 id を変えると別々の実行時 mod を同時に持てる
// 生成   : mint(x)（int / ll / ull / __int128、負数可）, mint::raw(x)（0 ≤ x < mod のとき剰余を取らずに作る）
// 演算   : + - * / += -= *= /= ++ -- 単項- == != << >>（整数との混在可: a * 2, 1 - a, a == 0）
// val()  : 値（long long）  mod() : 法
// inv()  : 拡張ユークリッド。合成数 mod でも gcd(x, mod) = 1 なら求まる（それ以外は assert）
// pow(n) : n < 0 なら逆元の累乗
template <unsigned long long M, int ID = 0>
struct ModInt {
    using u32 = unsigned int;
    using u64 = unsigned long long;
    using u128 = unsigned __int128;
    static_assert(M < (1ULL << 63));
    static constexpr bool SMALL = M != 0 && M <= (1ULL << 31);
    using U = std::conditional_t<SMALL, u32, u64>;
    inline static u64 dmod = 998244353, dinv = ~0ULL / 998244353 + 1;  // 実行時 mod と Barrett 用の ceil(2^64 / mod)
    U v = 0;

    static constexpr u64 mod() requires (M != 0) { return M; }
    static u64 mod() requires (M == 0) { return dmod; }
    static void set_mod(u64 m) requires (M == 0) { assert(1 <= m && m < (1ULL << 63)); dmod = m, dinv = ~0ULL / m + 1; }

    ModInt() = default;
    template <std::integral T>
    ModInt(T x) {
        if constexpr (std::is_signed_v<T>) {
            using S = std::conditional_t<(sizeof(T) > 8), __int128, long long>;
            S y = S(x) % S(mod());
            v = U(y < 0 ? y + S(mod()) : y);
        } else {
            using W = std::conditional_t<(sizeof(T) > 8), u128, u64>;
            v = U(W(x) % mod());
        }
    }
    static ModInt raw(u64 x) { ModInt r; r.v = U(x); return r; }
    long long val() const { return (long long)v; }

    ModInt& operator+=(const ModInt& r) { v += r.v; if (v >= mod()) v -= U(mod()); return *this; }
    ModInt& operator-=(const ModInt& r) { if (v < r.v) v += U(mod()); v -= r.v; return *this; }
    ModInt& operator*=(const ModInt& r) {
        if constexpr (M != 0) {  // コンパイル時 mod
            if constexpr (M < (1ULL << 32)) v = U((u64)v * r.v % M);
            else v = U((u128)v * r.v % M);
        } else if (dmod < (1ULL << 32)) {  // 実行時 mod（32bit 以内）：Barrett 還元
            u64 z = v * r.v, y = u64((u128)z * dinv >> 64) * dmod;
            v = z - y + (z < y ? dmod : 0);
        } else {
            v = U((u128)v * r.v % dmod);
        }
        return *this;
    }
    ModInt& operator/=(const ModInt& r) { return *this *= r.inv(); }
    ModInt& operator++() { if (++v == mod()) v = 0; return *this; }
    ModInt& operator--() { if (v == 0) v = U(mod()); --v; return *this; }
    ModInt operator++(int) { ModInt t = *this; ++*this; return t; }
    ModInt operator--(int) { ModInt t = *this; --*this; return t; }
    ModInt operator+() const { return *this; }
    ModInt operator-() const { return ModInt() - *this; }

    ModInt inv() const {
        long long a = (long long)v, b = (long long)mod(), x = 1, y = 0;
        while (b) { long long q = a / b; std::swap(a -= q * b, b); std::swap(x -= q * y, y); }
        assert(a == 1);  // gcd(val, mod) ≠ 1 なら逆元なし
        return ModInt(x);
    }
    ModInt pow(long long n) const {
        ModInt b = n < 0 ? inv() : *this, r = 1;
        for (u64 e = n < 0 ? 0 - (u64)n : (u64)n; e; e >>= 1, b *= b)
            if (e & 1) r *= b;
        return r;
    }

    friend ModInt operator+(ModInt a, const ModInt& b) { return a += b; }
    friend ModInt operator-(ModInt a, const ModInt& b) { return a -= b; }
    friend ModInt operator*(ModInt a, const ModInt& b) { return a *= b; }
    friend ModInt operator/(ModInt a, const ModInt& b) { return a /= b; }
    friend bool operator==(const ModInt&, const ModInt&) = default;
    friend std::ostream& operator<<(std::ostream& os, const ModInt& x) { return os << x.v; }
    friend std::istream& operator>>(std::istream& is, ModInt& x) { long long t; is >> t; x = ModInt(t); return is; }
};
using mint = ModInt<${1|998244353,1000000007|}>;
using dmint = ModInt<0>;  // dmint::set_mod(m);
$0
