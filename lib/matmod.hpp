// @name Matrix semiring: mod
// @prefix matmod
// @prefix matrix_pow
// @prefix matarr
// @description 行列用の半環: mod（コンパイル時 / 実行時）+ using Mat

template <class SR, int FixN> struct Matrix;  // 本体は matrix スニペット（前後どちらに置いてもよい）
template <long long MOD>  // MOD = 0 なら実行時 mod: ModRing<0>::set_mod(m)
struct ModRing {
	static_assert(0 <= MOD && MOD < (1LL << 62));
	using T = long long;
	static inline long long dyn = 0;
	static void set_mod(long long m) { assert(1 <= m && m < (1LL << 62)); dyn = m; }
	static long long mod() {
		if constexpr (MOD > 0) return MOD;
		else { assert(dyn > 0); return dyn; }  // set_mod 忘れ
	}
	static T zero() { return 0; }
	static T one() { return 1 % mod(); }
	static T norm(long long x) { const long long M = mod(); return 0 <= x && x < M ? x : (x % M + M) % M; }
	static T add(T a, T b) { return (a += b) >= mod() ? a - mod() : a; }
	static T sub(T a, T b) { return (a -= b) < 0 ? a + mod() : a; }
	static T mul(T a, T b) { const long long M = mod(); return M < (1LL << 31) ? a * b % M : (T)((__int128)a * b % M); }
	static T inv(T a) {  // gcd(a, mod) = 1 が前提
		long long b = mod(), u = 1, v = 0;
		while (b) { long long t = a / b; a -= t * b; swap(a, b); u -= t * v; swap(u, v); }
		assert(a == 1);
		return norm(u);
	}
};
// MOD = 0 なら先に ModRing<0>::set_mod(m)。det/inv/rank は MOD が素数のときのみ。GF(2) は ModRing<2>
using Mat = Matrix<ModRing<${1|998244353,1000000007,0|}>, ${2:0}>;
$0
