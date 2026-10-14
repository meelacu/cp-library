// @name Matrix semiring: tropical
// @prefix mattrop
// @description 行列用の半環: トロピカル（min-plus / max-plus）+ using TMat

template <class SR, int FixN> struct Matrix;  // 本体は matrix スニペット（前後どちらに置いてもよい）
template <bool MAX>  // false: min-plus（最短）, true: max-plus（最長）
struct Tropical {
	using T = long long;
	static constexpr T INF = (T)3e18;  // |x| >= INF は ±∞ 扱い
	static T zero() { return MAX ? -INF : INF; }  // 到達不能
	static T one() { return 0; }
	static T norm(T x) { return x >= INF ? INF : x <= -INF ? -INF : x; }
	static T add(T a, T b) { return MAX ? max(a, b) : min(a, b); }
	static T mul(T a, T b) { return a == zero() || b == zero() ? zero() : a + b; }  // 範囲外は行列側で norm
};
using MinPlus = Tropical<false>;
using MaxPlus = Tropical<true>;
// 零行列は全 INF（到達不能）。判定は x == MinPlus::INF / x == MaxPlus::zero()。K 辺以下は A.pow_sum(K + 1)
using TMat = Matrix<${1|MinPlus,MaxPlus|}, ${2:0}>;
$0
