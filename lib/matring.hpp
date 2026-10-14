// @name Matrix semiring: ring
// @prefix matring
// @description 行列用の半環: 型の演算そのまま（double / 自作 modint など）+ using RMat

template <class SR, int FixN> struct Matrix;  // 本体は matrix スニペット（前後どちらに置いてもよい）
template <class Val>  // Val の演算をそのまま使う（自作 modint, double, 素の long long など）
struct Ring {
	using T = Val;
	static T zero() { return T(0); }
	static T one() { return T(1); }
	static T norm(const T& x) { return x; }
	static T add(const T& a, const T& b) { return a + b; }
	static T sub(const T& a, const T& b) { return a - b; }
	static T mul(const T& a, const T& b) { return a * b; }
	static T inv(const T& a) { return T(1) / a; }
};
// double なら det/inv/rank は部分ピボット選択付き。自作 modint は Ring<mint>（== と / が要る）
using RMat = Matrix<Ring<${1|double,long double,mint,long long|}>, ${2:0}>;
$0
