// @name Matrix semiring: bool
// @prefix matbool
// @description 行列用の半環: 到達可能性（OR-AND）+ using BMat

template <class SR, int FixN> struct Matrix;  // 本体は matrix スニペット（前後どちらに置いてもよい）
struct OrAnd {  // 到達可能性（0 以外は true）
	using T = int;
	static T zero() { return 0; }
	static T one() { return 1; }
	static T norm(T x) { return x != 0; }
	static T add(T a, T b) { return a | b; }
	static T mul(T a, T b) { return a & b; }
};
// A.pow(K)[s][t]: ちょうど K 歩で行けるか, A.pow_sum(K + 1)[s][t]: K 歩以内で行けるか
using BMat = Matrix<OrAnd, ${1:0}>;
$0
