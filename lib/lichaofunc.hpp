// @name Li Chao Function
// @prefix lichaofunc
// @prefix LiChaoFunc
// @description Li Chao Tree に載せる関数族 F の原形 (評価関数だけ書く)

// Li Chao に載せる関数族。追加されるものをメンバに、クエリで変わるものを X にする
// 要件: 任意の 2 つの差 f - g の符号が X の順 (Cmp 順) に高々 1 回しか変わらないこと
struct ${1:Func} {
	using X = ${2:ll};  // クエリ点の型 (順序は less<X>。自作型なら operator< を定義)
	${3:ll} ${4:a, b};  // 追加されるもののパラメータ
	int id = -1;        // 復元に使わなければ削除してよい
	$3 operator()(const X& t) const { return ${5:a * t + b}; }
};
$0
