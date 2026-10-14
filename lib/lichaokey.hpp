// @name Li Chao Key
// @prefix lichaokey
// @prefix LiChaoKey
// @description Li Chao Tree のクエリ点型 X の原形 (順序を書く)

// Li Chao のクエリ点 X を自作する場合の原形 (既定は分数 p/q, q >= 0 で q = 0 は +INF)
// 順序は「差の符号変化が高々 1 回」になるように選ぶ。同値な点は関数の順序を変えないこと
struct ${1:Key} {
	${2:ll p, q};
	bool operator<(const $1& o) const { return ${3:p * o.q < o.p * q}; }
};
$0
