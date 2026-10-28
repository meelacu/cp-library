// @name SqrtDecomp Block skeleton
// @prefix sqrtblock
// @description 平方分割 Block 雛形（サマリ / lazy / whole / part を埋める）
 
struct ${1:MyBlock} {
	vector<${2:long long}> raw;   // 生値（lazy 抜き）。push() 後は現在値
	${3:// サマリ: ソート列 / 頻度表 / 部分和 など}
	${4:// lazy : 一様加算 / 一律代入 など（サマリと両立するか確認）}
 
	explicit $1(vector<$2> v) : raw(std::move(v)) { rebuild(); }
	void rebuild() {
		${5:/* raw からサマリ再構築 O(B)〜O(B log B) */}
	}
	void push() {
		${6:/* lazy を raw・サマリへ反映して消す O(B)（ブロック長の自動調整でも使うので必ず実装） */}
	}
 
	// whole 版 op(args) と part 版 op(s, t, args) を同名で並べる（[s, t) はブロック内ローカル。
	// ブロックを丸ごと覆うときは必ず whole 版が呼ばれる）
	//   更新: whole は lazy へ O(1) / part は push → 素朴更新 → rebuild
	//   取得: whole はサマリ + lazy / part は raw + lazy を素朴に O(t - s)
	//   part でサマリを差分で直せるなら push / rebuild を省いて [s, t) だけ触る（O(B) → O(t - s)。端が軽いほど B を小さく取れて速い）
	//   1 点操作は sd.at(i, [&](auto& k, int j) { ... }) から（自動調整がその費用も測る）
	$0
};
 