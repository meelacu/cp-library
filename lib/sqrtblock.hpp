// @name SqrtDecomp Block skeleton
// @prefix sqrtblock
// @description 平方分割 Block 雛形（サマリ / lazy / whole / part を埋める）

struct ${1:MyBlock} {
	vector<long long> raw;   // 生値
	${2:// サマリ: ソート列 / 頻度表 / 部分和 など}
	${3:// lazy : 一様加算 / 一律代入 など（サマリと両立するか確認）}

	explicit $1(vector<long long> v) : raw(std::move(v)) { rebuild(); }
	void rebuild() { 
        ${4:/* raw からサマリ再構築 O(B)〜O(B log B) */}
    }
	void push()    { 
        ${5:/* lazy を raw・サマリへ反映して消す O(B) */} 
    }

	// ---- whole 系: 更新は lazy へ O(1) / クエリはサマリ+lazy で O(1)〜O(log B) ----
	$6
	// ---- part 系: 更新は push→素朴更新→rebuild / クエリは lazy 読み替えで O(B) ----
	$0
};
