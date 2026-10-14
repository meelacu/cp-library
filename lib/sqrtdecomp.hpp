// @name SqrtDecomp (framework)
// @prefix sqrtdecomp
// @description 平方分割 汎用フレームワーク（whole/part ディスパッチ + Block 差し替え）

// ===== 平方分割 汎用フレームワーク =====
// run(l, r, whole, part): 区間 [l, r) を分解してディスパッチ
//   whole(b)      : ブロック b 全体への操作（lazy / サマリで O(1)〜O(log B)）
//   part(b, s, t) : ブロック b 内ローカル [s, t) への操作（push→素朴更新→rebuild）
// Block の規約: Block(vector<T> raw) で構築。操作はラムダ注入なのでシグネチャ自由。
template <class Block>
struct SqrtDecomp {
	int n, B, nb;
	vector<Block> blk;
	template <class T>
	explicit SqrtDecomp(const vector<T>& a, int block_size = -1) : n((int)a.size()) {
		B  = block_size > 0 ? block_size : max(1, (int)llround(sqrt((double)max(n, 1))));
		nb = (n + B - 1) / B;
		blk.reserve(nb);
		for (int b = 0; b < nb; ++b)
			blk.emplace_back(vector<T>(a.begin() + b * B, a.begin() + min(n, (b + 1) * B)));
	}
	int block_of(int i) const { return i / B; }
	int local(int i)    const { return i % B; }
	template <class Whole, class Part>
	void run(int l, int r, Whole whole, Part part) {
		assert(0 <= l && l <= r && r <= n);
		if (l == r) return;
		int lb = l / B, rb = (r - 1) / B;
		if (lb == rb) { part(lb, l - lb * B, r - lb * B); return; }
		part(lb, l - lb * B, B);
		for (int b = lb + 1; b < rb; ++b) whole(b);
		part(rb, 0, r - rb * B);
	}
};
/* 使用例:
	SqrtDecomp<MyBlock> sd(a);              // B は省略で √N
	sd.run(l, r, [&](int b){ sd.blk[b].xxx_all(v); },
	             [&](int b, int s, int t){ sd.blk[b].xxx_part(s, t, v); });
	// クエリは外側の変数に畳み込む: long long res = 0; の res += ... を両ラムダで
*/
