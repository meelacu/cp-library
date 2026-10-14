// @name SqrtDecomp Block example (add + count<=x)
// @prefix sqrtblockex
// @description Block 実装例: 区間加算 & x 以下カウント

// 例: 区間加算 & 区間の「x 以下の要素数」（サマリ=ソート列, lazy=一様加算）
struct AddCountBlock {
	vector<long long> raw, srt;
	long long lz = 0;
	explicit AddCountBlock(vector<long long> v) : raw(std::move(v)) { rebuild(); }
	void rebuild() { srt = raw; sort(srt.begin(), srt.end()); }
	void push() {
		if (!lz) return;
		for (auto& x : raw) x += lz;
		for (auto& x : srt) x += lz;   // 一様加算なので順序不変
		lz = 0;
	}
	void add_all(long long x) { lz += x; }                                  // O(1)
	long long cnt_le_all(long long x) const {                               // O(log B)
		return upper_bound(srt.begin(), srt.end(), x - lz) - srt.begin();
	}
	void add_part(int s, int t, long long x) {                              // O(B log B)
		push();
		for (int i = s; i < t; ++i) raw[i] += x;
		rebuild();
	}
	long long cnt_le_part(int s, int t, long long x) const {                // O(B)
		long long c = 0;
		for (int i = s; i < t; ++i) c += (raw[i] + lz <= x);
		return c;
	}
};
