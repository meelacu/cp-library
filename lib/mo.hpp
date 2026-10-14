// @name Mo's algorithm
// @prefix mo
// @description Mo's algorithm: 区間 [l, r) クエリのオフライン処理 (構造体, main の外に置く)

// Mo's algorithm: 区間クエリ [l, r) をオフラインで処理する
//   Mo mo(n);  mo.insert(l, r);                // 0 <= l <= r <= n
//   mo.run(add, del, out);                     // 左右で更新が同じとき
//   mo.run(add_l, add_r, del_l, del_r, out);   // 左右で更新が異なるとき (転倒数など)
//   add(i) / del(i): 位置 i を区間に入れる / 外す,  out(q): q 番目のクエリの答えを記録
//   add / del の呼び出しは計 O(n√q) 回
struct Mo {
	int n;
	vector<int> L, R;
	explicit Mo(int n) : n(n) {}
	int insert(int l, int r) {
		assert(0 <= l && l <= r && r <= n);
		L.push_back(l), R.push_back(r);
		return (int)L.size() - 1;
	}
	// 処理順: 奇偶ブロック順 (幅 B は左端の移動 qB/3 と右端の移動 n^2/2B の和を最小化)
	vector<int> order() const {
		const int q = L.size();
		const int B = max(1, (int)(n / sqrt(max(1.0, 2.0 * q / 3))));
		vector<pair<unsigned long long, int>> key(q);
		for (int i = 0; i < q; i++) {
			const unsigned long long b = L[i] / B;
			key[i] = {b << 32 | (b & 1 ? n - R[i] : R[i]), i};
		}
		sort(key.begin(), key.end());
		vector<int> ord(q);
		for (int i = 0; i < q; i++) ord[i] = key[i].second;
		return ord;
	}
	template <class AL, class AR, class DL, class DR, class O>
	void run(AL&& add_l, AR&& add_r, DL&& del_l, DR&& del_r, O&& out) const {
		int l = 0, r = 0;
		for (const int i : order()) {
			const int ql = L[i], qr = R[i];
			while (l > ql) add_l(--l);
			while (r < qr) add_r(r++);
			while (l < ql) del_l(l++);
			while (r > qr) del_r(--r);
			out(i);
		}
	}
	template <class A, class D, class O>
	void run(A&& add, D&& del, O&& out) const { run(add, add, del, del, out); }
};
