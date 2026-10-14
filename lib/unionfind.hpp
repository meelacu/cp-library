// @name UnionFind
// @prefix unionfind
// @description UnionFind (union by size + iterative path compression, merge hook, count, groups)

struct UnionFind {
	int n, cnt;
	vector<int> p; // 根なら -size、そうでなければ親
	UnionFind(int n = 0) : n(n), cnt(n), p(n, -1) {}
	int find(int x) {
		assert(0 <= x && x < n);
		int r = x;
		while (p[r] >= 0) r = p[r];
		while (p[x] >= 0) { int nx = p[x]; p[x] = r; x = nx; }
		return r;
	}
	// 併合したら true。f(root, child) を「繋ぐ直前」に呼ぶ（成分ごとのデータ統合用。f 内で size(root)/size(child) は併合前の値）
	template <class F>
	bool unite(int x, int y, F f) {
		x = find(x), y = find(y);
		if (x == y) return false;
		if (-p[x] < -p[y]) swap(x, y);
		f(x, y);
		p[x] += p[y];
		p[y] = x;
		cnt--;
		return true;
	}
	bool unite(int x, int y) { return unite(x, y, [](int, int) {}); }
	bool same(int x, int y) { return find(x) == find(y); }
	int size(int x) { return -p[find(x)]; }
	int count() const { return cnt; } // 連結成分数
	vector<vector<int>> groups() {
		vector<vector<int>> g(n);
		for (int i = 0; i < n; i++) g[find(i)].push_back(i);
		g.erase(remove_if(g.begin(), g.end(), [](const vector<int>& v) { return v.empty(); }), g.end());
		return g;
	}
};
