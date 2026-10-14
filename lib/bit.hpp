// @name BIT
// @prefix bit
// @description Fenwick tree (0-indexed, half-open [l, r), O(n) build, lower_bound)

template <class T>
struct BIT {
	int n;
	vector<T> d;
	BIT(int n = 0) : n(n), d(n + 1, T()) {}
	// O(n) 構築
	BIT(const vector<T>& a) : n((int)a.size()), d(a.size() + 1, T()) {
		for (int i = 1; i <= n; i++) {
			d[i] += a[i - 1];
			int j = i + (i & -i);
			if (j <= n) d[j] += d[i];
		}
	}
	// a[i] += v
	void add(int i, T v) {
		assert(0 <= i && i < n);
		for (i++; i <= n; i += i & -i) d[i] += v;
	}
	// sum of [0, r)
	T sum(int r) const {
		assert(0 <= r && r <= n);
		T s = T();
		for (; r > 0; r -= r & -r) s += d[r];
		return s;
	}
	// sum of [l, r)
	T sum(int l, int r) const {
		assert(0 <= l && l <= r && r <= n);
		return sum(r) - sum(l);
	}
	T get(int i) const { return sum(i, i + 1); }
	void set(int i, T v) { add(i, v - get(i)); }
	// min i s.t. sum(0, i + 1) >= x (要 a[j] >= 0)。存在しなければ n
	int lower_bound(T x) const {
		if (x <= T()) return 0;
		int i = 0, k = 1;
		while (k << 1 <= n) k <<= 1;
		for (; k > 0; k >>= 1)
			if (i + k <= n && d[i + k] < x) i += k, x -= d[i];
		return i;
	}
};
