// @name SegTree
// @prefix segtree
// @description Segment tree (ACL-style, with max_right / min_left)

template <class T, auto op, auto e>
struct SegTree {
	int sz, n;
	vector<T> dat;

	SegTree() : SegTree(0) {}
	explicit SegTree(int n_) : SegTree(vector<T>(n_, e())) {}
	explicit SegTree(const vector<T>& v) : sz(v.size()) {
		n = 1;
		while (n < sz) n *= 2;
		dat.assign(2 * n, e());
		for (int i = 0; i < sz; i++) dat[i + n] = v[i];
		for (int i = n - 1; i >= 1; i--) update(i);
	}

	int size() const { return sz; }

	void set(int i, T x) {
		assert(0 <= i && i < sz);
		i += n;
		dat[i] = x;
		while (i >>= 1) update(i);
	}
	// a[i] <- op(a[i], x)
	void apply(int i, T x) { set(i, op(get(i), x)); }

	T get(int i) const {
		assert(0 <= i && i < sz);
		return dat[i + n];
	}
	T operator[](int i) const { return get(i); }

	// [l, r)
	T prod(int l, int r) const {
		assert(0 <= l && l <= r && r <= sz);
		T resL = e(), resR = e();
		l += n; r += n;
		while (l < r) {
			if (l & 1) resL = op(resL, dat[l++]);
			if (r & 1) resR = op(dat[--r], resR);
			l >>= 1; r >>= 1;
		}
		return op(resL, resR);
	}
	T all_prod() const { return dat[1]; }

	// f(prod(l, r)) が true となる最大の r （f(e()) == true, f は単調）
	template <class F> int max_right(int l, F f) const {
		assert(0 <= l && l <= sz);
		assert(f(e()));
		if (l == sz) return sz;
		l += n;
		T sm = e();
		do {
			while (l % 2 == 0) l >>= 1;
			if (!f(op(sm, dat[l]))) {
				while (l < n) {
					l = 2 * l;
					if (f(op(sm, dat[l]))) { sm = op(sm, dat[l]); l++; }
				}
				return l - n;
			}
			sm = op(sm, dat[l]); l++;
		} while ((l & -l) != l);
		return sz;
	}

	// f(prod(l, r)) が true となる最小の l
	template <class F> int min_left(int r, F f) const {
		assert(0 <= r && r <= sz);
		assert(f(e()));
		if (r == 0) return 0;
		r += n;
		T sm = e();
		do {
			r--;
			while (r > 1 && (r % 2)) r >>= 1;
			if (!f(op(dat[r], sm))) {
				while (r < n) {
					r = 2 * r + 1;
					if (f(op(dat[r], sm))) { sm = op(dat[r], sm); r--; }
				}
				return r + 1 - n;
			}
			sm = op(dat[r], sm);
		} while ((r & -r) != r);
		return 0;
	}

private:
	void update(int k) { dat[k] = op(dat[2 * k], dat[2 * k + 1]); }
};
