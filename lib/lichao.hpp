// @name Li Chao Tree
// @prefix lichao
// @prefix LiChaoTree
// @description Li Chao Tree (一般化版: 一回交差関数族 F を載せる, min/max, 線分, optional<F> で復元)

// Li Chao Tree（一般化版）— セグ木にモノイドを載せるのと同じ感覚で「関数族 F」を載せる
// F の要件:
//   using X = ...;            // 引数（クエリ点）の型。Cmp (既定 less<X>) で全順序が入ること
//   T operator()(X x) const;  // x での値。T は自動推論 (< で比較可能で numeric_limits を持つこと)
//   任意の f, g ∈ F について f(x) - g(x) が Cmp 順に見て高々 1 回しか 0 にならない (交点が高々 1 つ)
// MINIMIZE: true=最小値, false=最大値
// 使い方:
//   LiChaoTree<LiChaoLine<ll>> lct(xs);  // 直線 y = a*x + b は付属の LiChaoLine で。xs: クエリ点の一覧
//   lct.add({a, b, id});                 // 関数を追加 (id 省略可, 復元用)                O(log N)
//   lct.add_segment({a, b}, xl, xr);     // 閉区間 xl <= x <= xr でのみ有効                O(log^2 N)
//   lct.add_segment_idx({a, b}, l, r);   // 添字の半開区間 [l, r) 版
//   lct.query(x);                        // x での最小(最大)値。無ければ INF                O(log N)
//   lct.query_func(x);                   // 最適値を与える F を optional で。無ければ nullopt
//   lct.query_idx(i) / query_func_idx(i) // xs (ソート・重複除去後) の i 番目での版
template <class T>
struct LiChaoLine {
	using X = T;
	T a, b;
	int id = -1;
	T operator()(T x) const { return a * x + b; }
};

template <class F, bool MINIMIZE = true, class Cmp = less<typename F::X>>
struct LiChaoTree {
	using X = typename F::X;
	using T = decltype(declval<const F&>()(declval<X>()));
	static constexpr T INF = MINIMIZE ? numeric_limits<T>::max() : numeric_limits<T>::lowest();

	int n = 0, sz = 1;
	Cmp cmp;
	vector<X> xs;
	vector<F> seg;
	vector<char> used;

	LiChaoTree() = default;
	explicit LiChaoTree(vector<X> xs_, Cmp cmp_ = Cmp()) : cmp(cmp_), xs(std::move(xs_)) {
		sort(xs.begin(), xs.end(), cmp);
		auto eq = [&](const X& p, const X& q) { return !cmp(p, q) && !cmp(q, p); };
		xs.erase(unique(xs.begin(), xs.end(), eq), xs.end());
		n = (int)xs.size();
		while (sz < n) sz <<= 1;
		seg.resize(2 * sz);
		used.assign(2 * sz, 0);
	}

	static bool better(const T& p, const T& q) { return MINIMIZE ? p < q : q < p; }

	int index(const X& x) const {
		int i = lower_bound(xs.begin(), xs.end(), x, cmp) - xs.begin();
		assert(i < n && !cmp(x, xs[i]));  // x はコンストラクタで渡した点のいずれかでなければならない
		return i;
	}

	void add(F f) {
		if (n == 0) return;
		add_node(f, 1, 0, sz);
	}
	void add_segment_idx(F f, int l, int r) {
		l = max(l, 0), r = min(r, n);
		if (l >= r) return;
		int w = 1;
		for (int p = l + sz, q = r + sz; p < q; p >>= 1, q >>= 1, w <<= 1) {
			if (p & 1) add_node(f, p, p * w - sz, (p + 1) * w - sz), p++;
			if (q & 1) q--, add_node(f, q, q * w - sz, (q + 1) * w - sz);
		}
	}
	void add_segment(F f, X xl, X xr) {
		int l = lower_bound(xs.begin(), xs.end(), xl, cmp) - xs.begin();
		int r = upper_bound(xs.begin(), xs.end(), xr, cmp) - xs.begin();
		add_segment_idx(f, l, r);
	}

	// 内部用: i 番目の点で最適な関数を持つノード (無ければ 0)
	int best_node(int i) const {
		assert(0 <= i && i < n);
		T res = INF;
		int best = 0;
		for (int k = i + sz; k >= 1; k >>= 1)
			if (used[k]) {
				T v = seg[k](xs[i]);
				if (better(v, res)) res = v, best = k;
			}
		return best;
	}
	T query_idx(int i) const { int k = best_node(i); return k ? seg[k](xs[i]) : INF; }
	optional<F> query_func_idx(int i) const { int k = best_node(i); if (k) return seg[k]; return nullopt; }
	// 比較は xs の代表点で行い、値は渡された x で評価する (Cmp で同値な点は同じ順序を与えること)
	T query(const X& x) const { int k = best_node(index(x)); return k ? seg[k](x) : INF; }
	optional<F> query_func(const X& x) const { return query_func_idx(index(x)); }

	// 内部用: ノード k (担当添字 [l, r)) に関数 f を挿入
	void add_node(F f, int k, int l, int r) {
		while (true) {
			if (!used[k]) { seg[k] = f, used[k] = 1; return; }
			int m = (l + r) >> 1;
			const X &xl = xs[min(l, n - 1)], &xm = xs[min(m, n - 1)], &xr = xs[min(r - 1, n - 1)];
			F& g = seg[k];
			bool bl = better(f(xl), g(xl)), bm = better(f(xm), g(xm)), br = better(f(xr), g(xr));
			if (bl && br) { g = f; return; }   // f が区間全体で優位 → 置換して終了
			if (!bl && !br) return;            // f が区間全体で劣位 → 捨てて終了
			if (bm) swap(g, f), bl = !bl, br = !br;  // 中点で優位な方をノードに残す
			if (r - l == 1) return;
			if (bl) k = 2 * k, r = m;          // 敗者は片側でしか勝てない（一回交差性）
			else k = 2 * k + 1, l = m;
		}
	}
};
