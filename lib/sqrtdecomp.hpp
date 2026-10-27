// @name SqrtDecomp (framework)
// @prefix sqrtdecomp
// @description 平方分割 汎用フレームワーク（whole/part ディスパッチ + Block 差し替え + ブロック長の自動調整）
 
// ===== 平方分割 汎用フレームワーク =====
// SqrtDecomp<Block> sd(a);        ブロック長 B は実行中に自動調整（sd(a, B) なら B 固定）
//   sd.apply(l, r, f)             [l, r) に操作。f(k) がブロック k 全体、f(k, s, t) がブロック内ローカル [s, t)
//   sd.sum(l, r, f)               同じ分解で f の戻り値を足して返す
//   sd.fold(l, r, init, op, f)    op で左から順に畳む（min / max / 非可換も可）
//   sd.at(i, f)                   1 点: f(k, j)（k は i を含むブロック、j はブロック内の位置）
//   sd.to_vector()                全要素（各ブロックを push して raw を連結）
//   f の代わりに whole, part の 2 本を渡してもよい: sd.apply(l, r, whole, part)
//     whole(k) / part(k, s, t) の代わりに旧形式の whole(b) / part(b, s, t)（ブロック番号）も可。run は apply の別名
// Block の規約: Block(vector<T> raw) で構築 / raw に生値 / push() で lazy を raw へ反映して消す
//   whole 版 op(args) と part 版 op(s, t, args) を同名で用意すると、ラムダ 1 本で書ける:
//     sd.apply(l, r, [&](auto& k, auto... st) { k.add(st..., x); });
//     ll c = sd.sum(l, r, [&](auto& k, auto... st) { return k.cnt_le(st..., x); });
//   ブロックを丸ごと覆う端は whole に回るので、part は真に一部のときだけ呼ばれる
// 自動調整: whole と part に実際にかかった時間を抜き取りで測り、両者が釣り合うブロック長（2 の冪）を予測して、
//   組み直しの元が取れるときだけ push() 後の raw から Block を作り直す（教科書の √N は実測で 0.15〜3 倍ずれる）
//   - Block の状態が raw だけで決まらない（raw の外に要素ごとの値を持つ等）なら B を指定して固定する
//   - whole の計算量がデータ次第で O(1)〜O(B) に跳ねる償却型の Block は測定がぶれるので、B 固定の方が安全
//   - B は途中で変わるので、ブロックへの参照や block_of(i) の結果を操作をまたいで持ち越さない
template <class Block>
struct SqrtDecomp {
	int n, B, nb;
	vector<Block> blk;
 
	template <class T>
	explicit SqrtDecomp(const vector<T>& a, int block_size = -1) : n((int)a.size()) {
		tuning = block_size <= 0 && reblockable;
		lg0 = 31 - __builtin_clz(max(1, (int)sqrt((double)n)));
		if (block_size > 0) B = block_size;
		else if (tuning) B = 1 << clamp(31 - __builtin_clz(max(1, (int)sqrt((double)n))), LGMIN, lgmax());
		else B = max(1, (int)llround(sqrt((double)max(n, 1))));
		const auto t0 = clk();
		build(a);
		build_cost = double(clk() - t0);
	}
	int block_of(int i) const { return i / B; }
	int local(int i) const { return i % B; }
	int len(int b) const { return min(B, n - b * B); }
 
	template <class F> void apply(int l, int r, F&& f) { apply(l, r, f, f); }
	template <class W, class P> void apply(int l, int r, W&& whole, P&& part) {
		go(l, r, [&](int b) { call(whole, b); }, [&](int b, int s, int t) { call(part, b, s, t); });
	}
	template <class... F> void run(int l, int r, F&&... f) { apply(l, r, f...); }  // 別名（旧版の呼び方）
 
	template <class F> auto sum(int l, int r, F&& f) { return sum(l, r, f, f); }
	template <class W, class P> auto sum(int l, int r, W&& whole, P&& part) {
		decay_t<decltype(call(whole, 0))> res{};
		go(l, r, [&](int b) { res += call(whole, b); }, [&](int b, int s, int t) { res += call(part, b, s, t); });
		return res;
	}
	template <class S, class Op, class F> S fold(int l, int r, S init, Op op, F&& f) { return fold(l, r, init, op, f, f); }
	template <class S, class Op, class W, class P> S fold(int l, int r, S init, Op op, W&& whole, P&& part) {
		go(l, r, [&](int b) { init = op(init, call(whole, b)); },
		   [&](int b, int s, int t) { init = op(init, call(part, b, s, t)); });
		return init;
	}
	template <class F> decltype(auto) at(int i, F&& f) {
		assert(0 <= i && i < n);
		return f(blk[i / B], i % B);
	}
	auto to_vector() {
		vector<typename decltype(Block::raw)::value_type> v;
		v.reserve(n);
		for (auto& k : blk) k.push(), v.insert(v.end(), k.raw.begin(), k.raw.end());
		return v;
	}
 
private:
	static constexpr int LGMIN = 3;
	template <class K> static auto can_reblock(int) -> decltype(declval<K&>().push(), declval<K&>().raw.begin(), true_type{});
	template <class K> static false_type can_reblock(...);
	static constexpr bool reblockable = decltype(can_reblock<Block>(0))::value;
 
	// チューナ: 抜き取った操作の whole / part の時間と個数、ブロック長 2^g ごとの単価（丸ごと 1 個、端の 1 要素）
	bool tuning = false;
	unsigned ops = 0;
	long long NW = 0, NP = 0, sampled = 0, window = 256;
	unsigned long long Wc = 0, Pc = 0;
	double cw[32] = {}, cp[32] = {}, build_cost = 0, spent = 0;
	bool sw[32] = {}, sp[32] = {};
	int lg0 = 0;
 
	static unsigned long long clk() {
#if defined(__x86_64__) || defined(__i386__)
		return __builtin_ia32_rdtsc();
#else
		return chrono::steady_clock::now().time_since_epoch().count();
#endif
	}
	int lgmax() const { return max(LGMIN, 31 - __builtin_clz(max(n, 1))); }
 
	template <class T> void build(const vector<T>& a) {
		nb = (n + B - 1) / B;
		blk.clear();
		blk.reserve(nb);
		for (int b = 0; b < nb; b++) blk.emplace_back(vector<T>(a.begin() + b * B, a.begin() + min(n, (b + 1) * B)));
	}
	template <class F> decltype(auto) call(F& f, int b) {
		if constexpr (is_invocable_v<F&, Block&>) return f(blk[b]);
		else return f(b);
	}
	template <class F> decltype(auto) call(F& f, int b, int s, int t) {
		if constexpr (is_invocable_v<F&, Block&, int, int>) return f(blk[b], s, t);
		else return f(b, s, t);
	}
	// [l, r) を 端（部分）/ 中央（丸ごと）に分けて呼ぶ。丸ごと覆われた端は whole に回す
	template <class W, class P> void go(int l, int r, W&& whole, P&& part) {
		assert(0 <= l && l <= r && r <= n);
		if (l == r) return;
		const bool smp = tuning && (++ops & 3) == 0;
		const int lb = l / B, rb = (r - 1) / B;
		if (lb == rb) {
			edge(lb, l - lb * B, r - lb * B, whole, part, smp);
		} else {
			edge(lb, l - lb * B, len(lb), whole, part, smp);
			const auto t0 = smp ? clk() : 0;
			for (int b = lb + 1; b < rb; b++) whole(b);
			if (smp) Wc += clk() - t0, NW += rb - lb - 1;
			edge(rb, 0, r - rb * B, whole, part, smp);
		}
		if (smp) tick();
	}
	template <class W, class P> void edge(int b, int s, int t, W& whole, P& part, bool smp) {
		const auto t0 = smp ? clk() : 0;
		if (s == 0 && t == len(b)) {
			whole(b);
			if (smp) Wc += clk() - t0, NW++;
		} else {
			part(b, s, t);
			if (smp) Pc += clk() - t0, NP++;
		}
	}
	// 均衡則: いまの作業量（丸ごと NW 個・端 NP 個）を B' = 2^h で処理した時間を
	//   cw[h]·NW·(B/B') + cp[h]·NP·B'   と予測する（訪問済みの h は実測単価、未訪問はいまの単価で外挿）。
	// - 1.5 割以上速くなり、「この先もこれまでと同じだけ操作が来る」として節約が組み直し費用の 2 倍を超えるなら移る
	// - 予測で得しなくても、未訪問の隣（B を 2 倍 / 半分）は元が取れる範囲で 1 度だけ実測で試す
	// 暴走よけ: 丸ごと・端の片方しか観測できなかった窓では動かない（観測ゼロを費用ゼロと外挿しない）。
	//   B は √N の 1/8〜8 倍に制限し、判定間隔は最大約 8000 操作に抑えて傾向の変化にすぐ追従する
	void tick() {
		if (++sampled < window) return;
		const int g = __builtin_ctz(B);
		if (NW > 0) cw[g] = sw[g] ? (cw[g] + Wc / double(NW)) / 2 : Wc / double(NW), sw[g] = true;
		if (NP > 0) cp[g] = sp[g] ? (cp[g] + Pc / (double(NP) * B)) / 2 : Pc / (double(NP) * B), sp[g] = true;
		const int lo = max({LGMIN, g - 2, lg0 - 3}), hi = min({lgmax(), g + 2, lg0 + 3});
		double pr[32];
		for (int h = lo; h <= hi; h++)
			pr[h] = (sw[h] ? cw[h] : cw[g]) * NW * ldexp(1.0, g - h) + (sp[h] ? cp[h] : cp[g]) * NP * ldexp(1.0, h);
		int bg = g, probe = g;
		if (NW > 0 && NP > 0)
			for (int h = lo; h <= hi; h++) {
				if (pr[h] < pr[bg]) bg = h;
				if (abs(h - g) == 1 && !(sw[h] && sp[h]) && (probe == g || pr[h] < pr[probe])) probe = h;
			}
		spent += 4.0 * double(Wc + Pc);
		Wc = Pc = 0, NW = NP = 0, sampled = 0;
		window = min(window * 2, 2048LL);
		if (bg != g && pr[bg] < 0.85 * pr[g]) {
			if ((1 - pr[bg] / pr[g]) * spent >= 2 * build_cost) reblock(bg);
		} else if (probe != g && spent >= 8 * build_cost) {
			reblock(probe);
		}
	}
	void reblock(int h) {
		if constexpr (reblockable) {
			const auto t0 = clk();
			const auto v = to_vector();
			B = 1 << h;
			build(v);
			build_cost = double(clk() - t0);
			window = 256;
			reblocks++;
		}
	}
	int reblocks = 0;
};
/* 使用例（Block は sqrtblockex の AddCountBlock）:
	SqrtDecomp<AddCountBlock<>> sd(a);
	sd.apply(l, r, [&](auto& k, auto... st) { k.add(st..., x); });                    // [l, r) に x を加算
	ll c = sd.sum(l, r, [&](auto& k, auto... st) { return k.cnt_le(st..., x); });     // [l, r) で x 以下の個数
	ll v = sd.at(i, [&](auto& k, int j) { return k.get(j); });                         // a[i]
*/