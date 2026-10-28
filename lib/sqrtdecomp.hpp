// @name SqrtDecomp (framework)
// @prefix sqrtdecomp
// @description 平方分割 汎用フレームワーク（whole/part ディスパッチ + Block 差し替え + ブロック長の自動調整）
 
// ===== 平方分割 汎用フレームワーク =====
// SqrtDecomp<Block> sd(a);        ブロック長 B は実行中に自動調整（sd(a, B) なら B 固定）
//   sd.apply(l, r, f)             [l, r) に操作。f(k) がブロック k 全体、f(k, s, t) がブロック内ローカル [s, t)
//   sd.sum(l, r, f)               同じ分解で f の戻り値を足して返す
//   sd.fold(l, r, init, op, f)    op で左から順に畳む（min / max / 非可換も可）
//   sd.at(i, f)                   1 点: f(k, j)（k は i を含むブロック、j はブロック内の位置）。1 点更新もここから
//   sd.to_vector()                全要素（各ブロックを push して raw を連結）
//   f の代わりに whole, part の 2 本を渡してもよい: sd.apply(l, r, whole, part)
//     whole(k) / part(k, s, t) の代わりに旧形式の whole(b) / part(b, s, t)（ブロック番号）も可。run は apply の別名
// Block の規約: Block(vector<T> raw) で構築 / raw に生値 / push() で lazy を raw へ反映して消す
//   whole 版 op(args) と part 版 op(s, t, args) を同名で用意すると、ラムダ 1 本で書ける:
//     sd.apply(l, r, [&](auto& k, auto... st) { k.add(st..., x); });
//     ll c = sd.sum(l, r, [&](auto& k, auto... st) { return k.cnt_le(st..., x); });
//   ブロックを丸ごと覆う端は whole に回るので、part は真に一部のときだけ呼ばれる
// 自動調整: 操作を抜き取り、whole / part / at にかかった時間を呼び出し箇所（ラムダ）ごとに測って、
//   B を 2 倍 / 半分にした段を実測で比べながら速い方へ移る（山登り）。組み直しは push() 後の raw から作り直す
//   - 保険: 同じ B のまま測った whole の単価が 4 倍を超えて跳ね、次の窓もその水準に留まったら（データ次第で
//     O(1)〜O(B) に跳ねる償却型の Block で状態が急に変わった）、測定を信用せず B = √N/2 付近に固定して以後は調整しない
//   - Block の状態が raw だけで決まらない（raw の外に要素ごとの値を持つ等）なら B を指定して固定する
//   - B は途中で変わるので、ブロックへの参照や block_of(i) の結果を操作をまたいで持ち越さない
//     （at が返した参照は次の sd の操作まで有効）。コールバックの中から sd の操作を呼ばない
template <class Block>
struct SqrtDecomp {
	int n, B, nb;
	vector<Block> blk;
 
	template <class T>
	explicit SqrtDecomp(const vector<T>& a, int block_size = -1) : n((int)a.size()) {
		tuning = block_size <= 0 && reblockable;
		lg0 = lg(max(1, (int)sqrt((double)n)));
		if (block_size > 0) B = block_size;
		else if (tuning) B = 1 << clamp(lg0, LGMIN, lgmax());
		else B = max(1, (int)llround(sqrt((double)max(n, 1))));
		if (tuning) calibrate();
		const auto t0 = clk();
		build(a);
		build_cost = double(clk() - t0);
	}
	int block_of(int i) const { return sh >= 0 ? i >> sh : i / B; }
	int local(int i) const { return i - block_of(i) * B; }
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
		if (due) due = false, close();
		const int b = block_of(i), j = i - b * B;
		if (!tuning || (++ops & mask)) return f(blk[b], j);
		// 抜き取り: 1 点操作は端 1 個として測る（窓の締めは次の操作の頭で行い、返す参照を組み直しで壊さない）
		Win& w = win[site<decay_t<F>>()];
		const auto t0 = clk();
		if constexpr (is_void_v<decltype(f(blk[b], j))>) {
			f(blk[b], j);
			note(w, false, 1, t0), sampled++, due = full();
		} else {
			decltype(auto) res = f(blk[b], j);
			note(w, false, 1, t0), sampled++, due = full();
			return res;
		}
	}
	auto to_vector() {
		vector<typename decltype(Block::raw)::value_type> v;
		v.reserve(n);
		for (auto& k : blk) k.push(), v.insert(v.end(), k.raw.begin(), k.raw.end());
		return v;
	}
	bool tuning_now() const { return tuning; }  // 自動調整中か（B 指定・保険での固定なら false）
 
private:
	static constexpr int LGMIN = 3, MS = 8;  // MS: 呼び出し箇所（ラムダの型）ごとに統計を分ける数
	template <class K> static auto can_reblock(int) -> decltype(declval<K&>().push(), declval<K&>().raw.begin(), true_type{});
	template <class K> static false_type can_reblock(...);
	static constexpr bool reblockable = decltype(can_reblock<Block>(0))::value;
	static int lg(int x) { return 31 - __builtin_clz(x); }
	int sh = -1;  // B が 2 の冪なら log2(B)（割り算をシフトにする）
 
	template <class T> void build(const vector<T>& a) {
		nb = (n + B - 1) / B;
		sh = (B & (B - 1)) == 0 ? lg(B) : -1;
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
		if (due) due = false, close();
		const int lb = block_of(l), rb = block_of(r - 1), s = l - lb * B, t = r - rb * B;
		if (!tuning || (++ops & mask)) {
			if (lb == rb) {
				if (s == 0 && t == len(lb)) whole(lb); else part(lb, s, t);
				return;
			}
			if (s == 0) whole(lb); else part(lb, s, B);
			for (int b = lb + 1; b < rb; b++) whole(b);
			if (t == len(rb)) whole(rb); else part(rb, 0, t);
			return;
		}
		// 抜き取った操作: 端・中央・端の時間を続けて測る
		Win& w = win[site<W>()];
		auto t0 = clk();
		if (lb == rb) {
			const bool wh = s == 0 && t == len(lb);
			if (wh) whole(lb); else part(lb, s, t);
			note(w, wh, 1, t0);
		} else {
			if (s == 0) whole(lb); else part(lb, s, B);
			t0 = note(w, s == 0, 1, t0);
			if (rb - lb > 1) {
				for (int b = lb + 1; b < rb; b++) whole(b);
				t0 = note(w, true, rb - lb - 1, t0);
			}
			const bool wh = t == len(rb);
			if (wh) whole(rb); else part(rb, 0, t);
			note(w, wh, 1, t0);
		}
		tick();
	}
 
	// ---- 自動調整 ----
	// 窓: 抜き取った操作を呼び出し箇所 s ごとに集計（丸ごと nw 個に wc、端 np 個に pc、8 個以上まとめて測った丸ごと nm 個に wm）
	struct Win { unsigned long long wc = 0, pc = 0, wm = 0; long long nw = 0, np = 0, nm = 0; };
	bool tuning = false, frozen = false, due = false;  // due: 窓が満ちていて次の操作の頭で締める
	unsigned ops = 0, mask = 3;          // ops & mask == 0 の操作を測る（1 操作が軽いほど間引く）
	long long sampled = 0, window = 256, ovh = 0, nwin = 0, stp[32];
	double acc = 0, rate = 0;            // いまの窓の合計時間、前の窓の 1 操作あたり時間
	Win win[MS];
	// 段 h（B = 2^h）での実測の累計（箇所 s ごと）: 丸ごとの時間 / 個数、端の時間 / (個数·B)。単価はその比
	double tw[MS][32] = {}, nw[MS][32] = {}, tp[MS][32] = {}, np[MS][32] = {};
	long long th[32] = {}, kh[32] = {};  // 段 h を最後に測った窓の番号（0 = 未測定）と、続けて測った窓の数
	double last[MS], pend[MS], build_cost = 0, spent = 0, rprobe = 0;  // last / pend: 跳ねの検出用、rprobe: 試しの組み直しの費用
	int lg0 = 0, reblocks = 0, from = -1, to = -1;
	static inline int nsite = 0;
	template <class W> static int site() { static const int id = min(nsite++, MS - 1); return id; }
 
	static unsigned long long clk() {
#if defined(__x86_64__) || defined(__i386__)
		return __builtin_ia32_rdtsc();
#else
		return chrono::steady_clock::now().time_since_epoch().count();
#endif
	}
	void calibrate() {  // 計測そのものの費用
		long long m = LLONG_MAX;
		for (int i = 0; i < 64; i++) { const auto a = clk(), b = clk(); m = min(m, (long long)(b - a)); }
		ovh = m;
		for (auto& x : last) x = -1;
		for (auto& x : pend) x = -1;
		for (auto& x : stp) x = STALE;
	}
	unsigned long long note(Win& w, bool wh, int cnt, unsigned long long t0) {
		const auto t1 = clk();
		const long long d0 = (long long)(t1 - t0) - ovh;
		const unsigned long long d = d0 > 0 ? d0 : 0;
		acc += d;
		if (!wh) w.pc += d, w.np++;
		else {
			w.wc += d, w.nw += cnt;
			if (cnt >= 8) w.wm += d, w.nm += cnt;
		}
		return t1;
	}
	int lgmax() const { return max(LGMIN, lg(max(n, 1))); }
	bool fresh(int h) const { return th[h] && nwin - th[h] <= stp[h]; }
 
	// 均衡則: いまの窓の作業量（箇所 s ごとに丸ごと nw 個・端 np 個）を B' = 2^h で処理した時間を
	//   Σ_s (丸ごと単価)·nw·(B/B') + (端の 1 要素単価)·np·B'   と予測する（最近測った段は実測単価、他はいまの段の単価で外挿）
	// - 移る先は実測済みの段だけ: 2 割（1 窓しか測っていない段なら 3 割）以上速く、「この先もこれまでと同じだけ
	//   操作が来る」として節約が組み直し費用の 2 倍を超えるなら移る
	// - いまの段が実測済みの中で最良なら、未測定（か古い）段を外挿で得しそうなものから試す（丸ごと・端の両方が来て
	//   いて 4 割以上得しそうなら 2 段先、それ以外は隣）。外挿は whole の単価が B によらない前提なので、償却型の
	//   Block では当てにならない。だから移る前に必ず測る。試して悪かった段は次に試すまでの間隔を倍にし、
	//   試しの組み直しは最初の 2 回を除いてそれまでの時間の 1 割まで
	// - B は √N の 1/8〜8 倍。窓は 256〜2048 回の抜き取り（組み直し直後は費用の 2 倍程度の時間）で、
	//   1 操作の時間が前の窓の 4 倍を超えたら早めに閉じる。抜き取り間隔は計測の費用が全体の 1% 程度になるように
	static constexpr double FLUCT = 4.0;
	static constexpr long long STALE = 16;
	bool full() const { return sampled >= window || (sampled >= 64 && rate > 0 && acc > 4 * rate * sampled); }
	void tick() {
		if (++sampled, full()) close();
	}
	void close() {
		rate = acc / sampled;
		const double per = mask + 1;
		unsigned p = 4;
		while (p < 64 && p * rate < 15000) p *= 2;
		mask = p - 1;
		acc = 0, sampled = 0, nwin++;
		const int g = __builtin_ctz(B);
		if (!fresh(g)) {  // 古い測定は捨てて測り直す
			kh[g] = 0;
			for (int s = 0; s < MS; s++) tw[s][g] = nw[s][g] = tp[s][g] = np[s][g] = 0;
		}
		th[g] = nwin, kh[g]++;
		long long NW = 0, NP = 0;
		bool jump = false, pending = false;
		Win cur[MS];
		for (int s = 0; s < MS; s++) {
			const Win& w = cur[s] = win[s];
			win[s] = Win{};
			NW += w.nw, NP += w.np;
			spent += per * double(w.wc + w.pc);
			tw[s][g] += w.wc, nw[s][g] += w.nw, tp[s][g] += w.pc, np[s][g] += double(w.np) * B;
			if (w.nm >= 64) {
				// 跳ね: 前の窓と 4 倍を超えて違い、次の窓もその新しい水準に留まったら確定（1 窓だけの雑音は無視）
				const double x = w.wm / double(w.nm);
				auto near = [](double a, double b) { return max(a, b) <= 2 * min(a, b) + 4; };
				if (pend[s] >= 0) {
					if (near(x, pend[s]) && !near(x, last[s])) jump = true;
					pend[s] = -1, last[s] = x;
				} else if (last[s] >= 0 && max(x, last[s]) > FLUCT * min(x, last[s]) + 4) {
					pend[s] = x;
				} else {
					last[s] = x;
				}
			}
			pending |= pend[s] >= 0;
		}
		if (jump) return freeze();
		window = min(window * 2, 2048LL);
		if (NW + NP == 0 || nwin < 2 || pending) return;  // 最初の B は 2 窓測ってから。跳ねの確認待ちの間も動かない
		const int lo = max({LGMIN, g - 2, lg0 - 3}), hi = min({lgmax(), g + 2, lg0 + 3});
		double pr[32] = {};
		bool seen[32];
		for (int h = lo; h <= hi; h++) {
			seen[h] = fresh(h);
			for (int s = 0; s < MS; s++) {
				const Win& w = cur[s];
				const bool fw = seen[h] && nw[s][h] > 0, fp = seen[h] && np[s][h] > 0;
				if (w.nw) pr[h] += (fw ? tw[s][h] / nw[s][h] : tw[s][g] / nw[s][g]) * w.nw * ldexp(1.0, g - h), seen[h] &= fw;
				if (w.np) pr[h] += (fp ? tp[s][h] / np[s][h] : tp[s][g] / np[s][g]) * w.np * ldexp(1.0, h), seen[h] &= fp;
			}
		}
		int bg = g, probe = g;
		for (int h = lo; h <= hi; h++)
			if (seen[h] && pr[h] < pr[bg]) bg = h;
		if (bg == g) {
			for (int h = lo; h <= hi; h++) {
				const int m = (g + h) / 2;
				if (h == g || seen[h] || pr[h] >= pr[probe]) continue;
				if (abs(h - g) == 2 && (!(NW > 0 && NP > 0) || pr[h] > 0.6 * pr[g] || (seen[m] && pr[m] >= pr[g]))) continue;
				probe = h;
			}
			if (probe == g)
				for (int h : {g - 1, g + 1})
					if (lo <= h && h <= hi && !seen[h] && (probe == g || pr[h] < pr[probe])) probe = h;
		}
		if (bg != g) {
			if (pr[bg] < (kh[bg] >= 2 ? 0.8 : 0.7) * pr[g] && (1 - pr[bg] / pr[g]) * spent >= 2 * build_cost) {
				if (g == to && bg == from) stp[g] = min(stp[g] * 2, 96LL);  // 試した段が悪かった
				from = to = -1;
				reblock(bg);
			}
		} else if (probe != g && spent >= (pr[probe] < pr[g] ? 2 : 8) * build_cost && rprobe + build_cost <= 0.1 * spent + 2 * build_cost) {
			from = g, to = probe;
			rprobe += build_cost;
			reblock(probe);
		} else if (g == to) {
			from = to = -1;
		}
	}
	void freeze() {
		tuning = false, frozen = true;
		const int h = clamp((int)lround(log2(max(1.0, sqrt((double)n) / 2))), LGMIN, lgmax());
		if ((1 << h) != B) reblock(h);
	}
	void reblock(int h) {
		if constexpr (reblockable) {
			const auto t0 = clk();
			const auto v = to_vector();
			B = 1 << h;
			build(v);
			build_cost = double(clk() - t0);
			window = clamp<long long>((long long)(2 * build_cost / max(1.0, rate * (mask + 1))), 256, 8192);
			sampled = 0, acc = 0, due = false;
			reblocks++;
			for (auto& x : last) x = -1;
			for (auto& x : pend) x = -1;
			for (auto& w : win) w = Win{};
		}
	}
};
/* 使用例（Block は sqrtblockex の AddCountBlock）:
	SqrtDecomp<AddCountBlock<>> sd(a);
	sd.apply(l, r, [&](auto& k, auto... st) { k.add(st..., x); });                    // [l, r) に x を加算
	ll c = sd.sum(l, r, [&](auto& k, auto... st) { return k.cnt_le(st..., x); });     // [l, r) で x 以下の個数
	ll v = sd.at(i, [&](auto& k, int j) { return k.get(j); });                         // a[i]
	sd.at(i, [&](auto& k, int j) { k.set(j, x); });                                    // a[i] = x
*/
 