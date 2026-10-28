// @name SqrtDecomp Block example (add + count<=x)
// @prefix sqrtblockex
// @description Block 実装例: 区間加算 & x 以下カウント（多段索引 + O(B) 再構築の高速版）
 
// 例: 区間加算 & 区間の「x 以下の要素数」  AddCountBlock<T>（T 省略で long long。int に収まるなら int で約 1.5〜1.8 倍速い）
//   whole: add(x) O(1)  cnt_le(x) O(log_K B)        part: add(s, t, x) O(B)  cnt_le(s, t, x) O(t - s)
//   1 点: get(i)  set(i, x)（O(B)）
//   x 未満 = cnt_le(x - 1)  [lo, hi] 内 = cnt_le(hi) - cnt_le(lo - 1)  k 番目に小さい値 = cnt_le(v) >= k となる最小の v
// 仕組み（教科書版との差）:
//   - 丸ごとの cnt_le は二分探索ではなく、ソート列を K 個（64 バイト = キャッシュライン 1 本）おきに間引いた
//     多段索引を、上の段から K 個ずつまとめて比較して降りる（最上段はブロック本体に直に置く）
//   - 部分加算はソート列を「範囲内（+x）/ 範囲外」に安定分割してマージ（sort の O(B log B) → O(B)）。
//     マージは前半・後半をそれぞれ前後から、計 4 本の独立な依存鎖で分岐なしに進める
//   - 構築（ブロック長の自動調整での組み直しも）は 8 bit ずつの基数ソートで、全要素で同じ桁は飛ばす
//   - 値は numeric_limits<T> の max / min を番兵に使うので、途中経過も含めてその手前に収まること
#ifdef __AVX2__
#include <immintrin.h>
#endif
template <class T = long long>
struct alignas(64) AddCountBlock {
	static constexpr int K = 64 / sizeof(T);
	static constexpr T INF = numeric_limits<T>::max(), NEG = numeric_limits<T>::min();
	T top[K];           // 索引の最上段（K 個以下。余りは INF）
	T lz = 0;           // 一様加算の遅延
	int m = 0, dep = 0, al = 0, off[10] = {}, cnt[11] = {};  // 要素数、段数、各段の位置と要素数（段 0 = ソート列）
	vector<T> buf;      // [ソート列][1 段目]…[dep-1 段目]（各段 K の倍数に切り上げ、余りは INF）
	vector<T> raw;      // 生値（lz 抜き）
	vector<int> pos;    // ソート列の各要素の元の位置
 
	explicit AddCountBlock(vector<T> v) : m((int)v.size()), raw(std::move(v)), pos(m) {
		cnt[0] = m;
		while (cnt[dep] > K) cnt[dep + 1] = (cnt[dep] + K - 1) / K, dep++;
		int tot = 0;
		for (int d = 0; d < max(dep, 1); d++) off[d] = tot, tot += max(1, (cnt[d] + K - 1) / K) * K;
		buf.assign(tot + K, INF);
		al = int((64 - reinterpret_cast<uintptr_t>(buf.data()) % 64) % 64 / sizeof(T));  // 64 バイト境界に揃える
		sort_init();
		reindex();
	}
	void push() {
		if (!lz) return;
		for (auto& x : raw) x += lz;
		for (int i = 0; i < m; i++) lv(0)[i] += lz;
		lz = 0;
		reindex();
	}
 
	// ---- whole ----
	void add(T x) { lz += x; }
	int cnt_le(T x) const {
		const T y = x - lz;
		int c = le(top, y);
		for (int d = dep - 1; d >= 0; d--) {
			const int s = K * max(c - 1, 0);
			c = s + le(lv(d) + s, y);
		}
		return c;
	}
	// ---- part（[s, t) はブロック内ローカル）----
	void add(int s, int t, T x) {
		for (int i = s; i < t; i++) raw[i] += x;
		T* v = lv(0);
		int* q = pos.data();
		if constexpr (sizeof(T) == 4) {  // 値と位置を 64 bit に詰めて 1 回で比較・移動する
			static vector<long long> A, C;
			if ((int)A.size() < m + 2) A.resize(m + 2), C.resize(m + 2);
			long long *a = A.data(), *c = C.data();
			const unsigned long long dx = (unsigned long long)(long long)x << 32;
			int na = 0, nc = 0;
			for (int k = 0; k < m; k++) {  // 安定分割（分岐なし）
				const int p = q[k], in = unsigned(p - s) < unsigned(t - s);
				const long long key = (long long)((unsigned long long)(long long)v[k] << 32 | unsigned(p));
				a[na + 1] = (long long)((unsigned long long)key + dx), c[nc + 1] = key;
				na += in, nc += 1 - in;
			}
			merge4(a, c, na, nc, [&](int k, long long w) { v[k] = T(w >> 32), q[k] = int(w & 0xffffffff); });
		} else {
			static vector<pair<T, int>> A, C;
			if ((int)A.size() < m + 2) A.resize(m + 2), C.resize(m + 2);
			pair<T, int> *a = A.data(), *c = C.data();
			int na = 0, nc = 0;
			for (int k = 0; k < m; k++) {  // 安定分割（分岐なし）
				const int p = q[k], in = unsigned(p - s) < unsigned(t - s);
				a[na + 1] = {T(make_unsigned_t<T>(v[k]) + make_unsigned_t<T>(x)), p}, c[nc + 1] = {v[k], p};
				na += in, nc += 1 - in;
			}
			merge4(a, c, na, nc, [&](int k, const pair<T, int>& w) { v[k] = w.first, q[k] = w.second; });
		}
		reindex();
	}
	int cnt_le(int s, int t, T x) const {
		const T y = x - lz;
		int c = 0, i = s;
		for (; i + K <= t; i += K) c += le(&raw[i], y);
		for (; i < t; i++) c += raw[i] <= y;
		return c;
	}
	// ---- 1 点 ----
	T get(int i) const { return raw[i] + lz; }
	void set(int i, T x) { add(i, i + 1, x - get(i)); }
 
private:
	T* lv(int d) { return buf.data() + al + off[d]; }
	const T* lv(int d) const { return buf.data() + al + off[d]; }
	static int le(const T* p, T y) {  // p[0, K) のうち y 以下の個数
#ifdef __AVX2__
		unsigned g = 0;
		if constexpr (sizeof(T) == 4) {
			const __m256i yv = _mm256_set1_epi32(y);
			for (int j = 0; j < K; j += 8)
				g |= unsigned(_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_loadu_si256((const __m256i*)(p + j)), yv)))) << j;
		} else {
			const __m256i yv = _mm256_set1_epi64x(y);
			for (int j = 0; j < K; j += 4)
				g |= unsigned(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(_mm256_loadu_si256((const __m256i*)(p + j)), yv)))) << j;
		}
		return K - __builtin_popcount(g);
#else
		int c = 0;
		for (int j = 0; j < K; j++) c += p[j] <= y;
		return c;
#endif
	}
	// ソート列と pos を作る: 符号を反転した値の LSD 基数ソート（8 bit ずつ、全要素で同じ桁は飛ばす）。
	// 桁が多いのにブロックが小さいときは std::sort
	void sort_init() {
		using U = make_unsigned_t<T>;
		const U sb = U(1) << (8 * sizeof(T) - 1);
		U o = 0, n = ~U(0);
		for (int i = 0; i < m; i++) o |= U(raw[i]) ^ sb, n &= U(raw[i]) ^ sb;
		int passes = 0;
		for (int sh = 0; sh < 8 * (int)sizeof(T); sh += 8) passes += ((o ^ n) >> sh & 255) != 0;
		T* z = lv(0);
		if (passes * (2 * m + 256) > 3 * m * (32 - __builtin_clz(m | 1))) {
			vector<pair<T, int>> s(m);
			for (int i = 0; i < m; i++) s[i] = {raw[i], i};
			sort(s.begin(), s.end());
			for (int i = 0; i < m; i++) z[i] = s[i].first, pos[i] = s[i].second;
			return;
		}
		static vector<U> k0, k1;
		static vector<int> q1;
		if ((int)k0.size() < m) k0.resize(m), k1.resize(m), q1.resize(m);
		U *a = k0.data(), *b = k1.data();
		int *p = pos.data(), *q = q1.data();
		for (int i = 0; i < m; i++) a[i] = U(raw[i]) ^ sb, p[i] = i;
		for (int sh = 0; sh < 8 * (int)sizeof(T); sh += 8) {
			if (((o ^ n) >> sh & 255) == 0) continue;  // 全要素で同じ桁
			int c[257] = {};
			for (int i = 0; i < m; i++) c[(a[i] >> sh & 255) + 1]++;
			for (int d = 0; d < 256; d++) c[d + 1] += c[d];
			for (int i = 0; i < m; i++) {
				const int d = c[a[i] >> sh & 255]++;
				b[d] = a[i], q[d] = p[i];
			}
			swap(a, b), swap(p, q);
		}
		for (int i = 0; i < m; i++) z[i] = T(a[i] ^ sb);
		if (p != pos.data()) copy(p, p + m, pos.begin());
	}
	void reindex() {
		for (int d = 1; d < dep; d++)
			for (int j = 0; j < cnt[d]; j++) lv(d)[j] = lv(d - 1)[K * j];
		const T* lo = lv(max(dep - 1, 0));
		for (int j = 0; j < K; j++) top[j] = j < cnt[dep] ? lo[dep ? K * j : j] : INF;
	}
	// a[1..na], c[1..nc]（それぞれソート済み、同値は c が先）を併合し、k 番目を out(k, 要素) で渡す。
	// マージパスで前半 h 個と後半に分け、それぞれを前後から進める
	template <class R, class Out> void merge4(R* a, R* c, int na, int nc, Out&& out) {
		auto key = [](const R& r) { if constexpr (is_arithmetic_v<R>) return r; else return r.first; };
		if constexpr (is_arithmetic_v<R>) a[0] = c[0] = numeric_limits<R>::min(), a[na + 1] = c[nc + 1] = numeric_limits<R>::max();
		else a[0] = c[0] = {NEG, 0}, a[na + 1] = c[nc + 1] = {INF, 0};
		const int h = m / 2;
		int lo = max(0, h - nc), hi = min(h, na);
		while (lo < hi) {
			const int md = (lo + hi) / 2;
			if (key(a[md + 1]) < key(c[h - md])) lo = md + 1; else hi = md;
		}
		int i1 = 1, j1 = 1, i2 = lo, j2 = h - lo, i3 = lo + 1, j3 = h - lo + 1, i4 = na, j4 = nc;
		auto fw = [&](int& i, int& j, int k) { const int u = key(a[i]) < key(c[j]); out(k, u ? a[i] : c[j]), i += u, j += 1 - u; };
		auto bw = [&](int& i, int& j, int k) { const int u = key(a[i]) >= key(c[j]); out(k, u ? a[i] : c[j]), i -= u, j -= 1 - u; };
		const int L1 = h / 2, L2 = h - L1, L3 = (m - h) / 2, L4 = m - h - L3, mn = min(L1, L3);
		for (int k = 0; k < mn; k++) fw(i1, j1, k), bw(i2, j2, h - 1 - k), fw(i3, j3, h + k), bw(i4, j4, m - 1 - k);
		for (int k = mn; k < L1; k++) fw(i1, j1, k);
		for (int k = mn; k < L2; k++) bw(i2, j2, h - 1 - k);
		for (int k = mn; k < L3; k++) fw(i3, j3, h + k);
		for (int k = mn; k < L4; k++) bw(i4, j4, m - 1 - k);
	}
};
 