// @name CRT (Chinese Remainder Theorem)
// @prefix crt
// @description 中国剰余定理: lin_cong / crt（2式・n式）/ crt_mod（巨大 lcm を mod で）

// ===== CRT（中国剰余定理）=====
// 法は 1 以上の long long、剰余は long long 全域（負・法以上も可）。中間計算は溢れない
//   lin_cong(a, b, m)   : a*x ≡ b (mod m) → {x0, m/g}（解は x0 + k*m/g）。解なし {0, 0}
//   crt(r0, m0, r1, m1) : 2 式を併合 → {r, lcm}（0 <= r < lcm）。解なし {0, 0}
//   crt(r, m)           : n 式（法は互いに素でなくてよい）→ 同上（ACL と同じ規約）。n = 0 なら {0, 1}
//                         lcm >= 2^63 になると assert で落ちる → その可能性があるなら crt_mod
//   crt_mod(r, m, md)   : lcm が巨大でも可 → {x mod md, lcm mod md}。解なし {-1, -1}
//                         O(n log + k^2 log)。k は 64bit に収まる範囲で前から併合した後の式数（k = 2000 で 0.4 s 程度）
// 解は最小「非負」解。最小「正」の解が要るときは x == 0（⇔ 全 r_i ≡ 0）なら lcm に置き換える
namespace crt_lib {
using i64 = long long;
__extension__ typedef __int128 i128;

i64 safe_mod(i64 a, i64 m) { a %= m; return a < 0 ? a + m : a; }

// {g, x}: g = gcd(a, m), a*x ≡ g (mod m), 0 <= x < m/g
std::pair<i64, i64> inv_gcd(i64 a, i64 m) {
	a = safe_mod(a, m);
	if (a == 0) return {m, 0};
	i64 s = m, t = a, x0 = 0, x1 = 1;
	while (t) {
		i64 q = s / t;
		s -= t * q, x0 -= x1 * q;
		std::swap(s, t), std::swap(x0, x1);
	}
	if (x0 < 0) x0 += m / s;
	return {s, x0};
}

std::pair<i64, i64> lin_cong(i64 a, i64 b, i64 m) {
	assert(m >= 1);
	auto [g, x] = inv_gcd(a, m);
	if (b % g) return {0, 0};
	i64 u = m / g;
	return {(i64)((i128)safe_mod(b / g, u) * x % u), u};
}

std::pair<i64, i64> crt(i64 r0, i64 m0, i64 r1, i64 m1) {
	assert(m0 >= 1 && m1 >= 1);
	r0 = safe_mod(r0, m0), r1 = safe_mod(r1, m1);
	auto [k, u] = lin_cong(m0, r1 - r0, m1);  // r0 + m0*k ≡ r1 (mod m1)
	if (u == 0) return {0, 0};
	assert(m0 <= LLONG_MAX / u && "crt: lcm overflow");
	return {r0 + m0 * k, m0 * u};
}

// テンプレートなのは vector<int> 対応と、ACL（using namespace atcoder）との曖昧さ回避のため
template <class T = i64>
std::pair<i64, i64> crt(const std::vector<T>& r, const std::vector<T>& m) {
	assert(r.size() == m.size());
	i64 r0 = 0, m0 = 1;
	for (size_t i = 0; i < r.size() && m0; i++) std::tie(r0, m0) = crt(r0, m0, (i64)r[i], (i64)m[i]);
	return {r0, m0};
}

template <class T = i64>
std::pair<i64, i64> crt_mod(const std::vector<T>& rs, const std::vector<T>& ms, i64 md) {
	assert(rs.size() == ms.size() && md >= 1);
	// 1. lcm が 64bit に収まる限り前から併合して k 本にまとめる
	std::vector<i64> r, m;
	for (size_t i = 0; i < rs.size(); i++) {
		i64 ri = (i64)rs[i], mi = (i64)ms[i];
		assert(mi >= 1);
		if (!m.empty() && m.back() / std::gcd(m.back(), mi) <= LLONG_MAX / mi) {
			auto [x, l] = crt(r.back(), m.back(), ri, mi);
			if (!l) return {-1, -1};
			r.back() = x, m.back() = l;
		} else {
			r.push_back(safe_mod(ri, mi)), m.push_back(mi);
		}
	}
	int n = (int)r.size();
	// 2. 素因数ごとに指数の大きい側へ寄せて法を互いに素にする（解集合・lcm は不変）
	for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) {
		i64 g = std::gcd(m[i], m[j]);
		if (g == 1) continue;
		if ((r[i] - r[j]) % g) return {-1, -1};
		m[i] /= g, m[j] /= g;
		i64 gi = std::gcd(m[i], g), gj = g / gi;
		for (i64 h; (h = std::gcd(gi, gj)) != 1;) gi *= h, gj /= h;
		m[i] *= gi, m[j] *= gj, r[i] %= m[i], r[j] %= m[j];
	}
	// 3. Garner: x = Σ t_k Π_{j<k} m_j。cs[i], cf[i] は m_i（と md）での部分和・係数
	m.push_back(md);
	std::vector<i64> cs(n + 1, 0), cf(n + 1);
	for (int i = 0; i <= n; i++) cf[i] = 1 % m[i];
	for (int k = 0; k < n; k++) {
		if (m[k] == 1) continue;
		i64 t = (i64)((i128)safe_mod(r[k] - cs[k], m[k]) * inv_gcd(cf[k], m[k]).second % m[k]);
		for (int i = k + 1; i <= n; i++) {
			cs[i] = (i64)((cs[i] + (i128)t * cf[i]) % m[i]);
			cf[i] = (i64)((i128)cf[i] * m[k] % m[i]);
		}
	}
	return {cs[n], cf[n]};
}
}  // namespace crt_lib
using crt_lib::crt, crt_lib::lin_cong, crt_lib::crt_mod;
