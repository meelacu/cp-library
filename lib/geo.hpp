// @name geometry (integer)
// @prefix geo
// @description 整数幾何: Pt/cross/ccw/ccw5/normalize_dir/偏角ソート/線分交差/面積/凸判定/内外判定/凸包/直径/格子点(Pick)/最近点対

// 整数幾何テンプレート。座標は T（long long）、浮動小数点なし。使い方: using namespace geo;
// 末尾の 2: 距離系 (dist2/norm2/diameter2/closest_pair2) は 2 乗、面積系 (area2) は 2 倍
// |座標| ≤ 1e9 なら cross/dot/dist2 は long long に収まる。不安なら T を __int128 に（I/O はそのまま動く）
namespace geo {

using T = long long;

T iabs(T v) { return v < 0 ? -v : v; }
int sgn(T v) { return (v > 0) - (v < 0); }
T gcd_(T a, T b) { while (b) { a %= b; std::swap(a, b); } return a; }

// ---------------------------------------------------------------- 点・ベクトル
struct Pt {
	T x, y;
	Pt(T x = 0, T y = 0) : x(x), y(y) {}
	Pt operator+(const Pt& o) const { return {x + o.x, y + o.y}; }
	Pt operator-(const Pt& o) const { return {x - o.x, y - o.y}; }
	Pt operator*(T k) const { return {x * k, y * k}; }
	Pt operator-() const { return {-x, -y}; }
	Pt& operator+=(const Pt& o) { x += o.x; y += o.y; return *this; }
	Pt& operator-=(const Pt& o) { x -= o.x; y -= o.y; return *this; }
	bool operator==(const Pt& o) const { return x == o.x && y == o.y; }
	bool operator!=(const Pt& o) const { return !(*this == o); }
	bool operator<(const Pt& o) const { return x != o.x ? x < o.x : y < o.y; }  // x → y の辞書順
	T norm2() const { return x * x + y * y; }  // |p|^2
};
// I/O は long long 経由（T を __int128 にしても動くように）
std::istream& operator>>(std::istream& is, Pt& p) { long long x, y; is >> x >> y; p = Pt(x, y); return is; }
std::ostream& operator<<(std::ostream& os, const Pt& p) { return os << (long long)p.x << ' ' << (long long)p.y; }

T dot(const Pt& a, const Pt& b) { return a.x * b.x + a.y * b.y; }
T cross(const Pt& a, const Pt& b) { return a.x * b.y - a.y * b.x; }
// 3 点版: cross(b - a, c - a)。a→b→c が左折なら正・右折なら負・共線なら 0（= 三角形 abc の符号付き面積 ×2）
T cross(const Pt& a, const Pt& b, const Pt& c) { return cross(b - a, c - a); }
T dist2(const Pt& a, const Pt& b) { return (a - b).norm2(); }
Pt rot90(const Pt& p) { return {-p.y, p.x}; }  // 反時計回りに 90° 回転

// +1: 反時計回り, -1: 時計回り, 0: 共線
int ccw(const Pt& a, const Pt& b, const Pt& c) { return sgn(cross(a, b, c)); }

// AOJ CGL 互換の 5 値版（共線のケースを細かく区別する）
//   +1: a→b→c が反時計回り   -1: 時計回り
//   +2: c-a-b の順に一直線 (ONLINE_BACK)   -2: a-b-c の順 (ONLINE_FRONT)
//    0: c が線分 ab 上 (ON_SEGMENT)
int ccw5(const Pt& a, const Pt& b, const Pt& c) {
	Pt u = b - a, v = c - a;
	T cr = cross(u, v);
	if (cr > 0) return +1;
	if (cr < 0) return -1;
	if (dot(u, v) < 0) return +2;
	if (u.norm2() < v.norm2()) return -2;
	return 0;
}

// 方向ベクトルの正規化: gcd で割って既約にする。「同じ傾き・同じ直線」を map のキーで同一視する用
// undirected=true なら逆向きも同一視する（x > 0, または x == 0 かつ y > 0 に揃える）。零ベクトルはそのまま
Pt normalize_dir(Pt d, bool undirected = true) {
	T g = gcd_(iabs(d.x), iabs(d.y));
	if (g == 0) return d;
	d.x /= g, d.y /= g;
	if (undirected && (d.x < 0 || (d.x == 0 && d.y < 0))) d = -d;
	return d;
}

// ---------------------------------------------------------------- 偏角ソート
// 偏角 [0, 2π) の昇順に並べる比較関数（浮動小数点なし）。原点は最初に来る
// 使い方: sort(v.begin(), v.end(), arg_less);
//        点 c の周りに並べるなら sort(v.begin(), v.end(), [&](auto& a, auto& b){ return arg_less(a - c, b - c); });
int arg_half(const Pt& p) {
	if (p.x == 0 && p.y == 0) return -1;
	return (p.y > 0 || (p.y == 0 && p.x > 0)) ? 0 : 1;
}
bool arg_less(const Pt& a, const Pt& b) {
	int ha = arg_half(a), hb = arg_half(b);
	if (ha != hb) return ha < hb;
	return cross(a, b) > 0;
}

// ---------------------------------------------------------------- 線分
// 点 p が線分 ab 上にあるか（端点を含む）
bool on_segment(const Pt& a, const Pt& b, const Pt& p) { return ccw5(a, b, p) == 0; }
// 線分 ab と線分 cd が共有点を持つか（端点で接する・共線で重なる場合も true）
bool seg_intersect(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
	return ccw5(a, b, c) * ccw5(a, b, d) <= 0 && ccw5(c, d, a) * ccw5(c, d, b) <= 0;
}
// 線分 ab と線分 cd が内部で真に交差するか（端点接触・共線は false）
bool seg_cross(const Pt& a, const Pt& b, const Pt& c, const Pt& d) {
	return ccw(a, b, c) * ccw(a, b, d) < 0 && ccw(c, d, a) * ccw(c, d, b) < 0;
}

// ---------------------------------------------------------------- 多角形
using Polygon = std::vector<Pt>;

// 符号付き面積の 2 倍（反時計回りなら正）
T area2(const Polygon& p) {
	T s = 0;
	int n = p.size();
	for (int i = 0; i < n; i++) s += cross(p[i], p[(i + 1) % n]);
	return s;
}
// 時計回りなら反時計回りに直す
void make_ccw(Polygon& p) {
	if (area2(p) < 0) std::reverse(p.begin(), p.end());
}

// 凸判定（単純多角形前提。向きは問わない）。strict=true なら共線な 3 頂点も不可
bool is_convex(const Polygon& p, bool strict = false) {
	int n = p.size();
	bool pos = false, neg = false;
	for (int i = 0; i < n; i++) {
		int s = ccw(p[i], p[(i + 1) % n], p[(i + 2) % n]);
		if (s > 0) pos = true;
		else if (s < 0) neg = true;
		else if (strict) return false;
	}
	return !(pos && neg);
}

// 点 q と多角形の位置関係: 2=内部, 1=境界上, 0=外部。凸でなくてよい。O(N)
int contains(const Polygon& p, const Pt& q) {
	int n = p.size();
	bool in = false;
	for (int i = 0; i < n; i++) {
		Pt a = p[i] - q, b = p[(i + 1) % n] - q;
		if (a.y > b.y) std::swap(a, b);
		if (a.y <= 0 && 0 < b.y && cross(a, b) < 0) in = !in;
		if (cross(a, b) == 0 && dot(a, b) <= 0) return 1;
	}
	return in ? 2 : 0;
}

// 凸多角形版: 2=内部, 1=境界上, 0=外部。O(log N)
// 前提: 反時計回り・N ≥ 3・共線な頂点がない（怪しければ convex_hull(p, true) を通してから使う）
int contains_convex(const Polygon& p, const Pt& q) {
	int n = p.size();
	T c1 = cross(p[0], p[1], q), c2 = cross(p[0], p[n - 1], q);
	if (c1 < 0 || c2 > 0) return 0;
	// cross(p[0], p[lo], q) >= 0 を満たす最大の lo (1 ≤ lo ≤ n-2) を二分探索
	int lo = 1, hi = n - 1;
	while (hi - lo > 1) {
		int mid = (lo + hi) / 2;
		(cross(p[0], p[mid], q) >= 0 ? lo : hi) = mid;
	}
	T c3 = cross(p[lo], p[lo + 1], q);  // 三角形 p[0], p[lo], p[lo+1] の外側の辺で判定
	if (c3 < 0) return 0;
	if (c3 == 0 || c1 == 0 || c2 == 0) return 1;
	return 2;
}

// 凸包（Andrew's monotone chain）。反時計回り、(x, y) 最小の点から始まる。O(N log N)
// strict=true なら辺上の点を含めない。重複点は除去される
Polygon convex_hull(Polygon p, bool strict = true) {
	std::sort(p.begin(), p.end());
	p.erase(std::unique(p.begin(), p.end()), p.end());
	int n = p.size();
	if (n <= 2) return p;
	auto bad = [&](const Pt& a, const Pt& b, const Pt& c) {
		T cr = cross(a, b, c);
		return strict ? cr <= 0 : cr < 0;
	};
	Polygon h;
	for (int i = 0; i < n; i++) {  // 下側
		while (h.size() >= 2 && bad(h[h.size() - 2], h.back(), p[i])) h.pop_back();
		h.push_back(p[i]);
	}
	for (int i = n - 2, t = h.size() + 1; i >= 0; i--) {  // 上側
		while ((int)h.size() >= t && bad(h[h.size() - 2], h.back(), p[i])) h.pop_back();
		h.push_back(p[i]);
	}
	h.pop_back();
	if (!strict && area2(h) == 0) return p;  // 全点が共線のときの往復を防ぐ
	return h;
}

// 点集合の最遠点対（直径）。凸包 + 回転カリパス O(N log N)。凸でなくてよい、N ≥ 1
std::pair<Pt, Pt> farthest_pair(const std::vector<Pt>& pts) {
	Polygon p = convex_hull(pts, true);  // 共線頂点を落として反時計回りにする（カリパスの前提）
	int n = p.size();
	if (n <= 2) return {p.front(), p.back()};
	std::pair<Pt, Pt> res = {p[0], p[1]};
	T best = -1;
	// 辺 e_i に対し、e_j が反対向き（cross(e_j, e_i) >= 0）になるまで j を進める。i < j の対だけ見れば十分
	for (int i = 0, j = 1; i < j; i++) {
		for (;; j = (j + 1) % n) {
			T d = dist2(p[i], p[j]);
			if (d > best) best = d, res = {p[i], p[j]};
			if (cross(p[(j + 1) % n] - p[j], p[(i + 1) % n] - p[i]) >= 0) break;
		}
	}
	return res;
}
// 点集合の直径の 2 乗
T diameter2(const std::vector<Pt>& pts) {
	auto [a, b] = farthest_pair(pts);
	return dist2(a, b);
}

// ---------------------------------------------------------------- 格子点
// 多角形の境界上の格子点の個数
T lattice_boundary(const Polygon& p) {
	T b = 0;
	int n = p.size();
	for (int i = 0; i < n; i++) {
		Pt d = p[(i + 1) % n] - p[i];
		b += gcd_(iabs(d.x), iabs(d.y));
	}
	return b;
}
// 単純多角形の内部の格子点の個数（Pick の定理: S = I + B/2 - 1）
T lattice_interior(const Polygon& p) { return (iabs(area2(p)) - lattice_boundary(p) + 2) / 2; }

// ---------------------------------------------------------------- 最近点対
// 最近点対の距離の 2 乗。分割統治 O(N log N)。N ≥ 2
T closest_pair2(std::vector<Pt> p) {
	int n = p.size();
	std::sort(p.begin(), p.end());
	T best = dist2(p[0], p[1]);
	auto rec = [&](auto& self, int l, int r) -> std::vector<Pt> {  // [l, r) を y 順に並べて返す
		if (r - l <= 1) return std::vector<Pt>(p.begin() + l, p.begin() + r);
		int m = (l + r) / 2;
		T xm = p[m].x;
		auto L = self(self, l, m), R = self(self, m, r);
		std::vector<Pt> res;
		res.reserve(r - l);
		std::merge(L.begin(), L.end(), R.begin(), R.end(), std::back_inserter(res),
				   [](const Pt& a, const Pt& b) { return a.y < b.y; });
		std::vector<Pt> strip;
		for (const Pt& q : res) {
			if ((q.x - xm) * (q.x - xm) >= best) continue;
			for (int k = (int)strip.size() - 1; k >= 0 && (q.y - strip[k].y) * (q.y - strip[k].y) < best; k--)
				best = std::min(best, dist2(q, strip[k]));
			strip.push_back(q);
		}
		return res;
	};
	rec(rec, 0, n);
	return best;
}

}  using namespace geo;
