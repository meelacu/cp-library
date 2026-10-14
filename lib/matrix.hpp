// @name Matrix core
// @prefix matrix
// @description 行列の原型（半環ジェネリック、動的・固定サイズ、pow / pow_sum / det / inv / rank、セグ木対応）。半環は matmod / mattrop / matbool / matring

// ===== 行列の原型（半環は matmod / mattrop / matbool / matring で乗せる）=====
// Matrix<SR, FixN>: FixN = 0 で動的サイズ / FixN > 0 で FixN×FixN 固定（ヒープ確保なし）
// 半環 SR の要件: T, zero(), one(), add, mul, norm（生の値→正規形）。det/inv/rank は体のみで sub, inv も要る
// 演算時に各要素を norm で正規化するので要素には生の値を入れてよい。add/mul は範囲外を返してもよい（積和の最後に norm）
template <class SR, class = void> struct mat_has_mod : false_type {};
template <class SR> struct mat_has_mod<SR, void_t<decltype(SR::mod())>> : true_type {};
struct MatDynDim { int n = 0, m = 0; };
struct MatFixDim {};

template <class SR, int FixN = 0>
struct Matrix : conditional_t<FixN == 0, MatDynDim, MatFixDim> {
	using T = typename SR::T;
	template <class Elem, int Len> using Buf = conditional_t<FixN == 0, vector<Elem>, array<Elem, (FixN ? Len : 1)>>;
	Buf<T, FixN * FixN> d;

	Matrix() { if constexpr (FixN > 0) d.fill(SR::zero()); }  // 動的なら 0×0（* と + の単位元として振る舞う）
	explicit Matrix(int n_) : Matrix(n_, n_) {}
	Matrix(int n_, int m_) {
		if constexpr (FixN > 0) { assert(n_ == FixN && m_ == FixN); d.fill(SR::zero()); }
		else { assert(n_ >= 0 && m_ >= 0); this->n = n_, this->m = m_; d.assign((size_t)n_ * m_, SR::zero()); }
	}
	template <class Src> Matrix(const vector<vector<Src>>& v) : Matrix((int)v.size(), v.empty() ? 0 : (int)v[0].size()) {
		for (int i = 0; i < rows(); i++) {
			assert((int)v[i].size() == cols());
			for (int j = 0; j < cols(); j++) (*this)[i][j] = SR::norm(v[i][j]);
		}
	}
	Matrix(initializer_list<initializer_list<T>> l) : Matrix(vector<vector<T>>(l.begin(), l.end())) {}

	static Matrix identity(int n_ = FixN) {
		Matrix res(n_, n_);
		for (int i = 0; i < n_; i++) res[i][i] = SR::one();
		return res;
	}
	// セグ木用。op: 左の要素から順に適用（行ベクトル規約）, op_col: 列ベクトル規約用
	static Matrix e() { if constexpr (FixN > 0) return identity(); else return Matrix(); }
	static Matrix op(const Matrix& a, const Matrix& b) { return a * b; }
	static Matrix op_col(const Matrix& a, const Matrix& b) { return b * a; }

	int rows() const { if constexpr (FixN > 0) return FixN; else return this->n; }
	int cols() const { if constexpr (FixN > 0) return FixN; else return this->m; }
	bool empty() const { if constexpr (FixN > 0) return false; else return this->n == 0 && this->m == 0; }
	T* operator[](int i) { assert(0 <= i && i < rows()); return d.data() + (size_t)i * cols(); }
	const T* operator[](int i) const { assert(0 <= i && i < rows()); return d.data() + (size_t)i * cols(); }
	Matrix normed() const { Matrix C = *this; for (auto& x : C.d) x = SR::norm(x); return C; }

	Matrix operator*(const Matrix& B) const {
		if (empty()) return B.normed();
		if (B.empty()) return normed();
		assert(cols() == B.rows());
		if constexpr (mat_has_mod<SR>::value)
			if (SR::mod() < (1LL << 30)) return mul_mod(B);
		const int R = rows(), K = cols(), W = B.cols();
		const T* bp = B.d.data();
		Buf<T, FixN * FixN> tmp{};
		if (!all_of(B.d.begin(), B.d.end(), [](const T& x) { return SR::norm(x) == x; })) {
			if constexpr (FixN == 0) tmp.resize(B.d.size());
			for (size_t t = 0; t < B.d.size(); t++) tmp[t] = SR::norm(B.d[t]);
			bp = tmp.data();
		}
		Matrix C(R, W);
		for (int i = 0; i < R; i++) {
			const T* a = d.data() + (size_t)i * K;
			T* c = C.d.data() + (size_t)i * W;
			for (int k = 0; k < K; k++) {
				const T x = SR::norm(a[k]);
				if (x == SR::zero()) continue;
				const T* b = bp + (size_t)k * W;
				for (int j = 0; j < W; j++) c[j] = SR::add(c[j], SR::mul(x, b[j]));
			}
			for (int j = 0; j < W; j++) c[j] = SR::norm(c[j]);
		}
		return C;
	}
	Matrix operator+(const Matrix& B) const {
		if (empty()) return B.normed();
		if (B.empty()) return normed();
		assert(rows() == B.rows() && cols() == B.cols());
		Matrix C = *this;
		for (size_t t = 0; t < d.size(); t++) C.d[t] = SR::add(SR::norm(d[t]), SR::norm(B.d[t]));
		return C;
	}
	Matrix& operator*=(const Matrix& B) { return *this = *this * B; }
	Matrix& operator+=(const Matrix& B) { return *this = *this + B; }

	// 列ベクトル: (A v)[i] = Σ_j A[i][j] v[j]
	template <class Src> vector<T> operator*(const vector<Src>& v) const {
		const vector<T> nv = normed_vec(v);
		if (empty()) return nv;
		assert((int)v.size() == cols());
		vector<T> r(rows(), SR::zero());
		for (int i = 0; i < rows(); i++)
			for (int j = 0; j < cols(); j++) r[i] = SR::add(r[i], SR::mul(SR::norm((*this)[i][j]), nv[j]));
		return normed_vec(r);
	}
	// 行ベクトル: (v A)[j] = Σ_i v[i] A[i][j]
	template <class Src> friend vector<T> operator*(const vector<Src>& v, const Matrix& A) {
		if (A.empty()) return normed_vec(v);
		assert((int)v.size() == A.rows());
		vector<T> r(A.cols(), SR::zero());
		for (int i = 0; i < A.rows(); i++) {
			const T x = SR::norm(v[i]);
			if (x == SR::zero()) continue;
			for (int j = 0; j < A.cols(); j++) r[j] = SR::add(r[j], SR::mul(x, SR::norm(A[i][j])));
		}
		return normed_vec(r);
	}

	Matrix pow(long long p) const {
		assert(rows() == cols() && p >= 0);
		Matrix res = identity(rows()), base = *this;
		for (; p; p >>= 1) {
			if (p & 1) res *= base;
			if (p > 1) base *= base;
		}
		return res;
	}
	// I + A + A^2 + ... + A^(p-1)
	Matrix pow_sum(long long p) const {
		assert(rows() == cols() && p >= 0);
		if (p == 0) return Matrix(rows(), cols());
		Matrix sum = identity(rows()), pw = normed();  // sum = Σ_{i<q} A^i, pw = A^q（q = 1 から上位ビット順に伸ばす）
		for (int b = 62 - __builtin_clzll(p); b >= 0; b--) {
			sum += pw * sum, pw *= pw;
			if (p >> b & 1) sum += pw, pw *= *this;
		}
		return sum;
	}
	Matrix transpose() const {
		Matrix C(cols(), rows());
		for (int i = 0; i < rows(); i++)
			for (int j = 0; j < cols(); j++) C[j][i] = SR::norm((*this)[i][j]);
		return C;
	}
	bool operator==(const Matrix& B) const {
		if (rows() != B.rows() || cols() != B.cols()) return false;
		for (size_t t = 0; t < d.size(); t++)
			if (!(SR::norm(d[t]) == SR::norm(B.d[t]))) return false;
		return true;
	}
	bool operator!=(const Matrix& B) const { return !(*this == B); }
	friend ostream& operator<<(ostream& os, const Matrix& A) {
		for (int i = 0; i < A.rows(); i++)
			for (int j = 0; j < A.cols(); j++) os << SR::norm(A[i][j]) << (j + 1 == A.cols() ? "\n" : " ");
		return os;
	}

	// ---- 以下は体（素数 MOD の ModRing, Ring<double> など）専用 ----
	T det() const { assert(rows() == cols()); auto M = dyn(); return sweep(M, cols(), false).second; }
	int rank() const { auto M = dyn(); return sweep(M, cols(), false).first; }
	optional<Matrix> inv() const {  // 正則でなければ nullopt
		assert(rows() == cols());
		const int k = rows();
		Matrix<SR> M(k, 2 * k);
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) M[i][j] = SR::norm((*this)[i][j]);
			M[i][k + i] = SR::one();
		}
		if (sweep(M, k, true).first < k) return nullopt;
		Matrix res(k, k);
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k; j++) res[i][j] = M[i][k + j];
		return res;
	}

private:
	template <class Src> static vector<T> normed_vec(const vector<Src>& v) {
		vector<T> r(v.size());
		for (size_t i = 0; i < v.size(); i++) r[i] = SR::norm(v[i]);
		return r;
	}
	// mod < 2^30: 積和を 64bit にためて剰余は 16 項ごと
	Matrix mul_mod(const Matrix& B) const {
		const unsigned long long md = SR::mod();
		const int R = rows(), K = cols(), W = B.cols();
		Buf<unsigned, FixN * FixN> b{};
		Buf<unsigned long long, FixN> acc{};
		if constexpr (FixN == 0) b.resize(B.d.size()), acc.resize(W);
		for (size_t t = 0; t < B.d.size(); t++) b[t] = (unsigned)SR::norm(B.d[t]);
		Matrix C(R, W);
		for (int i = 0; i < R; i++) {
			fill(acc.begin(), acc.end(), 0ULL);
			const T* a = d.data() + (size_t)i * K;
			for (int k = 0; k < K; k++) {
				const unsigned long long x = (unsigned)SR::norm(a[k]);
				if (x) {
					const unsigned* bk = b.data() + (size_t)k * W;
					for (int j = 0; j < W; j++) acc[j] += x * bk[j];
				}
				if ((k & 15) == 15) for (auto& s : acc) s %= md;
			}
			T* c = C.d.data() + (size_t)i * W;
			for (int j = 0; j < W; j++) c[j] = (T)(acc[j] % md);
		}
		return C;
	}
	Matrix<SR> dyn() const {  // 正規化した動的サイズのコピー
		Matrix<SR> M(rows(), cols());
		for (int i = 0; i < rows(); i++)
			for (int j = 0; j < cols(); j++) M[i][j] = SR::norm((*this)[i][j]);
		return M;
	}
	// 先頭 lim 列で掃き出して {rank, det} を返す。full なら上側も消去
	static pair<int, T> sweep(Matrix<SR>& M, int lim, bool full) {
		const int R = M.rows(), W = M.cols();
		int r = 0;
		T det = SR::one();
		for (int c = 0; c < lim && r < R; c++) {
			int p = -1;
			if constexpr (is_floating_point_v<T>) {
				T best = 1e-12;
				for (int i = r; i < R; i++)
					if (abs(M[i][c]) > best) best = abs(M[i][c]), p = i;
			} else {
				for (int i = r; i < R && p < 0; i++)
					if (!(M[i][c] == SR::zero())) p = i;
			}
			if (p < 0) { det = SR::zero(); continue; }
			if (p != r) {
				for (int j = 0; j < W; j++) swap(M[p][j], M[r][j]);
				det = SR::sub(SR::zero(), det);
			}
			det = SR::mul(det, M[r][c]);
			const T iv = SR::inv(M[r][c]);
			for (int j = c; j < W; j++) M[r][j] = SR::mul(M[r][j], iv);
			for (int i = full ? 0 : r + 1; i < R; i++) {
				if (i == r || M[i][c] == SR::zero()) continue;
				const T f = M[i][c];
				for (int j = c; j < W; j++) M[i][j] = SR::sub(M[i][j], SR::mul(f, M[r][j]));
			}
			r++;
		}
		return {r, det};
	}
};

// ---- 使い方（Mat は半環スニペット側の using）----
// Mat A(n), B(n, m), F{{1, 1}, {1, 0}}, I = Mat::identity(n);  // 零行列で初期化
// A[i][j] = x;  // 生の値（負の数・MOD 以上・自前の INF）もそのまま入れてよい
// A * B, A + B, A *= B, A == B, A.pow(k), A.pow_sum(k)（= I + A + ... + A^(k-1)）, A.transpose(), cout << A
// v * A: 行ベクトル（A[from][to]）, A * v: 列ベクトル（A[to][from]）。v は vector
// A.det(), A.rank(), A.inv(): 体のみ。inv は非正則なら nullopt
// セグ木: SegTree<Mat, Mat::op, Mat::e> seg(v);  // prod(l, r) = A_l * ... * A_{r-1}（左から順に適用）
//         列ベクトル規約なら Mat::op_col。固定サイズ推奨。動的サイズの空区間 prod は 0×0（v * prod で使う）
// 自作半環: struct SR { using T = ...; static T zero(), one(), norm(T), add(T, T), mul(T, T); };
$0
