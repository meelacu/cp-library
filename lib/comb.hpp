// @name modint combination
// @prefix comb
// @description 二項係数 Comb<mint>：C/P/H/fact/inv_fact/inv。init 不要で表が自動で伸びる。巨大 n・小さい素数 mod・dmint 対応
 
// C(n, k) : nCk          P(n, k) : nPk = n!/(n-k)!          H(n, k) : nHk = C(n+k-1, k)（重複組合せ）
// fact(n) : n!（n ≥ mod なら 0）   inv_fact(n) : 1/n!（n < mod）   inv(n) : 1/n
// ・引数は long long。k < 0, k > n, n < 0 は 0（ただし H(0, 0) = 1）
// ・init もサイズ指定も不要。表は使った n まで自動で伸び、dmint で set_mod し直したら作り直す
// ・n < min(mod, LIM)          : 表を引いて O(1)（表づくりは合計 O(n)）
template <class mint>
struct Comb {
    static constexpr long long LIM = 10'000'000;  
    std::vector<mint> f{1}, fi{1};                // f[i] = i!, fi[i] = 1/i!
    long long md = 0;                             // 表を作ったときの mod
 
    static long long mod() { return (long long)mint::mod(); }
    void ensure(long long n) {  
        if (md == mod() && n < (long long)f.size()) return;
        grow(n);
    }
    void grow(long long n) {  
        if (md != mod()) md = mod(), f.assign(1, 1), fi.assign(1, 1);
        int sz = (int)f.size();
        if (n < sz) return;
        int m = (int)std::min(std::max(n + 1, std::min(2LL * sz, LIM)), md);
        f.resize(m), fi.resize(m);
        for (int i = sz; i < m; i++) f[i] = f[i - 1] * mint(i);
        fi[m - 1] = f[m - 1].inv();
        for (int i = m - 1; i > sz; i--) fi[i - 1] = fi[i] * mint(i);
    }
 
    mint fact(long long n) {
        assert(n >= 0);
        if (n >= mod()) return 0;
        ensure(n);
        return f[n];
    }
    mint inv_fact(long long n) {
        assert(0 <= n && n < mod());
        ensure(n);
        return fi[n];
    }
    mint inv(long long n) {
        if (n <= 0 || n >= std::min(mod(), LIM)) return mint(n).inv();
        ensure(n);
        return fi[n] * f[n - 1];
    }
 
    mint C(long long n, long long k) {
        if (k < 0 || k > n) return 0;
        if (n >= std::min(mod(), LIM)) return C_large(n, std::min(k, n - k));
        ensure(n);  // 先に伸ばす（式の中で伸ばすと、評価順しだいで伸ばす前の表を読む）
        return f[n] * fi[k] * fi[n - k];
    }
    mint C_large(long long n, long long k) {  // 表に乗らない n（k ≤ n - k）
        const long long p = mod();
        if (p < LIM) {  // Lucas：p 進の各桁 (n_i, k_i) について C(n_i, k_i) の積
            if (p == 1) return 0;  // mod 1 だと桁が減らずループが終わらない
            mint r = 1;
            for (; k > 0; n /= p, k /= p) {
                long long a = n % p, b = k % p;
                if (b > a) return 0;
                ensure(a);
                r *= f[a] * fi[b] * fi[a - b];
            }
            return r;
        }
        assert(k < p);  // k! が mod で割り切れると割れない
        mint num = 1, den = 1;
        for (long long i = 1; i <= k; i++) num *= mint(n - k + i), den *= mint(i);
        return num * den.inv();
    }
    mint P(long long n, long long k) {
        if (k < 0 || k > n) return 0;
        if (n >= std::min(mod(), LIM)) {
            if (k >= mod()) return 0;  // 連続する k 個の整数の積は mod の倍数を含む
            mint r = 1;
            for (long long i = 0; i < k; i++) r *= mint(n - i);
            return r;
        }
        ensure(n);
        return f[n] * fi[n - k];
    }
    mint H(long long n, long long k) {
        if (n < 0 || k < 0) return 0;
        return k == 0 ? mint(1) : C(n + k - 1, k);
    }
};
Comb<${1|mint,dmint|}> comb;
$0
