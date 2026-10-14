// @name SPF (Smallest Prime Factor)
// @prefix spf
// @description 最小素因数テーブルでO(log n)素因数分解

struct SPF {
	vector<int> spf;
	SPF(int n) : spf(n + 1) {
		iota(spf.begin(), spf.end(), 0);
		for (int i = 2; (long long)i * i <= n; i++) {
			if (spf[i] == i)
				for (int j = i * i; j <= n; j += i)
					if (spf[j] == j) spf[j] = i;
		}
	}
	// 素因数分解: O(log x), 戻り値 {素因数, 指数}
	vector<pair<int, int>> factorize(int x) {
		vector<pair<int, int>> res;
		while (x > 1) {
			int p = spf[x], e = 0;
			while (x % p == 0) x /= p, e++;
			res.emplace_back(p, e);
		}
		return res;
	}
	// 素因数列挙(重複あり): O(log x)
	vector<int> primes(int x) {
		vector<int> res;
		while (x > 1) { res.push_back(spf[x]); x /= spf[x]; }
		return res;
	}
};
$0
