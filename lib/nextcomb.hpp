// @name nCk 全探索 (index, next_permutation 風)
// @prefix nextcomb
// @description k個の添字の組合せを辞書順に列挙（next_permutation の nCk 版）

// nCk 全探索（添字版, n の制限なし）: idx[0] < idx[1] < ... < idx[k-1] を辞書順に列挙
// 使い方: vector<int> idx(k); iota(idx.begin(), idx.end(), 0);
//         do { ... } while (next_combination(idx, n));
// 1ステップ最悪 O(k)、全体almost O(nCk)
bool next_combination(vector<int>& idx, int n) {
	int k = idx.size();
	for (int i = k - 1; i >= 0; i--) {
		if (idx[i] < n - k + i) {
			idx[i]++;
			for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
			return true;
		}
	}
	return false;
}
$0
