// @name Inversion Number
// @prefix inversion
// @description 転倒数 O(n log n)(座圧不要・任意の比較可能な型)、要素ごとの転倒数、隣接swapで a を b に並べ替える最小回数(重複可)

// 転倒数（マージソートで数えるので座標圧縮は要らない）。main の外に置く
//   ll inversion(a)            i < j かつ a[i] > a[j] の組の数（同じ値どうしは数えない）。O(n log n)
//   ll inversion(a, cmp)       i < j かつ cmp(a[j], a[i]) の組の数（greater<>() なら a[i] < a[j] の組）
//   auto [L, R] = inversion_each(a);   要素ごとの内訳（vector<int> 2 本）。O(n log n)
//                              L[i] = i より左にある、a[i] より大きい要素の数
//                              R[i] = i より右にある、a[i] より小さい要素の数      （inversion(a) = L の総和 = R の総和）
//                              inversion_each(a, greater<>()) なら L = 左の小さい要素の数、R = 右の大きい要素の数
//   ll swap_distance(a, b)     隣接 swap で a を b に並べ替える最小回数（同じ値があってもよい）
//                              並べ替えても一致しない（要素の多重集合が違う）なら -1。O(n log n)
// 例: バブルソートの交換回数 = inversion(a)、文字列 s を t にする隣接 swap の最小回数 = swap_distance(s, t)
//     順列 p が辞書順で何番目か（0 始まり）= Σ R[i] * (n - 1 - i)!
template <class Vec, class Cmp = less<>>
long long inversion(Vec a, Cmp cmp = Cmp()) {
    int n = a.size();
    Vec buf(a);
    long long res = 0;
    for (int w = 1; w < n; w *= 2) {
        for (int l = 0; l + w < n; l += 2 * w) {
            int m = l + w, r = min(l + 2 * w, n), i = l, j = m, k = l;
            while (i < m && j < r) {
                if (cmp(a[j], a[i])) buf[k++] = a[j++], res += m - i;
                else buf[k++] = a[i++];
            }
            while (i < m) buf[k++] = a[i++];
            while (j < r) buf[k++] = a[j++];
            for (int t = l; t < r; t++) a[t] = buf[t];
        }
    }
    return res;
}

template <class Vec, class Cmp = less<>>
pair<vector<int>, vector<int>> inversion_each(const Vec& a, Cmp cmp = Cmp()) {
    int n = a.size();
    vector<int> id(n), buf(n), L(n), R(n);  // id: 添字を値の昇順に並べていく（マージソート）
    iota(id.begin(), id.end(), 0);
    for (int w = 1; w < n; w *= 2) {
        for (int l = 0; l + w < n; l += 2 * w) {
            int m = l + w, r = min(l + 2 * w, n), i = l, j = m, k = l;
            while (i < m || j < r) {
                if (j < r && (i == m || cmp(a[id[j]], a[id[i]]))) L[id[j]] += m - i, buf[k++] = id[j++];  // 左半分の残りは全部これより大きい
                else R[id[i]] += j - m, buf[k++] = id[i++];                                             // 右半分で先に出たものは全部これより小さい
            }
            for (int t = l; t < r; t++) id[t] = buf[t];
        }
    }
    return {L, R};
}

template <class Vec>
long long swap_distance(const Vec& a, const Vec& b) {
    int n = a.size();
    if ((int)b.size() != n) return -1;
    // 同じ値どうしは順番を保って対応させるのが最小（a の k 個目の x → b の k 個目の x）
    vector<int> ia(n), ib(n);
    iota(ia.begin(), ia.end(), 0);
    iota(ib.begin(), ib.end(), 0);
    stable_sort(ia.begin(), ia.end(), [&](int x, int y) { return a[x] < a[y]; });
    stable_sort(ib.begin(), ib.end(), [&](int x, int y) { return b[x] < b[y]; });
    vector<int> to(n);
    for (int k = 0; k < n; k++) {
        if (a[ia[k]] < b[ib[k]] || b[ib[k]] < a[ia[k]]) return -1;
        to[ia[k]] = ib[k];
    }
    return inversion(to);
}
$0
