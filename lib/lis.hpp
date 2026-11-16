// @name LIS (Longest Increasing Subsequence)
// @prefix lis
// @description 最長増加部分列 O(n log n): 復元(添字)・各位置で終わる/始まる最長の長さ。狭義/広義、比較関数で減少列も

// 最長増加部分列（LIS）。どれも O(n log n)。main の外に置く
//   a は vector でも string でもよい
//   vi lis(a)                    狭義増加（a[i] < a[j]）の最長部分列の添字（昇順）。長さは .size()
//   vi lis(a, false)             広義増加（a[i] <= a[j]）
//   vi lis(a, strict, cmp)       cmp で「小さい」を決める。greater<>() で減少列、pair の片方だけで比べるラムダなども可
//   vi lis_len_at(a, ...)        len[i] = a[i] で終わる最長の長さ（引数は lis と同じ）
//   vi lis_len_from(a, ...)      len[i] = a[i] から始まる最長の長さ
// 例: a[i] を含む最長増加部分列がある ⇔ lis_len_at(a)[i] + lis_len_from(a)[i] - 1 == lis(a).size()
//     箱の入れ子（縦横とも真に小さい）: 縦の昇順・横の降順に並べて、横で狭義の lis
namespace lis_internal {
// pos[i] = a[i] を末尾にしたときの長さ - 1、prv[i] = 直前の添字（なければ -1）
template <class Vec, class Cmp>
pair<vector<int>, vector<int>> run(const Vec& a, bool strict, Cmp cmp) {
    using T = typename Vec::value_type;
    int n = a.size();
    vector<T> tail;       // tail[k] = 長さ k+1 の増加部分列の末尾として最小のもの
    vector<int> tid;      // tail[k] の添字
    vector<int> pos(n), prv(n);
    for (int i = 0; i < n; i++) {
        auto it = strict ? lower_bound(tail.begin(), tail.end(), a[i], cmp) : upper_bound(tail.begin(), tail.end(), a[i], cmp);
        int k = it - tail.begin();
        if (it == tail.end()) tail.push_back(a[i]), tid.push_back(i);
        else *it = a[i], tid[k] = i;
        pos[i] = k;
        prv[i] = k ? tid[k - 1] : -1;
    }
    return {pos, prv};
}
}  // namespace lis_internal

template <class Vec, class Cmp = less<>>
vector<int> lis(const Vec& a, bool strict = true, Cmp cmp = Cmp()) {
    auto [pos, prv] = lis_internal::run(a, strict, cmp);
    int n = a.size(), last = -1;
    for (int i = 0; i < n; i++)
        if (last < 0 || pos[i] > pos[last]) last = i;
    vector<int> res;
    for (int i = last; i >= 0; i = prv[i]) res.push_back(i);
    reverse(res.begin(), res.end());
    return res;
}

template <class Vec, class Cmp = less<>>
vector<int> lis_len_at(const Vec& a, bool strict = true, Cmp cmp = Cmp()) {
    vector<int> len = lis_internal::run(a, strict, cmp).first;
    for (int& x : len) x++;
    return len;
}

template <class Vec, class Cmp = less<>>
vector<int> lis_len_from(const Vec& a, bool strict = true, Cmp cmp = Cmp()) {
    using T = typename Vec::value_type;
    vector<T> b(a.rbegin(), a.rend());  // 逆順にして「大きい方へ」を「小さい方へ」に読み替える
    vector<int> len = lis_len_at(b, strict, [&](const T& x, const T& y) { return cmp(y, x); });
    reverse(len.begin(), len.end());
    return len;
}
$0
