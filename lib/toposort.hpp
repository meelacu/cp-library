// @name Topological Sort
// @prefix toposort
// @description トポロジカルソート: 閉路判定・優先度つき(辞書順最小など)・一意性判定・閉路の復元

// 有向グラフのトポロジカルソート（Kahn 法）。main の外に置く
//   G[v] = v から出る辺のリスト。要素は 行き先(int) / pair・tuple・array（先頭が行き先）/ .to を持つ struct のどれでもよい
//
// vi ord = toposort(G);                  O(V + E)
//     辺 u → v があれば u が v より前に来る頂点の並び。
//     ord.size() == n なら閉路なし（DAG）。閉路があると n より短くなり、
//     「どの閉路からも到達できない頂点」だけがトポロジカル順に並ぶ
// vi ord = toposort(G, less<int>());     O((V + E) log V)
//     取り出せる頂点が複数あるとき、cmp(a, b) = true になる a を先に出す（cmp は狭義弱順序）。
//     less<int>() で辞書順最小、greater<int>() で辞書順最大。キーで比べるラムダも渡せる
// bool uq = toposort_unique(G, ord);     O(V + E)
//     トポロジカル順序がただ 1 通りか。ord は上のどちらの結果でもよい（閉路があれば false）
// vi cyc = find_cycle(G);                O(V + E)
//     閉路を 1 つ、辺の向きに沿った頂点列で返す（cyc[i] → cyc[i+1] → … → cyc[0]）。DAG なら空
//
// 例: 最長路（辺数）   vi dp(n); for (int v : ord) for (int to : G[v]) chmax(dp[to], dp[v] + 1);
//     各頂点から出るパスの数など「後ろから集める」DP は ord を逆順に回す
//     閉路に到達できる頂点（そこから無限に歩き続けられる頂点）= 逆辺グラフ RG で toposort(RG) に並ばなかった頂点
namespace toposort_internal {
template <class E>
int dest(const E& e) {
    if constexpr (is_integral_v<E>) return e;
    else if constexpr (requires { e.to; }) return e.to;
    else return get<0>(e);
}
template <class Graph>
vector<int> indegree(const Graph& G) {
    vector<int> deg(G.size());
    for (auto& es : G)
        for (auto& e : es) deg[dest(e)]++;
    return deg;
}
}  // namespace toposort_internal

template <class Graph>
vector<int> toposort(const Graph& G) {
    using toposort_internal::dest;
    int n = G.size();
    vector<int> deg = toposort_internal::indegree(G), ord;
    ord.reserve(n);
    for (int v = 0; v < n; v++)
        if (deg[v] == 0) ord.push_back(v);
    for (int i = 0; i < (int)ord.size(); i++)
        for (auto& e : G[ord[i]])
            if (--deg[dest(e)] == 0) ord.push_back(dest(e));
    return ord;
}

template <class Graph, class Cmp>
vector<int> toposort(const Graph& G, Cmp cmp) {
    using toposort_internal::dest;
    int n = G.size();
    vector<int> deg = toposort_internal::indegree(G), ord;
    ord.reserve(n);
    auto later = [&](int a, int b) { return cmp(b, a); };
    priority_queue<int, vector<int>, decltype(later)> pq(later);
    for (int v = 0; v < n; v++)
        if (deg[v] == 0) pq.push(v);
    while (!pq.empty()) {
        int v = pq.top();
        pq.pop();
        ord.push_back(v);
        for (auto& e : G[v])
            if (--deg[dest(e)] == 0) pq.push(dest(e));
    }
    return ord;
}

// DAG のトポロジカル順序が一意 ⇔ 順序の隣どうしがすべて辺で結ばれている（ハミルトン路になっている）
template <class Graph>
bool toposort_unique(const Graph& G, const vector<int>& ord) {
    using toposort_internal::dest;
    if (ord.size() != G.size()) return false;
    for (int i = 0; i + 1 < (int)ord.size(); i++) {
        bool ok = false;
        for (auto& e : G[ord[i]])
            if (dest(e) == ord[i + 1]) {
                ok = true;
                break;
            }
        if (!ok) return false;
    }
    return true;
}

template <class Graph>
vector<int> find_cycle(const Graph& G) {
    using toposort_internal::dest;
    int n = G.size();
    vector<int> ord = toposort(G);
    if ((int)ord.size() == n) return {};
    vector<char> alive(n, 1);  // トポロジカル順に並ばなかった頂点 = 閉路上か、閉路から到達できる頂点
    for (int v : ord) alive[v] = 0;
    vector<int> pre(n, -1);    // pre[v] = v に入る辺の始点のうち alive なもの（alive な v には必ずある）
    for (int u = 0; u < n; u++)
        if (alive[u])
            for (auto& e : G[u]) pre[dest(e)] = u;
    int v = find(alive.begin(), alive.end(), 1) - alive.begin();
    for (int i = 0; i < n; i++) v = pre[v];  // 逆向きに n 回たどれば必ず閉路の上にいる
    vector<int> cyc;
    for (int u = v;;) {
        cyc.push_back(u);
        u = pre[u];
        if (u == v) break;
    }
    reverse(cyc.begin(), cyc.end());
    return cyc;
}
$0
