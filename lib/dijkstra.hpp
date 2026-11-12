// @name Dijkstra's Algorithm
// @prefix dijkstra
// @description Dijkstra: 負辺のないグラフの単一始点最短路 (関数, main の外に置く)

// 負の辺がないグラフの単一始点最短経路 O((V + E) log V)
//   G[v] = {(行き先, 重み), ...} の隣接リスト（vector<vector<pair<int, ll>>> など）
//   戻り値 dist[v] : s から v への最短距離（到達不能は LINF）
//   vll dist = dijkstra(G, s);   関数なので、始点を変える・逆辺グラフを渡すなど何度でも呼べる
template <class Graph>
vll dijkstra(const Graph& G, int s) {
    vll dist(G.size(), LINF);
    dist[s] = 0;
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<pair<ll,int>>> pq;
    pq.push({0, s});

    while (!pq.empty()) {
        auto [d, v] = pq.top();
        pq.pop();
        if (dist[v] < d) continue;

        for (auto [next_v, cost] : G[v]) {
            if (dist[next_v] > dist[v] + cost) {
                dist[next_v] = dist[v] + cost;
                pq.push({dist[next_v], next_v});
            }
        }
    }
    return dist;
}
$0