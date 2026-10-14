// @name Dijkstra's Algorithm
// @prefix dijkstra
// @description Dijkstra's algorithm for shortest path

//負の辺がないグラフの単一始点最短経路、ダイクストラの待遇ストの部分
vll dist(n, LINF);
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
}$0
