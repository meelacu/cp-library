// @name Breadth-First Search (BFS)
// @prefix bfs
// @description Breadth-First Search (Graph traversal)

// BFS: 始点 s からの最短距離を求める
vi dist(n, -1);
queue<int> q;

dist[s] = 0;
q.push(s);

while (!q.empty()) {
    int v = q.front();
    q.pop();

    for (int next_v : G[v]) {
        if (dist[next_v] != -1) continue;
        dist[next_v] = dist[v] + 1;
        q.push(next_v);
    }
}$0
