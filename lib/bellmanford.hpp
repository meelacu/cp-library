// @name Bellman-Ford
// @prefix bellmanford
// @prefix bf
// @description Bellman-Ford (負辺/負閉路検出/-INF伝播/経路復元)

// ---------- Bellman-Ford ----------
// 負辺OK / 負閉路検出 / -INF伝播 / 経路復元
// build(s): O(V*E) を 2V 周 → O(V*E)
// 全体の負閉路検出は 仮想始点 n を作り n->v をコスト0で張る
struct BellmanFord {
	using ll = long long;
	static constexpr ll INF  = std::numeric_limits<ll>::max() / 4;
	static constexpr ll NINF = -INF;
	struct Edge { int from, to; ll cost; };

	int n;
	std::vector<Edge> es;
	std::vector<ll> dist;
	std::vector<int> pre;   // pre[v] = v へ入る辺の index
	bool neg_cycle = false; // s から到達可能な負閉路があるか

	explicit BellmanFord(int n_ = 0) : n(n_) {}
	void resize(int n_) { n = n_; }
	void add_edge(int u, int v, ll w) { es.push_back({u, v, w}); }
	void add_undirected(int u, int v, ll w) { add_edge(u, v, w); add_edge(v, u, w); }

	// 負閉路が無ければ true
	bool build(int s) {
		dist.assign(n, INF);
		pre.assign(n, -1);
		dist[s] = 0;
		neg_cycle = false;
		for (int it = 0; it < 2 * n; ++it) {
			bool upd = false;
			for (int i = 0; i < (int)es.size(); ++i) {
				const Edge& e = es[i];
				if (dist[e.from] == INF) continue;
				ll nd = (dist[e.from] == NINF ? NINF : dist[e.from] + e.cost);
				if (nd < dist[e.to]) {
					dist[e.to] = (it >= n - 1 ? NINF : nd);
					pre[e.to]  = i;
					upd = true;
					if (it >= n - 1) neg_cycle = true;
				}
			}
			if (!upd) break;
		}
		return !neg_cycle;
	}

	ll operator[](int v) const { return dist[v]; }
	bool unreachable(int v) const { return dist[v] == INF; }
	bool minus_inf(int v)  const { return dist[v] == NINF; } // 負閉路の影響下

	// s->t の最短路(頂点列)。到達不能 or -INF なら空
	std::vector<int> restore_path(int s, int t) const {
		if (dist[t] == INF || dist[t] == NINF) return {};
		std::vector<int> path;
		for (int v = t; v != s; v = es[pre[v]].from) path.push_back(v);
		path.push_back(s);
		std::reverse(path.begin(), path.end());
		return path;
	}
};
$0
