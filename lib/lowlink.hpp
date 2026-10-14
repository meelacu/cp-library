// @prefix lowlink
// @description low_link: 橋・関節点・二重辺連結成分 (非再帰, 多重辺OK)

// low_link: 橋・関節点・二重辺連結成分分解 (0-indexed / 非連結・多重辺・自己ループOK / 非再帰 / C++17)
//   low_link g(n); g.add_edge(u, v); g.build();
//   bridge        : 橋の辺id (昇順)      辺(u,v)が橋 <=> comp[u] != comp[v]
//   articulation  : 関節点 (昇順)        v が関節点 <=> cut_cnt[v] >= 2
//   cut_cnt[v]    : v を除去したとき v の属していた連結成分がいくつに分かれるか (孤立点は 0)
//   comp[v]       : 二重辺連結成分の id (0..num_comp-1)
//   bridge_tree() : 二重辺連結成分を縮約した森 (頂点は comp id, 辺は橋)
struct low_link {
	int n, m = 0;
	std::vector<std::vector<std::pair<int, int>>> g;  // (to, edge id)
	std::vector<std::pair<int, int>> edges;
	std::vector<int> ord, low, par, cut_cnt, comp, bridge, articulation;
	int num_comp = 0;

	low_link(int n) : n(n), g(n) {}

	int add_edge(int u, int v) {
		g[u].emplace_back(v, m);
		g[v].emplace_back(u, m);
		edges.emplace_back(u, v);
		return m++;
	}

	void build() {
		ord.assign(n, -1); low.assign(n, -1); par.assign(n, -1); cut_cnt.assign(n, 0);
		bridge.clear(); articulation.clear();
		std::vector<int> par_edge(n, -1), it(n, 0), pre, st;
		pre.reserve(n);
		int k = 0;
		for (int r = 0; r < n; r++) {
			if (ord[r] != -1) continue;
			ord[r] = low[r] = k++;
			pre.push_back(r); st.push_back(r);
			while (!st.empty()) {
				int v = st.back();
				if (it[v] < (int)g[v].size()) {
					auto [to, id] = g[v][it[v]++];
					if (id == par_edge[v]) continue;  // 親辺だけ除外 (多重辺は back edge 扱い)
					if (ord[to] == -1) {
						ord[to] = low[to] = k++;
						par[to] = v; par_edge[to] = id;
						pre.push_back(to); st.push_back(to);
					} else {
						low[v] = std::min(low[v], ord[to]);
					}
				} else {
					st.pop_back();
					int p = par[v];
					if (p == -1) continue;
					low[p] = std::min(low[p], low[v]);
					if (ord[p] <= low[v]) cut_cnt[p]++;
					if (ord[p] < low[v]) bridge.push_back(par_edge[v]);
				}
			}
		}
		for (int v = 0; v < n; v++) {
			if (par[v] != -1) cut_cnt[v]++;  // 非根は親側も1つの成分
			if (cut_cnt[v] >= 2) articulation.push_back(v);
		}
		std::sort(bridge.begin(), bridge.end());
		// 親辺が橋 <=> low[v] == ord[v] (根も含む)。行きがけ順に見るので親の comp は確定済み
		comp.assign(n, -1); num_comp = 0;
		for (int v : pre) comp[v] = (low[v] == ord[v]) ? num_comp++ : comp[par[v]];
	}

	std::vector<std::vector<int>> bridge_tree() const {
		std::vector<std::vector<int>> t(num_comp);
		for (int id : bridge) {
			auto [u, v] = edges[id];
			t[comp[u]].push_back(comp[v]); t[comp[v]].push_back(comp[u]);
		}
		return t;
	}
};
$0
