// @name Warshall-Floyd
// @prefix warshall-floyd
// @description 全点対最短路（負辺OK・負閉路検出可）

// d[i][j] : i->j の最短距離（辺なしは LINF, d[i][i] = 0）O(V^3)
void warshall_floyd(vector<vector<ll>>& d) {
	int n = d.size();
	for (int k = 0; k < n; k++)
		for (int i = 0; i < n; i++) {
			if (d[i][k] == LINF) continue;
			for (int j = 0; j < n; j++) {
				if (d[k][j] == LINF) continue;
				d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
			}
		}
}
$0
