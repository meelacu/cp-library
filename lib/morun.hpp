// @name Mo's algorithm (run)
// @prefix morun
// @description Mo's algorithm の呼び出し雛形 (main 内に置く)

Mo mo(${1:N});
rep(i, ${2:Q}) {
	int l, r;
	cin >> l >> r;
	mo.insert(l - 1, r);  // 1-indexed 閉区間 → [l, r)
}
${3:ll} cur = 0;
vector<$3> ans($2);
mo.run(
	[&](int i) { $4 },  // add
	[&](int i) { $5 },  // del
	[&](int q) { ans[q] = cur; });
rep(i, $2) cout << ans[i] << '\n';
$0
