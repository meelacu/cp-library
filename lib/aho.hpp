// @prefix aho
// @description Aho-Corasick (requires trie): build/go/hit/out/bfs/count_matches

// ===== Aho-Corasick (Trie の拡張) =====
// パターンを insert し終えてから build()
// build() 後は nxt が全遷移で埋まる (オートマトン化) ので
// Trie の find/count/walk/kth/count_less/lcp/build_order は使えない (cnt/term/par/pch/dep/str は可)
template <int K = 26, char B = 'a'>
struct AhoCorasick : Trie<K, B> {
    using T = Trie<K, B>;
    using T::nxt; using T::term; using T::idx;
    vector<int> fail;  // 失敗リンク (根は 0)
    vector<int> hit;   // hit[v] = fail 鎖上 (v 含む) の終端本数 = 状態 v で終わるパターン出現数
    vector<int> out;   // fail 鎖上で v より上の最も近い終端頂点 (-1: なし). マッチ列挙用
    vector<int> bfs;   // BFS 順 (根含む). 逆順に舐めると fail 方向へ集約できる
    void build() {
        int n = T::size();
        fail.assign(n, 0); hit.assign(n, 0); out.assign(n, -1);
        bfs.clear(); bfs.reserve(n); bfs.push_back(0);
        for (int i = 0; i < (int)bfs.size(); i++) {
            int v = bfs[i];
            hit[v] = term[v] + (v ? hit[fail[v]] : 0);
            if (v) out[v] = term[fail[v]] ? fail[v] : out[fail[v]];
            for (int c = 0; c < K; c++) {
                int u = nxt[v][c];
                if (u == -1) { nxt[v][c] = v ? nxt[fail[v]][c] : 0; continue; }
                fail[u] = v ? nxt[fail[v]][c] : 0;
                bfs.push_back(u);
            }
        }
    }
    int go(int v, char ch) const { return nxt[v][idx(ch)]; }  // 状態遷移
    long long count_matches(string_view s) const {  // s 中のパターン総出現数 (重複込み)
        int v = 0; long long r = 0;
        for (char ch : s) { v = go(v, ch); r += hit[v]; }
        return r;
    }
};
