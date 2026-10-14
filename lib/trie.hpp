// @prefix trie
// @description Trie (insert/erase with hook, find/count/count_prefix/lcp/walk, count_less/kth/str, build_order)

// ===== Trie =====
// K 種の文字を idx(c) = c - B で [0, K) に写す
//   小文字 <26,'a'> / 数字 <10,'0'> / 英大小 <58,'A'> ('A'..'z' を連続扱い)
// 頂点 0 が根。erase は cnt/term を減らすだけで頂点は残す
// メモリ: 頂点あたり約 4K+20 byte
template <int K = 26, char B = 'a'>
struct Trie {
    static constexpr int idx(char c) { return c - B; }
    static constexpr char chr(int c) { return (char)(B + c); }
    vector<array<int, K>> nxt;  // 子 (-1: なし)
    vector<int> cnt;   // v を接頭辞に持つ文字列の本数 (cnt[0] = 総数)
    vector<int> term;  // v で終わる文字列の本数
    vector<int> par;   // 親 (根は -1)
    vector<int> pch;   // 親から v への文字 (idx 済み)
    vector<int> dep;   // 深さ = 対応する文字列の長さ
    Trie() { new_node(-1, -1); }
    int size() const { return (int)nxt.size(); }
    void reserve(int n) {
        nxt.reserve(n); cnt.reserve(n); term.reserve(n);
        par.reserve(n); pch.reserve(n); dep.reserve(n);
    }
    int new_node(int p, int c) {
        nxt.push_back({}); nxt.back().fill(-1);
        cnt.push_back(0); term.push_back(0);
        par.push_back(p); pch.push_back(c); dep.push_back(p < 0 ? 0 : dep[p] + 1);
        return size() - 1;
    }
    int child(int v, char ch) const { return nxt[v][idx(ch)]; }
    int extend(int v, char ch) {  // 子がなければ作って返す
        int c = idx(ch);
        if (nxt[v][c] == -1) { int u = new_node(v, c); nxt[v][c] = u; }
        return nxt[v][c];
    }

    // --- 更新: 経路上の各頂点 v (根含む) で f(v, len) を呼ぶ (len = 根からの文字数) ---
    template <class F> int insert(string_view s, F f) {  // 戻り値: 終端頂点
        int v = 0, n = (int)s.size();
        for (int i = 0;; i++) {
            cnt[v]++;
            if (i == n) { term[v]++; f(v, i); break; }
            f(v, i);
            v = extend(v, s[i]);
        }
        return v;
    }
    int insert(string_view s) { return insert(s, [](int, int) {}); }
    template <class F> bool erase(string_view s, F f) {  // 1 本消す. 無ければ false
        int u = find(s);
        if (u == -1 || term[u] == 0) return false;
        int v = 0, n = (int)s.size();
        for (int i = 0;; i++) {
            cnt[v]--;
            if (i == n) { term[v]--; f(v, i); break; }
            f(v, i);
            v = nxt[v][idx(s[i])];
        }
        return true;
    }
    bool erase(string_view s) { return erase(s, [](int, int) {}); }

    // --- 検索 ---
    int find(string_view s) const {  // s の頂点 (無ければ -1). erase 済みでも頂点があれば返す
        int v = 0;
        for (char ch : s) { v = nxt[v][idx(ch)]; if (v == -1) return -1; }
        return v;
    }
    int count(string_view s) const { int v = find(s); return v < 0 ? 0 : term[v]; }
    int count_prefix(string_view s) const { int v = find(s); return v < 0 ? 0 : cnt[v]; }
    int lcp(string_view s) const {  // 集合内の文字列との最長共通接頭辞の長さ
        int v = 0, res = 0;
        for (char ch : s) { v = nxt[v][idx(ch)]; if (v == -1 || cnt[v] == 0) break; res++; }
        return res;
    }
    // s[start..] を根から辿り, 通った各頂点 v (根は含まない) で f(v, len) を呼ぶ. 戻り値: 辿れた文字数
    // 例 (単語分割 DP): t.walk(s, [&](int v, int len){ if (t.term[v]) dp[i + len] += dp[i]; }, i);
    template <class F> int walk(string_view s, F f, int start = 0) const {
        int v = 0;
        for (int i = start; i < (int)s.size(); i++) {
            v = nxt[v][idx(s[i])];
            if (v == -1) return i - start;
            f(v, i - start + 1);
        }
        return (int)s.size() - start;
    }

    // --- 順序 (重複込み. cnt/term に基づく) ---
    int count_less(string_view s) const {  // s より辞書順で小さい文字列の本数
        int v = 0, res = 0;
        for (char ch : s) {
            int c = idx(ch);
            res += term[v];
            for (int d = 0; d < c; d++) if (nxt[v][d] != -1) res += cnt[nxt[v][d]];
            v = nxt[v][c];
            if (v == -1) break;
        }
        return res;
    }
    int kth(int k) const {  // 辞書順 k 番目 (0-indexed) の文字列の終端頂点. 無ければ -1
        if (k < 0 || k >= cnt[0]) return -1;
        int v = 0;
        while (true) {
            if (k < term[v]) return v;
            k -= term[v];
            for (int c = 0; c < K; c++) {
                int u = nxt[v][c];
                if (u == -1) continue;
                if (k < cnt[u]) { v = u; break; }
                k -= cnt[u];
            }
        }
    }
    string str(int v) const {  // 頂点 v に対応する文字列
        string s(dep[v], ' ');
        for (int i = dep[v] - 1; i >= 0; i--) { s[i] = chr(pch[v]); v = par[v]; }
        return s;
    }

    // --- DFS 順 (子は文字順なので先行順 = 辞書順). v の部分木 = order[tin[v], tout[v]) ---
    // 頂点を増やしたら build_order() し直す
    vector<int> tin, tout, order;
    void build_order() {
        int n = size();
        tin.assign(n, 0); tout.assign(n, 0);
        order.clear(); order.reserve(n); order.push_back(0);
        vector<int> st = {0}, it(n, 0);
        while (!st.empty()) {
            int v = st.back();
            while (it[v] < K && nxt[v][it[v]] == -1) it[v]++;
            if (it[v] == K) { tout[v] = (int)order.size(); st.pop_back(); continue; }
            int u = nxt[v][it[v]++];
            tin[u] = (int)order.size(); order.push_back(u); st.push_back(u);
        }
    }
};
