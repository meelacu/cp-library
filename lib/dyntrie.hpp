// @name Trie (dynamic, integer labels)
// @prefix dyntrie
// @description Trie木(map版/ラベルが任意の整数)。子をラベル昇順に走査できる

// ===== DynTrie (ラベルが任意の整数 / アルファベットが巨大なとき) =====
struct DynTrie {
    vector<map<int,int>> ch;   // ch[v][c] = 子の頂点番号 (ラベル昇順に走査できる)
    vector<int> cnt, term;     // cnt: 通る本数, term: そこで終わる本数
    DynTrie() : ch(1), cnt(1,0), term(1,0) {}
    int size() const { return (int)ch.size(); }

    // 頂点 v からラベル c の子へ進む (無ければ作る)
    int extend(int v, int c) {
        auto it = ch[v].find(c);
        if (it != ch[v].end()) return it->second;
        int u = size();
        ch.emplace_back(); cnt.push_back(0); term.push_back(0);
        ch[v][c] = u;          // 再確保が起きるので emplace_back の後に代入
        return u;
    }

    int insert(const vector<int>& s) {
        int v = 0; cnt[v]++;
        for (int c : s) { v = extend(v, c); cnt[v]++; }
        term[v]++;
        return v;              // 終端の頂点番号
    }

    // s に対応する頂点番号 (存在しなければ -1)
    int node(const vector<int>& s) const {
        int v = 0;
        for (int c : s) {
            auto it = ch[v].find(c);
            if (it == ch[v].end()) return -1;
            v = it->second;
        }
        return v;
    }

    int count(const vector<int>& s) const { int v = node(s); return v < 0 ? 0 : term[v]; }
    int count_prefix(const vector<int>& s) const { int v = node(s); return v < 0 ? 0 : cnt[v]; }

    int lcp(const vector<int>& s) const {
        int v = 0, res = 0;
        for (int c : s) {
            auto it = ch[v].find(c);
            if (it == ch[v].end() || cnt[it->second] == 0) break;
            v = it->second; res++;
        }
        return res;
    }

    // s[start..] を辿り、登録済み列と一致するたび f(頂点番号, 一致長) を呼ぶ
    template <class F>
    void query(const vector<int>& s, F f, int start = 0) const {
        int v = 0;
        for (int i = start; i < (int)s.size(); i++) {
            auto it = ch[v].find(s[i]);
            if (it == ch[v].end()) return;
            v = it->second;
            if (term[v]) f(v, i - start + 1);
        }
    }
};
$0
