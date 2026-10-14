// @name Z-algorithm
// @prefix zalgo
// @description Z-algorithm: 各位置の前方一致長を O(n) で求める

// Z[i] = S[i..] と S[0..] の最長共通接頭辞長。全体 O(n)
vi z_algorithm(const string &s){
    int n = s.size();
    vi z(n);
    z[0] = n;
    int l = 0, r = 0;
    for(int i = 1; i < n; i++){
        if(i < r) z[i] = min(r - i, z[i - l]);
        while(i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
        if(i + z[i] > r){ l = i; r = i + z[i]; }
    }
    return z;
}
