// @name modint combination
// @prefix comb
// @description modint nCk with factorial precomputation

const long long MOD = ;
const int MAX_C = 2000005;
long long fact[MAX_C], inv_fact[MAX_C];

long long power(long long a, long long b, long long m) {
    long long res = 1;
    a %= m;
    while (b > 0) {
        if (b & 1) res = res * a % m;
        a = a * a % m;
        b >>= 1;
    }
    return res;
}

void init_comb() {
    fact[0] = 1;
    for (int i = 1; i < MAX_C; i++) fact[i] = fact[i-1] * i % MOD;
    inv_fact[MAX_C-1] = power(fact[MAX_C-1], MOD-2, MOD);
    for (int i = MAX_C-2; i >= 0; i--) {
        inv_fact[i] = inv_fact[i+1] * (i+1) % MOD;
    }
}

long long nCk(ll n, int k) {
    if (k < 0 || k > n) return 0;
    if(n<MAX_C) return fact[n] * inv_fact[k] % MOD * inv_fact[n-k] % MOD;
    else{
        ll re=1;
        for(ll i=0;i<k;i++) re=re*((n-i)%MOD)%MOD;
        return re*inv_fact[k]%MOD;
        }
}
$0
