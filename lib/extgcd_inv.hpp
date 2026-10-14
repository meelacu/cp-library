// @name extgcd inverse
// @prefix extgcd_inv
// @description 拡張ユークリッド互除法による mod 逆元

long long ext_gcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    long long x1, y1;
    long long g = ext_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

long long inv(long long n, long long p) {
    long long x, y;
    long long g = ext_gcd(n, p, x, y);
    assert(g == 1);
    return (x % p + p) % p;
}
