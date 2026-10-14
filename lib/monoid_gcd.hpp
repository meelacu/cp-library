// @name Monoid: Range GCD
// @prefix monoid_gcd
// @description GCD monoid for SegTree

long long op_gcd(long long a, long long b) { return __gcd(a, b); }
long long e_gcd() { return 0; }
// SegTree<long long, op_gcd, e_gcd> seg(n);
$0
