// @name Monoid: Range Min
// @prefix monoid_min
// @description Min monoid for SegTree

long long op_min(long long a, long long b) { return min(a, b); }
long long e_min() { return LLONG_MAX; }
// SegTree<long long, op_min, e_min> seg(n);
$0
