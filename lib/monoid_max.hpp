// @name Monoid: Range Max
// @prefix monoid_max
// @description Max monoid for SegTree

long long op_max(long long a, long long b) { return max(a, b); }
long long e_max() { return LLONG_MIN; }
// SegTree<long long, op_max, e_max> seg(n);
$0
