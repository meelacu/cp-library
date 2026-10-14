// @name Monoid: Range Sum
// @prefix monoid_sum
// @description Sum monoid for SegTree

long long op_sum(long long a, long long b) { return a + b; }
long long e_sum() { return 0; }
// SegTree<long long, op_sum, e_sum> seg(n);
$0
