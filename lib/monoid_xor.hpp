// @name Monoid: Range XOR
// @prefix monoid_xor
// @description XOR monoid for SegTree

long long op_xor(long long a, long long b) { return a ^ b; }
long long e_xor() { return 0; }
// SegTree<long long, op_xor, e_xor> seg(n);
$0
