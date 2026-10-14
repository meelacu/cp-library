// @name Monoid: 区間代入・区間min
// @prefix lazy_assign_min
// @description 区間代入・区間min

// 区間代入・区間min
const long long INF = 1LL << 60;
long long min_op(long long a, long long b) { return min(a, b); }
long long min_e() { return INF; }
const long long NONE = LLONG_MIN;
long long assign_min_mapping(long long f, long long x) { return f == NONE ? x : f; }
long long assign_composition(long long f, long long g) { return f == NONE ? g : f; }
long long assign_id() { return NONE; }
// LazySegTree<long long, min_op, min_e, long long, assign_min_mapping, assign_composition, assign_id> seg(v);
