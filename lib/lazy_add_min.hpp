// @name Monoid: 区間加算・区間min
// @prefix lazy_add_min
// @description 区間加算・区間min

// 区間加算・区間min  (max にするなら min→max, INF→-INF)
long long min_op(long long a, long long b) { return min(a, b); }
long long min_e() { return LINF; }
long long add_min_mapping(long long f, long long x) { return x == LINF ? LINF : x + f; }
long long add_composition(long long f, long long g) { return f + g; }
long long add_id() { return 0; }
// LazySegTree<long long, min_op, min_e, long long, add_min_mapping, add_composition, add_id> seg(v);
