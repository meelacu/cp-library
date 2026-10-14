// @name Monoid: 区間加算・区間和
// @prefix lazy_add_sum
// @description 区間加算・区間和

// 区間加算・区間和  初期化は {a[i], 1}
struct SumLen { long long sum; int len; };
SumLen sl_op(SumLen a, SumLen b) { return {a.sum + b.sum, a.len + b.len}; }
SumLen sl_e() { return {0, 0}; }
SumLen add_sum_mapping(long long f, SumLen x) { return {x.sum + f * x.len, x.len}; }
long long add_composition(long long f, long long g) { return f + g; }
long long add_id() { return 0; }
// LazySegTree<SumLen, sl_op, sl_e, long long, add_sum_mapping, add_composition, add_id> seg(v);
