// @name Monoid: 区間代入・区間和
// @prefix lazy_assign_sum
// @description 区間代入・区間和

// 区間代入・区間和  初期化は {a[i], 1}
struct SumLen { long long sum; int len; };
SumLen sl_op(SumLen a, SumLen b) { return {a.sum + b.sum, a.len + b.len}; }
SumLen sl_e() { return {0, 0}; }
const long long NONE = LLONG_MIN;
SumLen assign_sum_mapping(long long f, SumLen x) { return f == NONE ? x : SumLen{f * x.len, x.len}; }
long long assign_composition(long long f, long long g) { return f == NONE ? g : f; }
long long assign_id() { return NONE; }
// LazySegTree<SumLen, sl_op, sl_e, long long, assign_sum_mapping, assign_composition, assign_id> seg(v);
