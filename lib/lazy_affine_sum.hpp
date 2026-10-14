// @name Monoid: 区間アフィン・区間和(mod)
// @prefix lazy_affine_sum
// @description 区間アフィン・区間和(mod)

// 区間アフィン (x -> a*x+b)・区間和 mod  初期化は {a[i], 1}
const long long MOD = 998244353;
struct SumLen { long long sum; int len; };
SumLen sl_op(SumLen a, SumLen b) { return {(a.sum + b.sum) % MOD, a.len + b.len}; }
SumLen sl_e() { return {0, 0}; }
struct Affine { long long a, b; };
SumLen affine_mapping(Affine f, SumLen x) { return {(f.a * x.sum % MOD + f.b * x.len % MOD) % MOD, x.len}; }
Affine affine_composition(Affine f, Affine g) { return {f.a * g.a % MOD, (f.a * g.b % MOD + f.b) % MOD}; }
Affine affine_id() { return {1, 0}; }
// LazySegTree<SumLen, sl_op, sl_e, Affine, affine_mapping, affine_composition, affine_id> seg(v);
// seg.apply(l, r, {a, b});
