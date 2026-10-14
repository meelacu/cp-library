// @name Monoid: 区間加算・区間和/二乗和
// @prefix lazy_add_sqsum
// @description 区間加算・区間和/二乗和

// 区間加算・区間和/二乗和 mod  初期化は {a[i], a[i]*a[i]%MOD, 1} mod書き換え忘れるなよアホ
const long long MOD = 998244353;
struct SqSum { long long sum, sqsum; int len; };  // 和, 二乗和, 区間長
SqSum sqsum_op(SqSum a, SqSum b) { return {(a.sum+b.sum)%MOD, (a.sqsum+b.sqsum)%MOD, a.len+b.len}; }
SqSum sqsum_e() { return {0, 0, 0}; }
SqSum add_sqsum_mapping(long long f, SqSum x) {
    f %= MOD;
    long long nsqsum = (x.sqsum + 2*f%MOD*x.sum + f*f%MOD*x.len) % MOD;
    long long nsum = (x.sum + f*x.len) % MOD;
    return {nsum, nsqsum, x.len};
}
long long add_composition(long long f, long long g) { return (f + g) % MOD; }
long long add_id() { return 0; }
// LazySegTree<SqSum, sqsum_op, sqsum_e, long long, add_sqsum_mapping, add_composition, add_id> seg(v);
