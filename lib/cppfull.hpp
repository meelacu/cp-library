// @name Competitive Programming Full Template
// @prefix cppfull
// @description C++ Full template for Competitive Programming

#include <bits/stdc++.h>
#include <cassert>
using namespace std;

using ll = long long;
using ull = unsigned long long;
using cent =__int128_t;
using ucent =  __uint128_t;
using ld = long double;
using pii = pair<int, int>;
using pll = pair<ll,ll>;
using vb = vector<bool>;
using vvb = vector<vb>;
using vi = vector<int>;
using vll = vector<long long>;
using vvi = vector<vi>;
using vvll = vector<vll>;
using vpii = vector<pii>;
using vpll = vector<pll>;
using vs = vector<string>;
using vc = vector<char>;
using vvc = vector<vc>;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(),(a).rend()
#define uniq(v) sort(all(v)), (v).erase(unique(all(v)),(v).end())
#define lb(v,x) int(lower_bound(all(v),(x))-(v).begin())
#define ub(v,x) int(upper_bound(all(v),(x))-(v).begin())
#define rrep(i,n) for(int i = (int)(n) - 1; i>=0; i--)
#define rep1(i,n) for(int i=1; i <=(int)(n);i++)
#define reps(i,s,n) for(long long i=(long long)(s);i<(long long)(n); i++)
#define koz " "
#define tor "\n"
#define elso first
#define maso second
template<class T> using min_pq = priority_queue<T, vector<T>, greater<T>>;
template<class T, class U> bool chmin(T& a, const U& b){if(b<a){a=b;return true;}return false;}
template<class T, class U> bool chmax(T& a, const U& b){if(a<b){a=b;return true;}return false;}
template<typename T, typename U> T ceil(T x, U y){
    assert(y!=0);
    if constexpr (is_signed_v<T> && is_unsigned_v<U>) return ceil(x, (long long)y);
    return x / y + ((x % y != 0) && ((x < 0) == (y < 0)));
}
template<typename T, typename U> T floor(T x, U y){
    assert(y!=0);
    if constexpr (is_signed_v<T> && is_unsigned_v<U>) return floor(x, (long long)y);
    return x / y - ((x % y != 0) && ((x < 0) != (y < 0)));
}
template<class T> int popcnt(T x){return __builtin_popcountll((make_unsigned_t<T>)x);}

template<class T> int topbit(T x){return x == 0 ? -1 : 63 - __builtin_clzll((make_unsigned_t<T>)x);}

inline long long nxtcomb(long long s){
    if(s==0) return 1LL<<62;
    long long x=s&-s,y=s+x; return(((s&~y)/x)>>1)|y;
}

template<class T> vector<T> racs(size_t n, T v){return vector<T>(n, v);}
template<class... Ts> auto racs(size_t n, Ts... ts){return vector<decltype(racs(ts...))>(n, racs(ts...));}

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
ll random(ll l, ll r){return uniform_int_distribution<ll>(l, r)(rng);}  // [l, r]

ostream& operator<<(ostream& os, ucent u){
    string s; do{ s += char('0' + u % 10); u /= 10; }while(u);
    reverse(all(s)); return os << s;
}
ostream& operator<<(ostream& os, cent x){
    if(x < 0) os << '-';
    return os << (x < 0 ? -(ucent)x : (ucent)x);
}
istream& operator>>(istream& is, cent& x){ ll t; is >> t; x=t; return is;}
istream& operator>>(istream& is, ucent& x){ ull t; is >> t; x=t; return is;}

const int dx[]={1,0,-1,0,1,1,-1,-1};
const int dy[]={0,1,0,-1,1,-1,1,-1};

inline bool inside(ll i, ll j, ll H, ll W){return 0 <= i && i < H && 0 <= j && j < W;}

const int INF = 1045141919;
const ll LINF = 3643648101145141919;

void solve(){
    $0
}

int main() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(16);

    int _ = 1;
    //$1cin >> _;
    while(_--) solve();
    
}
