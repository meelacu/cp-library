// @name power & inv (mod)
// @prefix power
// @description Modular exponentiation and modular inverse (Fermat)

ll power(ll a, ll b){
	a %= mod; ll r = 1;
	while(b){ if(b & 1) r = r * a % mod; a = a * a % mod; b >>= 1; }
	return r;
}
ll inv(ll a){ return power(a, mod - 2); }
