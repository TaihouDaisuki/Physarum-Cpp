#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<algorithm>
#include<queue>
#include<set>
#include<vector>

#define For(i,l,r) for(int i=l;i<=r;++i)
#define getLink(k,x) for(int k=fst[x];k!=-1;k=pre[k])
#define add(x,y,c) to[++tot]=y,pre[tot]=fst[x],fst[x]=tot,rl[tot]=c
#define Min(a,b) ((a)<(b)?(a):(b))

using namespace std;

const int Maxn = 20010;
const int Maxm = 450010;
const int inf = 999999999 + 413;

int N, M;
int S, T, tot;
int fst[Maxn], pre[Maxm], to[Maxm], rl[Maxm];
int vis[Maxn], Q[Maxn], l, r;

inline int dfs(int u, int flow)
{
	if (u == T) return flow; int sf = 0;
	getLink(k, u)
	{
		int v = to[k]; if (!rl[k] || vis[v] != vis[u] + 1) continue;
		int tmp = dfs(v, Min(rl[k], flow));
		rl[k] -= tmp; rl[k ^ 1] += tmp; flow -= tmp; sf += tmp;
		if (!flow) break;
	}
	if (!sf) vis[u] = 0;
	return sf;
}
inline bool bfs()
{
	For(i, S, T) vis[i] = 0;
	l = 1; r = 2; Q[l] = S; vis[S] = 1;
	while (l < r)
	{
		int u = Q[l++];
		getLink(k, u)
		{
			int v = to[k]; if (vis[v] || !rl[k]) continue;
			vis[v] = vis[u] + 1; Q[r++] = v;
		}
	}
	return vis[T];
}
void work()
{
	int ans = 0;
	while (bfs()) ans += dfs(S, inf);
	printf("%d\n", ans);
}

inline int get()
{
	int v; char ch;
	while (!isdigit(ch = getchar())); v = ch - 48;
	while (isdigit(ch = getchar())) v = v * 10 + ch - 48;
	return v;
}
void init()
{
	int a, b, w;
	N = get(); M = get(); S = 1; T = N;
	tot = -1; For(i, S, T) fst[i] = -1;
	For(i, 1, M)
	{
		a = get(); b = get(); w = get();
		add(a, b, w); add(b, a, w);
	}
}

int main()
{
	init();
	work();

	system("pause");
	return 0;
}