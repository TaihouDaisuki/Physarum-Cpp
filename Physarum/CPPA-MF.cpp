#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<algorithm>
#include<queue>
#include<set>
#include<vector>

using namespace std;

const int Maxn = 110;
const int Maxm = 10010;
const int inf = 999999999 + 413;
const double eps = 1e-5;

int N, M, I0;
int S, T;
int C[Maxn][Maxn], L[Maxn][Maxn];
// var for Physarum Algorithm
double Q[Maxn][Maxn];
double D[Maxn][Maxn];
double p[Maxn]; // pressure
// var for solving p
double tI[Maxn];
double para[Maxn][Maxn], _para[Maxn][Maxn];
// var for ending
double preSumQ, sumQ;

// Gauss Algorithm solving Ab = C
void Gauss(const int n, double(*MatrixA)[Maxn], double* MatrixX, double* MatrixB)
{
	const double epx = 1e-6;

	for (int i = 1; i < n; i++)
	{
		int maxi = i;
		for (int j = i + 1; j <= n; j++)
			if (fabs(MatrixA[j][i]) > fabs(MatrixA[maxi][i]))
				maxi = j;
		if (maxi != i)
		{
			for (int j = i; j <= n; j++)
				swap(MatrixA[i][j], MatrixA[maxi][j]);
			swap(MatrixB[i], MatrixB[maxi]);
		}

		for (int j = i + 1; j <= n; j++)
		{
			if (fabs(MatrixA[j][i]) < 1e-6) continue;
			double tmp = MatrixA[i][i] / MatrixA[j][i];
			for (int k = i; k <= n; k++)
				MatrixA[j][k] = MatrixA[j][k] * tmp - MatrixA[i][k];
			MatrixB[j] = MatrixB[j] * tmp - MatrixB[i];
		}
	}

	for (int i = n; i; i--)
	{
		MatrixX[i] = MatrixB[i] / MatrixA[i][i];
		for (int j = 1; j < i; j++)
			MatrixB[j] -= MatrixX[i] * MatrixA[j][i];
	}
}

int CPPA_MF(const double k = 0.85) // k is the parameter of the capacity
{
	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
		{
			D[i][j] = (L[i][j] ? 0.5 : 0);
			Q[i][j] = 0;
		}
		p[i] = 0;
	}
	
	sumQ = 0;
	int count = 0;
	do
	{
		preSumQ = sumQ;
		/* Step one, caculate the value of p[i], using the D[i][j] from the last run */
		for (int i = 1; i <= N; ++i)
			for (int j = 1; j <= N; ++j)
				para[i][j] = 0;
		for (int i = 1; i <= N; ++i)
		{
			for (int j = 1; j <= N; ++j)
			{
				if (fabs(D[i][j]) <= eps)
					continue;
				para[i][i] += D[i][j] / L[i][j];
				para[i][j] = -D[i][j] / L[i][j];
			}

			if (i == S)
				tI[i] = I0;
			else if (i == T)
				tI[i] = -I0;
			else
				tI[i] = 0;
		}

		// set p[S] = 0, remove the para
		for (int i = 1; i <= N; ++i)
		{
			int _j = 0;
			for (int j = 1; j <= N; ++j)
			{
				if (j == T)
					continue;
				_para[i][++_j] = para[i][j];
			}
		}
		Gauss(N - 1, _para, p, tI);
		for (int i = N - 1; i >= T; --i)
			p[i + 1] = p[i];
		p[T] = 0;
		/*for (int i = 1; i <= N; ++i)
			printf("p[%d] = %.2lf, ", i, p[i]);
		puts("");*/

		/* Step two, caculate the value of Q[i][j], using the former p[i] */
		for (int i = 1; i <= N; ++i)
			for (int j = 1; j <= N; ++j)
				if (fabs(D[i][j]) > eps)
					Q[i][j] = (p[i] - p[j]) * D[i][j] / L[i][j];

		/*
		for (int i = 1; i <= N; ++i)
		{
			for (int j = 1; j <= N; ++j)
			{
				printf("Q[%d][%d] = %.4lf, ", i, j, Q[i][j]);
			}
			puts("");
		}
		puts("------------------------------------");
		*/

		/* Step three, caculate the value of D[i][j], for the next turn */
		for (int i = 1; i <= N; ++i)
			for (int j = 1; j <= N; ++j)
				// D[i][j] = (fabs(Q[i][j]) + D[i][j]) / 2; // without limit of edge ***
				if (fabs(D[i][j]) > eps)
					if (fabs(Q[i][j]) - k * C[i][j] <= eps)
						D[i][j] = (fabs(Q[i][j]) + D[i][j]) / 2;
					else
						D[i][j] = 1.0 * C[i][j] * L[i][j] / fabs(p[i] - p[j]);

		/*
		for (int i = 1; i <= N; ++i)
		{
			for (int j = 1; j <= N; ++j)
			{
				printf("D[%d][%d] = %.4lf, ", i, j, D[i][j]);
			}
			puts("");
		}
		puts("------------------------------------");
		*/

		sumQ = 0;
		for (int i = 1; i <= N; ++i)
			for (int j = 1; j <= N; ++j)
				if (Q[i][j] > -eps && L[i][j] == 1)
					sumQ += Q[i][j];

		++count;
	} while (fabs(preSumQ - sumQ) > eps * N);

	puts("===================================================");
	puts("=======================FINAL=======================");
	puts("===================================================");
	printf("Times: %d\n", count);

	for (int i = 1; i <= N; ++i)
	{
		for (int j = 1; j <= N; ++j)
		{
			printf("Q[%d][%d] = %.4lf, ", i, j, Q[i][j]);
		}
		puts("");
	}

	printf("MaxFlow = %.4lf\n", 1.0 * I0 - fabs(Q[S][T]));

	return 1; //******
}
void work()
{
	CPPA_MF();
	// int ans = CPPA();
	// printf("%d\n", ans);
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
	int a, b;
	N = get(); M = get(); 
	S = 1; T = N;
	int sumC = 0;
	for (int i = 1; i <= M; ++i)
	{
		a = get(); b = get();
		L[a][b] = L[b][a] = 1;
		sumC += (C[a][b] = C[b][a] = get());
	}

	// adding virtual path
	I0 = C[S][T] = C[T][S] = 100 * sumC;
	L[S][T] = L[T][S] = 100  * M;
}

int main()
{
	init();
	work();

	return 0;
}