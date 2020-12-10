#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<string>
#include<algorithm>
#include<queue>
#include<set>

using namespace std;

const int Maxn = 15;

void Gauss(const int n, double(*MatrixA)[Maxn], double* MatrixC, double* MatrixB)
{
	const double epx = 1e-6;

	for (int i = 1; i < n; i++)
	{
		int maxi = i;
		for (int j = i + 1; j <= n; j++)
			if (fabs(MatrixA[j][i]) > fabs(MatrixA[maxi][i])) maxi = j;
		if (maxi != i)
		{
			for (int j = i; j <= n; j++) swap(MatrixA[i][j], MatrixA[maxi][j]);
			swap(MatrixC[i], MatrixC[maxi]);
		}

		for (int j = i + 1; j <= n; j++)
		{
			if (fabs(MatrixA[j][i]) < 1e-6) continue;
			double tmp = MatrixA[i][i] / MatrixA[j][i];
			for (int k = i; k <= n; k++) MatrixA[j][k] = MatrixA[j][k] * tmp - MatrixA[i][k];
			MatrixC[j] = MatrixC[j] * tmp - MatrixC[i];
		}
	}

	for (int i = n; i; i--)
	{
		MatrixB[i] = MatrixC[i] / MatrixA[i][i];
		for (int j = 1; j < i; j++) MatrixC[j] -= MatrixB[i] * MatrixA[j][i];
	}
}

int main()
{
	int n;
	double a[Maxn][Maxn], c[Maxn];
	double ans[Maxn];

	scanf(" %d", &n);
	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= n; ++j)
			scanf(" %lf", &a[i][j]);
		scanf(" %lf", &c[i]);
	}

	Gauss(n, a, c, ans);

	for (int i = 1; i < n; i++) printf("%.3lf ", ans[i]);
	printf("%.3lf\n", ans[n]);

	system("pause");
	return 0;
}