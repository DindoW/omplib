#include "omplib.h"
#include <time.h>
#include <iostream>
#include <functional>

using namespace std;

int main()
{
	clock_t pt = 0, st = 0;
	int p = 4;

/*
//并行排序
	int n = 1000;
	int *a = new int[n];
	int *b = new int[n];
	int *c = new int[n];

	srand((unsigned)time(NULL));
	//	for (int i = 0; i < n; i++)
	//		a[i] = (double)rand() / (RAND_MAX / 10);
	for (int i = 0; i < n; i++)
		a[i] = i;
	random_shuffle(a, a + n);
	for (int i = 0; i < n; i++)
		b[i] = c[i] = a[i];

	cout << "A: ";
	//	for (int i = 0; i < n; i++)
	//		cout << a[i] << "  ";
	cout << endl;

	//*****串行排序时间
	{
		st -= clock();
		//		sort(b, b + n, greater<double>());
		sort(b, b + n);
		st += clock();
	}

	//******并行排序时间
	{
		pt -= clock();
		psort(c, n, p);
		pt += clock();
	}

	//*******验证排序结果
	cout << endl << "B: ";
	//	for (int i = 0; i < n; i++)
	//		cout << b[i] << "  ";
	cout << endl << "C: ";
	for (int i = 0; i < n; i++)
	{
		if (b[i] != (double)c[i])
		{
			cout << "e:";
			break;
		}
		//		cout << (double)c[i] << "  ";
	}
	cout << endl;

	delete[]a;
	delete[]b;
	delete[]c;
*/

/*
//并行矩阵乘法
	int m = 100, s = 200, n = 110;

	double **a = new double*[m];
	double **b = new double*[s];
	double **c = new double*[m];
	double **d = new double*[m];
	for (int i = 0; i < m; i++)
	{
		a[i] = new double[s];
		c[i] = new double[n];
		d[i] = new double[n];
	}
	for (int i = 0; i < s; i++)
		b[i] = new double[n];

	{//******初始化
		for (int i = 0; i < m; i++)
			for (int j = 0; j < s; j++)
				a[i][j] = 1;
		for (int i = 0; i < s; i++)
			for (int j = 0; j < n; j++)
				b[i][j] = i + j / 3;
	}

	{//*****串行排序时间
		st -= clock();
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
			{
				double tmp = 0;
				for (int k = 0; k < s; k++)
					tmp += a[i][k] * b[k][j];
				c[i][j] = tmp;
			}
		st += clock();
	}

	{//*****并行排序时间
		pt -= clock();
		pMatrixMult(a, b, d, m, s, n, p);
		pt += clock();
	}

	{//*****验证计算结果
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if (c[i][j] != d[i][j])
				{
					cout << "Error!" << endl;
					i = m;
					break;
				}
	}
*/



	cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
	cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
	cout << "speedup ratio:" << (double)st / pt << endl;

	system("pause");
	return 0;
}