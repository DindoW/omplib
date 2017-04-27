#include "omplib.h"
#include <time.h>
#include <iostream>
#include <functional>
#include <iomanip>

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
				if (fabs(c[i][j] - d[i][j]) <= EPSINON)
				{
					cout << "Error!" << endl;
					i = m;
					break;
				}
	}
*/

/*
//线性方程组求解
	int n = 100;
	double **a = new double*[n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n];
	double *b = new double[n];
	double **aa = new double*[n];
	for (int i = 0; i < n; i++)
		aa[i] = new double[n];
	double *bb = new double[n];
	double *x = new double[n];
	double *xx = new double[n];

	{//*******初始化
		srand((unsigned)time(NULL));
		//srand(4);
		for (int i = 0; i < n; i++)
		{
			bb[i] = b[i] = rand() % n;
			for (int j = 0; j < n; j++)
				aa[i][j] = a[i][j] = (double)(rand() % n) / 10 + rand() % n;
		}
#ifdef DEBUG
		cout << "A & B:" << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				cout << setprecision(3) << setw(8) << a[i][j];
			cout << setw(8) << b[i] << endl;
		}
#endif // DEBUG
	}

	{//*****串行排序时间
		st -= clock();
		if (sGaussElim(aa, bb, xx, n) == -1)
			cout << "Matrix rank error!" << endl;
		st += clock();
	}

	{//*****验证串行计算结果
#ifdef DEBUG
		cout << "Simulation A * XX = B :" << endl;
#endif // DEBUG
		for (int i = 0; i < n; i++)
		{
			double res = 0;
			for (int j = 0; j < n; j++)
			{
				res += xx[j] * a[i][j];
			}
			if (abs(res - b[i]) > EPSINON)
			{
				cout << " Serial Error!";
			}
#ifdef DEBUG
			cout << setiosflags(ios::fixed) << setprecision(3) << setw(8) << res;
#endif // DEBUG
		}
		cout << endl;
	}

	{//*****并行排序时间
		pt -= clock();
		if (pGaussElim(a, b, x, n) == -1)
			cout << "Matrix rank error!" << endl;
		pt += clock();
	}

	{//*****验证计算结果
		for (int i = 0; i < n; i++)
			if (fabs(x[i] - xx[i]) > EPSINON)
			{
				cout << " Parallel Error!" << endl;
				break;
			}
		cout << endl;
	}
*/


//快速傅里叶变换
	int n = 2 << 10;
	double *x = new double[n];
	CM *sres, *pres;

#ifdef DEBUG
	cout.setf(ios::fixed);
	cout.precision(3);
	cout.setf(ios::showpos);
#endif // DEBUG

	//*******初始化
	{
		for (int i = 0; i < n; i++)
			x[i] = (double)rand() / 10 + rand();
			//x[i] = i;
	}

	//*****串行排序时间
	{
		st -= clock();
		sres = sfft(x, n);
		st += clock();
	}

	//******并行排序时间
	{
		pt -= clock();
		pres = pfft(x, n);
		pt += clock();
	}

	//*******验证排序结果
	{
#ifdef DEBUG
		cout << "The s_result are as follows：" << endl;
		for (int i = 0; i < n; i++)
			cout << sres[i] << endl;
		cout << "The p_result are as follows：" << endl;
		for (int i = 0; i < n; i++)
			cout << pres[i] << endl;
#endif // DEBUG
		for (int i = 0; i < n; i++)
			if (fabs(sres[i].real() - pres[i].real()) > EPSINON || fabs(sres[i].imag() - pres[i].imag()) > EPSINON)
			{
				cout << "ERROR!" << i << endl;
				break;
			}
	}

	cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
	cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
	cout << "speedup ratio:" << (double)st / pt << endl;

	system("pause");
	return 0;
}