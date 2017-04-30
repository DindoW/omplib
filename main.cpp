#include <time.h>
#include <iostream>
#include <functional>
#include <iomanip>
#include "omplib.hpp"
#include "slib.hpp"

using namespace std;

#define PSORTDEBUG
#define PMATRIXTRANS
#define PMATRIXMULTY
#define PGAUSSELIM
#define PFFT


int main()
{
	clock_t pt = 0, st = 0;
	int p = 4;

	/**************psort测试**************/
#ifdef PSORTDEBUG
	{
		cout << "******************   psort   *****************" << endl;
		int n = 10000;
		int *a = new int[n];
		int *b = new int[n];
		int *c = new int[n];

		{	//*****初始化
			for (int i = 0; i < n; i++)
				a[i] = i;
			random_shuffle(a, a + n);
			for (int i = 0; i < n; i++)
				b[i] = c[i] = a[i];
		}

#ifdef DEBUG
		cout << "A: ";
		for (int i = 0; i < n; i++)
			cout << a[i] << "  ";
		cout << endl;
#endif // DEBUG
		
		{	//*****串行排序时间
			st -= clock();
			//sort(b, b + n, greater<double>());
			sort(b, b + n);
			st += clock();
		}

		{	//******并行排序时间
			pt -= clock();
			par::psort(c, n, p);
			pt += clock();
		}

		{	//*******验证排序结果
#ifdef DEBUG
			cout << endl << "B: ";
			for (int i = 0; i < n; i++)
				cout << b[i] << "  ";
			cout << endl << "C: ";
			for (int i = 0; i < n; i++)
				cout << c[i] << "  ";
#endif // DEBUG

			for (int i = 0; i < n; i++)
				if (b[i] != c[i])
				{
					cout << "Error : " << i << endl;
					break;
				}
		}
		delete[]a;
		delete[]b;
		delete[]c;
		cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
		cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
		cout << "speedup ratio:" << (double)st / pt << endl;
	}
#endif // PSORTDEBUG
	
	/**********pMatrixTrans测试***********/
#ifdef PMATRIXTRANS
	for (int test = 0; test < 5; test++)
	{
		std::cout << "****************  pMatrixTrans  **************" << endl;
		
		int m = 400, n = 400;

		int **a = new int*[m];
		for (int i = 0; i < m; i++)
			a[i] = new int[n];
		int **b, **c;

		{	//******初始化
			srand((unsigned int)time(NULL));
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					a[i][j] = rand();
		}

		{	//*****串行转置时间
			st = 0;
			st -= clock();
			b = slib::sMartrixTrans(a, m, n);
			st += clock();
		}

		{	//*****并行转置时间
			pt = 0;
			pt -= clock();
			c = par::pMatrixTrans(a, m, n, p);
			pt += clock();
		}

		{	//*****验证计算结果

#ifdef DEBUG
			cout.setf(ios::fixed);
			cout.precision(3);
			cout.setf(ios::showpos);
			cout << "A::::" << endl;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
					cout << setw(12) << a[i][j];
				cout << endl;
			}
			cout << "B::::" << endl;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					cout << setw(12) << b[i][j];
				cout << endl;
			}
			cout << "C::::" << endl;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
					cout << setw(12) << c[i][j];
				cout << endl;
			}
#endif // DEBUG

			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
				{
					if (fabs(b[j][i] - a[i][j]) > slib::EPSINON)
					{
						std::cout << "Serial Error!" << endl;
						i = m;
						break;
					}
					if (fabs(b[j][i] - c[j][i]) > par::EPSINON)
					{
						std::cout << "Parrallel Error!" << endl;
						i = m;
						break;
					}
				}
		}

		for (int i = 0; i < m; i++)
			delete[] a[i];		
		for (int i = 0; i < n; i++)
		{
			delete[] b[i];
			delete[] c[i];
		}
		delete[] a;
		delete[] b;
		delete[] c;
		std::cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
		std::cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
		std::cout << "speedup ratio:" << (double)st / pt << endl;
	}
#endif // PMATRIXTRANS

	/**********pMatrixMulty测试***********/
#ifdef PMATRIXMULTY
	{
		cout << "****************  pMatrixMult  ***************" << endl;
		int m = 2, s = 2, n = 2;

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

		{	//******初始化
			srand((unsigned int)time(NULL));
			for (int i = 0; i < m; i++)
				for (int j = 0; j < s; j++)
					a[i][j] = rand();
			for (int i = 0; i < s; i++)
				for (int j = 0; j < n; j++)
					b[i][j] = rand();
		}

		{	//*****串行排序时间
			st = 0;
			st -= clock();
			slib::sMartrixMult(a, b, c, m, s, n);
			st += clock();
		}

		{	//*****并行排序时间
			pt = 0;
			pt -= clock();
			par::pMatrixMult(a, b, d, m, s, n, p);
			pt += clock();
		}

		{	//*****验证计算结果

#ifdef DEBUG
			cout.setf(ios::fixed);
			cout.precision(3);
			cout.setf(ios::showpos);

			cout << "A::::" << endl;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < s; j++)
					cout << setw(12) << a[i][j];
				cout << endl;
			}
			cout << "B::::" << endl;
			for (int i = 0; i < s; i++)
			{
				for (int j = 0; j < n; j++)
					cout << setw(12) << b[i][j];
				cout << endl;
			}
			cout << "C::::" << endl;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
					cout << setw(12) << c[i][j];
				cout << endl;
			}
			cout << "D::::" << endl;
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < n; j++)
					cout << setw(12) << d[i][j];
				cout << endl;
			}
#endif // DEBUG

			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					if (fabs(c[i][j] - d[i][j]) > slib::EPSINON)
					{
						cout << "Error!" << endl;
						i = m;
						break;
					}
		}

		for (int i = 0; i < m; i++)
		{
			delete[] a[i];
			delete[] c[i];
			delete[] d[i];
		}
		for (int i = 0; i < s; i++)
			delete[] b[i];
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] d;
		cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
		cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
		cout << "speedup ratio:" << (double)st / pt << endl;
	}
#endif // PMATRIXMULTY
	
	/***********pGaussElim测试************/
#ifdef PGAUSSELIM
	{
		cout << "****************  pGaussElim  ****************" << endl;
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

		{	//*******初始化
			srand((unsigned)time(NULL));
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

		{	//*****串行排序时间
			st -= clock();
			xx = slib::sGaussElim(aa, bb, n);
			if (xx == nullptr)
				cout << "Matrix rank error!" << endl;
			st += clock();
		}

		{	//*****验证串行计算结果
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
				if (abs(res - b[i]) > slib::EPSINON)
				{
					cout << " Serial Error!";
				}
#ifdef DEBUG
				cout << setiosflags(ios::fixed) << setprecision(3) << setw(8) << res;
#endif // DEBUG
			}
		}

		{//*****并行排序时间
			pt -= clock();
			x = par::pGaussElim(a, b, n);
			if (x == nullptr)
				cout << "Matrix rank error!" << endl;
			pt += clock();
		}

		{//*****验证计算结果
			for (int i = 0; i < n; i++)
				if (fabs(x[i] - xx[i]) > par::EPSINON)
				{
					cout << " Parallel Error!" << endl;
					break;
				}
		}

		for (int i = 0; i < n; i++)
		{
			delete[] a[i];
			delete[] aa[i];
		}
		delete[] a;
		delete[] aa;
		delete[] b;
		delete[] bb;
		delete[] x;
		delete[] xx;
		cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
		cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
		cout << "speedup ratio:" << (double)st / pt << endl;
	}
#endif // PGAUSSELIM
	
	/**************pfft测试***************/
#ifdef PFFT
	{
		cout << "****************  pGaussElim  ****************" << endl;
		int n = 2 << 10;
		double *x = new double[n];
		par::CM *sres, *pres;

		//*******初始化
		{
			for (int i = 0; i < n; i++)
				x[i] = (double)rand() / 10 + rand();
			//x[i] = i;
		}

		//*****串行排序时间
		{
			st -= clock();
			sres = slib::sfft(x, n);
			st += clock();
		}

		//******并行排序时间
		{
			pt -= clock();
			pres = par::pfft(x, n);
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
				if (fabs(sres[i].real() - pres[i].real()) > par::EPSINON || fabs(sres[i].imag() - pres[i].imag()) > par::EPSINON)
				{
					cout << "ERROR!" << i << endl;
					break;
				}
		}
		delete[]x;
		delete[]sres;
		delete[]pres;
		cout << "serial time:" << (double)st / CLOCKS_PER_SEC << "s." << endl;
		cout << "parallel time:" << (double)pt / CLOCKS_PER_SEC << "s." << endl;
		cout << "speedup ratio:" << (double)st / pt << endl;
	}
#endif // PFFT

	std::system("pause");
	return 0;
}