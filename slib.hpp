#pragma once
#include <complex>
#include <iostream>
#include "omplib.hpp"

namespace slib
{
	const double EPSINON = 0.0001;
	const double PI = acos(-1.0);
	typedef std::complex<double> CM;

	template <typename T>
	T** sMartrixTrans(T **a, int m, int n)
	{
		T **transMx = new T*[n];
		for (int i = 0; i < n; i++)
			transMx[i] = new T[m];

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				transMx[j][i] = a[i][j];
			}
		}
		return transMx;
	}

	template <typename T1, typename T2>
	void sMartrixMult(T1 **a, T1 **b, T2**c, int m, int s, int n)
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				T2 tmp = 0;
				for (int k = 0; k < s; k++)
					tmp += a[i][k] * b[k][j];
				c[i][j] = tmp;
			}
		}
	}

	double* sGaussElim(double **a, double *b, int n)
	{
		int *marked = new int[n];
		int *pivot = new int[n];
		int picked, pos;
		double tmp;

		for (int i = 0; i < n; i++)
			marked[i] = 0;

		for (int i = 0; i < n; i++)
		{
			tmp = 0;
			for (int j = 0; j < n; j++)	//选绝对值最大者作主元
			{
				if (marked[j] == 0 && fabs(a[j][i]) > tmp)
				{
					tmp = fabs(a[j][i]);
					picked = j;
				}
			}
			if (fabs(tmp - 0) < EPSINON)
			{
				std::cout << "Rank of Matrix Error !" << std::endl;
				return nullptr;
			}
			marked[picked] = 1;
			pivot[picked] = i;

			for (int j = 0; j < n; j++)	//消元
			{
				if (picked != j)
				{
					tmp = a[j][i] / a[picked][i];
					int k;
					for (k = i; k < n; k++)
						a[j][k] -= a[picked][k] * tmp;
					b[j] -= b[picked] * tmp;
				}
			}
		}
		double *x = new double[n];
		for (int i = 0; i < n; i++)	//求解x
		{
			pos = pivot[i];
			x[pos] = b[i] / a[i][pos];
		}

#ifdef DEBUG
		std::cout << "after Gauss-Jordan Elimination A & B & pivot :" << std::endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(8) << a[i][j];
			std::cout << std::setprecision(3) << std::setw(8) << b[i] << std::setw(6) << pivot[i] << std::endl;
		}
		std::cout << "XX : " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << std::setw(8) << x[i];
		std::cout << std::endl;
#endif // DEBUG

		return x;
	}

	CM* sfft(const double *x, const int n)			//cormen迭代fft
	{
		if ((n & n - 1) != 0)	//n必须为2的幂
			return nullptr;

		CM *w = new CM[n];

		double theta = 2 * PI / n;			//用欧拉公式计算旋转因子 
		for (int i = 0; i < n; i++)
			w[i] = CM(cos(theta * i), sin(theta * i));

#ifdef DEBUG
		std::cout << "The sssssssssssssssssssssfft begin： " << std::endl;
		std::cout << "The w are as follows： " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << w[i] << std::endl;
#endif // DEBUG	

		CM *newx = new CM[n];				//利用按位与以及循环实现码位颠倒
		for (int i = 0; i < n; i++)
		{
			int k = i;
			int j = 0;
			for (int t = (int)(log(n) / log(2)); t > 0; t--)
			{
				j = j << 1;
				j |= (k & 1);
				k = k >> 1;
			}
			newx[j] = x[i];
		}

#ifdef DEBUG
		std::cout << "after bit reversed change the x[N] is: " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << newx[i] << "   ";
		std::cout << std::endl;
#endif // DEBUG

		for (int i = 0; i < log(n) / log(2); i++)		//一级蝶形运算 stage
		{
			int m = 2 << i;
			CM tmpw = w[n / m];
			CM z = 1;
			for (int j = 0; j < m / 2; j++)
			{
				for (int k = j; k < n; k += m)
				{
					CM v = z * newx[k + m / 2];
					CM u = newx[k];
					newx[k] = u + v;
					newx[k + m / 2] = u - v;
				}
				z *= tmpw;
			}
#ifdef DEBUG
			std::cout << "the *********** " << i << " ******** x[N] are as follows： " << std::endl;
			for (int i = 0; i < n; i++)
				std::cout << newx[i] << std::endl;
#endif // DEBUG
		}
		delete[]w;
		return newx;
	}
}