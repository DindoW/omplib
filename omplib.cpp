#include "omplib.h"

int sGaussElim(double **a, double *b, double *x, int n)
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
			return -1;
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

	return 0;
}

int pGaussElim(double **a, double *b, double *x, int n, int p)
{
	double max;
	int picked;
	int *marked = new int[n];
	for (int i = 0; i < n; i++)
		marked[i] = 0;
	int *pivot = new int[n];
	double *maxPivot = new double[p];

	omp_set_num_threads(p);
	for (int k = 0; k < n; k++)
	{
		max = 0;
		for (int i = 0; i < n; i++)		//选主元
		{
			if (marked[i] == 0 && fabs(a[i][k]) > max)
			{
				max = fabs(a[i][k]);
				picked = i;
			}
		}
		if (fabs(max - 0) < EPSINON)
		{
			delete[]marked;
			delete[]pivot;
			return -1;
		}
		marked[picked] = 1;
		pivot[picked] = k;

#pragma omp parallel for shared(a, b, max) schedule(dynamic)
		for (int i = 0; i < n; i++)	//并行消元
		{
			double tmp = a[i][k] / a[picked][k];
			if (i != picked)
			{
				for (int j = k; j < n; j++)
					a[i][j] -= tmp * a[picked][j];
				b[i] -= tmp * b[picked];
			}
		}
	}
	
#pragma omp parallel for shared(a, b, x, pivot) schedule(dynamic)
	for (int i = 0; i < n; i++)	//求解x
		x[pivot[i]] = b[i] / a[i][pivot[i]];


#ifdef DEBUG
	std::cout << std::endl << "after parallel Gauss-Jordan Elimination A & B & pivot :" << std::endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(8) << a[i][j];
		std::cout << std::setprecision(3) << std::setw(8) << b[i] << std::setw(6) << pivot[i] << std::endl;
	}
	std::cout << "X : " << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << std::setw(8) << x[i];
	std::cout << std::endl << "threads: " << p << std::endl;
#endif // DEBUG

	delete[]marked;
	delete[]pivot;
	return 0;
}