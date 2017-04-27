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

CM* pfft(double *x, const int n, int p)
{
	if ((n & n - 1) != 0)   //n必须为2的幂
		return nullptr;
	if ((n % p) != 0)		//p必须整除n
		return nullptr;

	CM **newx = new CM*[2];
	newx[0] = new CM[n];
	newx[1] = new CM[n];
	CM *src = newx[0];
	CM *dst = newx[1];
	CM *ww = new CM[n];
	double theta = 2 * PI / n;
	int gap = n / p;
#pragma omp parallel for shared(src, ww) schedule(dynamic)
	for (int i = 0; i < n; i++)     //************1*****位反存储
	{
		int k = i;
		int j = 0;
		for (int t = (int)(log(n) / log(2)); t > 0; t--)    //利用按位与以及循环实现码位颠倒  
		{
			j = j << 1;
			j |= (k & 1);
			k = k >> 1;
		}
		src[j] = x[i];
		ww[i] = CM(cos(theta * i), sin(theta * i));
	}

#ifdef DEBUG
	std::cout << "The pppppppppppppppppppppfft begin： " << std::endl;
	std::cout << "The w are as follows： " << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << ww[i] << std::endl;
	std::cout << "after change the newx[N] is: " << std::endl;
	for (int i = 0; i < n; i++)
		std::cout << src[i] << "   ";
	std::cout << std::endl;
#endif // DEBUG

#pragma omp parallel shared(src, dst)
	{
		int tid = omp_get_thread_num();		//************2*****无通信要求的迭代计算
		int startPos = tid * gap;
		for (int s = 1; s <= (int)(log(gap) / log(2)); s++)
		{
			int m = 1 << s;		//m = 2 ^ s;
			CM w = ww[n / m];
			CM z = 1;
			for (int j = 0; j < m / 2; j++)
			{
				for (int k = startPos + j; k < startPos + gap; k += m)
				{
					CM q = z * src[k + m / 2];
					CM r = src[k];
					src[k] = r + q;
					src[k + m / 2] = r - q;
				}
				z = z * w;
			}
		}
#pragma omp barrier
#ifdef DEBUG
#pragma omp single
		{
			std::cout << "after 2nd step the  ******** x[N] are as follows： " << std::endl;
			for (int i = 0; i < n; i++)
				std::cout << src[i] << std::endl;
		}
#endif // DEBUG
		for (int s = 1; s <= (int)(log(p) / log(2)); s++)		//************3*****p个核需要的通信次数
		{
			int m = 1 << (s + (int)(log(gap) / log(2)));
			CM w = ww[n / m];
#ifdef DEBUG
#pragma omp single
			std::cout << " m : " << m  << " w : " << n / m << std::endl;
#endif // DEBUG
			if (tid / (1 << (s - 1)) & 1)
			{
				CM z = ww[(startPos % (m / 2)) * n / m];
#ifdef DEBUG3
#pragma omp critical
				std::cout << "tid: "  << tid << " newx[j + 1] : " << src[tid * (n / p) + 1]
					<< "  z * w:  " << z * w << " newx * z : " << src[tid * (n / p) + 1] * z * w << std::endl;
#endif // DEBUG
				for (int j = startPos; j < startPos + gap; j++)
				{
					src[j] = z * src[j];
					z = z * w;
				}
			}
			int partnerPos = (tid ^ (1 << (s - 1))) * gap;
#pragma omp barrier
#ifdef DEBUG
#pragma omp single
			{
				std::cout << "the **** " << s << " **** x[N] * w[N] = newx ****** are as follows： " << std::endl;
				for (int i = 0; i < n; i++)
					std::cout << src[i] << std::endl;
			}
#endif // DEBUG
			if ((tid / (1 << (s - 1))) & 1)
			{
				for (int j = 0; j < gap; j++)
					dst[startPos + j] = src[partnerPos + j] - src[startPos + j];
			}
			else
			{
				for (int j = 0; j < gap; j++)
					dst[startPos + j] = src[startPos + j] + src[partnerPos + j];
			}
#pragma omp barrier
#pragma omp single
			{
				src = newx[s % 2];
				dst = newx[(s + 1) % 2];
			}
#ifdef DEBUG
#pragma omp single
			{
				std::cout << "the *********** " << s << " ******** x[N] are as follows： " << std::endl;
				for (int i = 0; i < n; i++)
					std::cout << src[i] << std::endl;
			}
#endif // DEBUG
		}
	}
	delete[]ww;
	delete[]dst;
	return src;
}
