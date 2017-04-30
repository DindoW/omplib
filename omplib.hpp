#pragma once
#include <omp.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <complex>

//#define DEBUG
#ifdef DEBUG	//*******输出调试信息
#include <iostream>
#include <iomanip>
#endif // DEBUG

namespace par
{
	const double EPSINON = 0.0001;
	const double PI = acos(-1.0);
	typedef std::complex<double> CM;

	template <typename T>
	void psort(T *arr, int n, int p = omp_get_max_threads(), int sortType = 1)
	{
		/*
		*(1)均匀划分: n个元素均匀地划分成p段，每台处理器有n/p个元素
		*(2)局部排序: 各处理器利用串行排序算法，排序n/p个数
		*(3)正则采样: 每台处理器各从自己的有序段中选取p个样本元素(各段的第1、n/p^2+1、2n/p^2+1、…个元素)
		*(4)样本排序: 用一台处理器将所有p^2个样本元素用串行排序算法排序之
		*(5)选择主元: 用一台处理器选取p-1个主元，并将其播送给其余处理器(第p+1、2p+1、…个元素)
		*(6)主元划分: 各处理器按主元将各自的有序段划分成p段
		*(7)全局交换: 各处理器将其辖段按段号交换到对应的处理器中
		*(8)归并排序: 各处理器使用归并排序将所接收的诸段施行排序
		*/
		assert(arr != NULL);
		assert(n % p == 0);
		assert(n > 0);

		int len = n / p;
		int square_p = p * p;
		int gap = n / square_p;

		T *tmp = new T[n];
		T *pivot_element = new T[p];		//储存主元,p-1个
		int **sgm_pvt_elm = new int*[p];	//储存主元划分的位置;
		for (int i = 0; i < p; i++)
			sgm_pvt_elm[i] = new int[p + 1];
		int *mergeIndx = new int[p + 1];	//全局交换归并排序每个处理器排序后的总长度
		mergeIndx[p] = 0;

		omp_set_num_threads(p);
#pragma omp parallel shared(arr, tmp, pivot_element, sgm_pvt_elm, mergeIndx)
		{
			int i, j, k;
			int tid = omp_get_thread_num();
			/***********2：局部排序*********/
			int startPos = tid * len;	//每段的起始位置
			std::sort(arr + startPos, arr + startPos + len);
#ifdef DEBUG	//**********测试输出
#pragma omp single
			for (i = 0; i < p; i++)
			{
				for (j = i * len; j < (i + 1) * len; j++)
					std::cout << setw(4) << arr[j];
			}
#endif // DEBUG

			/***********3：正则采样**********/
			k = tid * p;		//每段的样本的起始存储位置
			j = 0;
			for (i = startPos; i < startPos + len; i += gap)
			{
				tmp[k + j] = arr[i];		//每段选出P个样本
				j++;
			}
#ifdef DEBUG	//**********测试输出
#pragma omp single
			{
				for (i = 0; i < p; i++)
				{
					std::cout << endl << i << "tid: ";
					for (j = 0; j < p; j++)
						std::cout << setw(4) << tmp[i * p + j];
				}
			}
#endif // DEBUG
#pragma omp barrier
#pragma omp single
			{
				/*********4：样本排序**********/
				std::sort(tmp, tmp + square_p);
#ifdef DEBUG	//**********测试输出
				std::cout << endl << "sort of sample: " << endl;
				for (k = 0; k < square_p; k++)
					std::cout << setw(4) << tmp[k];
#endif // DEBUG
				/*********5：选择主元**********/
				for (i = p, j = 0; i < square_p; i += p)	//用一台处理器选取p-1个主元
					pivot_element[j++] = tmp[i];
#ifdef DEBUG	//**********测试输出
				std::cout << endl << "pivot_element: " << endl;
				for (j = 0; j < p - 1; j++)	//共有p-1个主元
					std::cout << setw(4) << pivot_element[j];
#endif // DEBUG
			}
			/*********6：主元划分**********/
			sgm_pvt_elm[tid][0] = startPos;	//各个处理器处理的首位置，为第一个位置的前一位
			sgm_pvt_elm[tid][p] = startPos + len;	//末尾位置
			for (j = 0, k = startPos; j < p - 1; j++)		//各处理器按p-1个主元将各自的有序段划分成p段
			{
				for (; k < startPos + len; k++)	//找到相应位置存储下标
				{
					if (arr[k] > pivot_element[j])
						break;
				}
				sgm_pvt_elm[tid][j + 1] = k;	//存储主元划分的位置
			}
#ifdef DEBUG	//**********测试输出
#pragma omp barrier
#pragma omp single
			{
				for (i = 0; i < p; i++)
				{
					std::cout << endl << "sgm_pvt_elm_pos_threads: " << i << endl;
					for (j = 0; j < p + 1; j++)
						std::cout << setw(4) << sgm_pvt_elm[i][j];
				}
			}
#endif // DEBUG		
			/************7全局交换**************/
			mergeIndx[tid] = 0;		//merge_pos[p] = 0; 申请空间时处理，最后一个元素为0
#pragma omp barrier
			for (i = 0; i < p; i++)	//每个处理器需要归并p段
			{	//各处理器将其辖段按段长度交换，计算出新的辖段长度，即每个线程归并排序需要的总长度
				mergeIndx[tid + 1] += (sgm_pvt_elm[i][tid + 1] - sgm_pvt_elm[i][tid]); //因为储存位置，计算长度不需要加1
			}
#pragma omp barrier
#pragma omp single
			{
#ifdef DEBUG	//*******测试输出
				std::cout << endl << "merge_length: " << endl;
				for (i = 0; i < p + 1; i++)
					std::cout << setw(4) << mergeIndx[i];
#endif // DEBUG
				for (i = 1; i <= p; i++)
					mergeIndx[i] = mergeIndx[i] + mergeIndx[i - 1];	//每个处理器归并结果存放区域	
#ifdef DEBUG	//*******测试输出
				std::cout << endl << "merge_sort_begin_at: " << endl;
				for (i = 0; i < p + 1; i++)
					std::cout << setw(4) << mergeIndx[i];
				std::cout << endl << "A: ";
				for (i = 0; i < n; i++)
					std::cout << setw(3) << arr[i];
#endif // DEBUG
			}
			/************8归并排序**************/
			T max;
			int maxIndex, ps;
			int *index = new int[p];
			int *endIndex = new int[p];

			for (i = 0; i < p; i++)
			{
				index[i] = sgm_pvt_elm[i][tid];
				endIndex[i] = sgm_pvt_elm[i][tid + 1] - 1;
			}

			max = arr[endIndex[0]];
			maxIndex = 0;
			for (i = 1; i < p; i++)
			{
				if (arr[endIndex[i]] >= max)
				{
					max = arr[endIndex[i]];
					maxIndex = i;
				}
			}	//得到maxNumber和maxIndex
#ifdef DEBUG	//**********测试输出
#pragma omp critical
			{
				std::cout << endl << "tid::" << tid;
				for (i = 0; i < p; i++)
				{
					std::cout << endl << "NO." << i << " from:" << index[i] << " to:" << endIndex[i] << endl;
					for (j = index[i]; j <= endIndex[i]; j++)
						std::cout << setw(4) << arr[j];
				}
				std::cout << endl << "max = " << max << " maxIndex = " << maxIndex;
			}
#endif // DEBUG

			for (i = mergeIndx[tid]; i < mergeIndx[tid + 1]; i++)	//p路归并存放位置
			{
				T m = max;
				ps = maxIndex;
				for (j = 0; j < p; j++)	//从p段中找出最小的
				{
					if (arr[index[j]] <= m && index[j] <= endIndex[j])
					{
						m = arr[index[j]];	//记录最小值
						ps = j;		//记录最小值所在的tid
					}
				}
				tmp[i] = m;
				index[ps]++;		//下次不比较该值
			}
			delete[] index;
			delete[] endIndex;
#pragma omp barrier
#ifdef DEBUG	//**********测试输出
#pragma omp single
			{
				std::cout << endl << "mergeIndx: ";
				for (i = 0; i <= p; i++)
					std::cout << setw(4) << mergeIndx[i];
			}
#pragma omp critical
			{
				std::cout << endl << "tid" << tid
					<< "  start:" << mergeIndx[tid] << "  end:" << mergeIndx[tid + 1] << endl;
				for (i = mergeIndx[tid]; i < mergeIndx[tid + 1]; i++)
					std::cout << setw(4) << tmp[i];
			}
#endif // DEBUG
			if (sortType)
			{
				for (i = startPos; i < startPos + len; i++)
					arr[i] = tmp[i];
			}
			else
			{
				for (i = startPos; i < startPos + len; i++)
					arr[n - 1 - i] = tmp[i];
			}
		}
		for (int i = 0; i < p; i++)
			delete[] sgm_pvt_elm[i];
		delete[] sgm_pvt_elm;
		delete[] pivot_element;
		delete[] mergeIndx;
		delete[] tmp;
	}

	template <typename T>
	T** pMatrixTrans(T **a, int m, int n, int p = omp_get_max_threads())
	{
		assert(a != NULL && m * n > 0);

		T **b = new T*[n];
		for (int i = 0; i < n; i++)
			b[i] = new T[m];

		omp_set_num_threads(p);
		if (m > n)
		{
#pragma omp parallel for shared(a, b) schedule(dynamic)
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					b[j][i] = a[i][j];
		}
		else
		{
#pragma omp parallel for shared(a, b) schedule(dynamic)
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
					b[i][j] = a[j][i];
		}
		return b;
	}

	template <typename T1, typename T2>
	void pMatrixMult(T1 **a, T1 **b, T2**c, int m, int s, int n, int p = omp_get_max_threads())
	{
		assert(a != NULL && b != NULL && c != NULL);
		assert(m > 0 && s > 0 && n > 0 && p > 0);

		//转置B矩阵，保证cache命中率
		T1 **tmpB;
		tmpB = pMatrixTrans(b, s, n);

		omp_set_num_threads(p);
		if (m > n)
		{
#pragma omp parallel for shared(a, tmpB, c) schedule(dynamic)
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
				{
					T2 temp = 0;
					for (int k = 0; k < s; k++)
						temp += a[i][k] * tmpB[j][k];
					c[i][j] = temp;
				}
		}
		else
		{
#pragma omp parallel for shared(a, tmpB, c) schedule(dynamic)
			for (int i = 0; i < n; i++)
				for (int j = 0; j < m; j++)
				{
					T2 temp = 0;
					for (int k = 0; k < s; k++)
						temp += tmpB[i][k] * a[j][k];
					c[j][i] = temp;
				}
		}

		for (int i = 0; i < n; i++)
			delete[]tmpB[i];
		delete[]tmpB;
	}

	template <typename T>
	double* pGaussElim(T **a, T *b, int n, unsigned int p = omp_get_max_threads())
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
				std::cout << "Rank of Matrix Error !" << std::endl;
				return nullptr;
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

		double *x = new double[n];

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
		return x;
	}

	template <typename T>
	CM* pfft(T *x, const int n, int p = omp_get_max_threads())
	{
		if ((n & n - 1) != 0)   //n必须为2的幂
		{
			std::cout << "Number of x[N] must be a power of 2 !" << std::endl;
			return nullptr;
		}
		while ((n % p) != 0)		//p必须整除n
			p--;

		CM **newx = new CM*[2];
		newx[0] = new CM[n];
		newx[1] = new CM[n];
		CM *src = newx[0];
		CM *dst = newx[1];
		CM *ww = new CM[n];

		double theta = 2 * PI / n;
		int gap = n / p;
#pragma omp parallel for shared(src, ww) schedule(dynamic)
		for (int i = 0; i < n; i++)     //************1*****位反存储及初始化
		{
			int k = i;
			int j = 0;
			for (int t = (int)(log(n) / log(2)); t > 0; t--)    //实现码位颠倒  
			{
				j = j << 1;
				j |= (k & 1);
				k = k >> 1;
			}
			src[j] = x[i];
			ww[i] = CM(cos(theta * i), sin(theta * i));
		}

#ifdef DEBUG
		std::cout << "The pfft begin： " << std::endl;
		std::cout << "The w[N] are as follows： " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << ww[i] << std::endl;
		std::cout << "After bit reversed change newx[N]: " << std::endl;
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
				std::cout << " m : " << m << " w : " << n / m << std::endl;
#endif // DEBUG

				if (tid / (1 << (s - 1)) & 1)
				{
					CM z = ww[(startPos % (m / 2)) * n / m];

#ifdef DEBUG
#pragma omp critical
					std::cout << "tid: " << tid << " newx[j + 1] : " << src[tid * (n / p) + 1]
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
}