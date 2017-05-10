#pragma once
#include <omp.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <complex>
#include <iostream>

namespace par
{
	const double EPSINON = 0.0001;
	const double PI = acos(-1.0);
	typedef std::complex<double> CM;

	template <typename T>
	void psort(T *arr, int n, int p = omp_get_max_threads(), int sortType = 1)
	{
		assert(arr != NULL);
		assert(n > 0);

		while (n % p || n / p < p)
			p--;

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
			/***********3：正则采样**********/
			k = tid * p;		//每段的样本的起始存储位置
			j = 0;
			for (i = startPos; i < startPos + len; i += gap)
			{
				tmp[k + j] = arr[i];		//每段选出P个样本
				j++;
			}
#pragma omp barrier
#pragma omp single
			{
				/*********4：样本排序**********/
				std::sort(tmp, tmp + square_p);
				/*********5：选择主元**********/
				for (i = p, j = 0; i < square_p; i += p)	//用一台处理器选取p-1个主元
					pivot_element[j++] = tmp[i];
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
				for (i = 1; i <= p; i++)
					mergeIndx[i] = mergeIndx[i] + mergeIndx[i - 1];	//每个处理器归并结果存放区域	
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

	template <typename T1, typename T2, typename T3>
	void pMatrixMult(T1 **a, T2 **b, T3**c, int m, int s, int n, int p = omp_get_max_threads())
	{
		assert(a != NULL && b != NULL && c != NULL);
		assert(m > 0 && s > 0 && n > 0 && p > 0);

		//转置B矩阵，保证cache命中率
		T2 **tmpB;
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
		assert(a != NULL && b != NULL);
		assert(n > 0 && p > 0);

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

		delete[]marked;
		delete[]pivot;
		return x;
	}

	template <typename T>
	CM* pfft(T *x, const int n, int p = omp_get_max_threads())
	{
		assert(x != NULL && n > 0);
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
		omp_set_num_threads(p);
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

			for (int s = 1; s <= (int)(log(p) / log(2)); s++)		//************3*****p个核需要的通信次数
			{
				int m = 1 << (s + (int)(log(gap) / log(2)));
				CM w = ww[n / m];
				if (tid / (1 << (s - 1)) & 1)
				{
					CM z = ww[(startPos % (m / 2)) * n / m];
					for (int j = startPos; j < startPos + gap; j++)
					{
						src[j] = z * src[j];
						z = z * w;
					}
				}
				int partnerPos = (tid ^ (1 << (s - 1))) * gap;
#pragma omp barrier
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
			}
		}
		delete[]ww;
		delete[]dst;
		return src;
	}
}