#pragma once
#include <omp.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <complex>

//#define DEBUG
#ifdef DEBUG	//*******���������Ϣ
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
		*(1)���Ȼ���: n��Ԫ�ؾ��ȵػ��ֳ�p�Σ�ÿ̨��������n/p��Ԫ��
		*(2)�ֲ�����: �����������ô��������㷨������n/p����
		*(3)�������: ÿ̨�����������Լ����������ѡȡp������Ԫ��(���εĵ�1��n/p^2+1��2n/p^2+1������Ԫ��)
		*(4)��������: ��һ̨������������p^2������Ԫ���ô��������㷨����֮
		*(5)ѡ����Ԫ: ��һ̨������ѡȡp-1����Ԫ�������䲥�͸����ദ����(��p+1��2p+1������Ԫ��)
		*(6)��Ԫ����: ������������Ԫ�����Ե�����λ��ֳ�p��
		*(7)ȫ�ֽ���: ������������Ͻ�ΰ��κŽ�������Ӧ�Ĵ�������
		*(8)�鲢����: ��������ʹ�ù鲢���������յ����ʩ������
		*/
		assert(arr != NULL);
		assert(n % p == 0);
		assert(n > 0);

		int len = n / p;
		int square_p = p * p;
		int gap = n / square_p;

		T *tmp = new T[n];
		T *pivot_element = new T[p];		//������Ԫ,p-1��
		int **sgm_pvt_elm = new int*[p];	//������Ԫ���ֵ�λ��;
		for (int i = 0; i < p; i++)
			sgm_pvt_elm[i] = new int[p + 1];
		int *mergeIndx = new int[p + 1];	//ȫ�ֽ����鲢����ÿ���������������ܳ���
		mergeIndx[p] = 0;

		omp_set_num_threads(p);
#pragma omp parallel shared(arr, tmp, pivot_element, sgm_pvt_elm, mergeIndx)
		{
			int i, j, k;
			int tid = omp_get_thread_num();
			/***********2���ֲ�����*********/
			int startPos = tid * len;	//ÿ�ε���ʼλ��
			std::sort(arr + startPos, arr + startPos + len);
#ifdef DEBUG	//**********�������
#pragma omp single
			for (i = 0; i < p; i++)
			{
				for (j = i * len; j < (i + 1) * len; j++)
					std::cout << setw(4) << arr[j];
			}
#endif // DEBUG

			/***********3���������**********/
			k = tid * p;		//ÿ�ε���������ʼ�洢λ��
			j = 0;
			for (i = startPos; i < startPos + len; i += gap)
			{
				tmp[k + j] = arr[i];		//ÿ��ѡ��P������
				j++;
			}
#ifdef DEBUG	//**********�������
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
				/*********4����������**********/
				std::sort(tmp, tmp + square_p);
#ifdef DEBUG	//**********�������
				std::cout << endl << "sort of sample: " << endl;
				for (k = 0; k < square_p; k++)
					std::cout << setw(4) << tmp[k];
#endif // DEBUG
				/*********5��ѡ����Ԫ**********/
				for (i = p, j = 0; i < square_p; i += p)	//��һ̨������ѡȡp-1����Ԫ
					pivot_element[j++] = tmp[i];
#ifdef DEBUG	//**********�������
				std::cout << endl << "pivot_element: " << endl;
				for (j = 0; j < p - 1; j++)	//����p-1����Ԫ
					std::cout << setw(4) << pivot_element[j];
#endif // DEBUG
			}
			/*********6����Ԫ����**********/
			sgm_pvt_elm[tid][0] = startPos;	//�����������������λ�ã�Ϊ��һ��λ�õ�ǰһλ
			sgm_pvt_elm[tid][p] = startPos + len;	//ĩβλ��
			for (j = 0, k = startPos; j < p - 1; j++)		//����������p-1����Ԫ�����Ե�����λ��ֳ�p��
			{
				for (; k < startPos + len; k++)	//�ҵ���Ӧλ�ô洢�±�
				{
					if (arr[k] > pivot_element[j])
						break;
				}
				sgm_pvt_elm[tid][j + 1] = k;	//�洢��Ԫ���ֵ�λ��
			}
#ifdef DEBUG	//**********�������
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
			/************7ȫ�ֽ���**************/
			mergeIndx[tid] = 0;		//merge_pos[p] = 0; ����ռ�ʱ�������һ��Ԫ��Ϊ0
#pragma omp barrier
			for (i = 0; i < p; i++)	//ÿ����������Ҫ�鲢p��
			{	//������������Ͻ�ΰ��γ��Ƚ�����������µ�Ͻ�γ��ȣ���ÿ���̹߳鲢������Ҫ���ܳ���
				mergeIndx[tid + 1] += (sgm_pvt_elm[i][tid + 1] - sgm_pvt_elm[i][tid]); //��Ϊ����λ�ã����㳤�Ȳ���Ҫ��1
			}
#pragma omp barrier
#pragma omp single
			{
#ifdef DEBUG	//*******�������
				std::cout << endl << "merge_length: " << endl;
				for (i = 0; i < p + 1; i++)
					std::cout << setw(4) << mergeIndx[i];
#endif // DEBUG
				for (i = 1; i <= p; i++)
					mergeIndx[i] = mergeIndx[i] + mergeIndx[i - 1];	//ÿ���������鲢����������	
#ifdef DEBUG	//*******�������
				std::cout << endl << "merge_sort_begin_at: " << endl;
				for (i = 0; i < p + 1; i++)
					std::cout << setw(4) << mergeIndx[i];
				std::cout << endl << "A: ";
				for (i = 0; i < n; i++)
					std::cout << setw(3) << arr[i];
#endif // DEBUG
			}
			/************8�鲢����**************/
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
			}	//�õ�maxNumber��maxIndex
#ifdef DEBUG	//**********�������
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

			for (i = mergeIndx[tid]; i < mergeIndx[tid + 1]; i++)	//p·�鲢���λ��
			{
				T m = max;
				ps = maxIndex;
				for (j = 0; j < p; j++)	//��p�����ҳ���С��
				{
					if (arr[index[j]] <= m && index[j] <= endIndex[j])
					{
						m = arr[index[j]];	//��¼��Сֵ
						ps = j;		//��¼��Сֵ���ڵ�tid
					}
				}
				tmp[i] = m;
				index[ps]++;		//�´β��Ƚϸ�ֵ
			}
			delete[] index;
			delete[] endIndex;
#pragma omp barrier
#ifdef DEBUG	//**********�������
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

		//ת��B���󣬱�֤cache������
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
			for (int i = 0; i < n; i++)		//ѡ��Ԫ
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
			for (int i = 0; i < n; i++)	//������Ԫ
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
		for (int i = 0; i < n; i++)	//���x
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
		if ((n & n - 1) != 0)   //n����Ϊ2����
		{
			std::cout << "Number of x[N] must be a power of 2 !" << std::endl;
			return nullptr;
		}
		while ((n % p) != 0)		//p��������n
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
		for (int i = 0; i < n; i++)     //************1*****λ���洢����ʼ��
		{
			int k = i;
			int j = 0;
			for (int t = (int)(log(n) / log(2)); t > 0; t--)    //ʵ����λ�ߵ�  
			{
				j = j << 1;
				j |= (k & 1);
				k = k >> 1;
			}
			src[j] = x[i];
			ww[i] = CM(cos(theta * i), sin(theta * i));
		}

#ifdef DEBUG
		std::cout << "The pfft begin�� " << std::endl;
		std::cout << "The w[N] are as follows�� " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << ww[i] << std::endl;
		std::cout << "After bit reversed change newx[N]: " << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << src[i] << "   ";
		std::cout << std::endl;
#endif // DEBUG

#pragma omp parallel shared(src, dst)
		{
			int tid = omp_get_thread_num();		//************2*****��ͨ��Ҫ��ĵ�������
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
				std::cout << "after 2nd step the  ******** x[N] are as follows�� " << std::endl;
				for (int i = 0; i < n; i++)
					std::cout << src[i] << std::endl;
			}
#endif // DEBUG

			for (int s = 1; s <= (int)(log(p) / log(2)); s++)		//************3*****p������Ҫ��ͨ�Ŵ���
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
					std::cout << "the **** " << s << " **** x[N] * w[N] = newx ****** are as follows�� " << std::endl;
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
					std::cout << "the *********** " << s << " ******** x[N] are as follows�� " << std::endl;
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