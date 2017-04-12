#include "omplib.h"
#include <time.h>
#include <iostream>
#include <functional>

using namespace std;

int main()
{
	clock_t pt = 0, st = 0;
	int n = 1000000, p = 4;

/*
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



	printf("serial time: %lf s.\n", (double)st / CLOCKS_PER_SEC);
	printf("parallel time: %lf s.\n", (double)pt / CLOCKS_PER_SEC);
	printf("speedup ratio: %lf \n", (double)st / pt);

	system("pause");
	return 0;
}