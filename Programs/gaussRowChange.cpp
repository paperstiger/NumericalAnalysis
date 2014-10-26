/* Dev by tiger at 2014/10/26
 * Subroutine to perform Gauss elimination method with row change
 */
#include <iostream>
#include "VectorMatrix.h"
using namespace std;
/* complex interface for Gauss elimination for sloving Ax=b
 * nA is dimension of A, nb is colomn of b
 * A is pointer to first element of matrix A, lda is leading column of A
 * b is pointer to first element of b, ldb is leading column of b
 */
int Gauss(int nA, int nb, double *A, int lda, double *b, int ldb)
{
	//establish an array to store order
	int *Order = new int[nA];
	if(Order == NULL)
	{
		cout << "Memory allocate error\n";
		return -1;
	}
	for(int i = 0;i < nA;i++)
		Order[i] = i;
	//nA - 1 operations needed to make A be an upper tri mat
	for(int i = 0;i < nA - 1;i++)
	{
		//find max col
		double *MaxRow = MaxAbs(A + lda*i + i, nA - i, lda);
		int maxrow = (MaxRow - (A + i))/lda;//C style
		//swap row maxrow and i, not in memory but in order
		if(i != maxrow)
			Swap(Order + i, Order + maxrow, 1);
		//row operation from Order[i+1] with Order[i]
		for(int j = i + 1;j < nA;j++)
		{
			double lji = A[lda*Order[j] + i]/A[lda*Order[i] + i];
			aXpY(-lji, A + lda*Order[i] + i, A + lda*Order[j] + i, nA - i);
			aXpY(-lji, b + ldb*Order[i], b + ldb*Order[j], nb);
		}
	}
	//it seems that U got, now back and got I
	for(int i = nA - 1;i >= 0;i--)
	{
		//operate row Order[i]
		Scale(b + ldb*Order[i], 1/A[lda*Order[i] + i], nb, 1);
		A[lda*Order[i] + i] = 1;
		for(int j = i - 1;j >= 0;j--)
		{
			aXpY(-A[lda*Order[j] + i], b + ldb*Order[i], b + ldb*Order[j], nb, 1, 1);
			A[lda*Order[j] + i] = 0.0;
		}
	}
	//Change something to make it right
	//Here Sacrifice space for time
	double *Saveb = new double[nA*nb];
	if(Saveb == NULL)
	{
		cout << "Memory allocate error\n";
		return -1;
	}
	//Perhaps Matrix Copy should be added later
	for(int i = 0;i < nA;i++)
	{
		Copy(Saveb + i*nb, b + i*ldb, nb, 1, 1);
	}
	//find each element order in Order
	for(int i = 0;i < nA;i++)
	{
		//int index = 0;
		//while(Order[index] != i)
		//	index++;
		Copy(b + i*ldb, Saveb + nb*Order[i], nb, 1, 1);
	}
	//Memory free
	delete[] Order;
	delete[] Saveb;
    return 1;
}
/* simple interface for Gauss
 * n is dimension of problem
 * A is pointer to first element of matrix A
 * b is pointer to first element of b
 * after this proc, A become unit matrix and b is the sol
 */
int Gauss(int n, double *A, double *b)
{
    return Gauss(n, 1, A, n, b, 1);
}
int main()
{
    double A[9] = {1,2,3,2,3,1,3,2,1};
    double b[9] = {6,7,10,6,8,11,6,9,12};
    int n = 3;
    int flag = Gauss(n, 2, A, n, b, n);
    //output b
    for(int i = 0;i < n*n;i++)
        cout << b[i] << endl;
    return 0;

}
