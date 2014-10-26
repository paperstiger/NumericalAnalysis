/* Dev by tiger at 2014/10/25
 * Program for LU decomposition using (2.6) and (2.7)
 */
#include "VectorMatrix.h"
#include <iostream>
using namespace std;
int LUDeComp(double *A, double *L, double *U, int n, int lda, int ldL, int ldU)
{
#define A(i, j) A[(i - 1)*lda + j - 1]
#define L(i, j) L[(i - 1)*ldL + j - 1]
#define U(i, j) U[(i - 1)*ldU + j - 1]
    for(int k = 1;k <= n;k++)
    {
        for(int j = k;j <= n;j++)
        {
            U(k, j) = A(k, j) - cDot(1.0, L + (k - 1)*ldL, U + j - 1, k - 1, 1, ldU);
        }
        L(k, k) = 1.0;
        for(int i = k + 1;i <= n;i++)
        {
            L(i, k) = (A(i, k) - cDot(1.0, L + (i - 1)*ldL, U + k - 1, k - 1, 1, ldU))/U(k, k);
        }
    }
    return 1;
#undef A
#undef L
#undef U
}
//New LU comp without use of L and U
int LUDeComp(double *A, int N, int lda)
{
#define A(i, j) A[(i - 1)*lda + j - 1]
	//Operate only N - 1 times
    for(int k = 1;k <= N - 1;k++)
    {
        //Row k no change, for col k
        for(int j = k + 1;j <= N;j++)
        {
        	double ljk = A(j, k)/A(k, k);
        	A(j, k) = ljk;
        	aXpY(-ljk, A + lda*(k - 1) + k, A + lda*(j - 1) + k, N - k, 1, 1);
        }
    }
#undef A
	return 1;
}
//With Row changing, PA = LU can be got
int LUDecomp(double *A, double *P, int N, int lda, int ldp)
{
	//Allocate an array to record row change
	int *Order = new int[N];
	if(Order == NULL)
	{
		cout << "Allocate error\n";
		return 0;
	}
	for(int i = 0;i < N;i++)
		Order[i] = i;
	//from row 1 to row N - 1
	for(int i = 0;i < N - 1;i++)
	{
		//search for maxabs
		double *MaxRow = MaxAbs(A + lda*i + i, N - i, lda);
		int maxrow = (MaxRow - (A + i))/lda;//C style
		//swap row maxrow and i in memory
		if(i != maxrow)
		{
			Swap(A + i*lda, A + maxrow*lda, N);
			Swap(Order + i, Order + maxrow, 1);
		}
		//Row i no change, for col i
        for(int j = i + 1;j < N;j++)
        {
        	double ljk = A[j*lda + i]/A[i*lda + i];
        	A[j*lda + i] = ljk;
        	aXpY(-ljk, A + lda*i + i + 1, A + lda*j + i + 1, N - i - 1, 1, 1);
        }
	}
	//Use Order to change P
	//set all to 0
	for(int i = 0;i < N;i++)
	{
		for(int j = 0;j < N;j++)
			P[i*ldp + j] = 0.0; 
	}
	for(int i = 0;i < N;i++)
		P[i*ldp + Order[i]] = 1.0;
	
	delete[] Order;
	return 1;
}
//Main function to test LUDeComp
int main()
{
    double D[20] = {6,2,1,-1,1,2,4,1,0,2,1,1,4,-1,9,-1,0,-1,3,9};
    double E[20];
    double F[9] = {1, 2, 3, 4, 5, 6, 7, 8, 0};
    for(int i = 0;i < 20;i++)
    	E[i] = D[i];
    double L[16], U[16], P[20];
    for(int i = 0;i < 16;i++)
    {
        L[i] = 0.0;
        U[i] = 0.0;
    }
    int n = 4, lda = 5, ldp = 5, nF = 3;
    cout << "LU with User Input L and U\n";
    LUDeComp(D, L, U, n, lda, n, n);
    Output(D, n, lda, lda);
    Output(L, n, n, n);
    Output(U, n, n, n);
    cout << "LU with no input\n";
    LUDeComp(D, n, lda);
    Output(D, n, n, lda);
    cout << "LU with P and RowMajor\n";
    LUDecomp(F, P, nF, nF, ldp);
    Output(F, nF, nF, nF);
    Output(P, nF, nF, ldp);
    return 0;
}
