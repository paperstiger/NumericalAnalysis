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
#undef A
#undef L
#undef U
}
//Main function to test LUDeComp
int main()
{
    double D[20] = {6,2,1,-1,1,2,4,1,0,2,1,1,4,-1,9,-1,0,-1,3,9};
    double L[16], U[16];
    for(int i = 0;i < 16;i++)
    {
        L[i] = 0.0;
        U[i] = 0.0;
    }
    int n = 4, lda = 5;
    LUDeComp(D, L, U, n, lda, n, n);
    Output(D, n, lda, lda);
    Output(L, n, n, n);
    Output(U, n, n, n);
    return 0;
}
