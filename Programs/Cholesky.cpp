/* Dev by tiger at 2014/11/6
 * function using Cholesky decomposition to solve Ax = b
 * which A is a symmetric matrix
 */
#include "VectorMatrix.h"
#include <iostream>
using namespace std;
using namespace TigerVecMat;
//A is stored in full mode
//which allow us to have space for computing
int Cholesky(int nA, int nb, double *A, int lda, double *b, int ldb)
{
    //A = LLt
    //First Row diffrent than others
    A[0] = sqrt(A[0]);
    Scale(A + lda, 1/A[0], nA - 1, lda);
    //two to last
    for(int i = 1;i < nA;i++)
    {
        double sumOfRow = cDot(1.0, A + lda*i, A + lda*i, i, 1, 1);
        A[lda*i + i] = sqrt(A[lda*i + i] - sumOfRow);
        for(int j = i + 1;j < nA;j++)
        {
            double sumOfTwoRow = cDot(1.0, A + lda*i, A + lda*j, i, 1, 1);
            A[lda*j + i] = (A[lda*j + i] - sumOfTwoRow)/A[lda*i + i];
        }
    }
    //Solve Ly=b
    for(int colb = 0;colb < nb;colb++)
    {
        for(int i = 0;i < nA;i++)
        {
            double sumOfLy = cDot(1.0, A + lda*i, b + colb, i, 1, ldb);
            b[ldb*i + colb] = (b[ldb*i + colb] - sumOfLy)/A[lda*i + i];
        }
    }
    //Solve Lt x = y
    for(int colb = 0;colb < nb;colb++)
    {
        //last row is diffrent to avoid out of bound of array
        b[(nA - 1)*ldb + colb] /= A[(nA - 1)*lda + nA - 1];
        for(int i = nA - 2;i >= 0;i--)
        {
            double sumOfLx = cDot(1.0, A + (i + 1)*lda + i, b + (i + 1)*ldb + colb, nA - i - 1, lda, ldb);
            b[i*ldb + colb] = (b[i*ldb + colb] - sumOfLx)/A[i*lda + i];
        }
    }
    return 1;
}
int main()
{
    int nA = 3, nb = 1, lda = 4, ldb = 2;
    double A[12] = {16, 4, 8, 1, 4, 5, -4, 2, 8, -4, 22, 1};
    double b[6] = {-4, 1, 3, 1, 10, 1};
    Cholesky(nA, nb, A, lda, b, ldb);
    cout << "A = " << endl;
    Output(A, nA, nA, lda);
    cout << "b = " << endl;
    Output(b, nA, nb, ldb);
    return 0;
}
