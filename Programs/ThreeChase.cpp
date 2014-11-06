/* Dev by tiger at 2014/11/6
 * Program for solving band3 matrix using Chasing method
 */
#include "VectorMatrix.h"
#include <iostream>
using namespace std;
using namespace TigerVecMat;
//Band stored matrix A which is stored in a 3*n length array
//a is stored in row(2:n) col(1) b stored in row(1:n) col(2) c stored in row (1:n-1) col(2)
int ChaseBand(int nA, int nb, double *A, int lda, double *b, int ldb)
{
    //LU decomposition store in A
    //first row no change of course
    for(int i = 1;i < nA;i++)
    {
        //calc l_i and set A[l_i + 1, 1]
        A[lda*i] = A[lda*i]/A[lda*(i - 1) + 1];
        //change U_i
        A[lda*i + 1] = A[lda*i + 1] - A[lda*i]*A[lda*(i - 1) + 2];
    }
    //solve Ly=d
    //first row no change
    for(int i = 1;i < nA;i++)
    {
        aXpY(-A[lda*i], b + ldb*(i - 1), b + ldb*i, nb, 1, 1);
    }
    //solve Ux=y
    //last row
    Scale(b + ldb*(nA - 1), 1/A[lda*(nA - 1) + 1], nb, 1);
    //last two to one
    for(int i = nA - 2;i >= 0;i--)
    {
        //yi -= ci*y(i+1)
        aXpY(-A[i*lda + 2], b + (i + 1)*ldb, b + i*ldb, ldb, 1, 1);
        //divide by ui
        Scale(b + i*ldb, 1/A[i*lda + 1], nb, 1);
    }
    //success
    return 1;
}
int main()
{
    int nA = 3, nb = 2, lda = 4, ldb = 3;
    double A[12] = {0, 1, 2, 0, 3, 4, 5, 1, 6, 7, 0, 1};
    double b[9] = {3, 5, 0, 12, 26, 3, 13, 33, 5};
    ChaseBand(nA, nb, A, lda, b, ldb);
    cout << "A = " << endl;
    Output(A, nA, nA, lda);
    cout << "b = " << endl;
    Output(b, nA, nb, ldb);
    return 0;
}
