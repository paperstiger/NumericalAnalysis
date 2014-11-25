/* Dev by tiger at 2014/11/25
 * subroutine to implement Jacobi Iter to solve Linear Equations
 */
#include "VectorMatrix.h"
#include <iostream>
using namespace std;
using namespace TigerVecMat;
const int MAXJACOBIITER = 100;
const int MAXGSITER = 100;
//subroutine for Jacobi Iter, Knowing x, calc Nextx
int JacobiNext(int N, double *A, int lda, double *b, int incb, double *x, int incx, double *NewX, int incNew)
{
    for(int i = 0;i < N;i++)
    {
        double tmp1 = cDot(1.0, A + i*lda, x, i, 1, incx);
        double tmp2 = cDot(1.0, A + i*lda + i + 1, x + i + 1, N - i - 1, 1, incx);
        NewX[i*incNew] = (b[i*incb] - tmp1 - tmp2)/A[i*lda + i];
    }
    return 1;
}
//subroutine for GS Iter, Knowing x, calc NewX
int GSNext(int N, double *A, int lda, double *b, int incb, double *x, int incx, double *NewX, int incNew, double omega = 1)
{
    if(1 == omega)
    {
        for(int i = 0;i < N;i++)
        {
            double tmp1 = cDot(1.0, A + i*lda, NewX, i, 1, incNew);
            double tmp2 = cDot(1.0, A + i*lda + i + 1, x + i + 1, N - i - 1, 1, incx);
            NewX[i*incNew] = (b[i*incb] - tmp1 - tmp2)/A[i*lda + i];
        }
    }
    else
    {
        for(int i = 0;i < N;i++)
        {
            double tmp1 = cDot(1.0, A + i*lda, NewX, i, 1, incNew);
            double tmp2 = cDot(1.0, A + i*lda + i + 1, x + i + 1, N - i - 1, 1, incx);
            NewX[i*incNew] = (1 - omega)*x[i*incx] + omega*(b[i*incb] - tmp1 - tmp2)/A[i*lda + i];
        }

    }
    return 1;
}
//Jacobi Iter to Solve Ax=b
//A (N*N) matrix, main element not equal to 0
//b(N) vector 
//x(N) stores initial guess when in and stores answer when out
//return Iter times, -1 if failed
int JacobiIter(int N, double *A, int lda, double *b, int incb, double *x, int incx, double tol)
{
    double *NewX = new double[N];
    double *Diff = new double[N];
    if(NULL == NewX || NULL == Diff)
    {
        cout << "Allocate error\n";
        return -1;
    }
    bool NowX = true;//true if x stores newest data
    int iter = 1;
    while(iter < MAXJACOBIITER)
    {
        if(NowX)
        {
            JacobiNext(N, A, lda, b, incb, x, incx, NewX, 1);
            NowX = false;
            VMinus(Diff, NewX, x, N, 1, incx);
            if(VNorm2(Diff, N, 1) < tol)
            {
                Copy(x, NewX, N, incx, 1);
                break;
            }
        }
        else
        {
            JacobiNext(N, A, lda, b, incb, NewX, 1, x, incx);
            NowX = true;
            VMinus(Diff, x, NewX, N, incx, 1);
            if(VNorm2(Diff, N, 1) < tol)
                break;
        }
        iter++;
    }


    
    delete[] NewX;
    delete[] Diff;
    if(iter == MAXJACOBIITER)
        return -1;
    return iter;
}
//GS Iter to Solve Ax=b
//A (N*N) matrix, main element not equal to 0
//b(N) vector 
//x(N) stores initial guess when in and stores answer when out
//return Iter times, -1 if failed
int GSIter(int N, double *A, int lda, double *b, int incb, double *x, int incx, double tol, double omega = 1)
{
    double *NewX = new double[N];
    double *Diff = new double[N];
    if(NULL == NewX || NULL == Diff)
    {
        cout << "Allocate error\n";
        return -1;
    }
    bool NowX = true;//true if x stores newest data
    int iter = 1;
    while(iter < MAXGSITER)
    {
        if(NowX)
        {
            GSNext(N, A, lda, b, incb, x, incx, NewX, 1, omega);
            NowX = false;
            VMinus(Diff, NewX, x, N, 1, incx);
            if(VNorm2(Diff, N, 1) < tol)
            {
                Copy(x, NewX, N, incx, 1);
                break;
            }
        }
        else
        {
            GSNext(N, A, lda, b, incb, NewX, 1, x, incx, omega);
            NowX = true;
            VMinus(Diff, x, NewX, N, incx, 1);
            if(VNorm2(Diff, N, 1) < tol)
                break;
        }
        iter++;
        if(8 == iter)//For Page 205
        {
            cout << "When itering 8 times\n";
            Output(NewX, N, 1, 1);
        }
    }
    delete[] NewX;
    delete[] Diff;
    if(iter == MAXGSITER)
        return -1;
    return iter;
}
//Test
int main()
{
    int N = 3;
    double A[12] = {4, 3, 0, 0, 3, 4, -1, 0, 0, -1, 4, 0};
    double b[3] = {24, 30, -24};
    double x[3] = {1, 1, 1};
    double tol = 1e-7;
    cout << "Jacobi Iter\n"; 
    int iter = JacobiIter(N, A, 4, b, 1, x, 1, tol);
    cout << "After itering " << iter << " times:\n";
    Output(x, N, 1, 1);
    Scale(x, 0.0, N, 1);
    cout << "GS Iter\n";
    iter = GSIter(N, A, 4, b, 1, x, 1, tol);
    cout << "After itering " << iter << " times:\n";
    Output(x, N, 1, 1);
    cout << "SOR using 1.25\n";
    Scale(x, 0.0, N, 1);
    double omega = 1.25;
    iter = GSIter(N, A, 4, b, 1, x, 1, tol, omega);
    cout << "After Itering " << iter << " times:" << endl;
    Output(x, 3, 1, 1);
    return 0;
}
