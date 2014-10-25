/* Dev by tiger at 2014/10/25
 * Header file for lots of Vector and Matrix operation
 */
#ifndef VECTORMATEIX_H
#define VECTORMATRIX_H
#include <iostream>
#include <iomanip>
using namespace std;
//y = ax + by + c
//a,b,c is num, x,y is vector
template<class T> void aXbYc(T a, T *x, T b, T *y, T c, int n, int incx = 1, int incy = 1)
{
    for(int i = 0;i < n;i++)
        y[i*incy] = a*x[i*incx] + b*y[i*incy] + c;
}
//N = c*dot(x,y)
template<class T> T cDot(T c, T *x, T *y, int n, int incx = 1, int incy = 1)
{
    if(n < 1)
        return 0;
    T tmp = 0;
    for(int i = 0;i < n;i++)
        tmp += c*x[i*incx]*y[i*incy];
    return tmp;
}
template<class T> void Output(T *A, int row, int col, int lda)
{
    for(int i = 0;i < row;i++)
    {
        for(int j = 0;j < col;j++)
        {
            cout << A[i*lda + j] << "\t";
        }
        cout << endl;
    }
}
#endif
