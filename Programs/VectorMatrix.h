/* Dev by tiger at 2014/10/25
 * Header file for lots of Vector and Matrix operation
 */
#ifndef VECTORMATEIX_H
#define VECTORMATRIX_H
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
//Level 1 subroutines
/* Operator for Vector*/
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
//sum of a vector
template<class T> T Sum(T *A, int dim, int incx = 1)
{
	T sum = 0;
	for(int i = 0;i < dim;i++)
		sum += A[incx*i];
	return sum;
}
//Swap of two vectors
template<class T> void Swap(T *X, T *Y, int N, int incx = 1, int incy = 1)
{
	//if equal to 1 ,use faster way??
	if(incx == 1 && incy == 1)
	{
		for(int i = 0;i < N;i++)
		{
			T tmp = X[i];
			X[i] = Y[i];
			Y[i] = tmp;
		}
		return;
	}
	//otherwise use ordinary way
	int ix = 0, iy = 0;
	for(int i = 0;i < N;i++)
	{
		T tmp = X[ix];
		X[ix] = Y[iy];
		Y[iy] = tmp;
		ix += incx;
		iy += incy
	}
	return;
}
//Scale by a num
template<class T> void Scale(T *A, T b, int N, int incx = 1)
{
	if(incx == 1)
	{
		for(int i = 0;i < N;i++)
			A[i] *= b;
		return;
	}
	//else use stupid way
	int ix = 0;
	for(int i = 0;i < N;i++)
	{
		A[ix] *= b;
		ix += incx;
	}
}
//Copy from X to Y
template<class T> void Copy(T *X, T *Y, int N, int incx = 1, int incy = 1)
{
	if(incx == 1 && incy == 1)
	{
		for(int i = 0;i < N;i++)
			Y[i] = X[i];
		return;
	}
	//other inc
	int ix = 0, iy = 0;
	for(int i = 0;i < N;i++)
	{
		Y[iy] = X[ix];
		ix += incx;
		iy += incy;
	}
}
//Norm2 of a vector
template<class T> T Vnorm2(T *A, int N, int incx = 1)
{
	T tmp = 0;
	if(incx == 1)
	{
		for(int i = 0;i < N;i++)
			tmp += A[i]*A[i];
		return sqrt(tmp);
	}
	int ix = 0;
	for(int i = 0;i < N;i++)
	{
		tmp += A[ix];
		ix += incx;
	}
	return sqrt(tmp);
}	
//sum of absolute value for a vector
template<class T> T Asum(T *A, int N, int incx = 1)
{
	T tmp = 0;
	for(int i = 0;i < N;i++)
		tmp += abs(A[i*incx]);
	return tmp;
}


//Output a matrix
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
