/* Dev by tiger at 2014/10/25
 * Header file for lots of Vector and Matrix operation
 */
#ifndef VECTORMATEIX_H
#define VECTORMATRIX_H
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
namespace TigerVecMat
{
//Level 1 subroutines
/* Operator for Vector*/
//y = ax + by + c
//a,b,c is num, x,y is vector
template<class T> void aXbYc(T a, T *x, T b, T *y, T c, int n, int incx = 1, int incy = 1)
{
	if(n < 1)
		return;
    for(int i = 0;i < n;i++)
        y[i*incy] = a*x[i*incx] + b*y[i*incy] + c;
}
//y = ax + y
//basic operation of matrix col or row
template<class T> void aXpY(T a, T *X, T *Y, int N, int incx = 1, int incy = 1)
{
	if(N < 1)
		return;
	if(incx == 1 && incy == 1)
	{
		for(int i = 0;i < N;i++)
			Y[i] += a*X[i];
		return;
	}
	int ix = 0, iy = 0;
	for(int i = 0;i < N;i++)
	{
		Y[iy] += a*X[ix];
		ix += incx;
		iy += incy;
	}
	return;
}
//A = B + C
template<class T> void VAdd(T *A, T *B, T *C, int N, int inca = 1, int incb = 1, int incc = 1)
{
    for(int i = 0;i < N;i++)
        A[i*inca] = B[i*inca] + C[i*inca];
}
//A = B - C
template<class T> void VMinus(T *A, T *B, T *C, int N, int inca = 1, int incb = 1, int incc = 1)
{
    for(int i = 0;i < N;i++)
        A[i*inca] = B[i*incb] - C[i*incc];
}
//N = c*dot(x,y)
template<class T> T cDot(T c, T *x, T *y, int n, int incx = 1, int incy = 1)
{
    if(n < 1)
        return 0;
    T tmp = 0;
    for(int i = 0;i < n;i++)
        tmp += x[i*incx]*y[i*incy];
    return c*tmp;
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
	if(N < 1)
		return;
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
		iy += incy;
	}
	return;
}
//Scale by a num
template<class T> void Scale(T *A, T b, int N, int incx = 1)
{
    if(0 == b)
    {
        int indx = 0;
        for(int i = 0;i < N;i++)
        {
            A[indx] = 0;
            indx += incx;
        }
    }
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
template<class T> void Copy(T *Y, T *X, int N, int incy = 1, int incx = 1)
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
template<class T> T VNorm2(T *A, int N, int incx = 1)
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
//NormInf of a vector
template<class T> T VNormInf(T *A, int N, int incx = 1)
{
    T result = 0;
    for(int i = 0;i < N;i++)
    {
        if(abs(A[i*incx]) > result)
            result = abs(A[i*incx]);
    }
    return result;
}
//NormOne of a vector
template<class T> T VNormOne(T *A, int N, int incx = 1)
{
    T result = 0;
    for(int i = 0;i < N;i++)
        result += abs(A[i*incx]);
    return result;
}
//sum of absolute value for a vector
template<class T> T Asum(T *A, int N, int incx = 1)
{
	T tmp = 0;
	for(int i = 0;i < N;i++)
		tmp += abs(A[i*incx]);
	return tmp;
}
//find max abs place 
template<class T> T* MaxAbs(T *A, int N, int incx = 1)
{
	if(N < 1)
		return NULL;
	T tmp = abs(A[0]);
	int index = 0;
	for(int i = 0;i < N;i++)
	{
		if(abs(A[i*incx]) > tmp)
		{
			tmp = abs(A[i*incx]);
			index = i;
		}
	}
	return A + index*incx;
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
//End of Level 1
//Start of Level 2, Matrix Vector Operation

//MatrixVectorMultiply
//y = alpha*A*x + beta*y
//Input:r:row,c:column,lda:leading dimension of A
template<class T>void MatVecMul(int r, int c, T alpha, T *A, int lda, T *x, int incx, double beta, T *y, int incy)
{
    //Try quick return
    if(0 == alpha)
    {
        Scale(y, beta, r, incy);
    }
    int indy = 0;
    for(int i = 0;i < r; i++)
    {
        T tmp = cDot(alpha, A + i*lda, x, c, 1, incy);
        y[indy] = beta*y[indy] + tmp;
        indy += incy;
    }
}
//Many other subroutines are ignored because they are suited for special matrix
//End of Level 2
//Start of Level 3, Matrix Matrix Operation
//C = alpha*A*B + beta*C
//A is M*K, B is K*N, C is M*N
//Input:M, N, K
template<class T>void MatMatMul(int M, int N, int K, T alpha, T *A, int lda, T *B, int ldb, T beta, T *C, int ldc)
{
    //Quick Return if possible
    if(0 == alpha && 1 == beta)
        return;
    if(0 == alpha)
    {
        if(0 == beta)
        {
            int indc = 0;
            for(int i = 0;i < M;i++)
            {
                for(int j = 0;j < N;j++)
                {
                    C[indc + j] = 0;
                }
                indc += ldc;
            }
        }
        else
        {
            int indc = 0;
            for(int i = 0;i < M;i++)
            {
                for(int j = 0;j < N;j++)
                {
                    C[indc + j] *= beta;
                }
                indc += ldc;
            }
        }
    }
    //if no quick return, start operation
    //Do Y=alpha*A*X + beta*Y N times
    for(int i = 0;i < N;i++)
    {
        MatVecMul(M, K, alpha, A, lda, B + i, ldb, beta, C + i, ldc);
    }
}
//MatrixMatrixMultiply
//End of namespace
}
#endif
