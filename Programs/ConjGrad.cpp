/* Dev by tiger at 2014/11/25
 * subroutine for Conjugate gradient algorithm
 */
#include "VectorMatrix.h"
#include <iostream>
using namespace TigerVecMat;
using namespace std;
const int MAXCGITER = 100;
//N:dimension of A
//x:In : initial guess Out: Final result
int ConjGrad(int N, double *A, int lda, double *b, int incb, double *x, int incx, double tol)
{
    double *r = new double[N];
    double *p = new double[N];
    double *Ap = new double[N];
    if(NULL == r || NULL == p)
        cout << "Allocate error\n";
    //init r(0) and p(0)
    //r=b-Ax r=-Ax+0*r, r=b+r
    //Scale(r, 0, N, 1);
    MatVecMul(N, N, -1.0, A, lda, x, incx, 0.0, r, 1);
    aXpY(1.0, b, r, N, incb, 1);
    //copy to p
    Copy(p, r, N, 1, 1);
    int iter = 1;
    while(iter < MAXCGITER)
    {
        //calc alphak
        double AlphaTop = cDot(1.0, r, r, N, 1, 1);
        //Why not use x to store
        MatVecMul(N, N, 1.0, A, lda, p, 1, 0.0, Ap, 1);
        double AlphaBelow = cDot(1.0, p, Ap, N, 1, 1);
        //if (p,Ap) < tol, converge already
        if(sqrt(AlphaBelow) < tol*tol)
            break;
        double alpha = AlphaTop/AlphaBelow;
        //Update x and r
        aXpY(alpha, p, x, N, 1, incx);
        double BetaBelow = cDot(1.0, r, r, N, 1);
        if(sqrt(BetaBelow) < tol)
            break;
        MatVecMul(N, N, -alpha, A, lda, p, 1, 1.0, r, 1);
        //Output(r, N, 1, 1);
        //calc beta
        double BetaTop = cDot(1.0, r, r, N, 1);
        double beta = BetaTop/BetaBelow;
        //update p
        aXbYc(1.0, r, beta, p, 0.0, N, 1, 1);
        iter++;
    }

    delete[] r;
    delete[] p;
    delete[] Ap;
    if(iter == MAXCGITER)
        return -1;
    return iter;
}
int main()
{
    
    int N = 4;
    double A[16] = {1, -0.30009, 0, -0.30898, -0.30009, 1, -0.46691, 0,
        0, -0.46691, 0, -0.27471, -0.30898, 0, -0.27471, 1};
    double b[4] = {5.32088, 6.07624, -8.80455, 2.676};
    /*
    int N = 2;
    double A[4] = {3, 1, 1, 2};
    double b[2] = {5, 5};
    */
    double x[4] = {0,0,0,0}, tol = 1e-6;
    int iter = ConjGrad(N, A, N, b, 1, x, 1, tol);
    cout << "After " << iter << " times:\n";
    Output(x, N, 1, 1);
    return 0;

}
