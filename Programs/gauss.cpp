/* Dev by tiger at 2014/10/21
 * Subroutine to perform Gauss elimination method
 */
#include <iostream>
using namespace std;
/* complex interface for Gauss elimination for sloving Ax=b
 * nA is dimension of A, nb is colomn of b
 * A is pointer to first element of matrix A, lda is leading column of A
 * b is pointer to first element of b, ldb is leading column of b
 */
int Gauss(int nA, int nb, double *A, int lda, double *b, int ldb)
{
#define A(i, j) A[(i - 1)*lda + j - 1]
#define b(i, j) b[(i - 1)*ldb + j - 1]
    //i is row to be multiply , perhaps some subroutines should be added to make code shorter
    for(int i = 1;i < nA;i++)
    {
        //j is the row to be add
        for(int j = i + 1;j <= nA;j++)
        {
            double Lij = A(j, i)/A(i, i);
            A(j, i) = 0.0;
            //k is the number (j to nA)
            for(int k = i + 1;k <= nA;k++)
                A(j, k) -= Lij*A(i, k);
            //change b
            for(int k = 1;k <= nb;k++)
                b(j, k) -= Lij*b(i, k);
        }
    }
    //U got , get x now, make A unit matrix and b should be the answer
    //from lowest to highest
    //i is the row to multiply
    for(int i = nA;i >= 1;i--)
    {
        //make A(i,i) to be one
        for(int k = 1;k <= nb;k++)
            b(i, k) /= A(i, i);
        A(i, i) = 1;
        //j is current row to make one
        for(int j = i - 1;j >= 1;j--)
        {
            //k is current number(nA to i - 1)
            for(int k = 1;k <= nb;k++)
                b(j, k) -= A(j, i)*b(i, k);
            A(j, i) = 0.0;
        }
    }
#undef A
#undef b
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
    int flag = Gauss(n, n, A, n, b, n);
    //output b
    for(int i = 0;i < n*n;i++)
        cout << b[i] << endl;
    double C[9] = {1, 2, 3, 2, 3, 1, 3, 2,1};
    double d[3] = {6,6,6};
    flag = Gauss(n, C, d);
    for(int i = 0;i < n;i++)
        cout << d[i] << endl;
    return 0;

}
