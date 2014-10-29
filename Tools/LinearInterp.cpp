/* Dev by tiget at 2014/10/29
 * Program for LinearInterp in an sorted array
 */
#include <iostream>
using namespace std;
//subroutime to return nearest position from an sorted array
//A is sorted array
//x is point to be search
//len is length of A
//dir = true if array is sorted in ascending order
//return value and value + 1 are two points for next interp
//return -1 if error happens:length too small, x less than min or larger than max of A
template<class T> int Nearest(const T *A, const T x, const int len, bool dir = true)
{
    if(len <= 1)
        return -1;
    if((true == dir) && (x < A[0] || x > A[len - 1]))
        return -1;
    if((false == dir) && (x > A[0] || x < A[len - 1]))
        return -1;
    int indexlow = 0, indexhigh = len - 1, index;
    while(indexhigh - indexlow > 1)
    {
        index = (indexlow + indexhigh)/2;
        if(A[index] == x)
            return index;
        if(A[index] < x)
        {
            if(true == dir)
                indexlow = index;
            else
                indexhigh = index;
        }
        else
        {
            if(true == dir)
                indexhigh = index;
            else
                indexlow = index;
        }
    }
}
//Interp using index A B and b
//A:Array to determine point of A
//B:Array to perform interp on
//x:point
//len:length of array A
double LinearInterp(const double *A, const double *B, const double x, const int len, bool dir = true)
{
    int index = Nearest(A, x, len, dir);
    if(-1 == index)
    {
        cout << "Linear Interpolation may be wrong" << endl;
        return -1;
    }
    if(A[index] == x)
        return B[index];
    double xlow = A[index], xupp = A[index + 1];
    double ylow = B[index], yupp = B[index + 1];
    double factor = (x - xlow)/(xupp - xlow);
    return ylow + factor*(yupp - ylow);
}
int main()
{
    const int len = 100;
    double testx[len], testy[len];
    for(int i = 0;i < len;i++)
    {
        testx[i] = 1 - (double)i/len;
        testy[i] = (double)i/len + 1;
    }
    double point[5] = {0, 0.99, 0.5, -0.3, 1.4};
    for(int i = 0;i < 5;i++)
        cout << point[i] << ":" << LinearInterp(testx, testy, point[i], len, false) << endl;
    return 0;
}
