#ifndef _ODE45_H_
#define _ODE45_H_
#include<math.h>
#include<iostream>
#include<stdio.h>

int ode45nw(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, const double* AbsTol, double RelTol=1.0e-3,
		  int NormControl=0, double MaxStep=-1, double InitialStep=-1, FILE* fid=NULL);

int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol=1.0e-3,
		  int NormControl=0, double MaxStep=-1, double InitialStep=-1, FILE* fid=NULL);
int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), double* x, const double* para,
		  double t0, double tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol,
		  int MaxIter, int NormControl=0, double MaxStep=-1, double InitialStep=-1, FILE* fid=NULL);
		  
		  
int ode45(int (*fun)(double t, const double* x, double* dx, const double* para), int (*eventfun)(double t, const double *x, const double *eventPara, double &fval, double &fgrad), double* x, const double* para,
 const double *eventPara, const double eventTol, int &eventFlag, 
		  double t0, double &tf, int n, int& NumPoint, double* work, const double* AbsTol, double RelTol,
		  int MaxIter, int NormControl, double MaxStep, double InitialStep, FILE* fid);
//符号函数
template<class T> inline int sign(const T & x)
{
	if(x>0)
		return 1;
	else if(x<0)
		return -1;
	else
		return 0;
}

//求最大值
template <class T>
inline T max(T x, T y) 
{ 
	return (x>y)?x:y;
}


//求最小值
template <class T>
inline T min(T x, T y) 
{
	return (x<y)?x:y;
}

#endif
